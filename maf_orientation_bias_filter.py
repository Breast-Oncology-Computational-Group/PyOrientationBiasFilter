#!/usr/bin/env python3
"""MAF front-end for the orientation-bias-variant-filter reimplementation.

Drop-in replacement for the compiled MATLAB ``orientationBiasFilter`` step: it
consumes the same ``<pairName>.OrientationBiasInfo.maf`` that the existing
pipeline python scripts prepare (run_local.sh steps 1-9) and writes the same
output files the downstream (step 11: tsvConcatFiles + add_judgement_column)
expects:

  * ``<output>``                          - filtered MAF (passing rows, cut == 0)
  * ``<output w/o .maf>.unfiltered.maf``   - all rows, with the computed columns
  * ``<output>.pass_count.txt``            - number passing  (no trailing newline)
  * ``<output>.reject_count.txt``          - number rejected (no trailing newline)

The heavy lifting is the faithful Java-port core (``OrientationBiasFilterer`` /
``ArtifactStatisticsScorer``). This module only maps MAF columns <-> the core's
model and re-emits the MAF, so the priority is consistency with the pipeline
(NOT numeric parity with the old MATLAB).

Computed columns appended (named in the pipeline's ``i_<stub>_*`` convention):
  * ``i_<stub>_F``        - fraction of alt reads in the artifact orientation (Java OBF/fob)
  * ``i_<stub>_mode``     - 1 if the variant is in the artifact mode or its complement, else 0
  * ``i_<stub>_p_value``  - orientation-bias artifact p-value (Java OBP)
  * ``i_<stub>_cut``      - 1 if cut as an artifact, else 0   (the column downstream keys on)

NOTE: MATLAB also emitted ``i_<stub>_q_value`` and ``i_<stub>_p_value_cutoff``.
The Java algorithm decides cuts by a Benjamini-Hochberg *count* (cut the top-N by
p-value) rather than a per-variant cutoff, so those two columns are intentionally
not produced. Everything downstream only needs ``i_<stub>_cut``.

Inputs mirror what run_local.sh passed to MATLAB: the reference/artifact alleles
are the *complement* bases read from the part1 files (``--reference-allele
$REF_ALLELE_COMP --artifact-allele $ARTIFACT_ALLELE_COMP``), and the numeric
parameters default to run_local's values (fdr/artifactThresholdRate=0.01,
biasQP1=30, biasQP2=1.5; BIAS_P is fixed at 0.96 in the core, matching pBias).

MEMORY MODEL (why this is streaming + key-joined)
--------------------------------------------------
A whole-genome OrientationBiasInfo.maf can be millions of rows x ~337 columns.
Parsing the entire file into a list-of-lists of strings and building one model
object per row would peak at tens of GB (the ~332 annotation columns are pure
passthrough - the algorithm only reads ~5 of them). Instead we do TWO streaming
passes and never hold the wide rows:

  Pass 1  read only the columns the core needs; build a transient VariantContext
          ONLY for the minority of artifact-mode rows, run the (unchanged) core
          per-variant annotation, and keep just small per-mode-row results plus
          the scalars the Benjamini-Hochberg cut needs. Then decide the cut set.
  Pass 2  re-read the file line by line; append the four computed columns to each
          raw line and write the filtered / unfiltered MAFs. O(1) rows in memory.

Rows are matched between the two passes by a VARIANT KEY (Chromosome,
Start_position, End_position, Reference_Allele, Tumor_Seq_Allele1,
Tumor_Seq_Allele2 - the same key the run_local.sh graft join uses), not by
position, so the result does not depend on the two passes iterating in lock-step.
Peak memory is proportional to the number of artifact-mode rows, not the file
size, so a multi-million-row whole-genome MAF runs in well under a couple of GB.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from typing import Dict, List, Optional, Tuple

# Allow running as a plain script (``python maf_orientation_bias_filter.py``)
# as the pipeline/WDL would: import sibling modules directly from this folder.
# We still prepend the parent dir to sys.path so the script works when invoked
# by path from outside the repository root.
_PARENT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _PARENT not in sys.path:
    sys.path.insert(0, _PARENT)

from artifact_statistics_scorer import ArtifactStatisticsScorer  # noqa: E402
from htsjdk_models import (  # noqa: E402
    Allele,
    GATKVCFConstants,
    Genotype,
    GenotypesContext,
    VariantContext,
)
from orientation_bias_filter_constants import OrientationBiasFilterConstants  # noqa: E402
from orientation_bias_filterer import OrientationBiasFilterer, _java_round  # noqa: E402
from orientation_bias_utils import OrientationBiasUtils  # noqa: E402
from transition import _COMPLEMENT, Transition  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger("maf_orientation_bias_filter")


# The identifying key for a variant row. Same columns the run_local.sh graft join
# uses (MAF cols 5,6,7,11,12,13). Used to match Pass-1 results to Pass-2 rows so
# we never depend on the two passes reading in the exact same order.
KEY_COLUMNS = [
    "Chromosome",
    "Start_position",
    "End_position",
    "Reference_Allele",
    "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2",
]

# The four columns appended to every row, in MATLAB's trailing order (cut last so
# it stays the final column, which downstream keys on).
_EXTRA_SUFFIXES = ["F", "mode", "p_value", "cut"]

# MAF I/O encoding. Annotation columns can carry non-UTF-8 bytes (e.g. accented
# characters from external annotation sources), so we read AND write with latin-1:
# it maps every byte 0x00-0xFF one-to-one, never raises on decode, and round-trips
# the passthrough columns byte-for-byte (our appended columns are pure ASCII).
_ENCODING = "latin-1"


# ---------------------------------------------------------------------------
# Small field helpers (a MAF here == leading '#' comment lines, one header line,
# then tab-separated data rows).
# ---------------------------------------------------------------------------
def _column_index(header: List[str]) -> Dict[str, int]:
    """Map column name -> position; keep the first occurrence of a duplicate name."""
    idx: Dict[str, int] = {}
    for i, name in enumerate(header):
        idx.setdefault(name, i)
    return idx


def _to_int(value: str) -> int:
    """Parse an int from a MAF cell, treating blanks / '.' as 0."""
    value = (value or "").strip()
    if value in ("", "."):
        return 0
    try:
        return int(float(value))   # tolerate "12.0"
    except ValueError:
        return 0


def _field(row: List[str], idx: Dict[str, int], name: str, default: str = "") -> str:
    """Safely fetch a named field from a row (default if column absent/short)."""
    i = idx.get(name)
    if i is None or i >= len(row):
        return default
    return row[i]


def _computed_value(value) -> str:
    """Render a FORMAT attribute for the MAF: blank for missing / EMPTY_ALLELE."""
    if value is None or value == "." or value == "":
        return ""
    return str(value)


def _variant_key(row: List[str], idx: Dict[str, int]) -> str:
    """Build the join key for a row from the identifying columns."""
    return "\t".join(_field(row, idx, c) for c in KEY_COLUMNS)


# ---------------------------------------------------------------------------
# MAF reading (streaming): header once, then data rows.
# ---------------------------------------------------------------------------
def _read_header(path: str) -> Tuple[List[str], List[str]]:
    """Return ``(comment_lines, header_fields)`` - stops right after the header."""
    comments: List[str] = []
    with open(path, encoding=_ENCODING, newline="") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("#"):
                comments.append(line)
            else:
                return comments, line.split("\t")
    raise ValueError(f"MAF {path} has no header line")


def _iter_data_rows(path: str):
    """Yield ``(raw_line_without_newline, fields)`` for each data row of the MAF."""
    header_seen = False
    with open(path, encoding=_ENCODING, newline="") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("#"):
                continue
            if not header_seen:
                header_seen = True          # first non-comment line is the header
                continue
            yield line, line.split("\t")


# ---------------------------------------------------------------------------
# MAF row  ->  core model (single row)
# ---------------------------------------------------------------------------
def build_variant_context(row: List[str], idx: Dict[str, int], sample_name: str) -> VariantContext:
    """Turn one MAF data row into a single-sample VariantContext for the core.

    The orientation-read counts (``i_t_ALT_F1R2`` / ``i_t_ALT_F2R1``) become the
    FORMAT fields the core reads (``ALT_F1R2`` / ``ALT_F2R1``); ``t_alt_count``
    becomes ``AD[1]``; ref/alt alleles come from ``Reference_Allele`` /
    ``Tumor_Seq_Allele2``.
    """
    ref = _field(row, idx, "Reference_Allele")
    alt = _field(row, idx, "Tumor_Seq_Allele2")
    contig = _field(row, idx, "Chromosome")
    start = _to_int(_field(row, idx, "Start_position"))
    end = _to_int(_field(row, idx, "End_position"))

    # Orientation-split alt-read counts -> the FORMAT fields the core reads.
    attributes = {
        GATKVCFConstants.OXOG_ALT_F1R2_KEY: _to_int(_field(row, idx, "i_t_ALT_F1R2")),
        GATKVCFConstants.OXOG_ALT_F2R1_KEY: _to_int(_field(row, idx, "i_t_ALT_F2R1")),
    }
    # AD = [ref_count, alt_count]; the core uses AD[1] as the total alt count.
    allele_depths = [_to_int(_field(row, idx, "t_ref_count")), _to_int(_field(row, idx, "t_alt_count"))]

    genotype = Genotype(
        sample_name,
        [Allele(ref, is_ref=True), Allele(alt)],
        attributes=attributes,
        ad=allele_depths,
    )
    return VariantContext(contig, start, end, GenotypesContext([genotype]))


def _snv_mode_transition(ref: str, alt: str, relevant_modes) -> Optional[Transition]:
    """Return the row's artifact-mode Transition if it is a relevant SNV, else None.

    Mirrors exactly the SNV gate inside
    ``OrientationBiasFilterer.annotate_variant_context_with_preprocessing_values``:
    a single-base ref, a single-base called non-'*' alt, both A/C/G/T. Only rows
    whose (ref->alt) transition is in ``relevant_modes`` (the artifact mode OR its
    complement) are artifact-mode rows; everything else gets the blank/0 columns.

    (The object path would raise on a non-ACGT single base here; we skip such rows
    instead - they cannot be one of the DNA artifact modes anyway.)
    """
    if len(ref) != 1 or len(alt) != 1:
        return None
    if alt in (".", "", "*"):           # uncalled / spanning-deletion allele
        return None
    ref_base = ref.upper()
    alt_base = alt.upper()
    if ref_base not in _COMPLEMENT or alt_base not in _COMPLEMENT:
        return None
    transition = Transition(ref_base, alt_base)
    return transition if transition in relevant_modes else None


# ---------------------------------------------------------------------------
# Pass 1: stream, annotate the artifact-mode rows, collect cut inputs.
# ---------------------------------------------------------------------------
class _Pass1Result:
    """The small state Pass 1 carries into the cut decision and Pass 2."""

    __slots__ = ("computed_by_key", "candidates", "unfiltered_genotype_count",
                 "pre_adapter_q", "n_data_rows")

    def __init__(self):
        # Artifact-mode rows only: key -> (F_string, p_value_string) for Pass 2.
        self.computed_by_key: Dict[str, Tuple[str, str]] = {}
        # One (p_value, order_index, key, transition) per artifact-mode candidate.
        self.candidates: List[Tuple[float, int, str, Transition]] = []
        # Denominator of the Benjamini-Hochberg calc (unfiltered non-ref genotypes).
        self.unfiltered_genotype_count = 0
        self.pre_adapter_q: Optional[float] = None
        self.n_data_rows = 0


def _run_pass1(
    input_maf: str,
    idx: Dict[str, int],
    stub: str,
    sample_name: str,
    relevant_transition: Transition,
    relevant_modes,
) -> _Pass1Result:
    """First streaming pass: compute per-variant values for artifact-mode rows."""
    result = _Pass1Result()
    q_column = f"i_{stub}_Q"
    relevant_without_complement = {relevant_transition}
    q_score_map: Dict[Transition, float] = {}

    fob_field = OrientationBiasFilterConstants.FOB
    p_field = OrientationBiasFilterConstants.P_ARTIFACT_FIELD_NAME

    for _raw, row in _iter_data_rows(input_maf):
        result.n_data_rows += 1

        # preAdapterQ is the same on every row (appended by run_local step 6); take
        # the first parseable value and build the score map for the annotation.
        if result.pre_adapter_q is None:
            raw_q = _field(row, idx, q_column).strip()
            if raw_q not in ("", "."):
                result.pre_adapter_q = float(raw_q)
                q_score_map = {relevant_transition: result.pre_adapter_q}

        ref = _field(row, idx, "Reference_Allele")
        alt = _field(row, idx, "Tumor_Seq_Allele2")

        # FDR denominator: unfiltered, non-ref genotypes (here: alt base != ref base;
        # inputs carry no genotype-level filter). Matches
        # OrientationBiasUtils.calculate_unfiltered_non_ref_genotype_count.
        if ref.upper() != alt.upper():
            result.unfiltered_genotype_count += 1

        mode_transition = _snv_mode_transition(ref, alt, relevant_modes)
        if mode_transition is None:
            continue                        # non-artifact-mode row -> blank columns in Pass 2

        if result.pre_adapter_q is None:
            raise ValueError(
                f"Hit an artifact-mode row before any '{q_column}' value was seen; "
                "cannot annotate without preAdapterQ."
            )

        # Build ONE transient VariantContext and run the unchanged core annotation,
        # so F / p-value are byte-identical to the object-graph implementation. The
        # object is discarded after we read the two values off it.
        vc = build_variant_context(row, idx, sample_name)
        annotated = OrientationBiasFilterer.annotate_variant_context_with_preprocessing_values(
            vc, relevant_without_complement, q_score_map
        )
        g = annotated.get_genotype(sample_name)

        f_string = _computed_value(g.get_extended_attribute(fob_field))
        p_string = _computed_value(g.get_extended_attribute(p_field))
        # p-value as the core reads it for the BH sort (empty -> 0.0).
        p_value = OrientationBiasUtils.get_genotype_double(g, p_field, 0.0)

        key = _variant_key(row, idx)
        if key in result.computed_by_key:
            raise ValueError(
                f"Duplicate variant key among artifact-mode rows: {key!r}. The join "
                "key must be unique; cannot map computed values back safely."
            )
        result.computed_by_key[key] = (f_string, p_string)
        result.candidates.append((p_value, len(result.candidates), key, mode_transition))

    return result


# ---------------------------------------------------------------------------
# Cut decision: Benjamini-Hochberg over the artifact-mode candidates.
# ---------------------------------------------------------------------------
def _decide_cuts(
    pass1: _Pass1Result,
    relevant_transition: Transition,
    fdr_threshold: float,
    bias_qp1: float,
    bias_qp2: float,
) -> set:
    """Return the set of variant keys to cut.

    Faithful, single-sample reimplementation of
    ``OrientationBiasFilterer.annotate_variant_contexts_with_filter_results``:
    sort candidates by descending p-value, use the (reused) BH threshold to get the
    total number to cut, split it across the mode and its complement in proportion
    to their candidate counts, scale each by the preAdapterQ suppression factor,
    then cut the top-N of each transition. The reused helpers guarantee the counts
    match the object path; the only difference is that ties at equal p-value are
    broken by file order (Pass-1 index) instead of the Java object hash, which can
    only change *which* equal-p variant is cut at a quota boundary, never how many.
    """
    complement_transition = relevant_transition.complement()
    if not pass1.candidates:
        return set()

    # Sort by descending p-value; deterministic tie-break by encounter order.
    ordered = sorted(pass1.candidates, key=lambda c: (-c[0], c[1]))

    # Candidate count per transition (mode / complement).
    transition_count: Dict[Transition, int] = {relevant_transition: 0, complement_transition: 0}
    for _p, _i, _key, transition in ordered:
        transition_count[transition] += 1
    all_transition_count = sum(transition_count.values())

    # Benjamini-Hochberg total-to-cut over [candidate p-values desc] + zero-padding
    # for the non-artifact genotypes (reuses the vetted core routine).
    scores = [c[0] for c in ordered]
    num_to_pad = pass1.unfiltered_genotype_count - len(scores)
    if num_to_pad > 0:
        scores = scores + [0.0] * num_to_pad
    total_num_to_cut = OrientationBiasFilterer.calculate_total_num_to_cut(
        fdr_threshold, pass1.unfiltered_genotype_count, scores
    )
    logger.info("Cutting (total) pre-preAdapterQ: %s", total_num_to_cut)

    # Split proportionally (integer arithmetic, as in the Java), then scale by the
    # preAdapterQ suppression factor. Both the mode and its complement suppress by
    # the mode's Q (the score map only carries the mode).
    suppression = ArtifactStatisticsScorer.calculate_suppression_factor_from_pre_adapter_q(
        pass1.pre_adapter_q, bias_qp1, bias_qp2
    )
    num_to_cut: Dict[Transition, int] = {}
    for transition in (relevant_transition, complement_transition):
        pre = 0 if all_transition_count == 0 else (total_num_to_cut * transition_count[transition]) // all_transition_count
        num_to_cut[transition] = _java_round(pre * suppression)
        logger.info("Cutting (%s) post-preAdapterQ: %s", transition, num_to_cut[transition])

    # Walk candidates in descending-p order, cutting the first N of each transition.
    cut_keys: set = set()
    cut_so_far: Dict[Transition, int] = {relevant_transition: 0, complement_transition: 0}
    for _p, _i, key, transition in ordered:
        if cut_so_far[transition] < num_to_cut[transition]:
            cut_keys.add(key)
            cut_so_far[transition] += 1
    return cut_keys


# ---------------------------------------------------------------------------
# Pass 2: stream again, append columns, write filtered / unfiltered MAFs.
# ---------------------------------------------------------------------------
def _run_pass2(
    input_maf: str,
    output_maf: str,
    unfiltered_maf: str,
    comments: List[str],
    header: List[str],
    idx: Dict[str, int],
    extra_col_names: List[str],
    pass1: _Pass1Result,
    cut_keys: set,
) -> Tuple[int, int]:
    """Second streaming pass. Returns ``(n_pass, n_total)``."""
    # MATLAB stamps a second comment line onto its output; mimic it for parity.
    out_comments = list(comments) + ["## Orientation Bias Filter "]
    header_line = "\t".join(header + extra_col_names) + "\n"
    expected_width = len(header)

    # Track that every artifact-mode key is consumed exactly once (catches a
    # duplicate/vanished key between the two passes).
    match_count: Dict[str, int] = {k: 0 for k in pass1.computed_by_key}

    n_pass = 0
    n_total = 0
    n_data_rows = 0
    width_mismatches = 0

    with open(unfiltered_maf, "w", encoding=_ENCODING, newline="") as unfiltered_out, \
            open(output_maf, "w", encoding=_ENCODING, newline="") as filtered_out:
        for comment in out_comments:
            unfiltered_out.write(comment + "\n")
            filtered_out.write(comment + "\n")
        unfiltered_out.write(header_line)
        filtered_out.write(header_line)

        for raw, row in _iter_data_rows(input_maf):
            n_data_rows += 1
            if len(row) != expected_width:
                width_mismatches += 1

            key = _variant_key(row, idx)
            computed = pass1.computed_by_key.get(key)
            if computed is not None:
                f_string, p_string = computed
                mode = "1"
                cut = "1" if key in cut_keys else "0"
                match_count[key] += 1
            else:
                f_string, mode, p_string, cut = "", "0", "", "0"

            out_line = raw + "\t" + "\t".join([f_string, mode, p_string, cut]) + "\n"
            unfiltered_out.write(out_line)
            n_total += 1
            if cut == "0":
                filtered_out.write(out_line)
                n_pass += 1

    # Guards: the two passes must agree on the row set.
    if n_data_rows != pass1.n_data_rows:
        raise ValueError(
            f"Row count changed between passes ({pass1.n_data_rows} then {n_data_rows}); "
            "the input MAF must not change mid-run."
        )
    bad_keys = [k for k, c in match_count.items() if c != 1]
    if bad_keys:
        raise ValueError(
            f"{len(bad_keys)} artifact-mode key(s) were not matched exactly once in "
            f"Pass 2 (e.g. {bad_keys[0]!r}); duplicate or missing variant keys."
        )
    if width_mismatches:
        logger.warning("%d row(s) did not have the expected %d columns.",
                       width_mismatches, expected_width)

    return n_pass, n_total


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def run(
    input_maf: str,
    output_maf: str,
    stub: str,
    reference_allele: str,
    artifact_allele: str,
    fdr_threshold: float,
    bias_qp1: float,
    bias_qp2: float,
    sample_name: str,
) -> None:
    comments, header = _read_header(input_maf)
    idx = _column_index(header)

    # Fail fast if the MAF is missing a column the algorithm truly needs.
    required = ["Chromosome", "Start_position", "End_position", "Reference_Allele",
                "Tumor_Seq_Allele2", "i_t_ALT_F1R2", "i_t_ALT_F2R1", "t_alt_count",
                f"i_{stub}_Q"]
    missing = [c for c in required if c not in idx]
    if missing:
        raise ValueError(f"OrientationBiasInfo MAF is missing required column(s): {missing}")

    # The artifact mode passed to MATLAB was (reference_allele -> artifact_allele);
    # the core internally also considers its reverse complement.
    relevant_transition = Transition.transition_of(reference_allele, artifact_allele)
    relevant_modes = {relevant_transition, relevant_transition.complement()}

    # ----- Pass 1: annotate artifact-mode rows, gather cut inputs -----
    pass1 = _run_pass1(input_maf, idx, stub, sample_name, relevant_transition, relevant_modes)
    logger.info(
        "Read %d variants from %s (%d artifact-mode candidates, %d unfiltered non-ref)",
        pass1.n_data_rows, input_maf, len(pass1.candidates), pass1.unfiltered_genotype_count,
    )
    if pass1.pre_adapter_q is None:
        raise ValueError(f"No usable value found in column 'i_{stub}_Q'")
    logger.info("Artifact mode %s (complement %s), preAdapterQ=%s, fdr=%s",
                relevant_transition, relevant_transition.complement(), pass1.pre_adapter_q, fdr_threshold)

    # ----- Cut decision (Benjamini-Hochberg over the candidates) -----
    cut_keys = _decide_cuts(pass1, relevant_transition, fdr_threshold, bias_qp1, bias_qp2)

    # ----- Pass 2: append columns and write the MAFs -----
    extra_col_names = [f"i_{stub}_{suffix}" for suffix in _EXTRA_SUFFIXES]

    base_no_ext = output_maf[:-4] if output_maf.endswith(".maf") else output_maf
    unfiltered_maf = base_no_ext + ".unfiltered.maf"
    pass_count_file = output_maf + ".pass_count.txt"
    reject_count_file = output_maf + ".reject_count.txt"

    n_pass, n_total = _run_pass2(
        input_maf, output_maf, unfiltered_maf, comments, header, idx,
        extra_col_names, pass1, cut_keys,
    )
    n_reject = n_total - n_pass

    # Counts, written with NO trailing newline to match MATLAB's format.
    with open(pass_count_file, "w") as handle:
        handle.write(str(n_pass))
    with open(reject_count_file, "w") as handle:
        handle.write(str(n_reject))

    logger.info("Wrote %s (%d passing) and %s (%d total); rejected=%d",
                output_maf, n_pass, unfiltered_maf, n_total, n_reject)


def parse_options(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the orientation-bias filter (Java-port reimplementation) on an "
                    "OrientationBiasInfo MAF; drop-in replacement for the MATLAB step."
    )
    parser.add_argument("-i", "--input", required=True, help="Input <pairName>.OrientationBiasInfo.maf")
    parser.add_argument("-o", "--output", required=True,
                        help="Output filtered MAF name, e.g. <pairName>.OrientationBiasFilter.maf "
                             "(the .unfiltered.maf and count files are derived from it)")
    parser.add_argument("-s", "--stub", required=True, help="Artifact stub, e.g. 'oxog' or 'ffpe'")
    parser.add_argument("--reference-allele", required=True,
                        help="Artifact-mode reference base (run_local passes $REF_ALLELE_COMP)")
    parser.add_argument("--artifact-allele", required=True,
                        help="Artifact-mode alt base (run_local passes $ARTIFACT_ALLELE_COMP)")
    parser.add_argument("--fdr-threshold", type=float, default=0.01,
                        help="Max FDR (MATLAB artifactThresholdRate; default 0.01)")
    parser.add_argument("--bias-qp1", type=float, default=30.0,
                        help="preAdapterQ sigmoid inflection point (run_local used 30; default 30)")
    parser.add_argument("--bias-qp2", type=float, default=1.5,
                        help="preAdapterQ sigmoid steepness (default 1.5)")
    parser.add_argument("--sample-name", default="TUMOR",
                        help="Sample name to attach to genotypes (internal; default 'TUMOR')")
    return parser.parse_args(argv)


def main(argv=None) -> None:
    args = parse_options(argv)
    run(
        input_maf=args.input,
        output_maf=args.output,
        stub=args.stub,
        reference_allele=args.reference_allele,
        artifact_allele=args.artifact_allele,
        fdr_threshold=args.fdr_threshold,
        bias_qp1=args.bias_qp1,
        bias_qp2=args.bias_qp2,
        sample_name=args.sample_name,
    )


if __name__ == "__main__":
    main()
