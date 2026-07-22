"""Translation of ``OrientationBiasFilterer.java``.

The heart of the filter. Two public entry points:
  * :meth:`annotate_variant_context_with_preprocessing_values` - annotate a single
    site's genotypes with the orientation-bias FORMAT fields (OBAM/OBP/OBF/...).
  * :meth:`annotate_variant_contexts_with_filter_results` - given all annotated
    variants, decide which to cut so the sample FDR stays below a threshold
    (Benjamini-Hochberg), scaled by the bam-level preAdapterQ suppression factor.
"""

from __future__ import annotations

import logging
import math
from typing import Dict, List, Optional

from artifact_statistics_scorer import ArtifactStatisticsScorer
from htsjdk_models import (
    Allele,
    GATKVCFConstants,
    Genotype,
    GenotypeBuilder,
    GenotypesContext,
    VariantContext,
    VariantContextBuilder,
    VCFConstants,
)
from orientation_bias_filter_constants import OrientationBiasFilterConstants
from orientation_bias_utils import OrientationBiasUtils
from transition import Transition

logger = logging.getLogger(__name__)


def _java_round(x: float) -> int:
    """Reproduce ``Math.round(double)``: floor(x + 0.5), and 0 for NaN.

    Python's built-in ``round`` uses banker's rounding, which differs on .5 ties,
    so we replicate Java's half-up behaviour to keep the cut counts identical.
    """
    if math.isnan(x):
        return 0
    return math.floor(x + 0.5)


def _java_boolean_string(value: bool) -> str:
    """``String.valueOf(boolean)`` -> "true"/"false" (lower-case, unlike Python)."""
    return "true" if value else "false"


class OrientationBiasFilterer:
    """Operations for the orientation bias filter."""

    # preAdapterQ value assigned to non-artifact modes (an effectively "clean" score).
    PRE_ADAPTER_METRIC_NOT_ARTIFACT_SCORE = 100.0

    # Mode of the binomial distribution of pair orientation for artifact alt reads.
    BIAS_P = 0.96

    def __init__(self):
        raise RuntimeError("OrientationBiasFilterer is not instantiable")

    # ==================================================================
    # Step 1: per-variant annotation
    # ==================================================================
    @staticmethod
    def annotate_variant_context_with_preprocessing_values(
        vc: VariantContext,
        relevant_transitions_without_complement,   # SortedSet[Transition]
        pre_adapter_q_score_map: Dict[Transition, float],
    ) -> VariantContext:
        """Add the FORMAT annotations used by the filter to each genotype of ``vc``.

        Java: ``annotateVariantContextWithPreprocessingValues``.
        """
        if vc is None or relevant_transitions_without_complement is None or pre_adapter_q_score_map is None:
            raise ValueError("Arguments cannot be null")

        relevant_transitions_complement = OrientationBiasUtils.create_reverse_complement_transitions(
            relevant_transitions_without_complement
        )

        vcb = VariantContextBuilder(vc)
        genotypes_context = vc.get_genotypes()

        new_genotypes: List[Genotype] = []
        for genotype in genotypes_context.iterate_in_sample_name_order():
            genotype_builder = GenotypeBuilder(genotype)
            # Default both artifact-mode flags to false; overwritten below for SNVs.
            genotype_builder.attribute(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_ARTIFACT_MODE, _java_boolean_string(False))
            genotype_builder.attribute(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_RC_ARTIFACT_MODE, _java_boolean_string(False))

            alleles = genotype.get_alleles()
            if genotype.get_ploidy() != 2:
                logger.warning(
                    "No action required:  This tool will skip non-diploid sites.  Saw GT: %s at %s",
                    genotype.get_genotype_string(), vc.to_string_without_genotypes(),
                )

            # Require exactly one reference allele of length 1 (i.e. a possible SNV site).
            ref_alleles = [a.get_base_string() for a in alleles if a.is_reference()]
            if len(ref_alleles) == 1 and len(ref_alleles[0]) == 1:
                ref_allele = ref_alleles[0][0]

                # Only the first alt allele is considered (that is where the F1R2/F2R1 counts live).
                allele = genotype.get_allele(1)
                if allele.is_called() and allele.is_non_reference() and not allele == Allele.SPAN_DEL \
                        and len(allele.get_base_string()) == 1:

                    genotype_mode = Transition.transition_of(ref_allele, allele.get_base_string()[0])
                    is_relevant_artifact = genotype_mode in relevant_transitions_without_complement
                    is_relevant_artifact_complement = genotype_mode in relevant_transitions_complement

                    genotype_builder.attribute(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_ARTIFACT_MODE, _java_boolean_string(is_relevant_artifact))
                    genotype_builder.attribute(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_RC_ARTIFACT_MODE, _java_boolean_string(is_relevant_artifact_complement))

                    genotype_builder.attribute(
                        OrientationBiasFilterConstants.PRE_ADAPTER_METRIC_FIELD_NAME,
                        pre_adapter_q_score_map.get(genotype_mode, OrientationBiasFilterer.PRE_ADAPTER_METRIC_NOT_ARTIFACT_SCORE),
                    )
                    genotype_builder.attribute(
                        OrientationBiasFilterConstants.PRE_ADAPTER_METRIC_RC_FIELD_NAME,
                        pre_adapter_q_score_map.get(genotype_mode.complement(), OrientationBiasFilterer.PRE_ADAPTER_METRIC_NOT_ARTIFACT_SCORE),
                    )

                    # FOB = fraction of alt reads in the artifact orientation (accounting for complement).
                    if is_relevant_artifact or is_relevant_artifact_complement:
                        alt_f1r2 = OrientationBiasUtils.get_genotype_integer(genotype, GATKVCFConstants.OXOG_ALT_F1R2_KEY, 0)
                        alt_f2r1 = OrientationBiasUtils.get_genotype_integer(genotype, GATKVCFConstants.OXOG_ALT_F2R1_KEY, 0)
                        # p-value uses the ORIENTATION-classified alt reads as the binomial n
                        # (n = ALT_F1R2 + ALT_F2R1) and the artifact-orientation count as k
                        # (ALT_F1R2 for the mode, ALT_F2R1 for its complement) -- matching the
                        # CGA MATLAB filter, verified against its i_<stub>_p_value column.
                        # NOTE: the original GATK port used n = AD[1] (t_alt_count) and
                        # k = round(FOB * AD[1]); that diverges from MATLAB whenever the total
                        # alt count differs from ALT_F1R2 + ALT_F2R1 (e.g. cdf(12,13) vs cdf(10,11)).
                        orientation_total = alt_f1r2 + alt_f2r1
                        artifact_count = alt_f1r2 if is_relevant_artifact else alt_f2r1
                        fob = OrientationBiasFilterer._calculate_fob(genotype, is_relevant_artifact)
                        if orientation_total > 0:
                            p_artifact = ArtifactStatisticsScorer.calculate_artifact_p_value(
                                orientation_total, artifact_count, OrientationBiasFilterer.BIAS_P
                            )
                        else:
                            # No orientation-classified alt reads (FOB is NaN): no artifact
                            # evidence, so score a near-zero p (not cut) rather than the
                            # degenerate cdf(0, 0, p) = 1 that would flag it as an artifact.
                            total_alt_allele_count = genotype.get_ad()[1] if genotype.has_ad() else 0
                            p_artifact = ArtifactStatisticsScorer.calculate_artifact_p_value(
                                total_alt_allele_count, 0, OrientationBiasFilterer.BIAS_P
                            )
                        genotype_builder.attribute(OrientationBiasFilterConstants.P_ARTIFACT_FIELD_NAME, p_artifact)
                        genotype_builder.attribute(OrientationBiasFilterConstants.FOB, fob)
                    else:
                        genotype_builder.attribute(OrientationBiasFilterConstants.P_ARTIFACT_FIELD_NAME, VCFConstants.EMPTY_ALLELE)
                        genotype_builder.attribute(OrientationBiasFilterConstants.FOB, VCFConstants.EMPTY_ALLELE)

            new_genotypes.append(genotype_builder.make())

        vcb.genotypes(new_genotypes)
        return vcb.make()

    @staticmethod
    def _calculate_fob(genotype: Genotype, is_relevant_artifact: bool) -> float:
        """Fraction of alt reads supporting orientation bias (Java: private ``calculateFob``)."""
        alt_f2r1 = OrientationBiasUtils.get_genotype_integer(genotype, GATKVCFConstants.OXOG_ALT_F2R1_KEY, 0)
        alt_f1r2 = OrientationBiasUtils.get_genotype_integer(genotype, GATKVCFConstants.OXOG_ALT_F1R2_KEY, 0)
        numerator = alt_f1r2 if is_relevant_artifact else alt_f2r1
        denominator = alt_f1r2 + alt_f2r1
        if denominator == 0:
            # Match Java double division: 0/0.0 -> NaN, positive/0.0 -> +Inf.
            return float("nan") if numerator == 0 else float("inf")
        return numerator / float(denominator)

    # ==================================================================
    # Step 2: FDR-controlled filtering
    # ==================================================================
    @staticmethod
    def annotate_variant_contexts_with_filter_results(
        fdr_threshold: float,
        relevant_transitions_without_complements,   # SortedSet[Transition]
        pre_adapter_q_annotated_variants: List[VariantContext],
        pre_adapter_q_score_map: Dict[Transition, float],
        bias_qp1: float = ArtifactStatisticsScorer.DEFAULT_BIASQP1,
        bias_qp2: float = ArtifactStatisticsScorer.DEFAULT_BIASQP2,
    ) -> List[VariantContext]:
        """Apply orientation-bias filtering keeping the sample FDR below the threshold.

        Java: ``annotateVariantContextsWithFilterResults``. The Java calls the
        one-arg suppression overload (fixed biasQP1=36 / biasQP2=1.5); here those
        shape parameters are exposed so the MAF pipeline can pass its own values
        (run_local used biasQP1=30). Defaults preserve the Java behaviour.
        """
        # Relevant modes AND their complements (we filter both).
        relevant_transitions = set()
        relevant_transitions.update(relevant_transitions_without_complements)
        relevant_transitions.update(OrientationBiasUtils.create_reverse_complement_transitions(relevant_transitions_without_complements))

        if len(pre_adapter_q_annotated_variants) == 0:
            logger.info("No samples found in this file.  NO FILTERING BEING DONE.")
            return pre_adapter_q_annotated_variants

        sample_names = pre_adapter_q_annotated_variants[0].get_sample_names_ordered_by_name()
        sample_name_to_variants = OrientationBiasFilterer.create_sample_to_genotype_variant_context_sorted_map(
            sample_names, pre_adapter_q_annotated_variants
        )

        # vc -> list of updated genotypes (accumulated across samples).
        new_genotypes: Dict[VariantContext, List[Genotype]] = {}

        for sample_name in sample_names:

            # Denominator of the FDR calc: unfiltered, non-ref genotypes for this sample.
            unfiltered_genotype_count = OrientationBiasUtils.calculate_unfiltered_non_ref_genotype_count(
                pre_adapter_q_annotated_variants, sample_name
            )

            # Ordered {genotype: vc} sorted by descending p_artifact (candidates to cut).
            genotypes_to_consider_for_filtering = sample_name_to_variants[sample_name]

            if len(genotypes_to_consider_for_filtering) == 0:
                logger.info("%s: Nothing to filter.", sample_name)
                continue

            # Count candidate genotypes per artifact mode (these can potentially be cut).
            transition_count = OrientationBiasFilterer._create_transition_count_map(
                relevant_transitions, genotypes_to_consider_for_filtering
            )

            # Number to cut per mode, before the preAdapterQ suppression adjustment.
            transition_num_to_cut = OrientationBiasFilterer._create_transition_to_num_cut_pre_pre_adapter_q(
                fdr_threshold, sample_name, unfiltered_genotype_count, genotypes_to_consider_for_filtering, transition_count
            )

            # Scale the per-mode cut count down by the preAdapterQ-based suppression factor.
            for transition in list(transition_num_to_cut.keys()):
                mode_or_reverse_complement = transition if transition in relevant_transitions_without_complements else transition.complement()
                suppression = ArtifactStatisticsScorer.calculate_suppression_factor_from_pre_adapter_q(
                    pre_adapter_q_score_map[mode_or_reverse_complement], bias_qp1, bias_qp2
                )
                transition_num_to_cut[transition] = _java_round(transition_num_to_cut[transition] * suppression)
                logger.info("%s: Cutting (%s) post-preAdapterQ: %s", sample_name, transition, transition_num_to_cut[transition])

            logger.info("%s: Adding orientation bias filter results to genotypes...", sample_name)

            # Walk candidates in sorted (descending p_artifact) order, cutting the first
            # N of each mode until the per-mode quota is reached.
            transition_cut_so_far: Dict[Transition, int] = {transition: 0 for transition in relevant_transitions}
            for genotype in genotypes_to_consider_for_filtering.keys():
                genotype_builder = GenotypeBuilder(genotype)
                transition = Transition.transition_of(
                    genotype.get_allele(0).get_base_string()[0], genotype.get_allele(1).get_base_string()[0]
                )

                if transition not in transition_num_to_cut:
                    logger.warning(
                        "Have to skip genotype: %s since it does not have the artifact mode in the first alt allele.  Total alleles: %s",
                        genotype, len(genotype.get_alleles()),
                    )
                else:
                    p_value = OrientationBiasUtils.get_genotype_double(genotype, OrientationBiasFilterConstants.P_ARTIFACT_FIELD_NAME, 0.0)
                    fraction_of_reads_supporting_orientation_bias = OrientationBiasUtils.get_genotype_double(genotype, OrientationBiasFilterConstants.FOB, 0.0)
                    if transition_cut_so_far[transition] < transition_num_to_cut[transition]:
                        updated_filter = OrientationBiasUtils.add_filter_to_genotype(
                            genotype.get_filters(), OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_CUT
                        )
                        genotype_builder.filter(updated_filter)
                        transition_cut_so_far[transition] += 1
                        logger.info("Cutting: %s %s %s p=%s Fob=%s", genotype.get_sample_name(), genotype.get_allele(0), genotype.get_allele(1), p_value, fraction_of_reads_supporting_orientation_bias)
                    else:
                        logger.info("Passing: %s %s %s p=%s Fob=%s", genotype.get_sample_name(), genotype.get_allele(0), genotype.get_allele(1), p_value, fraction_of_reads_supporting_orientation_bias)

                vc_for_genotype = genotypes_to_consider_for_filtering[genotype]
                new_genotypes.setdefault(vc_for_genotype, []).append(genotype_builder.make())

        # Rebuild the variant contexts with their updated genotypes.
        logger.info("Updating genotypes and creating final list of variants...")
        final_variants: List[VariantContext] = []
        for vc in pre_adapter_q_annotated_variants:
            if vc in new_genotypes:
                gcc = GenotypesContext.copy(vc.get_genotypes())
                new_genotypes_for_this_variant_context = new_genotypes[vc]
                for g in new_genotypes_for_this_variant_context:
                    gcc.replace(g)
                variant_context_builder = VariantContextBuilder(vc).genotypes(gcc)
                # If any genotype was cut, mark the whole site with the OB filter too.
                if any(
                    (g is not None) and (g.get_filters() is not None)
                    and (OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_CUT in g.get_filters())
                    for g in new_genotypes_for_this_variant_context
                ):
                    variant_context_builder.filter(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_CUT)
                final_variants.append(variant_context_builder.make())
            else:
                final_variants.append(vc)
        return final_variants

    @staticmethod
    def _create_transition_to_num_cut_pre_pre_adapter_q(
        fdr_thresh: float,
        sample_name: str,
        unfiltered_genotype_count: int,
        genotypes_to_consider_for_filtering: Dict[Genotype, VariantContext],
        transition_count: Dict[Transition, int],
    ) -> Dict[Transition, int]:
        """Split the total number-to-cut across modes proportional to their counts.

        Java: private ``createTransitionToNumCutPrePreAdapterQ``.
        """
        all_transition_count = sum(transition_count.values())
        total_num_to_cut = OrientationBiasFilterer._calculate_total_num_to_cut_from_genotypes(
            fdr_thresh, unfiltered_genotype_count, genotypes_to_consider_for_filtering
        )

        logger.info("%s: Cutting (total) pre-preAdapterQ: %s", sample_name, total_num_to_cut)

        transition_num_to_cut: Dict[Transition, int] = {transition: 0 for transition in transition_count}
        for transition in transition_num_to_cut:
            # NOTE: integer arithmetic, matching the Java (long * long / long).
            if all_transition_count == 0:
                transition_num_to_cut[transition] = 0
            else:
                transition_num_to_cut[transition] = (total_num_to_cut * transition_count[transition]) // all_transition_count
            logger.info("%s: Cutting (%s) pre-preAdapterQ: %s", sample_name, transition, transition_num_to_cut[transition])
        return transition_num_to_cut

    @staticmethod
    def _calculate_total_num_to_cut_from_genotypes(
        fdr_thresh: float,
        unfiltered_genotype_count: int,
        genotypes_to_consider_for_filtering: Dict[Genotype, VariantContext],
    ) -> int:
        """Java: private ``calculateTotalNumToCut(fdrThresh, count, SortedMap)`` overload.

        Pads the artifact p-values with zeros for the non-artifact genotypes so the
        Benjamini-Hochberg procedure controls FDR over ALL unfiltered variants.
        """
        p_artifact_scores = [
            OrientationBiasUtils.get_genotype_double(g, OrientationBiasFilterConstants.P_ARTIFACT_FIELD_NAME, 0.0)
            for g in genotypes_to_consider_for_filtering.keys()
        ]

        # Pad with zeros (p-value 0) for the non-artifact-mode SNVs.
        num_to_pad_zeroes = int(unfiltered_genotype_count) - len(p_artifact_scores)
        final_p_artifact_scores = list(p_artifact_scores) + [0.0] * num_to_pad_zeroes

        return OrientationBiasFilterer.calculate_total_num_to_cut(
            fdr_thresh, unfiltered_genotype_count, final_p_artifact_scores
        )

    @staticmethod
    def calculate_total_num_to_cut(
        fdr_threshold: float,
        unfiltered_genotype_count: int,
        p_artifact_scores_including_non_artifact: List[float],
    ) -> int:
        """Benjamini-Hochberg number to cut (Java: ``calculateTotalNumToCut`` list overload).

        ``p_artifact_scores_including_non_artifact`` must be sorted DESCENDING and
        include zeros for the non-artifact variants. Returns the first index ``i``
        whose score drops below the BH line ``fdrThreshold * (i+1) / n``; if none
        do, returns ``len - 1``.
        """
        # ParamUtils.isPositive
        if not fdr_threshold > 0:
            raise ValueError("FDR threshold must be positive and greater than zero.")

        # https://en.wikipedia.org/wiki/False_discovery_rate#Benjamini-Hochberg_procedure
        for i in range(len(p_artifact_scores_including_non_artifact)):
            if p_artifact_scores_including_non_artifact[i] < fdr_threshold * (i + 1) / unfiltered_genotype_count:
                return i
        return len(p_artifact_scores_including_non_artifact) - 1

    @staticmethod
    def _create_transition_count_map(
        relevant_transitions,                       # SortedSet[Transition]
        genotypes_to_consider_for_filtering: Dict[Genotype, VariantContext],
    ) -> Dict[Transition, int]:
        """Count candidate genotypes per artifact mode (Java: private ``createTransitionCountMap``)."""
        transition_count: Dict[Transition, int] = {transition: 0 for transition in relevant_transitions}
        for g in genotypes_to_consider_for_filtering.keys():
            for transition in relevant_transitions:
                if OrientationBiasUtils.is_genotype_in_transition(g, transition):
                    transition_count[transition] += 1
        return transition_count

    @staticmethod
    def create_sample_to_genotype_variant_context_sorted_map(
        sample_names: List[str],
        variants,                                   # Collection[VariantContext]
    ) -> Dict[str, "Dict[Genotype, VariantContext]"]:
        """Map each sample to a {genotype: vc} dict ordered by descending p_artifact.

        Java: ``createSampleToGenotypeVariantContextSortedMap`` (a TreeMap keyed by a
        comparator). Python dicts preserve insertion order, so we insert candidates
        already sorted to reproduce the ordered-map iteration behaviour.
        """
        def sort_key(pair):
            g = pair[0]
            # Negative p_artifact => descending; id(g) breaks ties (Java used hashCode()).
            return (-OrientationBiasUtils.get_genotype_double(g, OrientationBiasFilterConstants.P_ARTIFACT_FIELD_NAME, 0.0), id(g))

        sample_name_to_variants: Dict[str, Dict[Genotype, VariantContext]] = {}
        for sample_name in sample_names:
            candidates = []
            for vc in variants:
                for genotype in vc.get_genotypes(sample_name):
                    if OrientationBiasFilterer._is_filtering_candidate(genotype, vc):
                        candidates.append((genotype, vc))
            candidates.sort(key=sort_key)
            # Build the ordered {genotype: vc} map in sorted order.
            genotypes_to_consider_for_filtering: Dict[Genotype, VariantContext] = {}
            for genotype, vc in candidates:
                genotypes_to_consider_for_filtering[genotype] = vc
            sample_name_to_variants[sample_name] = genotypes_to_consider_for_filtering
        return sample_name_to_variants

    @staticmethod
    def _is_filtering_candidate(genotype: Genotype, vc: VariantContext) -> bool:
        """Whether this genotype can be considered for the OB filter.

        Java: private ``isFilteringCandidate``.
        """
        return (
            not vc.is_filtered()
            and not genotype.is_filtered()
            and (
                genotype.get_any_attribute(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_ARTIFACT_MODE) == _java_boolean_string(True)
                or genotype.get_any_attribute(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_RC_ARTIFACT_MODE) == _java_boolean_string(True)
            )
        )

    @staticmethod
    def create_vcf_header(input_vcf_header, command_line: str, transitions: List[str]):
        """Augment the VCF header with the FORMAT/FILTER lines the filter adds.

        Java: ``createVCFHeader``. htsjdk's ``VCFHeader`` / ``VCF*HeaderLine`` types
        are not reproduced here; instead a header line is represented as a plain
        dict and the "header" is the list of such lines. Swap in real htsjdk/pysam
        header objects if writing an actual VCF.
        """
        if input_vcf_header is None or transitions is None:
            raise ValueError("Arguments cannot be null")

        # Start from the existing header lines (assumed to be an iterable of dicts/str).
        header_lines = list(input_vcf_header)

        def fmt(field_id, number, value_type, description):
            return {"kind": "FORMAT", "ID": field_id, "Number": number, "Type": value_type, "Description": description}

        header_lines.append(fmt(OrientationBiasFilterConstants.PRE_ADAPTER_METRIC_FIELD_NAME, "A", "Float", "Measure (across entire bam file) of orientation bias for a given REF/ALT error."))
        header_lines.append(fmt(OrientationBiasFilterConstants.PRE_ADAPTER_METRIC_RC_FIELD_NAME, "A", "Float", "Measure (across entire bam file) of orientation bias for the complement of a given REF/ALT error."))
        header_lines.append(fmt(OrientationBiasFilterConstants.P_ARTIFACT_FIELD_NAME, "A", "Float", "Orientation bias p value for the given REF/ALT artifact or its complement."))
        header_lines.append(fmt(OrientationBiasFilterConstants.FOB, "A", "Float", "Fraction of alt reads indicating orientation bias error (taking into account artifact mode complement)."))
        header_lines.append(fmt(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_ARTIFACT_MODE, "A", "String", "Whether the variant can be one of the given REF/ALT artifact modes."))
        header_lines.append(fmt(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_RC_ARTIFACT_MODE, "A", "String", "Whether the variant can be one of the given REF/ALT artifact mode complements."))
        header_lines.append({"kind": "FORMAT", "ID": VCFConstants.GENOTYPE_FILTER_KEY, "Number": 1, "Type": "String", "Description": "Genotype-level filter"})
        header_lines.append({"kind": "FILTER", "ID": OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_CUT, "Description": "Orientation bias (in one of the specified artifact mode(s) or complement) seen in one or more samples."})
        header_lines.append({"kind": "SIMPLE", "ID": "orientation_bias_artifact_modes", "Value": "|".join(transitions), "Description": "The artifact modes that were used for orientation bias artifact filtering for this VCF"})
        header_lines.append({"kind": "OTHER", "key": "command", "value": command_line})
        return header_lines
