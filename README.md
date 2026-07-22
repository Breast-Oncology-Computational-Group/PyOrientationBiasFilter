# Orientation Bias Filter (Python)

A Python implementation of the orientation-bias variant filter, based on the GATK
**FilterByOrientationBias** function:
<https://gatk.broadinstitute.org/hc/en-us/articles/360036834571-FilterByOrientationBias-EXPERIMENTAL>.

The filter identifies and removes sequencing **orientation-bias artifacts** — most
commonly OxoG (8-oxoguanine, a `G>T` / `C>A` error) and FFPE deamination
(`C>T` / `G>A`). For each candidate SNV it looks at how the supporting alt reads
split between the two read orientations, scores how consistent that split is with
an artifact, and — controlling the false-discovery rate with a Benjamini-Hochberg
procedure and a per-sample preAdapterQ score — cuts the variants most likely to be
artifacts.

## What it does

The entry point reads an annotated MAF (one that already carries the per-variant
orientation read counts `i_t_ALT_F1R2` / `i_t_ALT_F2R1`, the alt-read count
`t_alt_count`, and a preAdapterQ value `i_<stub>_Q`), decides which rows are
orientation-bias artifacts, and writes the results back out as MAFs with four
appended columns:

| Column | Meaning |
|---|---|
| `i_<stub>_F` | fraction of alt reads in the artifact orientation |
| `i_<stub>_mode` | 1 if the variant is in the artifact mode (or its complement), else 0 |
| `i_<stub>_p_value` | orientation-bias artifact p-value |
| `i_<stub>_cut` | 1 if cut as an artifact, else 0 — the column downstream keys on |

`<stub>` is the artifact label (e.g. `oxog`, `ffpe`).

## Algorithm

The filter runs in two stages over a single tumor sample.

### Stage 1 — score each variant (`orientation_bias_filterer` / `artifact_statistics_scorer`)
For every SNV genotype:

1. **Artifact mode.** A substitution is `ref→alt`. You pass one artifact mode (e.g.
   `C>T` for FFPE, `G>T` for OxoG); the filter also considers its reverse
   complement. A variant is *in mode* if its `ref→alt` matches, *in complement* if
   it matches the reverse complement, otherwise it is not an artifact candidate.
   → `i_<stub>_mode` = 1 for in-mode-or-complement.
2. **FOB — fraction of orientation bias** (`i_<stub>_F`). Alt reads split by read-pair
   orientation into F1R2 and F2R1; the artifact concentrates in one of them. For an
   in-mode variant that orientation is F1R2, so
   `FOB = ALT_F1R2 / (ALT_F1R2 + ALT_F2R1)`; for the complement it is F2R1,
   `FOB = ALT_F2R1 / (…)`. (Both counts zero → NaN.)
3. **Artifact p-value** (`i_<stub>_p_value`). Model the artifact-oriented alt-read
   count as `Binomial(n, p = pBias = 0.96)` and take the cumulative probability
   (`scipy.stats.binom.cdf`) at the artifact-orientation count `k`, where:
   - `n = ALT_F1R2 + ALT_F2R1` — the alt reads that got a read-pair orientation, and
   - `k = ALT_F1R2` for an in-mode variant (or `ALT_F2R1` for the complement) — i.e.
     `k = FOB · n`.

   A real artifact has almost all alt reads in the artifact orientation (`k ≈ n`) →
   the count sits near the top of the distribution → **p ≈ 1**; a clean ~50/50
   variant → **p ≈ 0**. So a *higher* p-value is *more* artifact-like.

   > This matches the CGA MATLAB filter, verified against its `i_<stub>_p_value`
   > column. Note `n` is the **orientation-classified** alt total (`ALT_F1R2 +
   > ALT_F2R1`), **not** `t_alt_count` (`AD[1]`); the two often differ, and the
   > original GATK port used `AD[1]` with `k = round(FOB × AD[1])`, which is why its
   > p-values diverged from MATLAB (e.g. `cdf(10,11)` vs `cdf(12,13)`).

### Stage 2 — decide how many to cut (FDR + preAdapterQ suppression)
Across the sample's unfiltered non-ref genotypes:

1. **N** = the sample's total passing SNV calls = the FDR denominator. The
   artifact-mode variants are the candidates; every other variant contributes p = 0.
2. **Benjamini-Hochberg count.** Sort all p-values descending and walk them, cutting
   while `p[i] ≥ fdr · (i+1) / N`; stop at the first that falls below that line. The
   index reached is the number to cut. (Because "artifact" = high p, the cut is on the
   *high* tail — the reverse of a classic BH on small p-values.)
3. **Split by mode.** Divide that count between the mode and its complement in
   proportion to each one's candidate count (integer arithmetic).
4. **preAdapterQ suppression.** Scale each per-mode count by a sigmoid of the
   bam-level preAdapterQ: `f(Q) = 1 / (1 + exp(bias_qp2 · (Q − bias_qp1)))`. A high Q
   (bam shows little global artifact) → factor ≈ 0 → cut nothing; a low Q → factor ≈ 1
   → cut the full BH count. `bias_qp1` (inflection, default 30 here) and `bias_qp2`
   (steepness, 1.5) shape the curve.
5. **Cut.** For each mode, mark the top-`count` candidates by descending p-value as
   artifacts → `i_<stub>_cut = 1`.

### Parameters
`pBias = 0.96` (binomial mode, fixed in the core), `--fdr-threshold` (artifact FDR,
0.01), `--bias-qp1` / `--bias-qp2` (suppression sigmoid). preAdapterQ is taken from
the MAF's `i_<stub>_Q` column (see below).

### Worked example

One tumor sample with **100** unfiltered non-ref SNVs, filtering FFPE (mode `C>T`),
`i_ffpe_Q = 25`, `--fdr-threshold 0.01 --bias-qp1 30`. Three SNVs are `C>T` (the
candidates); the other 97 are other substitutions (not in mode → p treated as 0).

**Stage 1** (`C>T` is in-mode → artifact orientation is F1R2, so `FOB = F1R2/n`,
`k = F1R2`, `n = ALT_F1R2 + ALT_F2R1`; p = `cdf(k, n, 0.96)`):

| variant | ALT_F1R2 | ALT_F2R1 | n (=F1R2+F2R1) | `i_ffpe_F` | k | `i_ffpe_p_value` |
|---|---|---|---|---|---|---|
| v1 | 20 | 0 | 20 | 1.00 | 20 | 1.000 |
| v2 | 18 | 2 | 20 | 0.90 | 18 | 0.190 |
| v3 | 10 | 10 | 20 | 0.50 | 10 | 0.000 |

**Stage 2.** Candidate p-values sorted descending = `[1.000, 0.190, 0.000]`, padded to
N=100 with zeros.
- BH: cut while `p[i] ≥ 0.01·(i+1)/100` → v1 (`1.000 ≥ 0.0001` ✓), v2 (`0.190 ≥ 0.0002`
  ✓), v3 (`0.000 ≥ 0.0003` ✗, stop) ⇒ **2** to cut.
- Suppression `f(25) = 1/(1+exp(1.5·(25−30))) ≈ 0.999` → `round(2 × 0.999) = 2`.

Result: v1, v2 → `i_ffpe_cut = 1`; v3 → `0`. Note the suppression lever: on a clean bam
(e.g. `Q = 40`, `f ≈ 0.003`) the count becomes `round(2 × 0.003) = 0` — nothing is cut
even though v1 is a perfect-looking artifact.

## How to run

```bash
pip install -r requirements.txt        # scipy (one-time)

python maf_orientation_bias_filter.py \
  -i  <pairName>.OrientationBiasInfo.maf \
  -o  <pairName>.OrientationBiasFilter.maf \
  -s  oxog \                                # artifact stub
  --reference-allele G \                    # artifact-mode reference base
  --artifact-allele  T \                    # artifact-mode alt base
  --fdr-threshold 0.01 --bias-qp1 30 --bias-qp2 1.5
```

It writes:

- `<output>.maf` — filtered MAF (rows that pass, `i_<stub>_cut == 0`)
- `<output w/o .maf>.unfiltered.maf` — every row, with the computed columns
- `<output>.maf.pass_count.txt` / `<output>.maf.reject_count.txt` — the counts

### Memory

The filter is **streaming**: it reads the MAF in two passes and never holds the
whole (wide) file in memory. Only the small subset of artifact-mode rows is kept
between passes, and rows are matched between passes by a variant key (chromosome,
start, end, ref, and both tumor alleles) rather than by position. Peak memory is
proportional to the number of artifact-mode candidates, not to the file size, so
multi-million-row whole-genome MAFs run in well under a couple of GB.

## Files

### Entry point
- **`maf_orientation_bias_filter.py`** — the command-line program and MAF adapter.
  Streams the input MAF, builds the small per-variant model the filter needs,
  runs the two filter stages (per-variant annotation, then the FDR cut), appends
  the four `i_<stub>_*` columns, and writes the filtered / unfiltered MAFs and the
  pass / reject count files. This is the only file you invoke directly.

### Filter logic
- **`orientation_bias_filterer.py`** — the heart of the filter. Two operations:
  annotate a variant with its orientation-bias values (artifact-mode membership,
  artifact p-value, and the fraction of alt reads supporting the bias
  orientation), and apply the Benjamini-Hochberg cut across the sample (how many
  artifacts to cut per mode, scaled by the preAdapterQ suppression factor, and
  which genotypes/sites to mark).
- **`orientation_bias_utils.py`** — helper routines used by the filter: reading
  typed fields (float / int / string) off a genotype with default handling;
  testing whether a genotype falls in a given artifact mode or its reverse
  complement; building complement modes; and the per-sample unfiltered non-ref
  genotype count that is the FDR denominator.
- **`artifact_statistics_scorer.py`** — the statistics: the binomial **artifact
  p-value** (via `scipy.stats.binom`) and the **preAdapterQ suppression factor**,
  a sigmoid (shape parameters `bias_qp1`, `bias_qp2`) that scales down how many
  artifacts to cut when the sample's artifact signal is weak.

### Model and constants
- **`htsjdk_models.py`** — a lightweight in-memory variant model (`Allele`,
  `Genotype` / `GenotypeBuilder`, `VariantContext` / `VariantContextBuilder`,
  `GenotypesContext`) plus the VCF constant holders. This is the object model the
  adapter builds from each MAF row and the filter operates on. It is a simplified
  model with only the accessors the algorithm uses — not a full VCF parser.
- **`transition.py`** — `Transition`, a single-base substitution (reference base →
  alt base) that encodes an artifact "mode", with helpers to build it, take its
  reverse complement, and render / parse its name (e.g. `G>T`).
- **`orientation_bias_filter_constants.py`** — the field-name constants for the
  orientation-bias values written onto genotypes (the preAdapterQ metric and its
  complement, the artifact p-value, the fraction of orientation-bias reads, the
  artifact-mode / reverse-complement-mode flags, and the orientation-bias filter
  code).

### Support files
- **`__init__.py`** — package marker and module map.
- **`requirements.txt`** — the dependency: `scipy` (needed for the artifact
  p-value).

## preAdapterQ

preAdapterQ is read straight from the MAF's `i_<stub>_Q` column (computed earlier
in the pipeline), so this package does not recompute it from Picard metrics.
