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
