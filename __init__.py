"""Python implementation of the orientation-bias variant filter.

Based on the GATK FilterByOrientationBias function. The command-line entry point
is ``maf_orientation_bias_filter`` (reads a MAF, writes filtered / unfiltered
MAFs). See ``README.md`` for usage.

Modules:
  * maf_orientation_bias_filter        - command-line entry point / MAF adapter
  * orientation_bias_filterer          - per-variant annotation + Benjamini-Hochberg cut
  * orientation_bias_utils             - typed field getters, mode membership, counting
  * artifact_statistics_scorer         - artifact p-value + preAdapterQ suppression factor
  * orientation_bias_filter_constants  - VCF FORMAT/FILTER field-name constants
  * transition                         - single-base substitution (an artifact mode)
  * htsjdk_models                      - in-memory Allele / Genotype / VariantContext model
"""
