"""Translation of ``OrientationBiasFilterConstants.java``.

Just the string constants (VCF FORMAT keys / filter names) the filter reads and
writes. Kept as class attributes so call sites read the same as the Java
(``OrientationBiasFilterConstants.OBP`` etc.).
"""

from __future__ import annotations


class OrientationBiasFilterConstants:
    # preAdapterQ score for the artifact mode (bam-file-level orientation-bias level).
    PRE_ADAPTER_METRIC_FIELD_NAME = "OBQ"
    # preAdapterQ score for the reverse-complement artifact mode.
    PRE_ADAPTER_METRIC_RC_FIELD_NAME = "OBQRC"
    # per-variant p-value of being an orientation-bias artifact.
    P_ARTIFACT_FIELD_NAME = "OBP"
    # genotype/variant filter code applied when a variant is cut as an artifact.
    IS_ORIENTATION_BIAS_CUT = "orientation_bias"
    # whether the variant matches one of the requested artifact modes.
    IS_ORIENTATION_BIAS_ARTIFACT_MODE = "OBAM"
    # whether the variant matches the complement of a requested artifact mode.
    IS_ORIENTATION_BIAS_RC_ARTIFACT_MODE = "OBAMRC"
    # fraction of alt reads indicating orientation-bias error.
    FOB = "OBF"

    def __init__(self):
        # Java hides the constructor; mirror that (this class is a constant holder).
        raise RuntimeError("OrientationBiasFilterConstants is not instantiable")
