"""Translation of ``ArtifactStatisticsScorer.java``.

Two pieces of the orientation-bias statistics:
  * ``calculate_suppression_factor_from_pre_adapter_q`` - a sigmoid on preAdapterQ
    that scales down how many artifacts to cut when the bam-level artifact signal
    is weak (a high preAdapterQ means "few artifacts, cut nothing").
  * ``calculate_artifact_p_value`` - the binomial p-value that a variant is an
    orientation-bias artifact given its alt-read orientation split.
"""

from __future__ import annotations

import math

from scipy.stats import binom


class ArtifactStatisticsScorer:

    # Default sigmoid shape parameters (match the CGA MATLAB OxoG filter defaults).
    DEFAULT_BIASQP1 = 36    # inflection point (in preAdapterQ / Phred units)
    DEFAULT_BIASQP2 = 1.5   # steepness; 0 turns the sigmoid into a hard cutoff

    def __init__(self):
        raise RuntimeError("ArtifactStatisticsScorer is not instantiable")

    @staticmethod
    def calculate_suppression_factor_from_pre_adapter_q(
        pre_adapter_q: float,
        bias_qp1: float = DEFAULT_BIASQP1,
        bias_qp2: float = DEFAULT_BIASQP2,
    ) -> float:
        """Multiplier (0..1) applied to the number of artifacts to cut.

        Java: ``calculateSuppressionFactorFromPreAdapterQ`` (both the 3-arg form and
        the 1-arg convenience overload that supplies the DEFAULT_BIASQP* values -
        collapsed here into default arguments).

        Zero means "this sample probably has no artifacts, so cut nothing".
        """
        # ParamUtils.isPositive / isPositiveOrZero
        if not pre_adapter_q > 0:
            raise ValueError("preAdapter Q score must be positive and not zero.")
        if bias_qp1 < 0 or bias_qp2 < 0:
            raise ValueError("bias Q shape parameters must be positive.")

        # From MATLAB: fQ = 1 / (1 + exp(biasQP2 * (biasQ - biasQP1)))
        # where biasQ is the preAdapterQ score.
        if bias_qp2 == 0:
            # Sharp cutoff: cut nothing once preAdapterQ exceeds the inflection point.
            return 0.0 if pre_adapter_q > bias_qp1 else 1.0
        return 1.0 / (1.0 + math.exp(bias_qp2 * (pre_adapter_q - bias_qp1)))

    @staticmethod
    def calculate_artifact_p_value(
        total_alt_allele_count: int,
        artifact_alt_allele_count: int,
        bias_p: float,
    ) -> float:
        """p-value for the variant being an orientation-bias artifact.

        Java: ``calculateArtifactPValue``. Models the number of artifact-oriented
        alt reads as Binomial(n = total alt reads, p = bias_p) and returns the
        cumulative probability at ``artifact_alt_allele_count``.

        :param bias_p: believed bias p-value of the binomial (mode of the pair
                       orientation for artifact alt reads, e.g. 0.96).
        """
        if bias_p < 0:
            raise ValueError("bias parameter must be positive or zero.")
        if total_alt_allele_count < 0:
            raise ValueError("total alt allele count must be positive or zero.")
        if artifact_alt_allele_count < 0:
            raise ValueError("artifact supporting alt allele count must be positive or zero.")
        if total_alt_allele_count - artifact_alt_allele_count < 0:
            raise ValueError("Total alt count must be same or greater than the artifact alt count.")

        # commons-math BinomialDistribution(n, p).cumulativeProbability(k)
        # == scipy.stats.binom.cdf(k, n, p)
        return float(binom.cdf(artifact_alt_allele_count, total_alt_allele_count, bias_p))
