"""Translation of ``OrientationBiasUtils.java``.

Grab-bag of helpers used by the filterer: reading typed FORMAT fields off a
genotype, testing whether a genotype falls in a given artifact mode (Transition)
or its complement, building reverse-complement modes, and counting the unfiltered
non-ref genotypes that form the Benjamini-Hochberg FDR denominator.
"""

from __future__ import annotations

import logging
from typing import Collection, List, Optional

from htsjdk_models import Allele, Genotype, VariantContext, VCFConstants
from transition import Transition

logger = logging.getLogger(__name__)


class OrientationBiasUtils:

    def __init__(self):
        raise RuntimeError("OrientationBiasUtils is not instantiable")

    # ------------------------------------------------------------------
    # Typed FORMAT-field accessors
    # ------------------------------------------------------------------
    @staticmethod
    def get_genotype_double(g: Genotype, field_name: str, default_value: float) -> float:
        """FORMAT field as float, or ``default_value`` if empty/missing.

        Java: ``getGenotypeDouble``.
        """
        if g is None:
            raise ValueError("Genotype cannot be null")
        genotype_field_as_string = OrientationBiasUtils.get_genotype_string(g, field_name)
        if OrientationBiasUtils._is_vcf_genotype_field_empty(genotype_field_as_string):
            return default_value
        return float(genotype_field_as_string)

    @staticmethod
    def get_genotype_integer(g: Genotype, field_name: str, default_value: int) -> int:
        """FORMAT field as int, or ``default_value`` if empty/missing.

        Java: ``getGenotypeInteger``.
        """
        if g is None:
            raise ValueError("Genotype cannot be null")
        genotype_field_as_string = OrientationBiasUtils.get_genotype_string(g, field_name)
        if OrientationBiasUtils._is_vcf_genotype_field_empty(genotype_field_as_string):
            return default_value
        return int(genotype_field_as_string)

    @staticmethod
    def _is_vcf_genotype_field_empty(genotype_field_as_string: Optional[str]) -> bool:
        # Java: private isVcfGenotypeFieldEmpty
        return (genotype_field_as_string is None) or (genotype_field_as_string == VCFConstants.MISSING_VALUE_v4)

    @staticmethod
    def get_genotype_string(g: Genotype, field_name: str) -> str:
        """FORMAT field as a string, or "." if absent (Java: ``getGenotypeString``)."""
        if g is None:
            raise ValueError("Genotype cannot be null")
        genotype_extended_attribute = g.get_extended_attribute(field_name, VCFConstants.MISSING_VALUE_v4)
        return str(genotype_extended_attribute)

    # ------------------------------------------------------------------
    # Artifact-mode (Transition) membership tests
    # ------------------------------------------------------------------
    @staticmethod
    def is_genotype_in_transition(g: Genotype, transition: Transition) -> bool:
        """Whether any alt allele makes this genotype the given ref->alt mode.

        Complement is NOT considered here. Java: ``isGenotypeInTransition``.
        """
        if g is None:
            raise ValueError("Genotype cannot be null")
        if transition is None:
            raise ValueError("Artifact mode cannot be null")

        alleles = g.get_alleles()
        # allele[0] must match the mode's ref base, and some alt allele[i>=1] the call base.
        return any(
            g.get_allele(0).bases_match(Allele.create(transition.ref(), True))
            and g.get_allele(i).bases_match(Allele.create(transition.call()))
            for i in range(1, len(alleles))
        )

    @staticmethod
    def is_genotype_in_transition_with_complement(g: Genotype, transition: Transition) -> bool:
        """Like :meth:`is_genotype_in_transition` but also matches the complement.

        Java: ``isGenotypeInTransitionWithComplement``.
        """
        if g is None:
            raise ValueError("Genotype cannot be null")
        if transition is None:
            raise ValueError("Transition cannot be null")
        return (
            OrientationBiasUtils.is_genotype_in_transition(g, transition)
            or OrientationBiasUtils.is_genotype_in_transition(g, transition.complement())
        )

    @staticmethod
    def is_genotype_in_transitions_with_complement(g: Genotype, transitions: Collection[Transition]) -> bool:
        """Whether the genotype falls in ANY of the modes (or their complements).

        Java: ``isGenotypeInTransitionsWithComplement``.
        """
        if g is None:
            raise ValueError("Genotype cannot be null.")
        return any(
            OrientationBiasUtils.is_genotype_in_transition_with_complement(g, am) for am in transitions
        )

    @staticmethod
    def create_reverse_complement_transitions(transitions: Collection[Transition]) -> List[Transition]:
        """Complements of the given modes (Java: ``createReverseComplementTransitions``)."""
        if transitions is None:
            raise ValueError("Transitions cannot be null.")
        return [transition.complement() for transition in transitions]

    # ------------------------------------------------------------------
    # Per-sample counting helper (the FDR denominator)
    # ------------------------------------------------------------------
    @staticmethod
    def _get_genotype_stream(sample_name: str, variant_contexts: List[VariantContext]) -> List[Genotype]:
        # Java: private getGenotypeStream - this sample's genotype from each unfiltered site.
        genotypes = []
        for vc in variant_contexts:
            if not vc.is_filtered():
                g = vc.get_genotype(sample_name)
                if g is not None:
                    genotypes.append(g)
        return genotypes

    @staticmethod
    def calculate_unfiltered_non_ref_genotype_count(variants: List[VariantContext], sample_name: str) -> int:
        """Count unfiltered, non-ref/ref genotypes for a sample.

        Java: ``calculateUnfilteredNonRefGenotypeCount``. This is the denominator
        of the Benjamini-Hochberg FDR calculation.
        """
        if variants is None:
            raise ValueError("variants cannot be null")
        if sample_name is None:
            raise ValueError("sampleName cannot be null")
        return sum(
            1
            for g in OrientationBiasUtils._get_genotype_stream(sample_name, variants)
            if not g.is_filtered() and not g.get_allele(0).bases_match(g.get_allele(1))
        )

    # ------------------------------------------------------------------
    # Genotype filter-string helper
    # ------------------------------------------------------------------
    @staticmethod
    def add_filter_to_genotype(existing_filter_value: Optional[str], new_filter_to_add: str) -> str:
        """Append ``new_filter_to_add`` to a genotype filter string.

        Java: ``addFilterToGenotype``.
        """
        if new_filter_to_add is None:
            raise ValueError("newFilterToAdd cannot be null")

        if (existing_filter_value is None) or (len(existing_filter_value.strip()) == 0) \
                or (existing_filter_value == VCFConstants.UNFILTERED) \
                or (existing_filter_value == VCFConstants.PASSES_FILTERS_v4):
            # No meaningful existing filter -> just use the new one.
            return new_filter_to_add
        elif len(existing_filter_value) > 0:
            # Append with the standard ";" separator.
            return existing_filter_value + VCFConstants.FILTER_CODE_SEPARATOR + new_filter_to_add
        else:
            # Unreachable branch kept to mirror the Java (which warns here).
            appended_filter_string = existing_filter_value + VCFConstants.FILTER_CODE_SEPARATOR + new_filter_to_add
            logger.warning(
                "Existing genotype filter could be incorrect: %s ... Proceeding with %s ...",
                existing_filter_value, appended_filter_string,
            )
            return appended_filter_string
