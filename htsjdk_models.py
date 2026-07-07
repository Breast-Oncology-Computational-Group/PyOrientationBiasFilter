"""Simplified in-memory stand-ins for the htsjdk variant model used by the filter.

The Java code is written against htsjdk's ``Allele`` / ``Genotype`` /
``VariantContext`` classes (plus their immutable "builder" companions) and a few
VCF constant holders. Re-implementing htsjdk is out of scope, so this module
provides *just enough* of that object model - as plain Python classes - for the
translated algorithm to run and to read the same way as the original.

These are NOT a VCF parser. In a real deployment you would back these accessors
with pysam / cyvcf2 records; the translated logic only depends on the small
method surface reproduced here.
"""

from __future__ import annotations

import copy as _copy
from typing import Dict, List, Optional


# ---------------------------------------------------------------------------
# VCF constant holders (htsjdk VCFConstants / GATKVCFConstants)
# ---------------------------------------------------------------------------
class VCFConstants:
    MISSING_VALUE_v4 = "."          # the "." used for absent field values
    EMPTY_ALLELE = "."              # written to a genotype attribute when N/A
    PASSES_FILTERS_v4 = "PASS"      # the canonical "passed all filters" token
    UNFILTERED = "."                # "no filtering applied yet"
    FILTER_CODE_SEPARATOR = ";"     # separates multiple filter codes
    GENOTYPE_FILTER_KEY = "FT"      # FORMAT key holding per-genotype filters


class GATKVCFConstants:
    # FORMAT keys holding the F1R2 / F2R1 alt-read counts produced upstream.
    OXOG_ALT_F1R2_KEY = "ALT_F1R2"
    OXOG_ALT_F2R1_KEY = "ALT_F2R1"


# ---------------------------------------------------------------------------
# Allele
# ---------------------------------------------------------------------------
class Allele:
    """A single VCF allele (a base string plus a reference flag)."""

    # htsjdk exposes a shared spanning-deletion allele "*"; reproduced so the
    # filterer can compare against it (``allele.equals(Allele.SPAN_DEL)``).
    SPAN_DEL: "Allele"  # assigned just below the class definition

    def __init__(self, base_string: str, is_ref: bool = False):
        self._bases = str(base_string).upper()
        self._is_ref = bool(is_ref)

    @staticmethod
    def create(base_string, is_ref: bool = False) -> "Allele":
        """Java: ``Allele.create(bases, isRef)``. Accepts a str or an int byte."""
        if isinstance(base_string, int):        # a Java ``(byte) char`` value
            base_string = chr(base_string)
        return Allele(base_string, is_ref)

    def is_reference(self) -> bool:
        return self._is_ref

    def is_non_reference(self) -> bool:
        return not self._is_ref

    def is_called(self) -> bool:
        # An uncalled allele is the VCF no-call "."; anything else is called.
        return self._bases not in (".", "")

    def get_base_string(self) -> str:
        return self._bases

    def bases_match(self, other: "Allele") -> bool:
        """Compare only the bases, ignoring the reference flag (htsjdk semantics)."""
        return self._bases == other._bases

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Allele) and self._bases == other._bases and self._is_ref == other._is_ref

    def __hash__(self) -> int:
        return hash((self._bases, self._is_ref))

    def __str__(self) -> str:
        return self._bases


Allele.SPAN_DEL = Allele("*", is_ref=False)


# ---------------------------------------------------------------------------
# Genotype  (immutable-ish; edits go through GenotypeBuilder)
# ---------------------------------------------------------------------------
class Genotype:
    """One sample's genotype at one site.

    ``attributes`` holds the extended (FORMAT) fields, e.g. the OBP/OBQ/FOB
    values the filter writes. ``filters`` is the per-genotype filter string
    ("PASS"/None means unfiltered).
    """

    def __init__(
        self,
        sample_name: str,
        alleles: List[Allele],
        attributes: Optional[Dict[str, object]] = None,
        filters: Optional[str] = None,
        ad: Optional[List[int]] = None,
    ):
        self._sample_name = sample_name
        self._alleles = list(alleles)
        self._attributes = dict(attributes) if attributes else {}
        self._filters = filters
        self._ad = list(ad) if ad is not None else None

    # --- allele accessors ---
    def get_alleles(self) -> List[Allele]:
        return self._alleles

    def get_allele(self, i: int) -> Allele:
        return self._alleles[i]

    def get_ploidy(self) -> int:
        return len(self._alleles)

    # --- allele-depth accessors ---
    def has_ad(self) -> bool:
        return self._ad is not None

    def get_ad(self) -> List[int]:
        return self._ad

    # --- FORMAT / extended-attribute accessors ---
    def get_extended_attribute(self, field_name: str, default: object = None) -> object:
        return self._attributes.get(field_name, default)

    def get_any_attribute(self, field_name: str) -> object:
        return self._attributes.get(field_name)

    # --- filter / identity accessors ---
    def get_filters(self) -> Optional[str]:
        return self._filters

    def is_filtered(self) -> bool:
        # htsjdk: filtered iff a non-empty, non-PASS filter string is present.
        return self._filters not in (None, "", VCFConstants.PASSES_FILTERS_v4, VCFConstants.UNFILTERED)

    def get_sample_name(self) -> str:
        return self._sample_name

    def get_genotype_string(self) -> str:
        return "/".join(a.get_base_string() for a in self._alleles)

    # Object identity hash (matches the Java tie-break ``g.hashCode()``), so that
    # distinct genotype objects remain distinct keys in the sorted map.
    def __hash__(self) -> int:
        return id(self)

    def __eq__(self, other: object) -> bool:
        return self is other

    def __repr__(self) -> str:
        return f"Genotype({self._sample_name}, {self.get_genotype_string()})"


class GenotypeBuilder:
    """Mutable builder mirroring htsjdk ``GenotypeBuilder`` (copy -> edit -> make)."""

    def __init__(self, genotype: Genotype):
        # Start from a deep-ish copy of the source genotype's mutable state.
        self._sample_name = genotype._sample_name
        self._alleles = list(genotype._alleles)
        self._attributes = dict(genotype._attributes)
        self._filters = genotype._filters
        self._ad = list(genotype._ad) if genotype._ad is not None else None

    def attribute(self, key: str, value: object) -> "GenotypeBuilder":
        """Set a FORMAT/extended attribute; returns self for chaining."""
        self._attributes[key] = value
        return self

    def filter(self, value: str) -> "GenotypeBuilder":
        """Set the per-genotype filter string; returns self for chaining."""
        self._filters = value
        return self

    def make(self) -> Genotype:
        return Genotype(self._sample_name, self._alleles, self._attributes, self._filters, self._ad)


# ---------------------------------------------------------------------------
# GenotypesContext  (htsjdk's ordered collection of a site's genotypes)
# ---------------------------------------------------------------------------
class GenotypesContext:
    """Ordered collection of genotypes keyed by sample name."""

    def __init__(self, genotypes: Optional[List[Genotype]] = None):
        # Preserve insertion order; dict keeps sample -> genotype.
        self._by_sample: Dict[str, Genotype] = {}
        for g in (genotypes or []):
            self._by_sample[g.get_sample_name()] = g

    @staticmethod
    def copy(context: "GenotypesContext") -> "GenotypesContext":
        """Java: ``GenotypesContext.copy(vc.getGenotypes())`` - a mutable copy."""
        new = GenotypesContext()
        new._by_sample = dict(context._by_sample)
        return new

    def replace(self, genotype: Genotype) -> None:
        """Replace the genotype for this genotype's sample."""
        self._by_sample[genotype.get_sample_name()] = genotype

    def iterate_in_sample_name_order(self) -> List[Genotype]:
        """Genotypes ordered by sample name (Java: iterateInSampleNameOrder)."""
        return [self._by_sample[s] for s in sorted(self._by_sample)]

    def get(self, sample_name: str) -> Optional[Genotype]:
        return self._by_sample.get(sample_name)

    def __iter__(self):
        return iter(self._by_sample.values())


# ---------------------------------------------------------------------------
# VariantContext  (immutable-ish; edits go through VariantContextBuilder)
# ---------------------------------------------------------------------------
class VariantContext:
    """A single VCF record (site) with its per-sample genotypes."""

    def __init__(
        self,
        contig: str,
        start: int,
        end: int,
        genotypes: GenotypesContext,
        filters: Optional[str] = None,
    ):
        self._contig = contig
        self._start = start
        self._end = end
        self._genotypes = genotypes
        self._filters = filters

    # --- genotype accessors ---
    def get_genotypes(self, sample_name: Optional[str] = None):
        """No arg -> the full GenotypesContext; a sample name -> a list with that
        sample's genotype (mirrors htsjdk's two ``getGenotypes`` overloads)."""
        if sample_name is None:
            return self._genotypes
        g = self._genotypes.get(sample_name)
        return [g] if g is not None else []

    def get_genotype(self, sample_name: str) -> Optional[Genotype]:
        return self._genotypes.get(sample_name)

    def get_sample_names_ordered_by_name(self) -> List[str]:
        return sorted(self._genotypes._by_sample)

    # --- filter / locus accessors ---
    def is_filtered(self) -> bool:
        return self._filters not in (None, "", VCFConstants.PASSES_FILTERS_v4, VCFConstants.UNFILTERED)

    def get_contig(self) -> str:
        return self._contig

    def get_start(self) -> int:
        return self._start

    def get_end(self) -> int:
        return self._end

    def to_string_without_genotypes(self) -> str:
        return f"{self._contig}:{self._start}-{self._end}"

    # Identity hash so VariantContexts can be used as dict keys (Java relies on
    # object identity for the ``newGenotypes`` map).
    def __hash__(self) -> int:
        return id(self)

    def __eq__(self, other: object) -> bool:
        return self is other


class VariantContextBuilder:
    """Mutable builder mirroring htsjdk ``VariantContextBuilder``."""

    def __init__(self, variant_context: VariantContext):
        self._contig = variant_context._contig
        self._start = variant_context._start
        self._end = variant_context._end
        self._genotypes = variant_context._genotypes
        self._filters = variant_context._filters

    def genotypes(self, genotypes) -> "VariantContextBuilder":
        """Accept either a GenotypesContext or a list of Genotype (htsjdk allows both)."""
        if isinstance(genotypes, GenotypesContext):
            self._genotypes = genotypes
        else:
            self._genotypes = GenotypesContext(list(genotypes))
        return self

    def filter(self, value: str) -> "VariantContextBuilder":
        self._filters = value
        return self

    def make(self) -> VariantContext:
        return VariantContext(self._contig, self._start, self._end, self._genotypes, self._filters)
