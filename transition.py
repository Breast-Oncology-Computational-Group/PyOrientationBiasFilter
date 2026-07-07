"""Stand-in for ``picard ...analysis.artifacts.Transition``.

In the Java code ``Transition`` is an enum of the 12 single-base substitutions
(A>C, A>G, ... T>G). It encodes an artifact "mode" as an ordered pair
(reference base -> called/alt base). Here it is a small value class that mirrors
the handful of members the orientation-bias filter actually uses:

  * ``transition_of(ref, alt)``  -> build from two bases          (Java: transitionOf)
  * ``ref()`` / ``call()``       -> the reference / alt base       (Java: ref() / call())
  * ``complement()``             -> the reverse-complement mode     (Java: complement())
  * ``value_of(name)``           -> parse "GtoT" style names        (Java: valueOf)
  * ``str(transition)``          -> "G>T" style label               (Java: toString())

Instances are immutable and compared/hashed by their (ref, call) pair, so they
can be used as dict keys and placed in sorted sets exactly like the Java enum.
"""

from __future__ import annotations

from functools import total_ordering

# Watson-Crick complement of each DNA base. Used to build the reverse-complement
# artifact mode (e.g. the OxoG mode G>T has complement C>A on the opposite strand).
_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


@total_ordering
class Transition:
    """A single-base substitution reference_base -> alt/called_base."""

    def __init__(self, ref_base: str, call_base: str):
        # Normalise to an uppercase single character, matching Picard's char usage.
        ref_base = str(ref_base).upper()
        call_base = str(call_base).upper()
        if ref_base not in _COMPLEMENT or call_base not in _COMPLEMENT:
            raise ValueError(f"Transition bases must be A/C/G/T, got {ref_base}>{call_base}")
        self._ref = ref_base
        self._call = call_base

    # --- factory methods -------------------------------------------------

    @staticmethod
    def transition_of(ref_base: str, call_base: str) -> "Transition":
        """Java: ``Transition.transitionOf(char ref, char alt)``."""
        return Transition(ref_base, call_base)

    @staticmethod
    def value_of(name: str) -> "Transition":
        """Java: ``Transition.valueOf("GtoT")`` -> the G>T transition.

        The summary-table reader passes names of the form ``"<ref>to<call>"``
        (it rewrites the human-readable "G>T" into "GtoT" first).
        """
        ref_base, _, call_base = name.partition("to")
        return Transition(ref_base, call_base)

    # --- accessors -------------------------------------------------------

    def ref(self) -> str:
        """Reference base of the substitution (Java: ``ref()``)."""
        return self._ref

    def call(self) -> str:
        """Alt / called base of the substitution (Java: ``call()``)."""
        return self._call

    def complement(self) -> "Transition":
        """Reverse-complement artifact mode (Java: ``complement()``)."""
        return Transition(_COMPLEMENT[self._ref], _COMPLEMENT[self._call])

    # --- identity: value-based equality / ordering / hashing ------------

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Transition) and (self._ref, self._call) == (other._ref, other._call)

    def __lt__(self, other: "Transition") -> bool:
        # Provides a deterministic ordering so Transitions can live in sorted sets
        # (Java used TreeSet/TreeMap keyed by the enum's natural order).
        return (self._ref, self._call) < (other._ref, other._call)

    def __hash__(self) -> int:
        return hash((self._ref, self._call))

    def name(self) -> str:
        """Enum-style name, e.g. ``"GtoT"`` (Java: ``Transition.name()``)."""
        return f"{self._ref}to{self._call}"

    def __str__(self) -> str:
        # Java toString() renders the mode as "REF>CALL", e.g. "G>T".
        return f"{self._ref}>{self._call}"

    def __repr__(self) -> str:
        return f"Transition({self._ref!r}, {self._call!r})"
