from collections.abc import Iterable
from typing import Literal

from openff.toolkit.topology.molecule import Atom
from pydantic import Field

from openff.interchange.components.potentials import Collection


class ConstraintCollection(Collection):
    """Collection storing constraint potentials as produced by a SMIRNOFF force field."""

    type: Literal["Constraints"] = "Constraints"
    expression: Literal[""] = ""

    @classmethod
    def potential_parameters(cls) -> Iterable[str]:
        """Return a list of names of parameters included in each potential in this colletion."""
        return ("length",)


class BondCollection(Collection):
    """Collection storing bond potentials."""

    type: Literal["Bonds"] = "Bonds"
    expression: Literal["k/2*(r-length)**2"] = "k/2*(r-length)**2"

    @classmethod
    def potential_parameters(cls) -> Iterable[str]:
        """Return a list of names of parameters included in each potential in this colletion."""
        return "k", "length"

    @classmethod
    def valence_terms(cls, topology) -> list[tuple["Atom", ...]]:
        """Return all bonds in this topology."""
        return [tuple(b.atoms) for b in topology.bonds]


class AngleCollection(Collection):
    """Collection storing Angle potentials."""

    type: Literal["Angles"] = "Angles"
    expression: Literal["k/2*(theta-angle)**2"] = "k/2*(theta-angle)**2"

    @classmethod
    def potential_parameters(cls) -> Iterable[str]:
        """Return a list of names of parameters included in each potential in this colletion."""
        return "k", "angle"

    @classmethod
    def valence_terms(cls, topology):
        """Return all angles in this topology."""
        return [angle for angle in topology.angles]


class ProperTorsionCollection(Collection):
    """Handler storing periodic proper torsion potentials."""

    type: Literal["ProperTorsions"] = "ProperTorsions"
    expression: Literal[
        "k*(1+cos(periodicity*theta-phase))"
    ] = "k*(1+cos(periodicity*theta-phase))"

    @classmethod
    def supported_parameters(cls) -> Iterable[str]:
        """Return a list of supported parameter attribute names."""
        return "k", "periodicity", "phase"


class RyckaertBellemansTorsionCollection(Collection):
    """Handler storing Ryckaert-Bellemans torsion potentials."""

    type: Literal["RBTorsions"] = "RBTorsions"
    expression: str = Field(
        "c0 + "
        "c1 * (cos(phi - 180)) "
        "c2 * (cos(phi - 180)) ** 2 + "
        "c3 * (cos(phi - 180)) ** 3 + "
        "c4 * (cos(phi - 180)) ** 4 + "
        "c5 * (cos(phi - 180)) ** 5",
    )

    @classmethod
    def supported_parameters(cls) -> Iterable[str]:
        """Return a list of supported parameter attribute names."""
        return "c0", "c1", "c2", "c3", "c4", "c5"


class ImproperTorsionCollection(Collection):
    """Handler storing periodic improper torsion potentials."""

    type: Literal["ImproperTorsions"] = "ImproperTorsions"
    expression: Literal[
        "k*(1+cos(periodicity*theta-phase))"
    ] = "k*(1+cos(periodicity*theta-phase))"

    @classmethod
    def supported_parameters(cls) -> Iterable[str]:
        """Return a list of supported parameter attribute names."""
        return "k", "periodicity", "phase"
