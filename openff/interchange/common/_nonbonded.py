import abc
from collections import defaultdict
from collections.abc import Iterable
from typing import DefaultDict, Literal, Optional, Union

from openff.models.types import FloatQuantity
from openff.units import Quantity, unit
from pydantic import Field, PrivateAttr

from openff.interchange.components.potentials import Collection
from openff.interchange.constants import _PME
from openff.interchange.models import LibraryChargeTopologyKey, TopologyKey


class _NonbondedCollection(Collection, abc.ABC):  # noqa
    type: str = "nonbonded"

    cutoff: FloatQuantity["angstrom"] = Field(  # noqa
        unit.Quantity(9.0, unit.angstrom),
        description="The distance at which pairwise interactions are truncated",
    )

    scale_13: float = Field(
        0.0,
        description="The scaling factor applied to 1-3 interactions",
    )
    scale_14: float = Field(
        0.5,
        description="The scaling factor applied to 1-4 interactions",
    )
    scale_15: float = Field(
        1.0,
        description="The scaling factor applied to 1-5 interactions",
    )


class vdWCollection(_NonbondedCollection):
    """Handler storing vdW potentials."""

    type: Literal["vdW"] = "vdW"

    expression: Literal[
        "4*epsilon*((sigma/r)**12-(sigma/r)**6)"
    ] = "4*epsilon*((sigma/r)**12-(sigma/r)**6)"

    mixing_rule: str = Field(
        "lorentz-berthelot",
        description="The mixing rule (combination rule) used in computing pairwise vdW interactions",
    )

    switch_width: FloatQuantity["angstrom"] = Field(  # noqa
        unit.Quantity(1.0, unit.angstrom),
        description="The width over which the switching function is applied",
    )

    method: Literal["cutoff", "no-cutoff", "pme"] = Field("cutoff")

    @classmethod
    def default_parameter_values(cls) -> Iterable[float]:
        """Per-particle parameter values passed to Force.addParticle()."""
        return 1.0, 0.0

    @classmethod
    def potential_parameters(cls) -> Iterable[str]:
        """Return a list of names of parameters included in each potential in this colletion."""
        return "sigma", "epsilon"


class ElectrostaticsCollection(_NonbondedCollection):
    """Handler storing electrostatics interactions."""

    type: Literal["Electrostatics"] = "Electrostatics"

    expression: Literal["coul"] = "coul"

    periodic_potential: Literal[
        "Ewald3D-ConductingBoundary",
        "cutoff",
        "no-cutoff",
    ] = Field(_PME)
    nonperiodic_potential: Literal["Coulomb", "cutoff", "no-cutoff"] = Field("Coulomb")
    exception_potential: Literal["Coulomb"] = Field("Coulomb")

    _charges: dict[
        Union[TopologyKey, LibraryChargeTopologyKey],
        Quantity,
    ] = PrivateAttr(dict())
    _charges_cached_with_virtual_sites: Optional[bool] = PrivateAttr(None)

    @property
    def charges(self) -> dict[Union[TopologyKey, LibraryChargeTopologyKey], Quantity]:
        """Get the total partial charge on each atom, excluding virtual sites."""
        if self._charges is None or self._charges_cached_with_virtual_sites in (
            True,
            None,
        ):
            self._charges = self._get_charges(include_virtual_sites=False)
            self._charges_cached_with_virtual_sites = False

        return self._charges

    @property
    def charges_with_virtual_sites(
        self,
    ) -> dict[Union[TopologyKey, LibraryChargeTopologyKey], Quantity]:
        """Get the total partial charge on each atom, including virtual sites."""
        raise NotImplementedError()

    def _get_charges(
        self,
        include_virtual_sites: bool = False,
    ) -> dict[Union[TopologyKey, LibraryChargeTopologyKey], Quantity]:
        """Get the total partial charge on each atom or particle."""
        if include_virtual_sites:
            raise NotImplementedError()

        charges: DefaultDict[
            Union[TopologyKey, LibraryChargeTopologyKey],
            Quantity,
        ] = defaultdict(
            lambda: Quantity(0.0, unit.elementary_charge),
        )

        for topology_key, potential_key in self.key_map.items():
            potential = self.potentials[potential_key]

            charges[topology_key] = potential.parameters["charge"]

        return charges
