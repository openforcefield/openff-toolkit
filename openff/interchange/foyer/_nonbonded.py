from typing import Literal

from openff.models.types import FloatQuantity
from openff.toolkit.topology import Topology
from openff.units import Quantity, unit
from openff.utilities.utilities import has_package
from pydantic import Field, PrivateAttr

from openff.interchange.common._nonbonded import ElectrostaticsCollection, vdWCollection
from openff.interchange.components.potentials import Potential
from openff.interchange.foyer._base import _copy_params
from openff.interchange.models import PotentialKey, TopologyKey

if has_package("foyer"):
    from foyer.forcefield import Forcefield


class FoyerVDWHandler(vdWCollection):
    """Handler storing vdW potentials as produced by a Foyer force field."""

    force_field_key: str = "atoms"
    mixing_rule: Literal["lorentz-berthelot", "geometric"] = Field("geometric")
    method: Literal["cutoff", "pme", "no-cutoff"] = Field("cutoff")

    def store_matches(
        self,
        force_field: "Forcefield",
        topology: "Topology",
    ) -> None:
        """Populate self.key_map with key-val pairs of [TopologyKey, PotentialKey]."""
        from foyer.atomtyper import find_atomtypes
        from foyer.topology_graph import TopologyGraph

        top_graph = TopologyGraph.from_openff_topology(topology)

        type_map = find_atomtypes(top_graph, forcefield=force_field)
        for key, val in type_map.items():
            top_key = TopologyKey(atom_indices=(key,))
            self.key_map[top_key] = PotentialKey(id=val["atomtype"])

    def store_potentials(self, force_field: "Forcefield") -> None:
        """Extract specific force field potentials a Forcefield object."""
        for top_key in self.key_map:
            atom_params = force_field.get_parameters(
                self.force_field_key,
                key=self.key_map[top_key].id,
            )

            atom_params = _copy_params(
                atom_params,
                "charge",
                param_units={"epsilon": unit.kJ / unit.mol, "sigma": unit.nm},
            )

            self.potentials[self.key_map[top_key]] = Potential(parameters=atom_params)


class FoyerElectrostaticsHandler(ElectrostaticsCollection):
    """Handler storing electrostatics potentials as produced by a Foyer force field."""

    force_field_key: str = "atoms"
    cutoff: FloatQuantity["angstrom"] = 9.0 * unit.angstrom

    _charges: dict[TopologyKey, Quantity] = PrivateAttr(dict())  # type: ignore

    @property
    def charges(self) -> dict[TopologyKey, Quantity]:  # type: ignore[override]
        """Get the total partial charge on each atom, including virtual sites."""
        return self._charges

    @property
    def charges_with_virtual_sites(self):
        """Get the total partial charge on each atom, including virtual sites."""
        raise NotImplementedError("Foyer force fields do not support virtual sites.")

    def store_charges(
        self,
        atom_slots: dict[TopologyKey, PotentialKey],
        force_field: "Forcefield",
    ):
        """Look up fixed charges (a.k.a. library charges) from the force field and store them in self._charges."""
        for top_key, pot_key in atom_slots.items():
            foyer_params = force_field.get_parameters(self.force_field_key, pot_key.id)

            charge = Quantity(foyer_params["charge"], unit.elementary_charge)

            self._charges[top_key] = charge

            self.key_map[top_key] = pot_key
            self.potentials[pot_key] = Potential(
                parameters={"charge": charge},
            )
