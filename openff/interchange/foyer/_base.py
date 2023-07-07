from abc import abstractmethod
from copy import copy
from typing import Optional

from openff.toolkit import Topology
from openff.utilities.utilities import has_package

from openff.interchange.components.potentials import Collection, Potential
from openff.interchange.models import PotentialKey, TopologyKey

if has_package("foyer"):
    from foyer.exceptions import MissingForceError, MissingParametersError
    from foyer.forcefield import Forcefield


# Is this the safest way to achieve PotentialKey id separation?
POTENTIAL_KEY_SEPARATOR = "-"


def _copy_params(
    params: dict[str, float],
    *drop_keys: str,
    param_units: Optional[dict] = None,
) -> dict:
    """Copy parameters from a dictionary."""
    params_copy = copy(params)
    for drop_key in drop_keys:
        params_copy.pop(drop_key, None)
    if param_units:
        for unit_item, units in param_units.items():
            params_copy[unit_item] = params_copy[unit_item] * units
    return params_copy


def _get_potential_key_id(atom_slots: dict[TopologyKey, PotentialKey], idx):
    """From a dictionary of TopologyKey: PotentialKey, get the PotentialKey id."""
    top_key = TopologyKey(atom_indices=(idx,))
    return atom_slots[top_key].id


class FoyerConnectedAtomsHandler(Collection):
    """Base class for handlers storing valence potentials produced by a Foyer force field."""

    force_field_key: str
    connection_attribute: str = ""
    raise_on_missing_params: bool = True

    def store_matches(
        self,
        atom_slots: dict[TopologyKey, PotentialKey],
        topology: "Topology",
    ) -> None:
        """Populate self.key_map with key-val pairs of [TopologyKey, PotentialKey]."""
        for connection in getattr(topology, self.connection_attribute):
            try:
                atoms_iterable = connection.atoms
            except AttributeError:
                atoms_iterable = connection
            atom_indices = tuple(topology.atom_index(atom) for atom in atoms_iterable)

            top_key = TopologyKey(atom_indices=atom_indices)
            pot_key_ids = tuple(
                _get_potential_key_id(atom_slots, idx) for idx in atom_indices
            )

            self.key_map[top_key] = PotentialKey(
                id=POTENTIAL_KEY_SEPARATOR.join(pot_key_ids),
            )

    def store_potentials(self, force_field: "Forcefield") -> None:
        """Populate self.potentials with key-val pairs of [PotentialKey, Potential]."""
        for pot_key in self.key_map.values():
            try:
                params = force_field.get_parameters(
                    self.force_field_key,
                    key=pot_key.id.split(POTENTIAL_KEY_SEPARATOR),
                )
                params = self.get_params_with_units(params)
                self.potentials[pot_key] = Potential(parameters=params)
            except MissingForceError:
                # Here, we can safely assume that the ForceGenerator is Missing
                self.key_map = {}
                self.potentials = {}
                return
            except MissingParametersError as e:
                if self.raise_on_missing_params:
                    raise e
                else:
                    pass

    @abstractmethod
    def get_params_with_units(self, params):
        """Get the parameters of this handler, tagged with units."""
        raise NotImplementedError
