from openff.toolkit.topology import Topology
from openff.units import unit
from pydantic import Field

from openff.interchange.common._valence import (
    AngleCollection,
    BondCollection,
    ProperTorsionCollection,
    RyckaertBellemansTorsionCollection,
)
from openff.interchange.foyer._base import (
    POTENTIAL_KEY_SEPARATOR,
    FoyerConnectedAtomsHandler,
    _copy_params,
    _get_potential_key_id,
)
from openff.interchange.models import PotentialKey, TopologyKey


class FoyerHarmonicBondHandler(FoyerConnectedAtomsHandler, BondCollection):
    """Handler storing bond potentials as produced by a Foyer force field."""

    type = "Bonds"
    expression = "k/2*(r-length)**2"
    force_field_key = "harmonic_bonds"
    connection_attribute = "bonds"

    def get_params_with_units(self, params):
        """Get the parameters of this handler, tagged with units."""
        return _copy_params(
            params,
            param_units={"k": unit.kJ / unit.mol / unit.nm**2, "length": unit.nm},
        )

    def store_matches(
        self,
        atom_slots: dict[TopologyKey, PotentialKey],
        topology: "Topology",
    ) -> None:
        """Populate self.key_map with key-val pairs of [TopologyKey, PotentialKey]."""
        for bond in topology.bonds:
            atom_indices = (
                topology.atom_index(bond.atom1),
                topology.atom_index(bond.atom2),
            )
            top_key = TopologyKey(atom_indices=atom_indices)

            pot_key_ids = tuple(
                _get_potential_key_id(atom_slots, idx) for idx in atom_indices
            )

            self.key_map[top_key] = PotentialKey(
                id=POTENTIAL_KEY_SEPARATOR.join(pot_key_ids),
            )


class FoyerHarmonicAngleHandler(FoyerConnectedAtomsHandler, AngleCollection):
    """Handler storing angle potentials as produced by a Foyer force field."""

    type = "Angles"
    expression = "k/2*(theta-angle)**2"
    force_field_key: str = "harmonic_angles"
    connection_attribute: str = "angles"

    def get_params_with_units(self, params):
        """Get the parameters of this handler, tagged with units."""
        return _copy_params(
            {"k": params["k"], "angle": params["theta"]},
            param_units={
                "k": unit.kJ / unit.mol / unit.radian**2,
                "angle": unit.dimensionless,
            },
        )

    def store_matches(
        self,
        atom_slots: dict[TopologyKey, PotentialKey],
        topology: "Topology",
    ) -> None:
        """Populate self.key_map with key-val pairs of [TopologyKey, PotentialKey]."""
        for angle in topology.angles:
            atom_indices = tuple(topology.atom_index(atom) for atom in angle)
            top_key = TopologyKey(atom_indices=atom_indices)

            pot_key_ids = tuple(
                _get_potential_key_id(atom_slots, idx) for idx in atom_indices
            )

            self.key_map[top_key] = PotentialKey(
                id=POTENTIAL_KEY_SEPARATOR.join(pot_key_ids),
            )


class FoyerRBProperHandler(
    FoyerConnectedAtomsHandler,
    RyckaertBellemansTorsionCollection,
):
    """Handler storing Ryckaert-Bellemans proper torsion potentials as produced by a Foyer force field."""

    force_field_key = "rb_propers"
    type = "RBTorsions"
    expression: str = Field(
        "c0 + "
        "c1 * (cos(phi - 180)) "
        "c2 * (cos(phi - 180)) ** 2 + "
        "c3 * (cos(phi - 180)) ** 3 + "
        "c4 * (cos(phi - 180)) ** 4 + "
        "c5 * (cos(phi - 180)) ** 5",
    )
    connection_attribute: str = "propers"
    raise_on_missing_params: bool = False

    def get_params_with_units(self, params):
        """Get the parameters of this handler, tagged with units."""
        rb_params = params
        param_units = {k: unit.kJ / unit.mol for k in rb_params}
        return _copy_params(rb_params, param_units=param_units)

    def store_matches(
        self,
        atom_slots: dict[TopologyKey, PotentialKey],
        topology: "Topology",
    ) -> None:
        """Populate self.key_map with key-val pairs of [TopologyKey, PotentialKey]."""
        for proper in topology.propers:
            atom_indices = tuple(topology.atom_index(atom) for atom in proper)
            top_key = TopologyKey(atom_indices=atom_indices)

            pot_key_ids = tuple(
                _get_potential_key_id(atom_slots, idx) for idx in atom_indices
            )

            self.key_map[top_key] = PotentialKey(
                id=POTENTIAL_KEY_SEPARATOR.join(pot_key_ids),
            )


class FoyerRBImproperHandler(FoyerRBProperHandler):
    """Handler storing Ryckaert-Bellemans improper torsion potentials as produced by a Foyer force field."""

    type = "RBImpropers"  # type: ignore[assignment]
    connection_attribute: str = "impropers"


class FoyerPeriodicProperHandler(FoyerConnectedAtomsHandler, ProperTorsionCollection):
    """Handler storing periodic proper torsion potentials as produced by a Foyer force field."""

    force_field_key = "periodic_propers"
    connection_attribute: str = "propers"
    raise_on_missing_params: bool = False
    type = "ProperTorsions"
    expression = "k*(1+cos(periodicity*theta-phase))"

    def get_params_with_units(self, params):
        """Get the parameters of this handler, tagged with units."""
        return _copy_params(
            params,
            param_units={
                "k": unit.kJ / unit.mol / unit.nm**2,
                "phase": unit.dimensionless,
                "periodicity": unit.dimensionless,
            },
        )


class FoyerPeriodicImproperHandler(FoyerPeriodicProperHandler):
    """Handler storing periodic improper torsion potentials as produced by a Foyer force field."""

    type = "ImproperTorsions"  # type: ignore[assignment]
    connection_attribute: str = "impropers"


class _RBTorsionHandler(RyckaertBellemansTorsionCollection):
    pass
