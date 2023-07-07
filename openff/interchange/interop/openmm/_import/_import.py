from typing import TYPE_CHECKING, Optional

import openmm
import openmm.app

from openff.interchange._experimental import experimental
from openff.interchange.common._nonbonded import ElectrostaticsCollection, vdWCollection
from openff.interchange.common._valence import (
    AngleCollection,
    BondCollection,
    ConstraintCollection,
    ProperTorsionCollection,
)
from openff.interchange.exceptions import UnsupportedImportError

if TYPE_CHECKING:
    from openff.interchange import Interchange


@experimental
def from_openmm(
    topology: Optional["openmm.app.Topology"] = None,
    system: Optional[openmm.System] = None,
    positions=None,
    box_vectors=None,
) -> "Interchange":
    """Create an Interchange object from OpenMM data."""
    from openff.interchange import Interchange

    interchange = Interchange()

    if system:
        constraints = _convert_constraints(system)

        if constraints is not None:
            interchange.collections["Constraints"] = constraints

        for force in system.getForces():
            if isinstance(force, openmm.NonbondedForce):
                vdw, coul = _convert_nonbonded_force(force)
                interchange.collections["vdW"] = vdw
                interchange.collections["Electrostatics"] = coul
            elif isinstance(force, openmm.HarmonicBondForce):
                bonds = _convert_harmonic_bond_force(force)
                interchange.collections["Bonds"] = bonds
            elif isinstance(force, openmm.HarmonicAngleForce):
                angles = _convert_harmonic_angle_force(force)
                interchange.collections["Angles"] = angles
            elif isinstance(force, openmm.PeriodicTorsionForce):
                proper_torsions = _convert_periodic_torsion_force(force)
                interchange.collections["ProperTorsions"] = proper_torsions
            elif isinstance(force, openmm.CMMotionRemover):
                pass
            else:
                raise UnsupportedImportError(
                    f"Unsupported OpenMM Force type ({type(force)}) found.",
                )

    if topology is not None:
        from openff.interchange.components.toolkit import _simple_topology_from_openmm

        openff_topology = _simple_topology_from_openmm(topology)

        interchange.topology = openff_topology

    if positions is not None:
        interchange.positions = positions

    if box_vectors is not None:
        interchange.box = box_vectors

    return interchange


def _convert_constraints(
    system: openmm.System,
) -> Optional[ConstraintCollection]:
    from openff.units import unit

    from openff.interchange.components.potentials import Potential
    from openff.interchange.models import BondKey, PotentialKey

    if system.getNumConstraints() == 0:
        return None

    constraints = ConstraintCollection()

    # Map the unique distances (float, implicitly nanometer) to indices used for deduplication
    unique_distances: dict[float, int] = {
        distance: index
        for index, distance in enumerate(
            {
                system.getConstraintParameters(index)[2].value_in_unit(
                    openmm.unit.nanometer,
                )
                for index in range(system.getNumConstraints())
            },
        )
    }

    _keys: dict[float, PotentialKey] = dict()

    for distance, index in unique_distances.items():
        potential_key = PotentialKey(id=f"Constraint{index}")
        _keys[distance] = potential_key
        constraints.potentials[potential_key] = Potential(
            parameters={"distance": distance * unit.nanometer},
        )

    for index in range(system.getNumConstraints()):
        atom1, atom2, _distance = system.getConstraintParameters(index)

        distance = _distance.value_in_unit(openmm.unit.nanometer)

        constraints.key_map[BondKey(atom_indices=(atom1, atom2))] = _keys[distance]

    return constraints


def _convert_nonbonded_force(
    force: openmm.NonbondedForce,
) -> tuple["vdWCollection", "ElectrostaticsCollection"]:
    from openff.units.openmm import from_openmm as from_openmm_quantity

    from openff.interchange.common._nonbonded import (
        ElectrostaticsCollection,
        vdWCollection,
    )
    from openff.interchange.components.potentials import Potential
    from openff.interchange.models import PotentialKey, TopologyKey

    if force.getNonbondedMethod() != 4:
        raise UnsupportedImportError(
            "Importing from OpenMM only currently supported with `openmm.NonbondedForce.PME`.",
        )

    vdw = vdWCollection()
    electrostatics = ElectrostaticsCollection(version=0.4, scale_14=0.833333)

    n_parametrized_particles = force.getNumParticles()

    for idx in range(n_parametrized_particles):
        charge, sigma, epsilon = force.getParticleParameters(idx)
        top_key = TopologyKey(atom_indices=(idx,))
        pot_key = PotentialKey(id=f"{idx}")
        pot = Potential(
            parameters={
                "sigma": from_openmm_quantity(sigma),
                "epsilon": from_openmm_quantity(epsilon),
            },
        )
        vdw.key_map.update({top_key: pot_key})
        vdw.potentials.update({pot_key: pot})

        electrostatics.key_map.update({top_key: pot_key})
        electrostatics.potentials.update(
            {pot_key: Potential(parameters={"charge": from_openmm_quantity(charge)})},
        )

    if force.getNonbondedMethod() == 4:
        vdw.cutoff = force.getCutoffDistance()
    else:
        raise UnsupportedImportError(
            f"Parsing a non-bonded force of type {type(force)} with {force.getNonbondedMethod()} not yet supported.",
        )

    if force.getUseSwitchingFunction():
        vdw.switch_width = vdw.cutoff - from_openmm_quantity(
            force.getSwitchingDistance(),
        )
    else:
        vdw.switch_width = 0.0 * vdw.cutoff.units

    return vdw, electrostatics


def _convert_harmonic_bond_force(
    force: openmm.HarmonicBondForce,
) -> "BondCollection":
    from openff.units.openmm import from_openmm as from_openmm_quantity

    from openff.interchange.common._valence import BondCollection
    from openff.interchange.components.potentials import Potential
    from openff.interchange.models import BondKey, PotentialKey

    bonds = BondCollection()

    n_parametrized_bonds = force.getNumBonds()

    for idx in range(n_parametrized_bonds):
        atom1, atom2, length, k = force.getBondParameters(idx)
        top_key = BondKey(atom_indices=(atom1, atom2))
        pot_key = PotentialKey(id=f"{atom1}-{atom2}")
        pot = Potential(
            parameters={
                "length": from_openmm_quantity(length),
                "k": from_openmm_quantity(k),
            },
        )

        bonds.key_map.update({top_key: pot_key})
        bonds.potentials.update({pot_key: pot})

    return bonds


def _convert_harmonic_angle_force(
    force: openmm.HarmonicAngleForce,
) -> "AngleCollection":
    from openff.units.openmm import from_openmm as from_openmm_quantity

    from openff.interchange.common._valence import AngleCollection
    from openff.interchange.components.potentials import Potential
    from openff.interchange.models import AngleKey, PotentialKey

    angles = AngleCollection()

    n_parametrized_angles = force.getNumAngles()

    for idx in range(n_parametrized_angles):
        atom1, atom2, atom3, angle, k = force.getAngleParameters(idx)
        top_key = AngleKey(atom_indices=(atom1, atom2, atom3))
        pot_key = PotentialKey(id=f"{atom1}-{atom2}-{atom3}")
        pot = Potential(
            parameters={
                "angle": from_openmm_quantity(angle),
                "k": from_openmm_quantity(k),
            },
        )

        angles.key_map.update({top_key: pot_key})
        angles.potentials.update({pot_key: pot})

    return angles


def _convert_periodic_torsion_force(
    force: openmm.PeriodicTorsionForce,
) -> "ProperTorsionCollection":
    # TODO: Can impropers be separated out from a PeriodicTorsionForce?
    # Maybe by seeing if a quartet is in mol/top.propers or .impropers
    from openff.units import unit
    from openff.units.openmm import from_openmm as from_openmm_quantity

    from openff.interchange.common._valence import ProperTorsionCollection
    from openff.interchange.components.potentials import Potential
    from openff.interchange.models import PotentialKey, ProperTorsionKey

    proper_torsions = ProperTorsionCollection()

    n_parametrized_torsions = force.getNumTorsions()

    for idx in range(n_parametrized_torsions):
        atom1, atom2, atom3, atom4, per, phase, k = force.getTorsionParameters(idx)
        # TODO: Process layered torsions
        # TODO: Check if this torsion is an improper
        top_key = ProperTorsionKey(atom_indices=(atom1, atom2, atom3, atom4), mult=0)
        while top_key in proper_torsions.key_map:
            top_key.mult = top_key.mult + 1  # type: ignore[operator]

        pot_key = PotentialKey(id=f"{atom1}-{atom2}-{atom3}-{atom4}", mult=top_key.mult)
        pot = Potential(
            parameters={
                "periodicity": int(per) * unit.dimensionless,
                "phase": from_openmm_quantity(phase),
                "k": from_openmm_quantity(k),
                "idivf": 1 * unit.dimensionless,
            },
        )

        proper_torsions.key_map.update({top_key: pot_key})
        proper_torsions.potentials.update({pot_key: pot})

    return proper_torsions
