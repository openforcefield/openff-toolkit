"""Interfaces with OpenMM."""
from pathlib import Path
from typing import Union

import openmm

from openff.interchange.exceptions import PluginCompatibilityError
from openff.interchange.interop.openmm._import._import import from_openmm
from openff.interchange.interop.openmm._positions import to_openmm_positions
from openff.interchange.interop.openmm._topology import to_openmm_topology

__all__ = [
    "to_openmm",
    "to_openmm_topology",
    "to_openmm_positions",
    "from_openmm",
]

import openmm.app
from openff.toolkit.topology import Topology


def to_openmm(
    interchange,
    combine_nonbonded_forces: bool = False,
    add_constrained_forces: bool = False,
) -> openmm.System:
    """
    Convert an Interchange to an OpenmM System.

    Parameters
    ----------
    interchange : openff.interchange.Interchange
        An OpenFF Interchange object
    combine_nonbonded_forces : bool, default=False
        If True, an attempt will be made to combine all non-bonded interactions into a single openmm.NonbondedForce.
        If False, non-bonded interactions will be split across multiple forces.
    add_constrained_forces : bool, default=False,
        If True, add valence forces that might be overridden by constraints, i.e. call `addBond` or `addAngle`
        on a bond or angle that is fully constrained.

    Returns
    -------
    system : openmm.System
        The corresponding OpenMM System object

    """
    from openff.units import unit as off_unit

    from openff.interchange.interop.openmm._gbsa import _process_gbsa
    from openff.interchange.interop.openmm._nonbonded import _process_nonbonded_forces
    from openff.interchange.interop.openmm._valence import (
        _process_angle_forces,
        _process_bond_forces,
        _process_constraints,
        _process_improper_torsion_forces,
        _process_torsion_forces,
    )

    for collection in interchange.collections.values():
        if collection.is_plugin:
            try:
                collection.check_openmm_requirements(combine_nonbonded_forces)
            except AssertionError as error:
                raise PluginCompatibilityError(
                    f"Collection of type {type(collection)} failed a compatibility check.",
                ) from error

    system = openmm.System()

    if interchange.box is not None:
        box = interchange.box.m_as(off_unit.nanometer)
        system.setDefaultPeriodicBoxVectors(*box)

    particle_map = _process_nonbonded_forces(
        interchange,
        system,
        combine_nonbonded_forces=combine_nonbonded_forces,
    )

    constrained_pairs = _process_constraints(interchange, system, particle_map)

    _process_torsion_forces(interchange, system, particle_map)
    _process_improper_torsion_forces(interchange, system, particle_map)
    _process_angle_forces(
        interchange,
        system,
        add_constrained_forces=add_constrained_forces,
        constrained_pairs=constrained_pairs,
        particle_map=particle_map,
    )
    _process_bond_forces(
        interchange,
        system,
        add_constrained_forces=add_constrained_forces,
        constrained_pairs=constrained_pairs,
        particle_map=particle_map,
    )

    _process_gbsa(
        interchange,
        system,
    )

    for collection in interchange.collections.values():
        if collection.is_plugin:
            try:
                collection.modify_openmm_forces(
                    interchange,
                    system,
                    add_constrained_forces=add_constrained_forces,
                    constrained_pairs=constrained_pairs,
                    particle_map=particle_map,
                )
            except NotImplementedError:
                continue

    return system


def _to_pdb(file_path: Union[Path, str], topology: "Topology", positions):
    from openff.units.openmm import to_openmm

    openmm_topology = topology.to_openmm(ensure_unique_atom_names=False)

    positions = to_openmm(positions)

    with open(file_path, "w") as outfile:
        openmm.app.PDBFile.writeFile(openmm_topology, positions, outfile)
