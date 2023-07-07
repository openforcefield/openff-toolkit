"""
Helper functions for exporting positions to OpenMM.
"""

import openmm.unit
from openff.units import unit

from openff.interchange import Interchange
from openff.interchange.exceptions import MissingPositionsError


def to_openmm_positions(
    interchange: "Interchange",
    include_virtual_sites: bool = True,
) -> "openmm.unit.Quantity":
    """Generate an array of positions of all particles, optionally including virtual sites."""
    from collections import defaultdict

    import numpy

    if interchange.positions is None:
        raise MissingPositionsError(
            f"Positions are required, found {interchange.positions=}.",
        )

    if "VirtualSites" not in interchange.collections:
        return interchange.positions.to_openmm()
    elif len(interchange["VirtualSites"].key_map) == 0:
        return interchange.positions.to_openmm()

    topology = interchange.topology

    if include_virtual_sites:
        from openff.interchange.interop._virtual_sites import (
            _virtual_site_parent_molecule_mapping,
        )

        virtual_site_molecule_map = _virtual_site_parent_molecule_mapping(interchange)

        molecule_virtual_site_map = defaultdict(list)

        for virtual_site, molecule_index in virtual_site_molecule_map.items():
            molecule_virtual_site_map[molecule_index].append(virtual_site)

    particle_positions = unit.Quantity(
        numpy.empty(shape=(0, 3)),
        unit.nanometer,
    )

    for molecule in topology.molecules:
        molecule_index = topology.molecule_index(molecule)

        atom_indices = [topology.atom_index(atom) for atom in molecule.atoms]
        this_molecule_atom_positions = interchange.positions[atom_indices, :]

        n_virtual_sites_in_this_molecule: int = len(
            molecule_virtual_site_map[molecule_index],
        )
        this_molecule_virtual_site_positions = unit.Quantity(
            numpy.zeros((n_virtual_sites_in_this_molecule, 3)),
            unit.nanometer,
        )
        particle_positions = numpy.concatenate(
            [
                particle_positions,
                this_molecule_atom_positions,
                this_molecule_virtual_site_positions,
            ],
        )

    return particle_positions.to_openmm()
