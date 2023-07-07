from openff.toolkit.topology._mm_molecule import _SimpleMolecule
from openff.toolkit.topology.topology import Topology

from openff.interchange.interop.gromacs._import._import import GROMACSSystem


def _create_topology_from_system(system: GROMACSSystem) -> Topology:
    pass

    topology = Topology()

    for molecule_name, gromacs_molecule in system.molecule_types.items():
        molecule = _SimpleMolecule()
        molecule.name = molecule_name

        for atom in gromacs_molecule.atoms:
            atom_type = atom.atom_type
            atomic_number = system.atom_types[atom_type].atomic_number

            molecule.add_atom(atomic_number=atomic_number)

        for bond in gromacs_molecule.bonds:
            molecule.add_bond(bond.atom1 - 1, bond.atom2 - 1)

        for settle in gromacs_molecule.settles:
            # Assume this is a water, ordered as O H H
            assert len(gromacs_molecule.atoms) == 3
            assert molecule.atoms[0].atomic_number == 8
            assert molecule.atoms[1].atomic_number == 1
            assert molecule.atoms[2].atomic_number == 1

            # How to handle H O H? The documentation is unclear about what the first
            # value index must refer to - oxygen or just whatever the first atom is?
            #  > The [ settles ] directive defines the first atom of the water molecule.
            #   https://manual.gromacs.org/current/reference-manual/topologies/constraint-algorithm-section.html#constraintalg
            assert settle.first_atom == 1

            molecule.add_bond(molecule.atom(0), molecule.atom(1))
            molecule.add_bond(molecule.atom(0), molecule.atom(2))

        for _ in range(system.molecules[molecule_name]):
            topology.add_molecule(molecule)

        # TODO: Add residue metadata

    return topology
