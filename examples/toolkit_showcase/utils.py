"""Utility methods for the toolkit_showcase example"""

from simtk import openmm


def off_topology_to_omm(off_topology):
    """Convert openff.toolkit.topology.Topology to simtk.openmm.app.topology.Topology"""
    omm_topology = openmm.app.topology.Topology()
    chain = omm_topology.addChain()
    res = omm_topology.addResidue("LIG", chain)
    atoms = {}

    for atom in off_topology.topology_atoms:
        atom_name = (
            atom.atom.name
            if atom.atom.name
            else str(atom.element.symbol) + str(atom.topology_atom_index)
        )

        atoms[atom.topology_atom_index] = omm_topology.addAtom(
            atom_name, atom.element, res
        )

    for bond in off_topology.topology_bonds:
        (a, b) = bond.atoms
        omm_topology.addBond(
            atoms[a.topology_atom_index],
            atoms[b.topology_atom_index],
            order=bond.bond_order,
        )

    return omm_topology
