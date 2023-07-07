"""Utilities for processing and interfacing with the OpenFF Toolkit."""
from typing import Union

import networkx as nx
import numpy as np
from openff.toolkit import Molecule, Topology
from openff.toolkit.topology._mm_molecule import _SimpleMolecule
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.utils.collections import ValidatedList
from openff.units import Quantity
from openmm.app import Topology as OpenMMTopology


def _get_num_h_bonds(topology: "Topology") -> int:
    """Get the number of (covalent) bonds containing a hydrogen atom."""
    n_bonds_containing_hydrogen = 0

    for bond in topology.bonds:
        if 1 in (bond.atom1.atomic_number, bond.atom2.atomic_number):
            n_bonds_containing_hydrogen += 1

    return n_bonds_containing_hydrogen


def _get_14_pairs(topology_or_molecule: Union["Topology", "Molecule"]):
    """Generate tuples of atom pairs, including symmetric duplicates."""
    # TODO: A replacement of Topology.nth_degree_neighbors in the toolkit
    #       may implement this in the future.
    for bond in topology_or_molecule.bonds:
        atom_i = bond.atom1
        atom_j = bond.atom2
        for atom_i_partner in atom_i.bonded_atoms:
            for atom_j_partner in atom_j.bonded_atoms:
                if atom_i_partner == atom_j_partner:
                    continue

                if atom_i_partner in atom_j_partner.bonded_atoms:
                    continue

                if atom_j_partner in atom_i_partner.bonded_atoms:
                    continue

                if {*atom_i_partner.bonded_atoms}.intersection(
                    {*atom_j_partner.bonded_atoms},
                ):
                    continue

                else:
                    atom_i_partner_index = topology_or_molecule.atom_index(
                        atom_i_partner,
                    )
                    atom_j_partner_index = topology_or_molecule.atom_index(
                        atom_j_partner,
                    )
                    if atom_i_partner_index > atom_j_partner_index:
                        yield (atom_j_partner, atom_i_partner)
                    else:
                        yield (atom_i_partner, atom_j_partner)


def _validated_list_to_array(validated_list: "ValidatedList") -> "Quantity":
    from openff.units import unit

    unit_ = validated_list[0].units
    return unit.Quantity(np.asarray([val.m for val in validated_list]), unit_)


def _combine_topologies(topology1: Topology, topology2: Topology) -> Topology:
    topology1_ = Topology(other=topology1)
    topology2_ = Topology(other=topology2)

    for molecule in topology2_.molecules:
        topology1_.add_molecule(molecule)

    return topology1_


def _check_electrostatics_handlers(force_field: "ForceField") -> bool:
    """
    Return whether or not this ForceField should have an Electrostatics tag.
    """
    # Manually-curated list of names of ParameterHandler classes that are expected
    # to assign/modify partial charges
    partial_charge_handlers = [
        "LibraryCharges",
        "ToolkitAM1BCC",
        "ChargeIncrementModel",
        "VirtualSites",
        "GBSA",
    ]

    # A more robust solution would probably take place late in the parameterization
    # process, but this solution should cover a vast majority of cases with minimal
    # complexity. Most notably this will *not* behave well with
    #   * parameter handler plugins
    #   * future additions to -built-in handlers
    #   * handlers that _could_ assign partial charges but happen to not assign
    #       any for some particular topology

    for parameter_handler_name in force_field.registered_parameter_handlers:
        if parameter_handler_name in partial_charge_handlers:
            return True

    return False


def _simple_topology_from_openmm(openmm_topology: "OpenMMTopology") -> Topology:
    """Convert an OpenMM Topology into an OpenFF Topology consisting **only** of so-called `_SimpleMolecule`s."""
    # TODO: Residue metadata
    # TODO: Splice in fully-defined OpenFF `Molecule`s?
    graph = nx.Graph()

    for atom in openmm_topology.atoms():
        graph.add_node(
            atom.index,
            atomic_number=atom.element.atomic_number,
        )

    for bond in openmm_topology.bonds():
        graph.add_edge(
            bond.atom1.index,
            bond.atom2.index,
        )

    return _simple_topology_from_graph(graph)


def _simple_topology_from_graph(graph: nx.Graph) -> Topology:
    topology = Topology()

    for component in nx.connected_components(graph):
        subgraph = graph.subgraph(component)
        topology.add_molecule(_SimpleMolecule._from_subgraph(subgraph))

    return topology


# This is to re-implement:
#   https://github.com/openforcefield/openff-toolkit/blob/60014820e6a333bed04e8bf5181d177da066da4d/
#   openff/toolkit/typing/engines/smirnoff/parameters.py#L2509-L2515
# It doesn't seem ideal to assume that matching SMILES === isomorphism?
class _HashedMolecule(Molecule):
    def __hash__(self):
        return hash(self.to_smiles())


def _assert_all_isomorphic(molecule_list: list[Molecule]) -> bool:
    hashed_molecules = {_HashedMolecule(molecule) for molecule in molecule_list}

    return len(hashed_molecules) == len(molecule_list)
