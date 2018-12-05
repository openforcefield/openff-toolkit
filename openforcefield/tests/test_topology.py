#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Tests for Topology

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import pickle
from functools import partial
from unittest import TestCase

from openforcefield import utils, topology

import pytest
import copy
import numpy as np
from simtk import unit
from openforcefield.utils import BASIC_CHEMINFORMATICS_TOOLKITS, RDKIT_AVAILABLE, OPENEYE_AVAILABLE, AMBERTOOLS_AVAILABLE, RDKitToolkitWrapper, OpenEyeToolkitWrapper, AmberToolsToolkitWrapper
from openforcefield.tests.utils import get_data_filename
from openforcefield.topology import Topology, TopologyAtom, TopologyBond, TopologyMolecule, TopologyVirtualSite
from openforcefield.topology import Molecule

#=============================================================================================
# TESTS
#=============================================================================================


# IF we've done our jobs right, it shouldn't matter which toolkit the tests for Topology run using (both's behaviors
# should be indistinguishable)
def test_cheminformatics_toolkit_is_installed():
    """Ensure that at least one supported cheminformatics toolkit is installed."""
    if not(RDKIT_AVAILABLE) and not(OPENEYE_AVAILABLE):
        msg = 'No supported cheminformatics toolkits are installed. Please install a supported toolkit:\n'
        msg += str(BASIC_CHEMINFORMATICS_TOOLKITS)
        raise Exception(msg)

class TestTopology(TestCase):
    from openforcefield.topology import Topology

    def setUp(self):
        self.empty_molecule = Molecule()
        self.ethane_from_smiles = Molecule.from_smiles('CC')
        self.ethene_from_smiles = Molecule.from_smiles('C=C')
        self.propane_from_smiles = Molecule.from_smiles('CCC')

        filename = get_data_filename('molecules/toluene.sdf')
        self.toluene_from_sdf = Molecule.from_file(filename)
        filename = get_data_filename('molecules/toluene_charged.mol2')
        # TODO: This will require openeye to load
        self.toluene_from_charged_mol2 = Molecule.from_file(filename)
        self.charged_methylamine_from_smiles = Molecule.from_smiles('[H]C([H])([H])[N+]([H])([H])[H]')

        molecule = Molecule.from_smiles('CC')
        carbons = [atom for atom in molecule.atoms if atom.atomic_number == 6]
        c0_hydrogens = [atom for atom in carbons[0].bonded_atoms if atom.atomic_number == 1]
        molecule.add_bond_charge_virtual_site((carbons[0], carbons[1]),
                                              0.1*unit.angstrom,
                                              charge_increments=[0.1, 0.05] * unit.elementary_charge
                                             )
        molecule.add_monovalent_lone_pair_virtual_site((c0_hydrogens[0], carbons[0], carbons[1]),
                                                       0.2*unit.angstrom,
                                                       20*unit.degree,
                                                       25*unit.degree,
                                                       charge_increments=[0.01, 0.02, 0.03] * unit.elementary_charge
                                                       )
        self.ethane_from_smiles_w_vsites = Molecule(molecule)

        # Make a propane with virtual sites
        molecule = Molecule.from_smiles('CCC')
        carbons = [atom for atom in molecule.atoms if atom.atomic_number == 6]
        c0_hydrogens = [atom for atom in carbons[0].bonded_atoms if atom.atomic_number == 1]
        molecule.add_bond_charge_virtual_site((carbons[0], carbons[1]),
                                              0.1*unit.angstrom,
                                              charge_increments=[0.1, 0.05] * unit.elementary_charge
                                             )
        molecule.add_monovalent_lone_pair_virtual_site((c0_hydrogens[0], carbons[0], carbons[1]),
                                                       0.2*unit.angstrom,
                                                       20*unit.degree,
                                                       25*unit.degree,
                                                       charge_increments=[0.01, 0.02, 0.03] * unit.elementary_charge
                                                       )
        self.propane_from_smiles_w_vsites = Molecule(molecule)

    def test_empty(self):
        """Test creation of empty topology"""
        topology = Topology()
        assert topology.n_reference_molecules == 0
        assert topology.n_molecules == 0
        assert topology.n_atoms == 0
        assert topology.n_bonds == 0
        assert topology.n_particles == 0
        assert topology.n_virtual_sites == 0
        assert topology.box_vectors is None
        assert len(topology.constrained_atom_pairs.items()) == 0
        assert topology.is_periodic == False

    def test_box_vectors(self):
        """Test the getter and setter for box_vectors"""
        topology = Topology()
        good_box_vectors = unit.Quantity(np.array([10,20,30]), unit.angstrom)
        bad_box_vectors = np.array([10,20,30]) # They're bad because they're unitless
        assert topology.box_vectors is None

        with self.assertRaises(Exception) as context:
            topology.box_vectors = bad_box_vectors
        assert topology.box_vectors is None

        topology.box_vectors = good_box_vectors
        assert (topology.box_vectors == good_box_vectors).all()


    def test_from_smiles(self):
        """Test creation of a openforcefield Topology object from a SMILES string"""
        topology = Topology.from_molecules(self.ethane_from_smiles)

        assert topology.n_reference_molecules == 1
        assert topology.n_molecules == 1
        assert topology.n_atoms == 8
        assert topology.n_bonds == 7
        assert topology.n_particles == 8
        assert topology.n_virtual_sites == 0
        assert topology.box_vectors is None
        assert len(topology.constrained_atom_pairs.items()) == 0
        assert topology.is_periodic == False

        topology.add_molecule(self.ethane_from_smiles)

        assert topology.n_reference_molecules == 1
        assert topology.n_molecules == 2
        assert topology.n_atoms == 16
        assert topology.n_bonds == 14
        assert topology.n_particles == 16
        assert topology.n_virtual_sites == 0
        assert topology.box_vectors is None
        assert len(topology.constrained_atom_pairs.items()) == 0
        assert topology.is_periodic == False

    def test_from_smiles_unique_mols(self):
        """Test the addition of two different molecules to a topology"""
        topology = Topology.from_molecules([self.ethane_from_smiles, self.propane_from_smiles])
        assert topology.n_molecules == 2
        assert topology.n_reference_molecules == 2


    def test_n_topology_atoms(self):
        """Test n_atoms function"""
        topology = Topology()
        assert topology.n_atoms == 0
        assert topology.n_bonds == 0
        topology.add_molecule(self.ethane_from_smiles)
        assert topology.n_atoms == 8
        assert topology.n_bonds == 7

    def test_get_atom(self):
            """
            """
            topology = Topology()
            topology.add_molecule(self.ethane_from_smiles)
            with self.assertRaises(Exception) as context:
                topology_atom = topology.atom(-1)

            # Make sure we get 2 carbons and 8 hydrogens
            n_carbons = 0
            n_hydrogens = 0
            for index in range(8):
                if topology.atom(index).atomic_number == 6:
                    n_carbons += 1
                if topology.atom(index).atomic_number == 1:
                    n_hydrogens += 1
            assert n_carbons == 2
            assert n_hydrogens == 6

            with self.assertRaises(Exception) as context:
                topology_atom = topology.atom(8)

    def test_get_bond(self):
            """
            """
            topology = Topology()
            topology.add_molecule(self.ethane_from_smiles)
            topology.add_molecule(self.ethene_from_smiles)
            with self.assertRaises(Exception) as context:
                topology_atom = topology.bond(-1)

            n_single_bonds = 0
            n_double_bonds = 0
            n_ch_bonds = 0
            n_cc_bonds = 0
            for index in range(12): # 7 from ethane, 5 from ethene
                topology_bond = topology.bond(index)
                if topology_bond.bond_order == 1:
                    n_single_bonds += 1
                if topology_bond.bond_order == 2:
                    n_double_bonds += 1
                n_bond_carbons = 0
                n_bond_hydrogens = 0
                for atom in topology_bond.atoms:
                    if atom.atomic_number == 6:
                        n_bond_carbons += 1
                    if atom.atomic_number == 1:
                        n_bond_hydrogens += 1
                if n_bond_carbons == 2:
                    n_cc_bonds += 1
                if n_bond_carbons == 1 and n_bond_hydrogens == 1:
                    n_ch_bonds += 1

            assert n_single_bonds == 11
            assert n_double_bonds == 1
            assert n_cc_bonds == 2
            assert n_ch_bonds == 10

            with self.assertRaises(Exception) as context:
                topology_bond = topology.bond(12)


    def test_get_virtual_site(self):
        """
        """
        topology = Topology()
        topology.add_molecule(self.ethane_from_smiles_w_vsites)
        assert topology.n_virtual_sites == 2
        topology.add_molecule(self.propane_from_smiles_w_vsites)
        assert topology.n_virtual_sites == 4
        with self.assertRaises(Exception) as context:
            topology_vsite = topology.virtual_site(-1)
        with self.assertRaises(Exception) as context:
            topology_vsite = topology.virtual_site(4)
        topology_vsite1 = topology.virtual_site(0)
        topology_vsite2 = topology.virtual_site(1)
        topology_vsite3 = topology.virtual_site(2)
        topology_vsite4 = topology.virtual_site(3)
        assert topology_vsite1.type == "BondChargeVirtualSite"
        assert topology_vsite2.type == "MonovalentLonePairVirtualSite"
        assert topology_vsite3.type == "BondChargeVirtualSite"
        assert topology_vsite4.type == "MonovalentLonePairVirtualSite"

        n_equal_atoms = 0
        for topology_atom in topology.atoms:
            for vsite in topology.virtual_sites:
                for vsite_atom in vsite.atoms:
                    if topology_atom == vsite_atom:
                        n_equal_atoms += 1

        # There are four virtual sites -- Two BondCharges with 2 atoms, and two MonovalentLonePairs with 3 atoms
        assert n_equal_atoms == 10





    # test_create_from_file_without_partial_charges
    # test_create_from_file_with_partial_charges
    # test_get_fractional_bond_order
    # test_two_of_same_molecule
    # test_two_different_molecules
    # test_redundant_molecule
    # test_to_from_dict
    # test_get_atom
    # test_get_bond
    # test_get_virtual_site
    # test_get_molecule
    # test_get_topology_atom
    # test_get_topology_bond
    # test_get_topology_virtual_site
    # test_get_topology_molecule


    #def test_from_openmm(self):
    #    """Test creation of an openforcefield Topology object from an OpenMM Topology and component molecules"""
    #    pdbfile = app.PDBFile(get_packmol_pdbfile('cyclohexane_ethanol_0.4_0.6.pdb'))
    #    molecules = [ Molecule.from_file(get_monomer_mol2file(name)) for name in ('ethanol', 'cyclohexane') ]
    #    topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)
