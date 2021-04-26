#!/usr/bin/env python

# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Tests for Topology

"""

# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

from unittest import TestCase

import numpy as np
import pytest
from simtk import unit
from simtk.openmm.app import element

from openff.toolkit.tests.test_forcefield import (
    create_cyclohexane,
    create_ethanol,
    create_reversed_ethanol,
)
from openff.toolkit.tests.utils import (
    get_data_file_path,
    requires_openeye,
    requires_pkg,
    requires_rdkit,
)
from openff.toolkit.topology import (
    DuplicateUniqueMoleculeError,
    ImproperDict,
    InvalidBoxVectorsError,
    InvalidPeriodicityError,
    MissingUniqueMoleculesError,
    Molecule,
    Topology,
    ValenceDict,
)
from openff.toolkit.utils import (
    BASIC_CHEMINFORMATICS_TOOLKITS,
    OPENEYE_AVAILABLE,
    RDKIT_AVAILABLE,
    OpenEyeToolkitWrapper,
    RDKitToolkitWrapper,
)

# =============================================================================================
# UTILITY FUNCTIONS
# =============================================================================================


def assert_tuple_of_atoms_equal(
    atom_tuples1, atom_tuples2, transformed_dict_cls=ValenceDict
):
    """Check that two lists of atoms are the same.

    The function compares that the parent molecules are isomorphic and
    that the molecule index is the same.
    """
    assert len(atom_tuples1) == len(atom_tuples2)

    # They are atoms of isomorphic molecules. We assume here that all
    # atoms in the same list of tuples belong to the same molecule so
    # that we can perform the check only once.
    molecule1 = atom_tuples1[0][0]._molecule
    molecule2 = atom_tuples2[0][0]._molecule
    assert molecule1 == molecule2
    for atom_tuple in atom_tuples1:
        for a in atom_tuple:
            assert a._molecule is molecule1
    for atom_tuple in atom_tuples2:
        for a in atom_tuple:
            assert a._molecule is molecule2

    # All atoms are equal. Use ValenceDict for this
    atom_indices = []
    for atom_tuples in [atom_tuples1, atom_tuples2]:
        valence_dict = transformed_dict_cls()
        for atom_tuple in atom_tuples:
            key = tuple(a.molecule_atom_index for a in atom_tuple)
            valence_dict[key] = atom_tuple
        atom_indices.append(valence_dict)
    assert set(atom_indices[0]) == set(atom_indices[1])


# =============================================================================================
# TESTS
# =============================================================================================

# IF we've done our jobs right, it shouldn't matter which toolkit the tests for Topology run using (both's behaviors
# should be indistinguishable)
def test_cheminformatics_toolkit_is_installed():
    """Ensure that at least one supported cheminformatics toolkit is installed."""
    if not (RDKIT_AVAILABLE) and not (OPENEYE_AVAILABLE):
        msg = "No supported cheminformatics toolkits are installed. Please install a supported toolkit:\n"
        msg += str(BASIC_CHEMINFORMATICS_TOOLKITS)
        raise Exception(msg)


class TestTopology(TestCase):
    def setUp(self):
        self.empty_molecule = Molecule()
        self.ethane_from_smiles = Molecule.from_smiles("CC")
        self.ethene_from_smiles = Molecule.from_smiles("C=C")
        self.propane_from_smiles = Molecule.from_smiles("CCC")

        filename = get_data_file_path("molecules/toluene.sdf")
        self.toluene_from_sdf = Molecule.from_file(filename)
        if OpenEyeToolkitWrapper.is_available():
            filename = get_data_file_path("molecules/toluene_charged.mol2")
            # TODO: This will require openeye to load
            self.toluene_from_charged_mol2 = Molecule.from_file(filename)
        self.charged_methylamine_from_smiles = Molecule.from_smiles(
            "[H]C([H])([H])[N+]([H])([H])[H]"
        )

        molecule = Molecule.from_smiles("CC")
        carbons = [atom for atom in molecule.atoms if atom.atomic_number == 6]
        c0_hydrogens = [
            atom for atom in carbons[0].bonded_atoms if atom.atomic_number == 1
        ]
        molecule.add_bond_charge_virtual_site(
            (carbons[0], carbons[1]),
            0.1 * unit.angstrom,
            charge_increments=[0.1, 0.05] * unit.elementary_charge,
        )
        molecule.add_monovalent_lone_pair_virtual_site(
            (c0_hydrogens[0], carbons[0], carbons[1]),
            0.2 * unit.angstrom,
            20 * unit.degree,
            25 * unit.degree,
            charge_increments=[0.01, 0.02, 0.03] * unit.elementary_charge,
        )
        self.ethane_from_smiles_w_vsites = Molecule(molecule)

        # Make a propane with virtual sites
        molecule = Molecule.from_smiles("CCC")
        carbons = [atom for atom in molecule.atoms if atom.atomic_number == 6]
        c0_hydrogens = [
            atom for atom in carbons[0].bonded_atoms if atom.atomic_number == 1
        ]
        # This will add *two* particles (symmetric=True), *one* virtual site
        molecule.add_monovalent_lone_pair_virtual_site(
            (c0_hydrogens[0], carbons[0], carbons[1]),
            0.2 * unit.angstrom,
            20 * unit.degree,
            25 * unit.degree,
            charge_increments=[0.01, 0.02, 0.03] * unit.elementary_charge,
            symmetric=True,
        )
        # This will add *one* particle (symmetric=False), *one* virtual site
        molecule.add_bond_charge_virtual_site(
            (carbons[0], carbons[1]),
            0.1 * unit.angstrom,
            charge_increments=[0.1, 0.05] * unit.elementary_charge,
            symmetric=False,
        )
        self.propane_from_smiles_w_vsites = Molecule(molecule)

        # Make a TIP5 water
        molecule = Molecule.from_smiles("[H][O][H]")
        O1 = [atom for atom in molecule.atoms if atom.atomic_number == 8][0]
        H1, H2 = [atom for atom in O1.bonded_atoms if atom.atomic_number == 1]
        molecule.add_divalent_lone_pair_virtual_site(
            (H1, O1, H2),
            0.7 * unit.angstrom,
            54.71384225 * unit.degree,
            charge_increments=[0.1205, 0.00, 0.1205] * unit.elementary_charge,
            symmetric=True,
        )
        self.tip5_water = Molecule(molecule)

    def test_empty(self):
        """Test creation of empty topology"""
        topology = Topology()
        assert topology.n_reference_molecules == 0
        assert topology.n_topology_molecules == 0
        assert topology.n_topology_atoms == 0
        assert topology.n_topology_bonds == 0
        assert topology.n_topology_particles == 0
        assert topology.n_topology_virtual_sites == 0
        assert topology.box_vectors is None
        assert not topology.is_periodic
        assert len(topology.constrained_atom_pairs.items()) == 0

    def test_box_vectors(self):
        """Test the getter and setter for box_vectors"""
        topology = Topology()
        good_box_vectors = unit.Quantity(np.eye(3) * 20 * unit.angstrom)
        one_dim_vectors = unit.Quantity(np.ones(3) * 20 * unit.angstrom)
        list_vectors = [20, 20, 20] * unit.angstrom
        bad_shape_vectors = unit.Quantity(np.ones(2) * 20 * unit.angstrom)
        bad_units_vectors = unit.Quantity(np.ones(3) * 20 * unit.year)
        unitless_vectors = np.array([10, 20, 30])
        assert topology.box_vectors is None

        for bad_vectors in [
            bad_shape_vectors,
            bad_units_vectors,
            unitless_vectors,
        ]:
            with self.assertRaises(InvalidBoxVectorsError):
                topology.box_vectors = bad_vectors
            assert topology.box_vectors is None

        for good_vectors in [good_box_vectors, one_dim_vectors, list_vectors]:
            topology.box_vectors = good_vectors
            assert (topology.box_vectors == good_vectors * np.eye(3)).all()

    def test_is_periodic(self):
        """Test the getter and setter for is_periodic"""
        vacuum_top = Topology()
        assert vacuum_top.is_periodic is False

        with pytest.raises(InvalidPeriodicityError):
            vacuum_top.is_periodic = True

        solvent_box = Topology()
        solvent_box.box_vectors = np.eye(3) * 4 * unit.nanometer
        assert solvent_box.is_periodic is True

        with pytest.raises(InvalidPeriodicityError):
            solvent_box.is_periodic = False
        solvent_box.box_vectors = None
        assert solvent_box.is_periodic is False

    def test_from_smiles(self):
        """Test creation of a OpenFF Topology object from a SMILES string"""
        topology = Topology.from_molecules(self.ethane_from_smiles)

        assert topology.n_reference_molecules == 1
        assert topology.n_topology_molecules == 1
        assert topology.n_topology_atoms == 8
        assert topology.n_topology_bonds == 7
        assert topology.n_topology_particles == 8
        assert topology.n_topology_virtual_sites == 0
        assert topology.box_vectors is None
        assert len(topology.constrained_atom_pairs.items()) == 0

        topology.add_molecule(self.ethane_from_smiles)

        assert topology.n_reference_molecules == 1
        assert topology.n_topology_molecules == 2
        assert topology.n_topology_atoms == 16
        assert topology.n_topology_bonds == 14
        assert topology.n_topology_particles == 16
        assert topology.n_topology_virtual_sites == 0
        assert topology.box_vectors is None
        assert len(topology.constrained_atom_pairs.items()) == 0

    def test_from_smiles_unique_mols(self):
        """Test the addition of two different molecules to a topology"""
        topology = Topology.from_molecules(
            [self.ethane_from_smiles, self.propane_from_smiles]
        )
        assert topology.n_topology_molecules == 2
        assert topology.n_reference_molecules == 2

    def test_n_topology_atoms(self):
        """Test n_atoms function"""
        topology = Topology()
        assert topology.n_topology_atoms == 0
        assert topology.n_topology_bonds == 0
        topology.add_molecule(self.ethane_from_smiles)
        assert topology.n_topology_atoms == 8
        assert topology.n_topology_bonds == 7

    def test_n_topology_atoms_with_vsites(self):
        """Test n_atoms function when vsites are present"""
        topology = Topology()
        assert topology.n_topology_atoms == 0
        assert topology.n_topology_bonds == 0
        topology.add_molecule(self.ethane_from_smiles_w_vsites)
        assert topology.n_topology_atoms == 8
        assert topology.n_topology_bonds == 7

        topology = Topology()
        assert topology.n_topology_atoms == 0
        assert topology.n_topology_bonds == 0
        topology.add_molecule(self.tip5_water)
        assert topology.n_topology_atoms == 3
        assert topology.n_topology_bonds == 2

    def test_n_topology_virtual_sites(self):
        """Test n_atoms function"""
        topology = Topology()
        assert topology.n_topology_virtual_sites == 0
        topology.add_molecule(self.ethane_from_smiles_w_vsites)
        assert topology.n_topology_virtual_sites == 2

        topology = Topology()
        assert topology.n_topology_virtual_sites == 0
        topology.add_molecule(self.tip5_water)
        assert topology.n_topology_virtual_sites == 1
        assert topology.n_topology_particles == 5

    def test_get_atom(self):
        """Test Topology.atom function (atom lookup from index)"""
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

    def test_topology_atom_element(self):
        """Test getters of TopologyAtom element and atomic number"""
        topology = Topology()
        topology.add_molecule(self.toluene_from_sdf)

        first_element = topology.atom(0).element
        eighth_element = topology.atom(7).element

        # Check if types/instances are correct
        assert isinstance(first_element, element.Element)
        assert isinstance(eighth_element, element.Element)

        # Make sure first is a carbon element and eighth is a hydrogen element
        assert first_element == element.carbon
        assert eighth_element == element.hydrogen

    def test_get_bond(self):
        """Test Topology.bond function (bond lookup from index)"""
        topology = Topology()
        topology.add_molecule(self.ethane_from_smiles)
        topology.add_molecule(self.ethene_from_smiles)
        with self.assertRaises(Exception) as context:
            topology_atom = topology.bond(-1)

        n_single_bonds = 0
        n_double_bonds = 0
        n_ch_bonds = 0
        n_cc_bonds = 0
        for index in range(12):  # 7 from ethane, 5 from ethene
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
        """Test Topology.virtual_site function (get virtual site from index)"""
        topology = Topology()
        topology.add_molecule(self.ethane_from_smiles_w_vsites)
        assert topology.n_topology_virtual_sites == 2
        topology.add_molecule(self.propane_from_smiles_w_vsites)
        assert topology.n_topology_virtual_sites == 4
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
        assert topology_vsite3.type == "MonovalentLonePairVirtualSite"
        assert topology_vsite4.type == "BondChargeVirtualSite"

        n_equal_atoms = 0
        for topology_atom in topology.topology_atoms:
            for vsite in topology.topology_virtual_sites:
                for vsite_atom in vsite.atoms:
                    if topology_atom == vsite_atom:
                        n_equal_atoms += 1

        # There are four virtual sites -- Two BondCharges with 2 atoms, and two MonovalentLonePairs with 3 atoms
        assert n_equal_atoms == 10

    def test_topology_particles_virtualsites_indexed_last(self):
        """
        Test to ensure that virtualsites are strictly indexed after all atoms
        in topology.particles
        """
        from openff.toolkit.topology import TopologyAtom, TopologyVirtualParticle

        topology = Topology()
        topology.add_molecule(self.ethane_from_smiles_w_vsites)
        topology.add_molecule(self.propane_from_smiles_w_vsites)

        # Iterate through all TopologyParticles, ensuring that all atoms appear
        # before all virtualsides
        reading_atoms = True
        for particle in topology.topology_particles:
            if reading_atoms:
                if isinstance(particle, TopologyAtom):
                    pass
                else:
                    reading_atoms = False
            elif not (reading_atoms):
                assert isinstance(particle, TopologyVirtualParticle)

    def test_topology_virtual_site_particle_start_index(self):

        topology = Topology()
        topology.add_molecule(self.propane_from_smiles_w_vsites)
        assert topology.virtual_site(0).topology_virtual_particle_start_index == 11
        assert topology.virtual_site(1).topology_virtual_particle_start_index == 13

    def test_topology_virtual_site_n_particles(self):
        """
        Test if the virtual sites report the correct number of particles
        """
        topology = Topology()
        topology.add_molecule(self.propane_from_smiles_w_vsites)
        assert topology.virtual_site(0).n_particles == 2
        assert topology.virtual_site(1).n_particles == 1

        topology = Topology()
        topology.add_molecule(self.tip5_water)
        assert topology.virtual_site(0).n_particles == 2

    def test_topology_virtualsites_atom_indexing(self):
        """
        Add multiple instances of the same molecule, but in a different
        order, and ensure that virtual site particles are indexed correctly
        """
        topology = Topology()

        topology.add_molecule(create_ethanol())
        topology.add_molecule(create_ethanol())
        topology.add_molecule(create_reversed_ethanol())

        # Add a virtualsite to the reference ethanol
        for ref_mol in topology.reference_molecules:
            atoms = [ref_mol.atoms[i] for i in [0, 1]]
            ref_mol._add_bond_charge_virtual_site(
                atoms,
                0.5 * unit.angstrom,
            )

        virtual_site_topology_atom_indices = [(0, 1), (9, 10), (26, 25)]
        for top_vs, expected_indices in zip(
            topology.topology_virtual_sites, virtual_site_topology_atom_indices
        ):
            assert (
                tuple([at.topology_particle_index for at in top_vs.atoms])
                == expected_indices
            )
            assert top_vs.atom(0).topology_particle_index == expected_indices[0]
            assert top_vs.atom(1).topology_particle_index == expected_indices[1]

    def test_is_bonded(self):
        """Test Topology.virtual_site function (get virtual site from index)"""
        topology = Topology()
        topology.add_molecule(self.propane_from_smiles_w_vsites)
        topology.assert_bonded(0, 1)
        topology.assert_bonded(1, 0)
        topology.assert_bonded(1, 2)
        # C-H bond
        topology.assert_bonded(0, 4)
        with self.assertRaises(Exception) as context:
            topology.assert_bonded(0, 2)

    def test_angles(self):
        """Topology.angles should return image angles of all topology molecules."""
        molecule1 = self.ethane_from_smiles
        molecule2 = self.propane_from_smiles

        # Create topology.
        topology = Topology()
        topology.add_molecule(molecule1)
        topology.add_molecule(molecule1)
        topology.add_molecule(molecule2)

        # The topology should have the correct number of angles.
        topology_angles = list(topology.angles)
        assert len(topology_angles) == topology.n_angles
        assert topology.n_angles == 2 * molecule1.n_angles + molecule2.n_angles

        # Check that the topology angles are the correct ones.
        mol_angle_atoms1 = list(molecule1.angles)
        mol_angle_atoms2 = list(molecule2.angles)
        top_angle_atoms1 = [
            tuple(a._atom for a in atoms)
            for atoms in topology_angles[: molecule1.n_angles]
        ]
        top_angle_atoms2 = [
            tuple(a._atom for a in atoms)
            for atoms in topology_angles[molecule1.n_angles : 2 * molecule1.n_angles]
        ]
        top_angle_atoms3 = [
            tuple(a._atom for a in atoms)
            for atoms in topology_angles[2 * molecule1.n_angles :]
        ]

        assert_tuple_of_atoms_equal(top_angle_atoms1, mol_angle_atoms1)
        assert_tuple_of_atoms_equal(top_angle_atoms2, mol_angle_atoms1)
        assert_tuple_of_atoms_equal(top_angle_atoms3, mol_angle_atoms2)

    def test_propers(self):
        """Topology.propers should return image propers torsions of all topology molecules."""
        molecule1 = self.ethane_from_smiles
        molecule2 = self.propane_from_smiles

        # Create topology.
        topology = Topology()
        topology.add_molecule(molecule1)
        topology.add_molecule(molecule1)
        topology.add_molecule(molecule2)

        # The topology should have the correct number of propers.
        topology_propers = list(topology.propers)
        assert len(topology_propers) == topology.n_propers
        assert topology.n_propers == 2 * molecule1.n_propers + molecule2.n_propers

        # Check that the topology propers are the correct ones.
        mol_proper_atoms1 = list(molecule1.propers)
        mol_proper_atoms2 = list(molecule2.propers)
        top_proper_atoms1 = [
            tuple(a._atom for a in atoms)
            for atoms in topology_propers[: molecule1.n_propers]
        ]
        top_proper_atoms2 = [
            tuple(a._atom for a in atoms)
            for atoms in topology_propers[molecule1.n_propers : 2 * molecule1.n_propers]
        ]
        top_proper_atoms3 = [
            tuple(a._atom for a in atoms)
            for atoms in topology_propers[2 * molecule1.n_propers :]
        ]

        assert_tuple_of_atoms_equal(top_proper_atoms1, mol_proper_atoms1)
        assert_tuple_of_atoms_equal(top_proper_atoms2, mol_proper_atoms1)
        assert_tuple_of_atoms_equal(top_proper_atoms3, mol_proper_atoms2)

    def test_impropers(self):
        """Topology.impropers should return image impropers torsions of all topology molecules."""
        molecule1 = self.ethane_from_smiles
        molecule2 = self.propane_from_smiles

        # Create topology.
        topology = Topology()
        topology.add_molecule(molecule1)
        topology.add_molecule(molecule1)
        topology.add_molecule(molecule2)

        # The topology should have the correct number of impropers.
        topology_impropers = list(topology.impropers)
        assert len(topology_impropers) == topology.n_impropers
        assert topology.n_impropers == 2 * molecule1.n_impropers + molecule2.n_impropers

        # Check that the topology impropers are the correct ones.
        mol_improper_atoms1 = list(molecule1.impropers)
        mol_improper_atoms2 = list(molecule2.impropers)
        top_improper_atoms1 = [
            tuple(a._atom for a in atoms)
            for atoms in topology_impropers[: molecule1.n_impropers]
        ]
        top_improper_atoms2 = [
            tuple(a._atom for a in atoms)
            for atoms in topology_impropers[
                molecule1.n_impropers : 2 * molecule1.n_impropers
            ]
        ]
        top_improper_atoms3 = [
            tuple(a._atom for a in atoms)
            for atoms in topology_impropers[2 * molecule1.n_impropers :]
        ]

        assert_tuple_of_atoms_equal(
            top_improper_atoms1, mol_improper_atoms1, transformed_dict_cls=ImproperDict
        )
        assert_tuple_of_atoms_equal(
            top_improper_atoms2, mol_improper_atoms1, transformed_dict_cls=ImproperDict
        )
        assert_tuple_of_atoms_equal(
            top_improper_atoms3, mol_improper_atoms2, transformed_dict_cls=ImproperDict
        )

    def test_pruned_impropers(self):
        """Test {smirnoff|amber}_impropers from the Topology API"""
        top = Topology.from_molecules(
            [Molecule.from_smiles(smi) for smi in ["N", "C=C"]]
        )

        assert len([*top.smirnoff_impropers]) == 18
        assert len([*top.amber_impropers]) == 18

        # Order not guaranteed, so cannot zip and compare directly
        for smirnoff_imp in top.smirnoff_impropers:
            # Convert SMIRNOFF-style improper into AMBER-style
            mod_imp = (
                smirnoff_imp[1],
                smirnoff_imp[0],
                smirnoff_imp[2],
                smirnoff_imp[3],
            )
            assert mod_imp in top.amber_impropers

    # test_get_fractional_bond_order
    # test_two_of_same_molecule
    # test_two_different_molecules
    # test_to_from_dict
    # test_get_molecule
    # test_get_topology_atom
    # test_get_topology_bond
    # test_get_topology_virtual_site
    # test_get_topology_molecule
    # test_is_bonded
    # TODO: Test serialization

    def test_from_openmm(self):
        """Test creation of an OpenFF Topology object from an OpenMM Topology and component molecules"""
        from simtk.openmm import app

        pdbfile = app.PDBFile(
            get_data_file_path("systems/packmol_boxes/cyclohexane_ethanol_0.4_0.6.pdb")
        )

        with pytest.raises(
            MissingUniqueMoleculesError, match="requires a list of Molecule objects"
        ):
            Topology.from_openmm(pdbfile.topology)

        molecules = [create_ethanol(), create_cyclohexane()]

        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)
        assert topology.n_reference_molecules == 2
        assert topology.n_topology_molecules == 239

    def test_from_openmm_missing_reference(self):
        """Test creation of an OpenFF Topology object from an OpenMM Topology when missing a unique molecule"""
        from simtk.openmm import app

        pdbfile = app.PDBFile(
            get_data_file_path("systems/packmol_boxes/cyclohexane_ethanol_0.4_0.6.pdb")
        )

        molecules = [create_ethanol()]
        with pytest.raises(
            ValueError, match="No match found for molecule C6H12"
        ) as excinfo:
            topology = Topology.from_openmm(
                pdbfile.topology, unique_molecules=molecules
            )

    def test_from_openmm_missing_conect(self):
        """
        Test creation of an OpenFF Topology object from an OpenMM Topology
        when the origin PDB lacks CONECT records
        """
        from simtk.openmm import app

        pdbfile = app.PDBFile(
            get_data_file_path("systems/test_systems/1_ethanol_no_conect.pdb")
        )

        molecules = []
        molecules.append(Molecule.from_smiles("CCO"))
        with pytest.raises(
            ValueError,
            match="No match found for molecule C. This would be a "
            "very unusual molecule to try and parameterize, "
            "and it is likely that the data source it was "
            "read from does not contain connectivity "
            "information. If this molecule is coming from "
            "PDB, please ensure that the file contains CONECT "
            "records.",
        ) as excinfo:
            topology = Topology.from_openmm(
                pdbfile.topology, unique_molecules=molecules
            )

    def test_to_from_openmm(self):
        """Test a round-trip OpenFF -> OpenMM -> OpenFF Topology."""
        from simtk.openmm.app import Aromatic

        # Create OpenFF topology with 1 ethanol and 2 benzenes.
        ethanol = Molecule.from_smiles("CCO")
        benzene = Molecule.from_smiles("c1ccccc1")
        off_topology = Topology.from_molecules(molecules=[ethanol, benzene, benzene])

        # Convert to OpenMM Topology.
        omm_topology = off_topology.to_openmm()

        # Check that bond orders are preserved.
        n_double_bonds = sum([b.order == 2 for b in omm_topology.bonds()])
        n_aromatic_bonds = sum([b.type is Aromatic for b in omm_topology.bonds()])
        assert n_double_bonds == 6
        assert n_aromatic_bonds == 12

        # Check that there is one residue for each molecule.
        assert omm_topology.getNumResidues() == 3
        assert omm_topology.getNumChains() == 3

        # Convert back to OpenFF Topology.
        off_topology_copy = Topology.from_openmm(
            omm_topology, unique_molecules=[ethanol, benzene]
        )

        # The round-trip OpenFF Topology is identical to the original.
        # The reference molecules are the same.
        assert (
            off_topology.n_reference_molecules
            == off_topology_copy.n_reference_molecules
        )
        reference_molecules_copy = list(off_topology_copy.reference_molecules)
        for ref_mol_idx, ref_mol in enumerate(off_topology.reference_molecules):
            assert ref_mol == reference_molecules_copy[ref_mol_idx]

        # The number of topology molecules is the same.
        assert (
            off_topology.n_topology_molecules == off_topology_copy.n_topology_molecules
        )

        # Check atoms.
        assert off_topology.n_topology_atoms == off_topology_copy.n_topology_atoms
        for atom_idx, atom in enumerate(off_topology.topology_atoms):
            atom_copy = off_topology_copy.atom(atom_idx)
            assert atom.atomic_number == atom_copy.atomic_number

        # Check bonds.
        for bond_idx, bond in enumerate(off_topology.topology_bonds):
            bond_copy = off_topology_copy.bond(bond_idx)
            bond_atoms = [a.atomic_number for a in bond.atoms]
            bond_atoms_copy = [a.atomic_number for a in bond_copy.atoms]
            assert bond_atoms == bond_atoms_copy
            assert bond.bond_order == bond_copy.bond_order
            assert bond.bond.is_aromatic == bond_copy.bond.is_aromatic

    @requires_pkg("mdtraj")
    def test_from_mdtraj(self):
        """Test construction of an OpenFF Topology from an MDTraj Topology object"""
        import mdtraj as md

        pdb_path = get_data_file_path(
            "systems/test_systems/1_cyclohexane_1_ethanol.pdb"
        )
        trj = md.load(pdb_path)

        with pytest.raises(
            MissingUniqueMoleculesError, match="requires a list of Molecule objects"
        ):
            Topology.from_mdtraj(trj.top)

        unique_molecules = [
            Molecule.from_smiles(mol_name) for mol_name in ["C1CCCCC1", "CCO"]
        ]
        top = Topology.from_mdtraj(trj.top, unique_molecules=unique_molecules)

        assert top.n_topology_molecules == 2
        assert top.n_topology_bonds == 26

    @requires_rdkit
    def test_to_file_vsites(self):
        """
        Checks that Topology.to_file() doesn't write vsites
        """
        from tempfile import NamedTemporaryFile

        from openff.toolkit.topology import Molecule, Topology

        mol = Molecule.from_pdb_and_smiles(
            get_data_file_path("systems/test_systems/1_ethanol.pdb"), "CCO"
        )
        carbons = [atom for atom in mol.atoms if atom.atomic_number == 6]
        positions = mol.conformers[0]
        mol.add_bond_charge_virtual_site(
            (carbons[0], carbons[1]),
            0.1 * unit.angstrom,
            charge_increments=[0.1, 0.05] * unit.elementary_charge,
        )
        topology = Topology()
        topology.add_molecule(mol)
        count = 0
        # The file should be printed out with 9 atoms and 0 virtualsites, so we check to ensure that thtere are only 9 HETATM entries
        with NamedTemporaryFile(suffix=".pdb") as iofile:
            topology.to_file(iofile.name, positions)
            data = open(iofile.name).readlines()
            for line in data:
                if line.startswith("HETATM"):
                    count = count + 1
        assert count == 9

    @requires_rdkit
    def test_to_file_units_check(self):
        """
        Checks whether writing pdb with unitless positions, Angstrom positions,
        nanometer positions, result in the same output
        """
        from tempfile import NamedTemporaryFile

        from simtk.unit import nanometer

        from openff.toolkit.topology import Molecule, Topology

        topology = Topology()
        mol = Molecule.from_pdb_and_smiles(
            get_data_file_path("systems/test_systems/1_ethanol.pdb"), "CCO"
        )
        topology.add_molecule(mol)
        positions_angstrom = mol.conformers[0]
        count = 1
        # Write the molecule to PDB and ensure that the X coordinate of the first atom is 10.172
        with NamedTemporaryFile(suffix=".pdb") as iofile:
            topology.to_file(iofile.name, positions_angstrom)
            data = open(iofile.name).readlines()
            for line in data:
                if line.startswith("HETATM") and count == 1:
                    count = count + 1
                    coord = line.split()[-6]
        assert coord == "10.172"

        # Do the same check, but feed in equivalent positions measured in nanometers and ensure the PDB is still the same
        count = 1
        coord = None
        with NamedTemporaryFile(suffix=".pdb") as iofile:
            positions_nanometer = positions_angstrom.in_units_of(nanometer)
            topology.to_file(iofile.name, positions_nanometer)
            data = open(iofile.name).readlines()
            for line in data:
                if line.startswith("HETATM") and count == 1:
                    count = count + 1
                    coord = line.split()[-6]
        assert coord == "10.172"

        count = 1
        coord = "abc"
        with NamedTemporaryFile(suffix=".pdb") as iofile:
            positions_unitless = positions_angstrom._value
            topology.to_file(iofile.name, positions_unitless)
            data = open(iofile.name).readlines()
            for line in data:
                if line.startswith("HETATM") and count == 1:
                    count = count + 1
                    coord = line.split()[-6]
        assert coord == "10.172"

    @requires_rdkit
    def test_to_file_fileformat_lettercase(self):
        """
        Checks if fileformat specifier is indpendent of upper/lowercase
        """
        from tempfile import NamedTemporaryFile

        from openff.toolkit.topology import Molecule, Topology

        topology = Topology()
        mol = Molecule.from_pdb_and_smiles(
            get_data_file_path("systems/test_systems/1_ethanol.pdb"), "CCO"
        )
        topology.add_molecule(mol)
        positions = mol.conformers[0]
        count = 1
        with NamedTemporaryFile(suffix=".pdb") as iofile:
            topology.to_file(iofile.name, positions, file_format="pDb")
            data = open(iofile.name).readlines()
            for line in data:
                if line.startswith("HETATM") and count == 1:
                    count = count + 1
                    coord = line.split()[-6]
        assert coord == "10.172"

    @requires_rdkit
    def test_to_file_fileformat_invalid(self):
        """
        Checks for invalid file format
        """
        from openff.toolkit.topology import Molecule, Topology

        topology = Topology()
        mol = Molecule.from_pdb_and_smiles(
            get_data_file_path("systems/test_systems/1_ethanol.pdb"), "CCO"
        )
        topology.add_molecule(mol)
        positions = mol.conformers[0]
        fname = "ethanol_file.pdb"
        with pytest.raises(NotImplementedError):
            topology.to_file(fname, positions, file_format="AbC")

    def test_to_file_no_molecules(self):
        """
        Checks if Topology.to_file() writes a file with no topology and no coordinates
        """
        from tempfile import NamedTemporaryFile

        from openff.toolkit.topology import Topology

        topology = Topology()
        lines = []
        with NamedTemporaryFile(suffix=".pdb") as iofile:
            topology.to_file(iofile.name, [])
            data = open(iofile.name).readlines()
            for line in data:
                lines.append(line.split())
        assert lines[1] == ["END"]

    @requires_rdkit
    def test_to_file_multi_molecule_different_order(self):
        """
        Checks for the following if Topology.to_write maintains the order of atoms
         for the same molecule with different indexing
        """
        from tempfile import NamedTemporaryFile

        from openff.toolkit.tests.test_forcefield import (
            create_ethanol,
            create_reversed_ethanol,
        )
        from openff.toolkit.topology import Molecule, Topology

        topology = Topology()
        topology.add_molecule(create_ethanol())
        topology.add_molecule(create_reversed_ethanol())
        mol = Molecule.from_pdb_and_smiles(
            get_data_file_path("systems/test_systems/1_ethanol.pdb"), "CCO"
        )
        positions = mol.conformers[0]
        # Make up coordinates for the second ethanol by translating the first by 10 angstroms
        # (note that this will still be a gibberish conformation, since the atom order in the second molecule is different)
        positions = np.concatenate([positions, positions + 10.0 * unit.angstrom])
        element_order = []

        with NamedTemporaryFile(suffix=".pdb") as iofile:
            topology.to_file(iofile.name, positions)
            data = open(iofile.name).readlines()
            for line in data:
                if line.startswith("HETATM"):
                    element_order.append(line.strip()[-1])
        assert element_order == [
            "C",
            "C",
            "O",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
            "O",
            "C",
            "C",
        ]

    @requires_openeye
    def test_from_openmm_duplicate_unique_mol(self):
        """Check that a DuplicateUniqueMoleculeError is raised if we try to pass in two indistinguishably unique mols"""
        from simtk.openmm import app

        pdbfile = app.PDBFile(
            get_data_file_path("systems/packmol_boxes/cyclohexane_ethanol_0.4_0.6.pdb")
        )
        molecules = [
            Molecule.from_file(get_data_file_path(name))
            for name in (
                "molecules/ethanol.mol2",
                "molecules/ethanol_reordered.mol2",
                "molecules/cyclohexane.mol2",
            )
        ]
        with self.assertRaises(DuplicateUniqueMoleculeError) as context:
            topology = Topology.from_openmm(
                pdbfile.topology, unique_molecules=molecules
            )

    @pytest.mark.skip
    def test_from_openmm_distinguish_using_stereochemistry(self):
        """Test creation of an OpenFF Topology object from an OpenMM topology with stereoisomers"""
        # From Jeff: The graph representation created from OMM molecules during the matching
        # process doesn't encode stereochemistry.
        raise NotImplementedError

    def test_to_openmm_assign_unique_atom_names(self):
        """
        Ensure that OFF topologies with no pre-existing atom names have unique
        atom names applied when being converted to openmm
        """
        # Create OpenFF topology with 1 ethanol and 2 benzenes.
        ethanol = Molecule.from_smiles("CCO")
        benzene = Molecule.from_smiles("c1ccccc1")
        off_topology = Topology.from_molecules(molecules=[ethanol, benzene, benzene])
        omm_topology = off_topology.to_openmm()
        atom_names = set()
        for atom in omm_topology.atoms():
            atom_names.add(atom.name)
        # There should be 6 unique Cs, 6 unique Hs, and 1 unique O, for a total of 13 unique atom names
        assert len(atom_names) == 13

    def test_to_openmm_assign_some_unique_atom_names(self):
        """
        Ensure that OFF topologies with some pre-existing atom names have unique
        atom names applied to the other atoms when being converted to openmm
        """
        # Create OpenFF topology with 1 ethanol and 2 benzenes.
        ethanol = Molecule.from_smiles("CCO")
        for atom in ethanol.atoms:
            atom.name = f"AT{atom.molecule_atom_index}"
        benzene = Molecule.from_smiles("c1ccccc1")
        off_topology = Topology.from_molecules(molecules=[ethanol, benzene, benzene])
        omm_topology = off_topology.to_openmm()
        atom_names = set()
        for atom in omm_topology.atoms():
            atom_names.add(atom.name)
        # There should be 9 "ATOM#"-labeled atoms, 6 unique Cs, and 6 unique Hs,
        # for a total of 21 unique atom names
        assert len(atom_names) == 21

    def test_to_openmm_assign_unique_atom_names_some_duplicates(self):
        """
        Ensure that OFF topologies where some molecules have invalid/duplicate
        atom names have unique atom names applied while the other molecules are unaffected.
        """
        # Create OpenFF topology with 1 ethanol and 2 benzenes.
        ethanol = Molecule.from_smiles("CCO")

        # Assign duplicate atom names in ethanol (two AT0s)
        ethanol_atom_names_with_duplicates = [f"AT{i}" for i in range(ethanol.n_atoms)]
        ethanol_atom_names_with_duplicates[1] = "AT0"
        for atom, atom_name in zip(ethanol.atoms, ethanol_atom_names_with_duplicates):
            atom.name = atom_name

        # Assign unique atom names in benzene
        benzene = Molecule.from_smiles("c1ccccc1")
        benzene_atom_names = [f"AT{i}" for i in range(benzene.n_atoms)]
        for atom, atom_name in zip(benzene.atoms, benzene_atom_names):
            atom.name = atom_name

        off_topology = Topology.from_molecules(molecules=[ethanol, benzene, benzene])
        omm_topology = off_topology.to_openmm()
        atom_names = set()
        for atom in omm_topology.atoms():
            atom_names.add(atom.name)

        # There should be  12 "AT#"-labeled atoms (from benzene), 2 unique Cs,
        # 1 unique O, and 6 unique Hs, for a total of 21 unique atom names
        assert len(atom_names) == 21

    def test_to_openmm_do_not_assign_unique_atom_names(self):
        """
        Test disabling unique atom name assignment in Topology.to_openmm
        """
        # Create OpenFF topology with 1 ethanol and 2 benzenes.
        ethanol = Molecule.from_smiles("CCO")
        for atom in ethanol.atoms:
            atom.name = "eth_test"
        benzene = Molecule.from_smiles("c1ccccc1")
        benzene.atoms[0].name = "bzn_test"
        off_topology = Topology.from_molecules(molecules=[ethanol, benzene, benzene])
        omm_topology = off_topology.to_openmm(ensure_unique_atom_names=False)
        atom_names = set()
        for atom in omm_topology.atoms():
            atom_names.add(atom.name)
        # There should be 9 atom named "eth_test", 1 atom named "bzn_test",
        # and 12 atoms named "", for a total of 3 unique atom names
        assert len(atom_names) == 3

    @requires_openeye
    def test_chemical_environments_matches_OE(self):
        """Test Topology.chemical_environment_matches"""
        from simtk.openmm import app

        toolkit_wrapper = OpenEyeToolkitWrapper()
        pdbfile = app.PDBFile(
            get_data_file_path("systems/packmol_boxes/cyclohexane_ethanol_0.4_0.6.pdb")
        )
        # toolkit_wrapper = RDKitToolkitWrapper()
        molecules = [
            Molecule.from_file(get_data_file_path(name))
            for name in ("molecules/ethanol.mol2", "molecules/cyclohexane.mol2")
        ]
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)
        # Test for substructure match
        matches = topology.chemical_environment_matches(
            "[C:1]-[C:2]-[O:3]", toolkit_registry=toolkit_wrapper
        )
        assert len(matches) == 143
        assert matches[0].topology_atom_indices == (1728, 1729, 1730)
        # Test for whole-molecule match
        matches = topology.chemical_environment_matches(
            "[H][C:1]([H])([H])-[C:2]([H])([H])-[O:3][H]",
            toolkit_registry=toolkit_wrapper,
        )
        assert (
            len(matches) == 1716
        )  # 143 * 12 (there are 12 possible hydrogen mappings)
        assert matches[0].topology_atom_indices == (1728, 1729, 1730)
        # Search for a substructure that isn't there
        matches = topology.chemical_environment_matches(
            "[C][C:1]-[C:2]-[O:3]", toolkit_registry=toolkit_wrapper
        )
        assert len(matches) == 0

    @requires_rdkit
    def test_chemical_environments_matches_RDK(self):
        """Test Topology.chemical_environment_matches"""
        from simtk.openmm import app

        toolkit_wrapper = RDKitToolkitWrapper()
        pdbfile = app.PDBFile(
            get_data_file_path("systems/packmol_boxes/cyclohexane_ethanol_0.4_0.6.pdb")
        )
        # toolkit_wrapper = RDKitToolkitWrapper()
        # molecules = [Molecule.from_file(get_data_file_path(name)) for name in ('molecules/ethanol.mol2',
        #                                                                      'molecules/cyclohexane.mol2')]
        molecules = []
        molecules.append(Molecule.from_smiles("CCO"))
        molecules.append(Molecule.from_smiles("C1CCCCC1"))
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)
        # Count CCO matches
        matches = topology.chemical_environment_matches(
            "[C:1]-[C:2]-[O:3]", toolkit_registry=toolkit_wrapper
        )
        assert len(matches) == 143
        assert matches[0].topology_atom_indices == (1728, 1729, 1730)
        matches = topology.chemical_environment_matches(
            "[H][C:1]([H])([H])-[C:2]([H])([H])-[O:3][H]",
            toolkit_registry=toolkit_wrapper,
        )
        assert (
            len(matches) == 1716
        )  # 143 * 12 (there are 12 possible hydrogen mappings)
        assert matches[0].topology_atom_indices == (1728, 1729, 1730)
        # Search for a substructure that isn't there
        matches = topology.chemical_environment_matches(
            "[C][C:1]-[C:2]-[O:3]", toolkit_registry=toolkit_wrapper
        )
        assert len(matches) == 0
