#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Tests for molecular topology representations

At least one supported cheminformatics toolkit must be installed to run these tests.
Only the tests applicable to that toolkit will be run.

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import pickle
from functools import partial
from unittest import TestCase
from numpy.testing import assert_almost_equal
from tempfile import NamedTemporaryFile
from simtk import unit
import pytest
import numpy as np
from openforcefield import utils, topology
from openforcefield.topology.molecule import FrozenMolecule, Molecule, Atom, Bond, BondChargeVirtualSite, MonovalentLonePairVirtualSite, DivalentLonePairVirtualSite, TrivalentLonePairVirtualSite, ALLOWED_CHARGE_MODELS, ALLOWED_FRACTIONAL_BOND_ORDER_MODELS
from openforcefield.utils import get_data_filename
# TODO: Will the ToolkitWrapper allow us to pare that down?
#from openforcefield.utils import RDKIT_UNAVAILABLE, OPENEYE_UNAVAILABLE, SUPPORTED_TOOLKITS, TOOLKIT_PRECEDENCE, SUPPORTED_FILE_FORMATS
from openforcefield.utils.toolkits import ToolkitWrapper, OpenEyeToolkitWrapper, RDKitToolkitWrapper, AmberToolsToolkitWrapper, ToolkitRegistry, BASIC_CHEMINFORMATICS_TOOLKITS


#=============================================================================================
# TESTS
#=============================================================================================

# TODO: Generalize this to instead catch ToolkitWrapper exceptions in case the toolkit capability is unavailable
_OPENEYE_UNAVAILABLE_MESSAGE = 'requires the OpenEye toolkit'
_RDKIT_UNAVAILABLE_MESSAGE = 'requires RDKit'

# TODO: Add tests comparing RDKit and OpenEye aromaticity perception

def assert_molecule_is_equal(molecule1, molecule2, msg):
    """Compare whether two Molecule objects are equal

    Parameters
    ----------
    molecule1, molecule2 : openforcefield.topology.Molecule
        Molecules to be compared
    msg : str
        Message to include if molecules fail to match.

    """
    if not(molecule1.is_isomorphic(molecule2)):
        raise Exception(msg)


# Skipping this test -- The cheminformatics toolkit test is run inside of toolkits.py
@pytest.mark.skip
def test_cheminformatics_toolkit_is_installed():
    """Ensure that at least one supported cheminformatics toolkit is installed."""
    if not(any(tk.toolkit_is_available for tk in toolkits.BASIC_CHEMINFORMATICS_TOOLKITS)):

        msg = 'No supported cheminformatics toolkits are installed. Please install a supported toolkit:\n'
        for tk in basic_cheminf_toolkits:
            msg += '{} : {}\n'.format(tk._toolkit_name, tk._toolkit_installation_instructions)

        raise Exception(msg)

class TestMolecule(TestCase):
    @classmethod
    def setUpClass(cls):
        """
        This will have us just load the test molecules once, as the setUp function runs for each test
        """
        super(TestMolecule, cls).setUpClass()
        #filename = get_data_filename('molecules/zinc-subset-tripos.mol2.gz')
        #filename = get_data_filename('molecules/DrugBank_tripos.mol2')
        filename = get_data_filename('molecules/MiniDrugBank_tripos.mol2')
        molecules = Molecule.from_file(filename, exception_if_undefined_stereo=False)
        molecules = [mol for mol in molecules if not(mol is None)]
        cls.test_molecules = molecules


    def setUp(self):
        # TODO: Serialize the offmols instead so that we can run this test without toolkits
        import copy
        #self.molecules = pickle.load('zinc-subset-offmols.pkl')

        self.molecules = copy.deepcopy(TestMolecule.test_molecules)

    # TODO: Test getstate/setstate
    # TODO: Test {to_from}_{dict|yaml|toml|json|bson|messagepack|pickle}
    def test_pickle(self):
        """Test pickling"""
        serialized = pickle.dumps(self.molecules)
        molecules_copy = pickle.loads(serialized)
        for (molecule, molecule_copy) in zip(self.molecules, molecules_copy):
            assert molecule == molecule_copy

    def test_to_from_dict(self):
        """Test to_dict and from_dict functions"""
        for molecule in self.molecules:
            serialized = molecule.to_dict()
            molecule_copy = Molecule.from_dict(serialized)
            assert molecule == molecule_copy

    def test_create_atom(self):
        """Test Atom creation"""
        # Create a non-aromatic carbon atom
        atom = Atom(6, 0, False)
        # Create a chiral carbon atom
        atom = Atom(6, 0, False, stereochemistry='R', name='CT')

    def test_atom_properties(self):
        """Test Atom.element"""
        from simtk.openmm.app import element
        formal_charge = 0
        is_aromatic = False
        # Attempt to create all elements supported by OpenMM
        elements = [getattr(element, name) for name in dir(element) if (type(getattr(element, name)) == element.Element)]
        # The above runs into a problem with deuterium (fails name assertion)
        elements.remove(element.deuterium)
        for this_element in elements:
            atom = Atom(this_element.atomic_number, formal_charge, is_aromatic, name=this_element.name)
            assert atom.atomic_number == this_element.atomic_number
            assert atom.element == this_element
            assert atom.mass == this_element.mass
            assert atom.formal_charge == formal_charge
            assert atom.is_aromatic == is_aromatic
            assert atom.name == this_element.name

    def test_create_molecule(self):
        """Test creation of molecule by adding molecules and bonds"""
        # Define a methane molecule
        molecule = Molecule()
        molecule.name = 'methane'
        C = molecule.add_atom(6, 0, False)
        H1 = molecule.add_atom(1, 0, False)
        H2 = molecule.add_atom(1, 0, False)
        H3 = molecule.add_atom(1, 0, False)
        H4 = molecule.add_atom(1, 0, False)
        molecule.add_bond(C, H1, 1, False)
        molecule.add_bond(C, H2, 1, False)
        molecule.add_bond(C, H3, 1, False)
        molecule.add_bond(C, H4, 1, False)



    def test_add_conformers(self):
        """Test addition of conformers to a molecule"""
        import numpy as np
        from simtk import unit
        # Define a methane molecule
        molecule = Molecule()
        molecule.name = 'methane'
        C = molecule.add_atom(6, 0, False)
        H1 = molecule.add_atom(1, 0, False)
        H2 = molecule.add_atom(1, 0, False)
        H3 = molecule.add_atom(1, 0, False)
        H4 = molecule.add_atom(1, 0, False)
        molecule.add_bond(C, H1, 1, False)
        molecule.add_bond(C, H2, 1, False)
        molecule.add_bond(C, H3, 1, False)
        molecule.add_bond(C, H4, 1, False)

        assert molecule.n_conformers == 0
        # Add a conformer that should work
        conf1 = unit.Quantity(np.array([[ 1., 2.,3.] ,[4. ,5. ,6.],[7., 8., 9.],
                                        [10.,11.,12.],[13.,14.,15]]),
                              unit.angstrom)
        molecule.add_conformer(conf1)
        assert molecule.n_conformers == 1

        conf2 = unit.Quantity(np.array([[101., 102. ,103.], [104. ,105. ,106.], [107., 108., 109.],
                                        [110.,111.,112.],   [113.,114.,115]]),
                              unit.angstrom)
        molecule.add_conformer(conf2)
        assert molecule.n_conformers == 2

        # Add conformers with too few coordinates
        conf_missing_z = unit.Quantity(np.array([[101., 102. ,103.], [104. ,105. ,106.], [107., 108., 109.],
                                        [110.,111.,112.],   [113.,114.]]),
                                        unit.angstrom)
        with self.assertRaises(Exception) as context:
            molecule.add_conformer(conf_missing_z)

        conf_too_few_atoms = unit.Quantity(np.array([[101., 102. ,103.], [104. ,105. ,106.], [107., 108., 109.],
                                                     [110.,111.,112.]]),
                                                     unit.angstrom)
        with self.assertRaises(Exception) as context:
            molecule.add_conformer(conf_too_few_atoms)


        # Add a conformer with too many coordinates
        conf_too_many_atoms = unit.Quantity(np.array([[101., 102., 103.], [104., 105., 106.], [107., 108., 109.],
                                                      [110., 111., 112.], [113., 114., 115.], [116., 117., 118.]]),
                                            unit.angstrom)
        with self.assertRaises(Exception) as context:
            molecule.add_conformer(conf_too_many_atoms)

        # Add a conformer with no coordinates
        conf_no_coordinates = unit.Quantity(np.array([]),
                                            unit.angstrom)
        with self.assertRaises(Exception) as context:
            molecule.add_conformer(conf_no_coordinates)

        # Add a conforer with units of nanometers
        conf3 = unit.Quantity(np.array([[ 1., 2.,3.] ,[4. ,5. ,6.],[7., 8., 9.],
                                        [10.,11.,12.],[13.,14.,15]]),
                              unit.nanometer)
        molecule.add_conformer(conf3)
        assert molecule.n_conformers == 3
        assert molecule.conformers[2][0][0] == 10. * unit.angstrom

        # Add a conformer with units of nanometers
        conf_nonsense_units = unit.Quantity(np.array([[ 1., 2.,3.] ,[4. ,5. ,6.],[7., 8., 9.],
                                        [10.,11.,12.],[13.,14.,15]]),
                              unit.joule)
        with self.assertRaises(Exception) as context:
            molecule.add_conformer(conf_nonsense_units)

        # Add a conformer with no units
        conf_unitless = np.array([[ 1., 2.,3.] ,[4. ,5. ,6.],[7., 8., 9.],
                                  [10.,11.,12.],[13.,14.,15]])
        with self.assertRaises(Exception) as context:
            molecule.add_conformer(conf_unitless)


    def test_create_empty(self):
        """Test creation of an empty Molecule"""
        molecule = Molecule()

    def test_molecule_name(self):
        """Test creation of an empty Molecule"""
        molecule1 = Molecule()
        molecule1.name = None

        molecule2 = Molecule()
        molecule2.name = ''
        assert molecule1.name == molecule2.name

    def test_create_copy(self):
        """Test creation of a Molecule from another Molecule"""
        for molecule in self.molecules:
            molecule_copy = Molecule(molecule)
            assert molecule_copy == molecule

    @pytest.mark.skipif( not OpenEyeToolkitWrapper.toolkit_is_available(), reason='OpenEye Toolkit not available')
    def test_create_openeye(self):
        """Test creation of a molecule from an OpenEye oemol"""
        known_failures = ['ZINC05964684', 'ZINC05885163', 'ZINC05543156', 'ZINC17211981',
                          'ZINC17312986', 'ZINC06424847', 'ZINC04963126']
        for molecule in self.molecules:
            if molecule.name in known_failures:
                continue
            oemol = molecule.to_openeye()
            molecule_copy = Molecule(oemol)
            assert molecule == molecule_copy


    @pytest.mark.skipif( not RDKitToolkitWrapper.toolkit_is_available(), reason='RDKit Toolkit not available')
    def test_create_rdkit(self):
        """Test creation of a molecule from an RDKit rdmol"""
        # Using ZINC test set
        #known_failures = ['ZINC17060065', 'ZINC16448882', 'ZINC15772239','ZINC11539132',
        #                  'ZINC05975187', 'ZINC17111082', 'ZINC00265517']
        # Using DrugBank test set
        known_failures = ['DrugBank_349', 'DrugBank_1420', 'DrugBank_1671']
        failures = []
        fail_smileses = []
        for molecule in self.molecules:
            if molecule.name in known_failures:
                continue
            rdmol = molecule.to_rdkit()
            molecule_copy = Molecule(rdmol)
            if not(molecule == molecule_copy):
                failures.append(molecule.name)
                fail_smileses.append((molecule.to_smiles(), molecule_copy.to_smiles()))
            #assert molecule == molecule_copy

        print(failures)
        for name, (smi1, smi2) in zip(failures, fail_smileses):
            print(name)
            print(smi1)
            print(smi2)


    #These should be toolkit tests
    @pytest.mark.skip
    def test_create_from_file(self):
        """Test creation of a molecule from a filename or file-like object"""
        # TODO: Expand test to both openeye and rdkit toolkits
        filename = get_data_filename('molecules/toluene.mol2')
        molecule1 = Molecule(filename)
        with open(filename, 'r') as infile:
            molecule2 = Molecule(infile, file_format='MOL2')
        assert molecule1 == molecule2
        import gzip
        with gzip.GzipFile(filename + '.gz', 'r') as infile:
            molecule3 = Molecule(infile, file_format='MOL2')
        assert molecule3 == molecule1

        # Ensure that attempting to initialize a single Molecule from a file containing multiple molecules raises a ValueError
        with pytest.raises(ValueError) as exc_info:
            filename = get_data_filename('molecules/zinc-subset-tripos.mol2.gz')
            molecule = Molecule(filename)

    def test_create_from_serialized(self):
        """Test creation from serialized molecule"""
        for molecule in self.molecules:
            serialized_molecule = molecule.__getstate__()
            molecule_copy = Molecule(serialized_molecule)
            assert molecule == molecule_copy

    # TODO: Takes too long right now -- performance improvements might help
    @pytest.mark.skip
    def test_equality(self):
        """Test equality operator"""
        nmolecules = len(self.molecules)
        for i in range(nmolecules):
            for j in range(i, nmolecules):
                assert (self.molecules[i] == self.molecules[j]) == (i == j)


    def test_to_networkx(self):
        """Test generation of NetworkX graphs"""
        for molecule in self.molecules:
            graph = molecule.to_networkx()

    def test_add_atoms_and_bonds(self):
        """Test the creation of a molecule from the addition of atoms and bonds"""
        for molecule in self.molecules:
            molecule_copy = Molecule()
            for atom in molecule.atoms:
                molecule_copy.add_atom(atom.atomic_number, atom.formal_charge, atom.is_aromatic, stereochemistry=atom.stereochemistry)
            for bond in molecule.bonds:
                molecule_copy.add_bond(bond.atom1_index, bond.atom2_index, bond.bond_order, bond.is_aromatic,
                                       stereochemistry=bond.stereochemistry,
                                       fractional_bond_order=bond.fractional_bond_order)
            # Try to add the final bond twice, which should raise an Exception
            with self.assertRaises(Exception) as context:
                molecule_copy.add_bond(bond.atom1_index, bond.atom2_index, bond.bond_order, bond.is_aromatic,
                                       stereochemistry=bond.stereochemistry,
                                       fractional_bond_order=bond.fractional_bond_order)

            assert molecule == molecule_copy


    def test_add_virtual_site_units(self):
        """
        Tests the unit type checking of the VirtualSite base class
        """
        
        # TODO: Should these be using BondChargeVirtualSite, or should we just call the base class (which does the unit checks) directly?
        
        # Prepare values for unit checks
        distance_unitless = 0.4
        sigma_unitless = 0.1
        rmin_half_unitless = 0.2
        epsilon_unitless = 0.3
        charge_increments_unitless = [0.1, 0.2, 0.3, 0.4]
        distance = distance_unitless * unit.angstrom
        sigma = sigma_unitless * unit.angstrom
        rmin_half = rmin_half_unitless * unit.angstrom
        epsilon = epsilon_unitless * (unit.kilojoule / unit.mole)
        charge_increments = [i * unit.elementary_charge for i in charge_increments_unitless]

        for molecule in self.molecules:
            atom1 = molecule.atoms[0]
            atom2 = molecule.atoms[1]
            atom3 = molecule.atoms[2]
            atom4 = molecule.atoms[3]

            # Try to feed in unitless sigma
            with self.assertRaises(Exception) as context:
                molecule.add_bond_charge_virtual_site([atom1, atom2, atom3], distance, epsilon=epsilon, sigma=sigma_unitless)
            
            # Try to feed in unitless rmin_half
            with self.assertRaises(Exception) as context:
                molecule.add_bond_charge_virtual_site([atom1, atom2, atom3], distance, epsilon=epsilon, rmin_half=rmin_half_unitless)

            # Try to feed in unitless epsilon
            with self.assertRaises(Exception) as context:
                molecule.add_bond_charge_virtual_site([atom1, atom2, atom3], distance, epsilon=epsilon_unitless, sigma=sigma, rmin_half=rmin_half)

            # Try to feed in unitless charges
            with self.assertRaises(Exception) as context:
                molecule.add_bond_charge_virtual_site([atom1, atom2, atom3, atom4], distance, charge_incrtements=charge_increments_unitless)
            

            # We shouldn't be able to give both rmin_half and sigma VdW parameters.
            with self.assertRaises(Exception) as context:
                molecule.add_bond_charge_virtual_site([atom1, atom2, atom3], distance, epsilon=epsilon, sigma=sigma, rmin_half=rmin_half)
                
            # Try creating virtual site from sigma+epsilon
            vsite1_index = molecule.add_bond_charge_virtual_site([atom1, atom2, atom3], distance, epsilon=epsilon, sigma=sigma)
            # Try creating virutal site from rmin_half+epsilon
            vsite2_index = molecule.add_bond_charge_virtual_site([atom1, atom2, atom3], distance, epsilon=epsilon, rmin_half=rmin_half)

            # TODO: Test the @property getters for sigma, epsilon, and rmin_half
            
            # We should have to give as many charge increments as atoms (len(charge_increments) = 4
            with self.assertRaises(Exception) as context:
                molecule.add_bond_charge_virtual_site([atom1, atom2, atom3], distance, charge_increments=charge_increments)
                
            vsite3_index = molecule.add_bond_charge_virtual_site([atom1, atom2, atom3, atom4], distance, charge_increments=charge_increments)
            
            
    def test_add_bond_charge_virtual_site(self):
        """Test the addition of a BondChargeVirtualSite to a molecule.
           Also tests many of the inputs of the parent VirtualSite class
        """
        for molecule in self.molecules:
            atom1 = molecule.atoms[0]
            atom2 = molecule.atoms[1]
            atom3 = molecule.atoms[2]
            atom4 = molecule.atoms[3]
            
            # Prepare values for unit checks
            distance_unitless = 0.4
            distance = distance_unitless * unit.angstrom

            
            # Try to feed in a unitless distance
            with self.assertRaises(AssertionError) as context:
                vsite1_index = molecule.add_bond_charge_virtual_site([atom1, atom2, atom3], distance_unitless)
                

            vsite1_index = molecule.add_bond_charge_virtual_site([atom1, atom2, atom3], distance)
            vsite1 = molecule.virtual_sites[vsite1_index]
            assert atom1 in vsite1.atoms
            assert atom2 in vsite1.atoms
            assert atom3 in vsite1.atoms
            assert vsite1 in atom1.virtual_sites
            assert vsite1 in atom2.virtual_sites
            assert vsite1 in atom3.virtual_sites
            assert vsite1.distance == distance

            vsite2_index = molecule.add_bond_charge_virtual_site([atom1, atom2, atom3],
                                                                 distance,
                                                                 sigma=0.1*unit.angstrom,
                                                                 epsilon=1.0*unit.kilojoule_per_mole,
                                                                 charge_increments=unit.Quantity(np.array([0.1, 0.2, 0.3]), unit.elementary_charge)
                                                                 )
            vsite2 = molecule.virtual_sites[vsite2_index]

            # test serialization
            molecule_dict = molecule.to_dict()
            molecule2 = Molecule.from_dict(molecule_dict)
            assert molecule.to_dict() == molecule2.to_dict()
            assert 0

    # TODO: Make a test for to_dict and from_dict for VirtualSites (even though they're currently just unloaded using
    #      (for example) Molecule._add_bond_virtual_site functions
    def test_add_monovalent_lone_pair_virtual_site(self):
        """Test addition of a MonovalentLonePairVirtualSite to the Molecule"""
        for molecule in self.molecules:
            atom1 = molecule.atoms[0]
            atom2 = molecule.atoms[1]
            atom3 = molecule.atoms[2]
            atom4 = molecule.atoms[3]

            # Prepare values for unit checks
            distance_unitless = 0.3
            out_of_plane_angle_unitless = 30
            in_plane_angle_unitless = 0.2
            distance = distance_unitless * unit.angstrom
            out_of_plane_angle = out_of_plane_angle_unitless * unit.degree
            in_plane_angle = in_plane_angle_unitless * unit.radian

            # Try passing in a unitless distance
            with self.assertRaises(AssertionError) as context:
                vsite1_index = molecule.add_monovalent_lone_pair_virtual_site([atom1, atom2], distance_unitless, out_of_plane_angle, in_plane_angle)

            # Try passing in a unitless out_of_plane_angle
            with self.assertRaises(AssertionError) as context:
                vsite1_index = molecule.add_monovalent_lone_pair_virtual_site([atom1, atom2], distance, out_of_plane_angle_unitless, in_plane_angle)
                
            # Try passing in a unitless in_plane_angle
            with self.assertRaises(AssertionError) as context:
                vsite1_index = molecule.add_monovalent_lone_pair_virtual_site([atom1, atom2], distance, out_of_plane_angle, in_plane_angle_unitless)
                
            # Try giving two atoms
            with self.assertRaises(AssertionError) as context:
                vsite1_index = molecule.add_monovalent_lone_pair_virtual_site([atom1, atom2], distance, out_of_plane_angle, in_plane_angle)
                
            # Successfully make a virtual site
            vsite1_index = molecule.add_monovalent_lone_pair_virtual_site([atom1, atom2, atom3], distance, out_of_plane_angle, in_plane_angle)
            # TODO: Check if we get the same values back out from the @properties
            molecule_dict = molecule.to_dict()
            molecule2 = Molecule.from_dict(molecule_dict)
            assert molecule.to_dict() == molecule2.to_dict()

        
    def test_add_divalent_lone_pair_virtual_site(self):
        """Test addition of a DivalentLonePairVirtualSite to the Molecule"""
        for molecule in self.molecules:
            atom1 = molecule.atoms[0]
            atom2 = molecule.atoms[1]
            atom3 = molecule.atoms[2]
            atom4 = molecule.atoms[3]
            distance = 0.3 * unit.angstrom
            out_of_plane_angle = 30 * unit.degree
            in_plane_angle = 0.2 * unit.radian
            vsite1_index = molecule.add_divalent_lone_pair_virtual_site([atom1, atom2, atom3], distance, out_of_plane_angle, in_plane_angle)
            with self.assertRaises(AssertionError) as context:
                vsite1_index = molecule.add_divalent_lone_pair_virtual_site([atom1, atom2], distance, out_of_plane_angle, in_plane_angle)
            molecule_dict = molecule.to_dict()
            molecule2 = Molecule.from_dict(molecule_dict)
            assert molecule.to_dict() == molecule2.to_dict()

    def test_add_trivalent_lone_pair_virtual_site(self):
        """Test addition of a TrivalentLonePairVirtualSite to the Molecule"""
        for molecule in self.molecules:
            atom1 = molecule.atoms[0]
            atom2 = molecule.atoms[1]
            atom3 = molecule.atoms[2]
            atom4 = molecule.atoms[3]
            distance = 0.3 * unit.angstrom
            out_of_plane_angle = 30 * unit.degree
            in_plane_angle = 0.2 * unit.radian
            vsite1_index = molecule.add_trivalent_lone_pair_virtual_site([atom1, atom2, atom3, atom4], distance, out_of_plane_angle, in_plane_angle)
            # Test for assertion when giving too few atoms
            with self.assertRaises(AssertionError) as context:
                vsite1_index = molecule.add_trivalent_lone_pair_virtual_site([atom1, atom2, atom3], distance, out_of_plane_angle, in_plane_angle)
            molecule_dict = molecule.to_dict()
            molecule2 = Molecule.from_dict(molecule_dict)
            assert molecule.to_dict() == molecule2.to_dict()


    def test_n_particles(self):
        """Test n_particles property"""
        for molecule in self.molecules:
            n_particles = sum([1 for particle in molecule.particles])
            assert n_particles == molecule.n_particles

    def test_n_atoms(self):
        """Test n_atoms property"""
        for molecule in self.molecules:
            n_atoms = sum([1 for atom in molecule.atoms])
            assert n_atoms == molecule.n_atoms

    def test_n_virtual_sites(self):
        """Test n_virtual_sites property"""
        for molecule in self.molecules:
            n_virtual_sites = sum([1 for virtual_site in molecule.virtual_sites])
            assert n_virtual_sites == molecule.n_virtual_sites

    def test_n_bonds(self):
        """Test n_bonds property"""
        for molecule in self.molecules:
            n_bonds = sum([1 for bond in molecule.bonds])
            assert n_bonds == molecule.n_bonds

    def test_angles(self):
        """Test angles property"""
        for molecule in self.molecules:
            for angle in molecule.angles:
                assert angle[0].is_bonded_to(angle[1])
                assert angle[1].is_bonded_to(angle[2])

    def test_propers(self):
        """Test propers property"""
        for molecule in self.molecules:
            for proper in molecule.propers:
                assert proper[0].is_bonded_to(proper[1])
                assert proper[1].is_bonded_to(proper[2])
                assert proper[2].is_bonded_to(proper[3])
                assert not(proper[0].is_bonded_to(proper[2]))
                assert not(proper[0].is_bonded_to(proper[3]))
                assert not(proper[1].is_bonded_to(proper[3]))

    def test_impropers(self):
        """Test impropers property"""
        for molecule in self.molecules:
            for improper in molecule.impropers:
                assert improper[0].is_bonded_to(improper[1])
                assert improper[1].is_bonded_to(improper[2])
                assert improper[2].is_bonded_to(improper[3])
                
                assert ((improper[0].is_bonded_to(improper[2])) or
                        (improper[0].is_bonded_to(improper[3])) or
                        (improper[1].is_bonded_to(improper[3])))
                


    def test_torsions(self):
        """Test torsions property"""
        for molecule in self.molecules:
            assert frozenset(molecule.torsions) == frozenset(set(molecule.propers) | set(molecule.impropers))
            assert len(molecule.propers & molecule.impropers) == 0

    def test_total_charge(self):
        """Test total charge"""
        for molecule in self.molecules:
            total_charge = sum([atom.formal_charge for atom in molecule.atoms])
            assert total_charge == molecule.total_charge

    def test_chemical_environment_matches(self):
        """Test chemical environment matches"""
        # TODO: Move this to test_toolkits, test all available toolkits
        # Create chiral molecule
        from simtk.openmm.app import element
        molecule = Molecule()
        atom_C = molecule.add_atom(element.carbon.atomic_number, 0, False, stereochemistry='R', name='C')
        atom_H = molecule.add_atom(element.hydrogen.atomic_number, 0, False, name='H')
        atom_Cl = molecule.add_atom(element.chlorine.atomic_number, 0, False, name='Cl')
        atom_Br = molecule.add_atom(element.bromine.atomic_number, 0, False, name='Br')
        atom_F = molecule.add_atom(element.fluorine.atomic_number, 0, False, name='F')
        molecule.add_bond(atom_C, atom_H, 1, False)
        molecule.add_bond(atom_C, atom_Cl, 1, False)
        molecule.add_bond(atom_C, atom_Br, 1, False)
        molecule.add_bond(atom_C, atom_F, 1, False)
        # Test known cases
        matches = molecule.chemical_environment_matches('[#6:1]')
        assert len(matches) == 1 # there should be a unique match, so one atom tuple is returned
        assert len(matches[0]) == 1 # it should have one tagged atom
        assert set(matches[0]) == set([atom_C])
        matches = molecule.chemical_environment_matches('[#6:1]~[#1:2]')
        assert len(matches) == 1 # there should be a unique match, so one atom tuple is returned
        assert len(matches[0]) == 2 # it should have two tagged atoms
        assert set(matches[0]) == set([atom_C, atom_H])
        matches = molecule.chemical_environment_matches('[Cl:1]-[C:2]-[H:3]')
        assert len(matches) == 1 # there should be a unique match, so one atom tuple is returned
        assert len(matches[0]) == 3 # it should have three tagged atoms
        assert set(matches[0]) == set([atom_Cl, atom_C, atom_H])
        matches = molecule.chemical_environment_matches('[#6:1]~[*:2]')
        assert len(matches) == 4 # there should be four matches
        for match in matches:
            assert len(match) == 2 # each match should have two tagged atoms

    def test_name(self):
        """Test name property"""
        name = 'benzene'
        #molecule = Molecule(name=name)
        molecule = Molecule()
        molecule.name = name
        assert molecule.name == name

    # TODO: This should be a toolkit test
    #@pytest.mark.skipif(OPENEYE_UNAVAILABLE, reason=_OPENEYE_UNAVAILABLE_MESSAGE)
    @OpenEyeToolkitWrapper.requires_toolkit()
    def test_iupac_roundtrip(self):
        """Test IUPAC conversion"""
        for molecule in self.molecules:
            iupac = molecule.to_iupac()
            molecule_copy = Molecule.from_iupac(iupac)
            assert molecule == molecule_copy

    def test_topology_roundtrip(self):
        """Test Topology round-trip"""
        for molecule in self.molecules:
            topology = molecule.to_topology()
            molecule_copy = Molecule.from_topology(topology)
            assert molecule == molecule_copy

    # TODO: This should be a toolkit test
    def test_file_roundtrip(self):
        """Test to/from file"""
        import os
        # TODO: Test all file capabilities; the current test is minimal
        # TODO: This makes no sense as implemented (don't know which toolkit it uses for what). Make this separate tests in test_toolkits
        for molecule in self.molecules:
            # Write and read mol2 file
            with NamedTemporaryFile(suffix='.mol2', delete=False) as iofile:
                molecule.to_file(iofile.name, 'MOL2')
                molecule2 = Molecule.from_file(iofile.name)
                assert molecule == molecule2
                # TODO: Test to make sure properties are preserved?
                os.unlink(iofile.name)
            # Write and read SDF file
            with NamedTemporaryFile(suffix='.sdf', delete=False) as iofile:
                molecule.to_file(iofile.name, 'SDF')
                molecule2 = Molecule.from_file(iofile.name)
                assert molecule == molecule2
                # TODO: Test to make sure properties are preserved?
                os.unlink(iofile.name)
            # Write and read PDB file
            with NamedTemporaryFile(suffix='.pdb', delete=False) as iofile:
                molecule.to_file(iofile.name, 'PDB')
                # NOTE: We can't read pdb files and expect chemical information to be preserved
                os.unlink(iofile.name)

    #@pytest.mark.skipif(OPENEYE_UNAVAILABLE, reason=_OPENEYE_UNAVAILABLE_MESSAGE)
    @OpenEyeToolkitWrapper.requires_toolkit()
    def test_oemol_roundtrip(self):
        """Test creation of Molecule object from OpenEye OEMol
        """
        for molecule in self.molecules:
            oemol = molecule.to_openeye()
            molecule2 = Molecule.from_openeye(oemol)
            assert_molecule_is_equal(molecule, molecule2, "Molecule.to_openeye()/from_openeye() round trip failed")
            molecule3 = Molecule(oemol)
            assert_molecule_is_equal(molecule, molecule3, "Molecule(oemol) constructor failed")

    def test_compute_partial_charges(self):
        """Test computation/retrieval of partial charges"""
        # TODO: Test only one molecule for speed?
        # TODO: Do we need to deepcopy each molecule, or is setUp called separately for each test method?
        from simtk import unit
        import numpy as np
        # Test a single toolkit at a time
        # Removed  ['amber', 'amberff94'] from OE list, as those won't find the residue types they're expecting
        toolkit_to_charge_method = {OpenEyeToolkitWrapper:['mmff', 'mmff94', 'am1bcc', 'am1bccnosymspt', 'am1bccelf10'],
                                   AmberToolsToolkitWrapper:['bcc', 'gas', 'mul']}

        manual_skips = []

        manual_skips.append('ZINC1564378') # Warning: OEMMFF94Charges: assigning OEMMFFAtomTypes failed on mol .
        manual_skips.append('ZINC00265517') # Warning: OEMMFF94Charges: assigning OEMMFFAtomTypes failed on mol .

        for toolkit in list(toolkit_to_charge_method.keys()):
            toolkit_registry = ToolkitRegistry(toolkit_precedence=[toolkit])
            for charge_model in toolkit_to_charge_method[toolkit]:
                c = 0
                for molecule in self.molecules[:1]: # Just test first molecule to save time
                    c += 1
                    if molecule.name in manual_skips:  # Manual skips, hopefully rare
                        continue
                    molecule.compute_partial_charges(charge_model=charge_model, toolkit_registry=toolkit_registry)
                    charges1 = molecule._partial_charges
                    # Make sure everything isn't 0s
                    assert (abs(charges1 / unit.elementary_charge) > 0.01).any()
                    # Check total charge
                    charges_sum_unitless = charges1.sum() / unit.elementary_charge
                    #if abs(charges_sum_unitless - float(molecule.total_charge)) > 0.0001:
                    #    print('c {}  molecule {}    charge_sum {}     molecule.total_charge {}'.format(c, molecule.name,
                    #                                                                                   charges_sum_unitless,
                    #                                                                                   molecule.total_charge))
                    # assert_almost_equal(charges_sum_unitless, molecule.total_charge, decimal=4)

                    # Call should be faster second time due to caching
                    # TODO: Implement caching
                    molecule.compute_partial_charges(charge_model=charge_model, toolkit_registry=toolkit_registry)
                    charges2 = molecule._partial_charges
                    assert (np.allclose(charges1, charges2, atol=0.002))

    #
    #@pytest.mark.skipif(OPENEYE_UNAVAILABLE, reason=_OPENEYE_UNAVAILABLE_MESSAGE)    
    @OpenEyeToolkitWrapper.requires_toolkit()
    def test_assign_fractional_bond_orders(self):
        """Test assignment of fractional bond orders
        """
        # TODO: Test only one molecule for speed?
        # TODO: Do we need to deepcopy each molecule, or is setUp called separately for each test method?

        toolkit_to_bondorder_method = {OpenEyeToolkitWrapper:['am1','pm3']}
        for toolkit in list(toolkit_to_bondorder_method.keys()):
            toolkit_registry = ToolkitRegistry(toolkit_precedence=[toolkit])
            for charge_model in toolkit_to_bondorder_method[toolkit]:
                for molecule in self.molecules[:5]: # Just test first five molecules for speed
                    molecule.compute_wiberg_bond_orders(charge_model=charge_model, toolkit_registry=toolkit_registry)
                    fbo1 = [bond.fractional_bond_order for bond in molecule.bonds]
                    # Call should be faster the second time due to caching
                    molecule.compute_wiberg_bond_orders(charge_model=charge_model, toolkit_registry=toolkit_registry)
                    fbo2 = [bond.fractional_bond_order for bond in molecule.bonds]
                    assert fbo1 == fbo2
