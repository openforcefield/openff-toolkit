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

from openforcefield import utils, topology
from openforcefield.topology.molecule import FrozenMolecule, Molecule, Atom, Bond, BondChargeVirtualSite, MonovalentLonePairVirtualSite, DivalentLonePairVirtualSite, TrivalentLonePairVirtualSite, ALLOWED_CHARGE_MODELS, ALLOWED_FRACTIONAL_BONDORDER_MODELS
from openforcefield.utils import get_data_filename
# TODO: Will the ToolkitWrapper allow us to pare that down?
#from openforcefield.utils import RDKIT_UNAVAILABLE, OPENEYE_UNAVAILABLE, SUPPORTED_TOOLKITS, TOOLKIT_PRECEDENCE, SUPPORTED_FILE_FORMATS
from openforcefield.utils.toolkits import ToolkitWrapper, OpenEyeToolkitWrapper, RDKitToolkitWrapper, AmberToolsToolkitWrapper, ToolkitRegistry


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
    # TODO:
    pass

def test_cheminformatics_toolkit_is_installed():
    """Ensure that at least one supported cheminformatics toolkit is installed."""
    # Can't just call all_subclasses here -- AmberToolsToolkitWrapper won't suffice for this test
    if not(RDKitToolkitWrapper.toolkit_is_available) and not(OpenEyeToolkitWrapper.toolkit_is_available):

        msg = 'No supported cheminformatics toolkits are installed. Please install a supported toolkit:\n'
        msg += 'OpeneEyeToolkitWrapper or RDKitToolkitWrapper'
        raise Exception(msg)

class TestMolecule(TestCase):
    from openforcefield.topology import Atom, Bond, Molecule

    def setUp(self):
        # TODO: Serialize the offmols instead so that we can run this test without toolkits
        #self.molecules = pickle.load('zinc-subset-offmols.pkl')
        filename = get_data_filename('molecules/zinc-subset-tripos.mol2.gz')
        self.molecules = Molecule.from_file(filename)

    # TODO: Test getstate/setstate
    # TODO: Test {to_from}_{dict|yaml|toml|json|bson|messagepack|pickle}
    def test_pickle(self):
        """Test pickling"""
        serialized = pickle.dumps(self.molecules)
        molecules_copy = pickle.loads(serialized)
        for (molecule, molecule_copy) in zip(self.molecules, molecules_copy):
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
            atom = Atom(this_element.atomic_number, formal_charge, is_aromatic)
            assert atom.atomic_number == this_element.atomic_number
            assert atom.element == this_element
            assert atom.mass == this_element.mass
            assert atom.formal_charge == formal_charge
            assert atom.is_aromatic == is_aromatic

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

    def test_create_empty(self):
        """Test creation of an empty Molecule"""
        molecule = Molecule()

    def test_create_copy(self):
        """Test creation of a Molecule from another Molecule"""
        for molecule in self.molecules:
            molecule_copy = Molecule(molecule)
            assert molecule_copy == molecule

    def test_create_openeye(self):
        """Test creation of a molecule from an OpenEye oemol"""
        for molecule in self.molecules:
            oemol = molecule.to_openeye()
            molecule_copy = Molecule(oemol)
            assert molecule == molecule_copy

    def test_create_rdkit(self):
        """Test creation of a molecule from an RDKit rdmol"""
        for molecule in self.molecules:
            rdmol = molecule.to_rdkit()
            molecule_copy = Molecule(rdmol)
            assert molecule == molecule_copy

    def test_create_from_file(self):
        """Test creation of a molecule from a filename or file-like object"""
        # TODO: Expand test to both openeye and rdkit toolkits
        filename = get_data_filename('molecules/toluene.mol2')
        molecule1 = Molecule(filename)
        with open(filename, 'r') as infile:
            molecule2 = Molecule(infile)
        assert molecule1 == molecule2
        import gzip
        with gzip.GzipFile(filename + '.gz', 'r') as infile:
            molecule3 = Molecule(infile)
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

    #def test_equality(self):
    #    """Test equality operator"""
    #    nmolecules = len(self.molecules)
    #    for i in range(nmolecules):
    #        for j in range(i, nmolecules):
    #            assert (self.molecules[i] == self.molecules[j]) == (i == j)

    def test_smiles_round_trip(self):
        """Test SMILES round-trip"""
        for molecule in self.molecules:
            smiles = molecule.to_smiles()
            molecule_copy = Molecule.from_smiles(smiles)
            assert molecule == molecule_copy

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
                molecule_copy.add_bond(bond.atom1_index, bond.atom2_index, bond.bond_order, bond.is_aromatic, stereochemistry=bond.stereochemistry)
            assert molecule == molecule_copy
            
    def test_add_bond_charge_virtual_site(self):
        """Test the addition of a BondChargeVirtualSite to a molecule.
        Also tests many of the input tests of the parent VirtualSite class"""
        # TODO: Add test for units in VdW and electrostatic parameters
        for molecule in self.molecules:
            atom1 = molecule.atoms[0]
            atom2 = molecule.atoms[1]
            atom3 = molecule.atoms[2]
            atom4 = molecule.atoms[3]
            # Try to feed in a unitless distance
            distance = 0.3
            with self.assertRaises(AssertionError) as context:
                vsite1_index = molecule.add_bond_charge_virtual_site([atom1, atom2, atom3], distance)

            distance = 0.3 * unit.angstrom
            vsite1_index = molecule.add_bond_charge_virtual_site([atom1, atom2, atom3], distance)
            vsite1 = molecule.virtual_sites[vsite1_index]
            assert atom1 in vsite1.atoms
            assert atom2 in vsite1.atoms
            assert atom3 in vsite1.atoms
            assert vsite1 in atom1.virtual_sites
            assert vsite1 in atom2.virtual_sites
            assert vsite1 in atom3.virtual_sites
            assert vsite1.distance == distance

            # We shouldn't be able to give both rmin_half and sigma VdW parameters.
            with self.assertRaises(Exception) as context:
                molecule.add_bond_charge_virtual_site([atom1, atom2, atom3], distance, epsilon=0.1, sigma=0.1, rmin_half=0.1)
            vsite2_index = molecule.add_bond_charge_virtual_site([atom1, atom2, atom3], distance, epsilon=0.1, sigma=0.1)
            vsite3_index = molecule.add_bond_charge_virtual_site([atom1, atom2, atom3], distance, epsilon=0.1, rmin_half=0.1)

            # We should have to give as many charge increments as atoms
            with self.assertRaises(Exception) as context:
                molecule.add_bond_charge_virtual_site([atom1, atom2, atom3], distance, charge_increments=[0.1, 0.2, 0.25, 0.3], epsilon=0.1, rmin_half=0.1)
            vsite4_index = molecule.add_bond_charge_virtual_site([atom1, atom2, atom3], distance, charge_increments=[0.1, 0.2, 0.25])
            
    def test_add_monovalent_lone_pair_virtual_site(self):
        """Test addition of a MonovalentLonePairVirtualSite to the Molecule"""
        for molecule in self.molecules:
            atom1 = molecule.atoms[0]
            atom2 = molecule.atoms[1]
            atom3 = molecule.atoms[2]
            atom4 = molecule.atoms[3]
            distance = 0.3 * unit.angstrom
            out_of_plane_angle = 30 * unit.degree
            in_plane_angle = 0.2 * unit.radian
            vsite1_index = molecule.add_monovalent_lone_pair_virtual_site([atom1, atom2, atom3], distance, out_of_plane_angle, in_plane_angle)
            with self.assertRaises(AssertionError) as context:
                vsite1_index = molecule.add_monovalent_lone_pair_virtual_site([atom1, atom2], distance, out_of_plane_angle, in_plane_angle)

            out_of_plane_angle = 30 
            with self.assertRaises(AssertionError) as context:
                vsite1_index = molecule.add_monovalent_lone_pair_virtual_site([atom1, atom2], distance, out_of_plane_angle, in_plane_angle)
        
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

    def test_file_roundtrip(self):
        """Test to/from file"""
        # TODO: Test all file capabilities; the current test is minimal
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

    #@pytest.mark.skipif(RDKIT_UNAVAILABLE, reason=_RDKIT_UNAVAILABLE_MESSAGE)
    @RDKitToolkitWrapper.requires_toolkit()
    def test_rdkit_roundtrip(self):
        for molecule in self.molecules:
            rdmol = molecule.to_rdkit()
            molecule2 = Molecule.from_rdkit(rdmol)
            assert_molecule_is_equal(molecule, molecule2, "Molecule.to_rdkit()/from_rdkit() round trip failed")
            molecule3 = Molecule(rdmol)
            assert_molecule_is_equal(molecule, molecule3, "Molecule(rdmol) constructor failed")

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

    def test_get_partial_charges(self):
        """Test computation/retrieval of partial charges"""
        # TODO: Test only one molecule for speed?
        # TODO: Do we need to deepcopy each molecule, or is setUp called separately for each test method?

        # Test a single toolkit at a time
        old_toolkit_precedence = TOOLKIT_PRECEDENCE
        for toolkit in list(TOOLKIT_PRECEDENCE):
            for charge_model in ALLOWED_CHARGE_MODELS:
                for molecule in self.molecules:
                    charges1 = molecule.get_partial_charges(method=charge_model)
                    # Check total charge
                    assert_almost_equal(charges1.sum(), molecule.total_charge)

                    # Call should be faster second time due to caching
                    charges2 = molecule.get_partial_charges(method=charge_model)
                    assert charges1 == charges2


        # Restore toolkit precedence
        TOOLKIT_PRECEDENCE = old_toolkit_precedence

    #@pytest.mark.skipif(OPENEYE_UNAVAILABLE, reason=_OPENEYE_UNAVAILABLE_MESSAGE)    
    @OpenEyeToolkitWrapper.requires_toolkit()
    def test_assign_fractional_bond_orders(self):
        """Test assignment of fractional bond orders
        """
        # TODO: Test only one molecule for speed?
        # TODO: Do we need to deepcopy each molecule, or is setUp called separately for each test method?

        for molecule in self.molecules:
            for method in ALLOWED_FRACTIONAL_BONDORDER_MODELS:
                fbo1 = molecule.assign_fractional_bond_orders(method=method)
                # Call should be faster the second time due to caching
                fbo2 = molecule.assign_fractional_bond_orders(method=method)
                assert fbo1 == fbo2
