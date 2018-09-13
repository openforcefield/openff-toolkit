#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Tests for Molecule

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

import pytest

from openforcefield import utils, topology
from openforcefield.topology import Molecule, ALLOWED_CHARGE_MODELS, ALLOWED_FRACTIONAL_BONDORDER_MODELS
from openforcefield.utils import get_data_filename
from openforcefield.utils import RDKIT_UNAVAILABLE, OPENEYE_UNAVAILABLE, SUPPORTED_TOOLKITS, TOOLKIT_PRECEDENCE, SUPPORTED_FILE_FORMATS

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
    if RDKIT_UNAVAILABLE and OPENEYE_UNAVAILABLE:
        msg = 'No supported cheminformatics toolkits are installed. Please install a supported toolkit:\n'
        msg += str(SUPPORTED_TOOLKITS)
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
        for element in elements:
            atom = Atom(element.atomic_number, formal_charge, is_aromatic)
            assert atom.atomic_number == element.atomic_number
            assert atom.element == element
            assert atom.mass == element.mass
            assert atom.formal_charge == formal_charge
            assert atom.is_aromatic = is_aromatic

    def test_create_molecule(self):
        """Test creation of molecule by adding molecules and bonds"""
        # Define a methane molecule
        molecule = Molecule(name='methane')
        C = molecule.add_atom(6, 0, False)
        H1 = molecule.add_atom(1, 0, False)
        H2 = molecule.add_atom(1, 0, False)
        H3 = molecule.add_atom(1, 0, False)
        H4 = molecule.add_atom(1, 0, False)
        molecule.add_bond(C, H1, False, 1)
        molecule.add_bond(C, H2, False, 1)
        molecule.add_bond(C, H3, False, 1)
        molecule.add_bond(C, H4, False, 1)

    def test_create_empty(self):
        """Test creation of an empty Molecule"""
        molecule = Molecule()

    def test_create_copy(self):
        """Test creation of a Molecule from another Molecule"""
        for molecule in self.molecules:
            molecule_copy = Molecule(molecule)

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

    def test_create_file(self):
        """Test creation of a molecule from a filename or file-like object"""
        filename = get_data_filename('molecules/toluene.mol2')
        molecule1 = Molecule(filename)
        with open(filename, 'r') as infile:
            molecule2 = Molecule(infile)
        assert molecule1 == molecule2
        import gzip
        with gzip.GzipFile(filename + '.gz', 'r') as infile:
            molecule3 == Molecule(infile)
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

    def test_equality(self):
        """Test equality operator"""
        nmolecules = len(self.molecules)
        for i in range(nmolecules):
            for j in range(i, nmolecules):
                assert (self.molecules[i] == self.molecules[j]) == (i == j)

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
                molecule_copy.add_atom(atomic_number=atom.atomic_number, formal_charge=atom.formal_charge, is_aromatic=atom.is_aromatic, stereochemistry=atom.stereochemistry)
            for bond in molecule.bonds:
                molecule_copy.add_bond(bond.atom1_index, bond.atom2_index, is_aromatic=bond.is_aromatic, order=bond.order, stereochemistry=bond.stereochemistry)
            assert molecule == molecule_copy

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
                assert angle[0].bonded_to(angle[1])
                assert angle[1].bonded_to(angle[1])

    def test_propers(self):
        """Test propers property"""
        for molecule in self.molecules:
            for proper in molecule.propers:
                assert angle[0].bonded_to(angle[1])
                assert angle[1].bonded_to(angle[2])
                assert angle[2].bonded_to(angle[3])

    def test_impropers(self):
        """Test impropers property"""
        for molecule in self.molecules:
            for improper in molecule.impropers:
                # TODO: Check improper bonds
                pass

    def test_torsions(self):
        """Test torsions property"""
        for molecule in self.molecules:
            assert frozenset(molecule.torsions) == frozenset(set(molecule.propers) + set(molecule.impropers))

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
        molecule.add_bond(atom_C, atom_H)
        molecule.add_bond(atom_C, atom_Cl)
        molecule.add_bond(atom_C, atom_Br)
        molecule.add_bond(atom_C, atom_F)
        # Test known cases
        matches = molecule.chemical_environment_matches('[#6:1]')
        assert len(matches) == 1 # there should be a unique match, so one atom tuple is returned
        assert len(matches[0]) == 1 # it should have one tagged atom
        assert set(matches[0]) = set(atom_C)
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
        molecule = Molecule(name=name)
        assert molecule.name == name

    @pytest.mark.skipif(OPENEYE_UNAVAILABLE, reason=_OPENEYE_UNAVAILABLE_MESSAGE)
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
            with NamedTemporaryFile(suffix='.mol2', delete=False) as file:
                molecule.to_file(file.name, format='MOL2')
                molecule2 = Molecule.from_file(file.name)
                assert molecule == molecule2
                # TODO: Test to make sure properties are preserved?
                os.unlink(file.name)
            # Write and read SDF file
            with NamedTemporaryFile(suffix='.sdf', delete=False) as file:
                molecule.to_file(file.name, format='SDF')
                molecule2 = Molecule.from_file(file.name)
                assert molecule == molecule2
                # TODO: Test to make sure properties are preserved?
                os.unlink(file.name)
            # Write and read SDF file
            with NamedTemporaryFile(suffix='.pdb', delete=False) as file:
                molecule.to_file(file.name, format='PDB')
                # NOTE: We can't read mol2 files and expect chemical information to be preserved
                os.unlink(file.name)

    @pytest.mark.skipif(RDKIT_UNAVAILABLE, reason=_RDKIT_UNAVAILABLE_MESSAGE)
    def test_rdkit_roundtrip(self):
        for molecule in self.molecules:
            rdmol = molecule.to_rdkit()
            molecule2 = Molecule.from_rdmol(rdmol)
            assert_molecule_is_equal(molecule, molecule2, "Molecule.to_rdmol()/from_rdmol() round trip failed")
            molecule3 = Molecule(rdmol)
            assert_molecule_is_equal(molecule, molecule3, "Molecule(rdmol) constructor failed")

    @pytest.mark.skipif(OPENEYE_UNAVAILABLE, reason=_OPENEYE_UNAVAILABLE_MESSAGE)
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

    @pytest.mark.skipif(OPENEYE_UNAVAILABLE, reason=_OPENEYE_UNAVAILABLE_MESSAGE)
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
