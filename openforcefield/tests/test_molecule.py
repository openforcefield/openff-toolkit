#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Tests for molecular topology representations

At least one supported cheminformatics toolkit must be installed to run these tests.
Only the tests applicable to that toolkit will be run.

TODO:
- Add tests comparing RDKit and OpenEye aromaticity perception
- Right now, the test database of TestMolecule is read from mol2, requiring the OE
  toolkit. Find a different test set that RDKit can read, or make a database of
  serialized OFFMols.

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import copy
import os
import pickle
from tempfile import NamedTemporaryFile

import numpy as np
import pytest
from simtk import unit

from openforcefield.topology.molecule import Molecule, Atom
from openforcefield.utils import get_data_file_path
# TODO: Will the ToolkitWrapper allow us to pare that down?
from openforcefield.utils.toolkits import OpenEyeToolkitWrapper, RDKitToolkitWrapper, AmberToolsToolkitWrapper, ToolkitRegistry


#=============================================================================================
# TEST UTILITIES
#=============================================================================================

requires_openeye = pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(),
                                      reason='Test requires OE toolkit')
requires_rdkit = pytest.mark.skipif(not RDKitToolkitWrapper.is_available(),
                                    reason='Test requires RDKit')


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
        raise AssertionError(msg)


def is_four_memebered_ring_torsion(torsion):
    """Check that three atoms in the given torsion form a four-membered ring."""
    # Push a copy of the first and second atom in the end to make the code simpler.
    torsion = list(torsion) + [torsion[0], torsion[1]]

    is_four_membered_ring = True
    for i in range(4):
        # The atom is bonded to the next one.
        is_four_membered_ring &= torsion[i].is_bonded_to(torsion[i+1])
        # The atom is not bonded to the atom on its diagonal.
        is_four_membered_ring &= not torsion[i].is_bonded_to(torsion[i+2])

    return is_four_membered_ring


def is_three_memebered_ring_torsion(torsion):
    """Check that three atoms in the given torsion form a three-membered ring.

    In order to be 4 atoms with a three-membered ring, there must be
    1) A central atom connected to all other atoms.
    2) An atom outside the ring connected exclusively to the central atom.
    3) Two atoms in the ring connected to the central atom and to each other.

    """
    # A set of atom indices for the atoms in the torsion.
    torsion_atom_indices = set(a.molecule_atom_index for a in torsion)

    # Collect all the bonds involving exclusively atoms in the torsion.
    bonds_by_atom_idx = {i: set() for i in torsion_atom_indices}
    for atom in torsion:
        for bond in atom.bonds:
            # Consider the bond only if both atoms are in the torsion.
            if (bond.atom1_index in torsion_atom_indices and
                        bond.atom2_index in torsion_atom_indices):
                bonds_by_atom_idx[bond.atom1_index].add(bond.atom2_index)
                bonds_by_atom_idx[bond.atom2_index].add(bond.atom1_index)

    # Find the central atom, which is connected to all other atoms.
    atom_indices = [i for i in torsion_atom_indices if len(bonds_by_atom_idx[i]) == 3]
    if len(atom_indices) != 1:
        return False
    central_atom_idx = atom_indices[0]

    # Find the atom outside the ring.
    atom_indices = [i for i in torsion_atom_indices if len(bonds_by_atom_idx[i]) == 1]
    if len(atom_indices) != 1 or central_atom_idx not in bonds_by_atom_idx[atom_indices[0]]:
        return False
    outside_atom_idx = atom_indices[0]

    # Check that the remaining two atoms are non-central atoms in the membered ring.
    atom1, atom2 = [i for i in torsion_atom_indices if i not in [central_atom_idx, outside_atom_idx]]
    # The two atoms are bonded to each other.
    if atom2 not in bonds_by_atom_idx[atom1] or atom1 not in bonds_by_atom_idx[atom2]:
        return False
    # Check that they are both bonded to the central atom and none other.
    for atom_idx in [atom1, atom2]:
        if (central_atom_idx not in bonds_by_atom_idx[atom_idx] or
                    len(bonds_by_atom_idx[atom_idx]) != 2):
            return False

    # This is a torsion including a three-membered ring.
    return True


#=============================================================================================
# FIXTURES
#=============================================================================================

def mini_drug_bank(xfail_mols=None, wip_mols=None):
    """Load the full MiniDrugBank into Molecule objects.

    Parameters
    ----------
    xfail_mols : Dict[str, str or None]
        Dictionary mapping the molecule names that are allowed to
        failed to the failure reason.
    wip_mols : Dict[str, str or None]
        Dictionary mapping the molecule names that are work in progress
        to the failure reason.

    """
    # If we have already loaded the data set, return the cached one.
    if mini_drug_bank.molecules is not None:
        molecules = mini_drug_bank.molecules
    else:
        # Load the dataset.
        file_path = get_data_file_path('molecules/MiniDrugBank_tripos.mol2')
        try:
            # We need OpenEye to parse the molecules, but pytest execute this
            # whether or not the test class is skipped so if OE is not available
            # we just return an empty list of test cases as a workaround.
            molecules = Molecule.from_file(file_path, allow_undefined_stereo=True)
        except NotImplementedError as e:
            assert 'No toolkits in registry can read file' in str(e)
            mini_drug_bank.molecules = []
            return []
        else:
            mini_drug_bank.molecules = molecules

    # Check if we need to mark anything.
    if xfail_mols is None and wip_mols is None:
        return molecules

    # Handle mutable default.
    if xfail_mols is None:
        xfail_mols = {}
    if wip_mols is None:
        wip_mols = {}
    # There should be no molecule in both dictionaries.
    assert len(set(xfail_mols).intersection(set(wip_mols))) == 0

    # Don't modify the cached molecules.
    molecules = copy.deepcopy(molecules)
    for i, mol in enumerate(molecules):
        if mol.name in xfail_mols:
            marker = pytest.mark.xfail(reason=xfail_mols[mol.name])
        elif mol.name in wip_mols:
            marker = pytest.mark.wip(reason=wip_mols[mol.name])
        else:
            marker = None

        if marker is not None:
            molecules[i] = pytest.param(mol, marks=marker)

    return molecules

# Use a "static" variable as a workaround as fixtures cannot be
# used inside pytest.mark.parametrize (see issue #349 in pytest).
mini_drug_bank.molecules = None


#=============================================================================================
# TESTS
#=============================================================================================

class TestAtom:
    """Test Atom class."""

    def test_atom_constructor(self):
        """Test Atom creation"""
        # Create a non-aromatic carbon atom
        atom1 = Atom(6, 0, False)
        assert atom1.atomic_number == 6
        assert atom1.formal_charge == 0

        # Create a chiral carbon atom
        atom2 = Atom(6, 0, False, stereochemistry='R', name='CT')
        assert atom1.stereochemistry != atom2.stereochemistry

    def test_atom_properties(self):
        """Test that atom properties are correctly populated and gettable"""
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


@requires_openeye
class TestMolecule:
    """Test Molecule class."""

    # TODO: Test getstate/setstate
    # TODO: Test {to_from}_{dict|yaml|toml|json|bson|messagepack|pickle}

    @pytest.mark.parametrize('molecule', mini_drug_bank())
    def test_pickle_serialization(self, molecule):
        """Test pickling of a molecule object."""
        serialized = pickle.dumps(molecule)
        molecule_copy = pickle.loads(serialized)
        assert molecule == molecule_copy

    # ----------------------------------------------------
    # Test Molecule constructors and conversion utilities.
    # ----------------------------------------------------

    def test_create_empty(self):
        """Test empty constructor."""
        molecule = Molecule()
        assert len(molecule.atoms) == 0
        assert len(molecule.bonds) == 0

    @pytest.mark.parametrize('molecule', mini_drug_bank())
    def test_create_copy(self, molecule):
        """Test copy constructor."""
        molecule_copy = Molecule(molecule)
        assert molecule_copy == molecule

    # TODO: Should there be an equivalent toolkit test and leave this as an integration test?
    @pytest.mark.slow
    def test_create_from_file(self):
        """Test standard constructor taking a filename or file-like object."""
        # TODO: Expand test to both openeye and rdkit toolkits
        filename = get_data_file_path('molecules/toluene.mol2')

        molecule1 = Molecule(filename, allow_undefined_stereo=True)
        with open(filename, 'r') as infile:
            molecule2 = Molecule(infile, file_format='MOL2', allow_undefined_stereo=True)
        assert molecule1 == molecule2

        import gzip
        with gzip.GzipFile(filename + '.gz', 'r') as infile:
            molecule3 = Molecule(infile, file_format='MOL2', allow_undefined_stereo=True)
        assert molecule3 == molecule1

        # Ensure that attempting to initialize a single Molecule from a file
        # containing multiple molecules raises a ValueError
        with pytest.raises(ValueError) as exc_info:
            filename = get_data_file_path('molecules/zinc-subset-tripos.mol2.gz')
            molecule = Molecule(filename, allow_undefined_stereo=True)

    @pytest.mark.parametrize('molecule', mini_drug_bank())
    def test_create_from_serialized(self, molecule):
        """Test standard constructor taking the output of __getstate__()."""
        serialized_molecule = molecule.__getstate__()
        molecule_copy = Molecule(serialized_molecule)
        assert molecule == molecule_copy

    @pytest.mark.parametrize('molecule', mini_drug_bank())
    def test_to_from_dict(self, molecule):
        """Test that conversion/creation of a molecule to and from a dict is consistent."""
        serialized = molecule.to_dict()
        molecule_copy = Molecule.from_dict(serialized)
        assert molecule == molecule_copy

    @pytest.mark.parametrize('molecule', mini_drug_bank())
    def test_to_networkx(self, molecule):
        """Test conversion to NetworkX graph."""
        graph = molecule.to_networkx()

    @requires_rdkit
    @pytest.mark.parametrize('molecule', mini_drug_bank())
    def test_to_from_rdkit(self, molecule):
        """Test that conversion/creation of a molecule to and from an RDKit rdmol is consistent.

        This tests creating an OpenFF Molecule from an RDKit Mol both
        through __init__() and from_rdkit(). However, __init__() doesn't
        have an allow_undefined_stereo argument yet, so in that case, we
        check for equality only for the from_rdkit() molecule.

        """
        import pickle
        from openforcefield.utils.toolkits import UndefinedStereochemistryError

        rdmol = molecule.to_rdkit()

        # The constructor should not change the molecule.
        rdmol_pickle = pickle.dumps(rdmol)

        # Check if this is a molecule with undefined stereo.
        molecule_copies = []
        try:
            molecule_copies.append(Molecule(rdmol))
        except UndefinedStereochemistryError:
            allow_undefined_stereo = True
        else:
            allow_undefined_stereo = False
        molecule_copies.append(Molecule.from_rdkit(
            rdmol, allow_undefined_stereo=allow_undefined_stereo))

        # Check that the roundtrip did not change anything in the OpenFF Molecule.
        for molecule_copy in molecule_copies:
            assert molecule == molecule_copy

        # Check that the constructor didn't modify rdmol.
        assert rdmol_pickle == pickle.dumps(rdmol)

    # TODO: Should there be an equivalent toolkit test and leave this as an integration test?
    @requires_openeye
    @pytest.mark.parametrize('molecule', mini_drug_bank(
        xfail_mols={
            'DrugBank_2397': 'OpenEye cannot generate a correct IUPAC name and raises a "Warning: Incorrect name:" or simply return "BLAH".',
            'DrugBank_2543': 'OpenEye cannot generate a correct IUPAC name and raises a "Warning: Incorrect name:" or simply return "BLAH".',
            'DrugBank_2642': 'OpenEye cannot generate a correct IUPAC name and raises a "Warning: Incorrect name:" or simply return "BLAH".',
        },
        wip_mols={
            'DrugBank_1212': 'the roundtrip generates molecules with very different IUPAC/SMILES!',
            'DrugBank_2210': 'the roundtrip generates molecules with very different IUPAC/SMILES!',
            'DrugBank_4584': 'the roundtrip generates molecules with very different IUPAC/SMILES!',

            'DrugBank_390': 'raises warning "Unable to make OFFMol from OEMol: OEMol has unspecified stereochemistry."',
            'DrugBank_810': 'raises warning "Unable to make OFFMol from OEMol: OEMol has unspecified stereochemistry."',
            'DrugBank_4316': 'raises warning "Unable to make OFFMol from OEMol: OEMol has unspecified stereochemistry."',
            'DrugBank_7124': 'raises warning "Unable to make OFFMol from OEMol: OEMol has unspecified stereochemistry."',

            'DrugBank_4346': 'raises warning "Failed to parse name:"',
        }
    ))
    def test_to_from_iupac(self, molecule):
        """Test that conversion/creation of a molecule to and from a IUPAC name is consistent."""
        from openforcefield.utils.toolkits import UndefinedStereochemistryError

        # All the molecules that raise UndefinedStereochemistryError in Molecule.from_iupac()
        undefined_stereo_mols = {'DrugBank_977', 'DrugBank_1634', 'DrugBank_1700', 'DrugBank_1962',
                                 'DrugBank_2148', 'DrugBank_2178', 'DrugBank_2186', 'DrugBank_2208',
                                 'DrugBank_2519', 'DrugBank_2538', 'DrugBank_2592', 'DrugBank_2651',
                                 'DrugBank_2987', 'DrugBank_3332', 'DrugBank_3502', 'DrugBank_3622',
                                 'DrugBank_3726', 'DrugBank_3844', 'DrugBank_3930', 'DrugBank_4161',
                                 'DrugBank_4162', 'DrugBank_4778', 'DrugBank_4593', 'DrugBank_4959',
                                 'DrugBank_5043', 'DrugBank_5076', 'DrugBank_5176', 'DrugBank_5418',
                                 'DrugBank_5737', 'DrugBank_5902', 'DrugBank_6304', 'DrugBank_6305',
                                 'DrugBank_6329', 'DrugBank_6355', 'DrugBank_6401', 'DrugBank_6509',
                                 'DrugBank_6531', 'DrugBank_6647',

                                 # These test cases are allowed to fail.
                                 'DrugBank_390', 'DrugBank_810', 'DrugBank_4316', 'DrugBank_4346',
                                 'DrugBank_7124'
                                 }
        undefined_stereo = molecule.name in undefined_stereo_mols

        iupac = molecule.to_iupac()

        if undefined_stereo:
            with pytest.raises(UndefinedStereochemistryError):
                Molecule.from_iupac(iupac)

        molecule_copy = Molecule.from_iupac(iupac, allow_undefined_stereo=undefined_stereo)
        assert molecule.is_isomorphic(molecule_copy, compare_atom_stereochemistry=not undefined_stereo)

    @pytest.mark.parametrize('molecule', mini_drug_bank())
    def test_to_from_topology(self, molecule):
        """Test that conversion/creation of a molecule to and from a Topology is consistent."""
        topology = molecule.to_topology()
        molecule_copy = Molecule.from_topology(topology)
        assert molecule == molecule_copy

    # TODO: Should there be an equivalent toolkit test and leave this as an integration test?
    @pytest.mark.parametrize('molecule', mini_drug_bank())
    @pytest.mark.parametrize('format', [
        'mol2',
        'sdf',
        pytest.param('pdb', marks=pytest.mark.wip(reason='Read from pdb has not bee implemented properly yet'))
    ])
    def test_to_from_file(self, molecule, format):
        """Test that conversion/creation of a molecule to and from a file is consistent."""
        from openforcefield.utils.toolkits import UndefinedStereochemistryError
        # TODO: Test all file capabilities; the current test is minimal

        # TODO: This is only for OE. Expand to both OE and RDKit toolkits.
        # Molecules that are known to raise UndefinedStereochemistryError.
        undefined_stereo_mols = {'DrugBank_1700', 'DrugBank_2987', 'DrugBank_3502', 'DrugBank_4161',
                                 'DrugBank_4162', 'DrugBank_6531'}
        undefined_stereo = molecule.name in undefined_stereo_mols

        # The file is automatically deleted outside the with-clause.
        with NamedTemporaryFile(suffix='.' + format) as iofile:
            # If this has undefined stereo, check that the exception is raised.
            extension = os.path.splitext(iofile.name)[1][1:]
            molecule.to_file(iofile.name, extension)
            if undefined_stereo:
                with pytest.raises(UndefinedStereochemistryError):
                    Molecule.from_file(iofile.name)
            molecule2 = Molecule.from_file(iofile.name, allow_undefined_stereo=undefined_stereo)
            assert molecule == molecule2
            # TODO: Test to make sure properties are preserved?
            # NOTE: We can't read pdb files and expect chemical information to be preserved

    @requires_openeye
    @pytest.mark.parametrize('molecule', mini_drug_bank())
    def test_to_from_oemol(self, molecule):
        """Test that conversion/creation of a molecule to and from a OEMol is consistent."""
        from openforcefield.utils.toolkits import UndefinedStereochemistryError

        # Known failures raise an UndefinedStereochemistryError, but
        # the round-trip SMILES representation with the OpenEyeToolkit
        # doesn't seem to be affected.

        # ZINC test set known failures.
        # known_failures = {'ZINC05964684', 'ZINC05885163', 'ZINC05543156', 'ZINC17211981',
        #                   'ZINC17312986', 'ZINC06424847', 'ZINC04963126'}

        # DrugBank test set known failures.
        undefined_stereo_mols = {'DrugBank_1634', 'DrugBank_1700', 'DrugBank_1962',
                                 'DrugBank_2519', 'DrugBank_2987', 'DrugBank_3502',
                                 'DrugBank_3930', 'DrugBank_4161', 'DrugBank_4162',
                                 'DrugBank_5043', 'DrugBank_5418', 'DrugBank_6531'}
        undefined_stereo = molecule.name in undefined_stereo_mols

        error_messages = [
            "Molecule(oemol) constructor failed",
            "Molecule.to_openeye()/from_openeye() round trip failed"
        ]

        toolkit_wrapper = OpenEyeToolkitWrapper()

        oemol = molecule.to_openeye()
        molecule_copies = []

        # If this is a known failure, check that it raises UndefinedStereochemistryError
        # and proceed with the test ignoring it.
        if undefined_stereo:
            with pytest.raises(UndefinedStereochemistryError):
                Molecule(oemol)
        else:
            molecule_copies.append(Molecule(oemol))
        molecule_copies.append(Molecule.from_openeye(oemol, allow_undefined_stereo=undefined_stereo))

        molecule_smiles = molecule.to_smiles(toolkit_registry=toolkit_wrapper)
        for molecule_copy, error_msg in zip(molecule_copies, error_messages):
            # Check that the original and the copied molecules have the same SMILES representation.
            molecule_copy_smiles = molecule_copy.to_smiles(toolkit_registry=toolkit_wrapper)
            assert molecule_smiles == molecule_copy_smiles

            # Check that the two topologies are isomorphic.
            assert_molecule_is_equal(molecule, molecule_copy, error_msg)

    # ----------------------------------------------------
    # Test properties.
    # ----------------------------------------------------

    def test_name(self):
        """Test Molecule name property"""
        molecule1 = Molecule()
        molecule1.name = None

        molecule2 = Molecule()
        molecule2.name = ''
        assert molecule1.name == molecule2.name

        name = 'benzene'
        molecule = Molecule()
        molecule.name = name
        assert molecule.name == name

    @pytest.mark.parametrize('molecule', mini_drug_bank())
    def test_n_particles(self, molecule):
        """Test n_particles property"""
        n_particles = sum([1 for particle in molecule.particles])
        assert n_particles == molecule.n_particles

    @pytest.mark.parametrize('molecule', mini_drug_bank())
    def test_n_atoms(self, molecule):
        """Test n_atoms property"""
        n_atoms = sum([1 for atom in molecule.atoms])
        assert n_atoms == molecule.n_atoms

    @pytest.mark.parametrize('molecule', mini_drug_bank())
    def test_n_virtual_sites(self, molecule):
        """Test n_virtual_sites property"""
        n_virtual_sites = sum([1 for virtual_site in molecule.virtual_sites])
        assert n_virtual_sites == molecule.n_virtual_sites

    @pytest.mark.parametrize('molecule', mini_drug_bank())
    def test_n_bonds(self, molecule):
        """Test n_bonds property"""
        n_bonds = sum([1 for bond in molecule.bonds])
        assert n_bonds == molecule.n_bonds

    @pytest.mark.parametrize('molecule', mini_drug_bank())
    def test_angles(self, molecule):
        """Test angles property"""
        for angle in molecule.angles:
            assert angle[0].is_bonded_to(angle[1])
            assert angle[1].is_bonded_to(angle[2])

    @pytest.mark.parametrize('molecule', mini_drug_bank())
    def test_propers(self, molecule):
        """Test propers property"""
        for proper in molecule.propers:
            # The bonds should be in order 0-1-2-3 unless the
            # atoms form a three- or four-membered ring.
            is_chain = proper[0].is_bonded_to(proper[1])
            is_chain &= proper[1].is_bonded_to(proper[2])
            is_chain &= proper[2].is_bonded_to(proper[3])
            is_chain &= not proper[0].is_bonded_to(proper[2])
            is_chain &= not proper[0].is_bonded_to(proper[3])
            is_chain &= not proper[1].is_bonded_to(proper[3])

            assert (is_chain or
                    is_three_memebered_ring_torsion(proper) or
                    is_four_memebered_ring_torsion(proper))

    @pytest.mark.parametrize('molecule', mini_drug_bank())
    def test_impropers(self, molecule):
        """Test impropers property"""
        for improper in molecule.impropers:
            assert improper[0].is_bonded_to(improper[1])
            assert improper[1].is_bonded_to(improper[2])
            assert improper[1].is_bonded_to(improper[3])

            # The non-central atoms can be connected only if
            # the improper atoms form a three-membered ring.
            is_not_cyclic = not((improper[0].is_bonded_to(improper[2])) or
                                (improper[0].is_bonded_to(improper[3])) or
                                (improper[2].is_bonded_to(improper[3])))
            assert is_not_cyclic or is_three_memebered_ring_torsion(improper)

    @pytest.mark.parametrize('molecule', mini_drug_bank())
    def test_torsions(self, molecule):
        """Test torsions property"""
        # molecule.torsions should be exactly equal to the union of propers and impropers.
        assert set(molecule.torsions) == set(molecule.propers) | set(molecule.impropers)

        # The intersection of molecule.propers and molecule.impropers should be largely null.
        # The only exception is for molecules containing 3-membered rings (e.g., DrugBank_5514).
        common_torsions = molecule.propers & molecule.impropers
        if len(common_torsions) > 0:
            for torsion in common_torsions:
                assert is_three_memebered_ring_torsion(torsion)

    @pytest.mark.parametrize('molecule', mini_drug_bank())
    def test_total_charge(self, molecule):
        """Test total charge"""
        total_charge = sum([atom.formal_charge for atom in molecule.atoms])
        assert total_charge == molecule.total_charge

    # ----------------------------------------------------
    # Test magic methods.
    # ----------------------------------------------------

    def test_equality(self):
        """Test equality operator"""
        molecules = mini_drug_bank()
        nmolecules = len(molecules)
        # TODO: Performance improvements should let us un-restrict this test
        for i in range(nmolecules):
            for j in range(i, min(i+3, nmolecules)):
                assert (molecules[i] == molecules[j]) == (i == j)

    # ----------------------
    # Test Molecule methods.
    # ----------------------

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
        with pytest.raises(Exception) as excinfo:
            molecule.add_conformer(conf_missing_z)

        conf_too_few_atoms = unit.Quantity(np.array([[101., 102. ,103.], [104. ,105. ,106.], [107., 108., 109.],
                                                     [110.,111.,112.]]),
                                                     unit.angstrom)
        with pytest.raises(Exception) as excinfo:
            molecule.add_conformer(conf_too_few_atoms)


        # Add a conformer with too many coordinates
        conf_too_many_atoms = unit.Quantity(np.array([[101., 102., 103.], [104., 105., 106.], [107., 108., 109.],
                                                      [110., 111., 112.], [113., 114., 115.], [116., 117., 118.]]),
                                            unit.angstrom)
        with pytest.raises(Exception) as excinfo:
            molecule.add_conformer(conf_too_many_atoms)

        # Add a conformer with no coordinates
        conf_no_coordinates = unit.Quantity(np.array([]),
                                            unit.angstrom)
        with pytest.raises(Exception) as excinfo:
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
        with pytest.raises(Exception) as excinfo:
            molecule.add_conformer(conf_nonsense_units)

        # Add a conformer with no units
        conf_unitless = np.array([[ 1., 2.,3.] ,[4. ,5. ,6.],[7., 8., 9.],
                                  [10.,11.,12.],[13.,14.,15]])
        with pytest.raises(Exception) as excinfo:
            molecule.add_conformer(conf_unitless)

    @pytest.mark.parametrize('molecule', mini_drug_bank())
    def test_add_atoms_and_bonds(self, molecule):
        """Test the creation of a molecule from the addition of atoms and bonds"""
        molecule_copy = Molecule()
        for atom in molecule.atoms:
            molecule_copy.add_atom(atom.atomic_number, atom.formal_charge, atom.is_aromatic, stereochemistry=atom.stereochemistry)
        for bond in molecule.bonds:
            molecule_copy.add_bond(bond.atom1_index, bond.atom2_index, bond.bond_order, bond.is_aromatic,
                                   stereochemistry=bond.stereochemistry,
                                   fractional_bond_order=bond.fractional_bond_order)
        # Try to add the final bond twice, which should raise an Exception
        with pytest.raises(Exception) as excinfo:
            molecule_copy.add_bond(bond.atom1_index, bond.atom2_index, bond.bond_order, bond.is_aromatic,
                                   stereochemistry=bond.stereochemistry,
                                   fractional_bond_order=bond.fractional_bond_order)

        assert molecule == molecule_copy

    @pytest.mark.parametrize('molecule', mini_drug_bank())
    def test_add_virtual_site_units(self, molecule):
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
        charge_increments = charge_increments_unitless * unit.elementary_charge

        # Do not modify the original molecule.
        molecule = copy.deepcopy(molecule)

        atom1 = molecule.atoms[0]
        atom2 = molecule.atoms[1]
        atom3 = molecule.atoms[2]
        atom4 = molecule.atoms[3]

        # Try to feed in unitless sigma
        with pytest.raises(Exception) as excinfo:
            molecule.add_bond_charge_virtual_site([atom1, atom2, atom3], distance, epsilon=epsilon, sigma=sigma_unitless)

        # Try to feed in unitless rmin_half
        with pytest.raises(Exception) as excinfo:
            molecule.add_bond_charge_virtual_site([atom1, atom2, atom3], distance, epsilon=epsilon, rmin_half=rmin_half_unitless)

        # Try to feed in unitless epsilon
        with pytest.raises(Exception) as excinfo:
            molecule.add_bond_charge_virtual_site([atom1, atom2, atom3], distance, epsilon=epsilon_unitless, sigma=sigma, rmin_half=rmin_half)

        # Try to feed in unitless charges
        with pytest.raises(Exception) as excinfo:
            molecule.add_bond_charge_virtual_site([atom1, atom2, atom3, atom4], distance, charge_incrtements=charge_increments_unitless)


        # We shouldn't be able to give both rmin_half and sigma VdW parameters.
        with pytest.raises(Exception) as excinfo:
            molecule.add_bond_charge_virtual_site([atom1, atom2, atom3], distance, epsilon=epsilon, sigma=sigma, rmin_half=rmin_half)

        # Try creating virtual site from sigma+epsilon
        vsite1_index = molecule.add_bond_charge_virtual_site([atom1, atom2, atom3], distance, epsilon=epsilon, sigma=sigma)
        # Try creating virutal site from rmin_half+epsilon
        vsite2_index = molecule.add_bond_charge_virtual_site([atom1, atom2, atom3], distance, epsilon=epsilon, rmin_half=rmin_half)

        # TODO: Test the @property getters for sigma, epsilon, and rmin_half

        # We should have to give as many charge increments as atoms (len(charge_increments)) = 4
        with pytest.raises(Exception) as excinfo:
            molecule.add_bond_charge_virtual_site([atom1, atom2, atom3], distance, charge_increments=charge_increments)

        vsite3_index = molecule.add_bond_charge_virtual_site([atom1, atom2, atom3, atom4], distance, charge_increments=charge_increments)
            
    @pytest.mark.parametrize('molecule', mini_drug_bank())
    def test_add_bond_charge_virtual_site(self, molecule):
        """Test the addition of a BondChargeVirtualSite to a molecule.
           Also tests many of the inputs of the parent VirtualSite class
        """
        # Do not modify the original molecule.
        molecule = copy.deepcopy(molecule)

        atom1 = molecule.atoms[0]
        atom2 = molecule.atoms[1]
        atom3 = molecule.atoms[2]
        atom4 = molecule.atoms[3]

        # Prepare values for unit checks
        distance_unitless = 0.4
        distance = distance_unitless * unit.angstrom


        # Try to feed in a unitless distance
        with pytest.raises(AssertionError) as excinfo:
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

        # Make an "everything bagel" virtual site
        vsite2_index = molecule.add_bond_charge_virtual_site([atom1, atom2, atom3],
                                                             distance,
                                                             sigma=0.1*unit.angstrom,
                                                             epsilon=1.0*unit.kilojoule_per_mole,
                                                             charge_increments=unit.Quantity(np.array([0.1, 0.2, 0.3]),
                                                                                             unit.elementary_charge)
                                                             )
        vsite2 = molecule.virtual_sites[vsite2_index]

        # test serialization
        molecule_dict = molecule.to_dict()
        molecule2 = Molecule.from_dict(molecule_dict)

        assert hash(molecule) == hash(molecule2)

    # TODO: Make a test for to_dict and from_dict for VirtualSites (even though they're currently just unloaded using
    #      (for example) Molecule._add_bond_virtual_site functions
    @pytest.mark.parametrize('molecule', mini_drug_bank())
    def test_add_monovalent_lone_pair_virtual_site(self, molecule):
        """Test addition of a MonovalentLonePairVirtualSite to the Molecule"""
        # Do not modify the original molecule.
        molecule = copy.deepcopy(molecule)

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
        with pytest.raises(AssertionError) as excinfo:
            vsite1_index = molecule.add_monovalent_lone_pair_virtual_site([atom1, atom2], distance_unitless, out_of_plane_angle, in_plane_angle)

        # Try passing in a unitless out_of_plane_angle
        with pytest.raises(AssertionError) as excinfo:
            vsite1_index = molecule.add_monovalent_lone_pair_virtual_site([atom1, atom2], distance, out_of_plane_angle_unitless, in_plane_angle)

        # Try passing in a unitless in_plane_angle
        with pytest.raises(AssertionError) as excinfo:
            vsite1_index = molecule.add_monovalent_lone_pair_virtual_site([atom1, atom2], distance, out_of_plane_angle, in_plane_angle_unitless)

        # Try giving two atoms
        with pytest.raises(AssertionError) as excinfo:
            vsite1_index = molecule.add_monovalent_lone_pair_virtual_site([atom1, atom2], distance, out_of_plane_angle, in_plane_angle)

        # Successfully make a virtual site
        vsite1_index = molecule.add_monovalent_lone_pair_virtual_site([atom1, atom2, atom3], distance, out_of_plane_angle, in_plane_angle)
        # TODO: Check if we get the same values back out from the @properties
        molecule_dict = molecule.to_dict()
        molecule2 = Molecule.from_dict(molecule_dict)
        assert molecule.to_dict() == molecule2.to_dict()

    @pytest.mark.parametrize('molecule', mini_drug_bank())
    def test_add_divalent_lone_pair_virtual_site(self, molecule):
        """Test addition of a DivalentLonePairVirtualSite to the Molecule"""
        # Do not modify the original molecule.
        molecule = copy.deepcopy(molecule)

        atom1 = molecule.atoms[0]
        atom2 = molecule.atoms[1]
        atom3 = molecule.atoms[2]
        atom4 = molecule.atoms[3]
        distance = 0.3 * unit.angstrom
        out_of_plane_angle = 30 * unit.degree
        in_plane_angle = 0.2 * unit.radian
        vsite1_index = molecule.add_divalent_lone_pair_virtual_site([atom1, atom2, atom3], distance, out_of_plane_angle, in_plane_angle)
        with pytest.raises(AssertionError) as excinfo:
            vsite1_index = molecule.add_divalent_lone_pair_virtual_site([atom1, atom2], distance, out_of_plane_angle, in_plane_angle)
        molecule_dict = molecule.to_dict()
        molecule2 = Molecule.from_dict(molecule_dict)
        assert molecule_dict == molecule2.to_dict()

    @pytest.mark.parametrize('molecule', mini_drug_bank())
    def test_add_trivalent_lone_pair_virtual_site(self, molecule):
        """Test addition of a TrivalentLonePairVirtualSite to the Molecule"""
        # Do not modify the original molecule.
        molecule = copy.deepcopy(molecule)

        atom1 = molecule.atoms[0]
        atom2 = molecule.atoms[1]
        atom3 = molecule.atoms[2]
        atom4 = molecule.atoms[3]
        distance = 0.3 * unit.angstrom
        out_of_plane_angle = 30 * unit.degree
        in_plane_angle = 0.2 * unit.radian
        vsite1_index = molecule.add_trivalent_lone_pair_virtual_site([atom1, atom2, atom3, atom4], distance, out_of_plane_angle, in_plane_angle)
        # Test for assertion when giving too few atoms
        with pytest.raises(AssertionError) as excinfo:
            vsite1_index = molecule.add_trivalent_lone_pair_virtual_site([atom1, atom2, atom3], distance, out_of_plane_angle, in_plane_angle)
        molecule_dict = molecule.to_dict()
        molecule2 = Molecule.from_dict(molecule_dict)
        assert molecule.to_dict() == molecule2.to_dict()

    @requires_openeye
    def test_chemical_environment_matches_OE(self):
        """Test chemical environment matches"""
        # TODO: Move this to test_toolkits, test all available toolkits
        # Create chiral molecule
        from simtk.openmm.app import element
        toolkit_wrapper = OpenEyeToolkitWrapper()
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
        matches = molecule.chemical_environment_matches('[#6:1]', toolkit_registry=toolkit_wrapper)
        assert len(matches) == 1 # there should be a unique match, so one atom tuple is returned
        assert len(matches[0]) == 1 # it should have one tagged atom
        assert set(matches[0]) == set([atom_C])
        matches = molecule.chemical_environment_matches('[#6:1]~[#1:2]', toolkit_registry=toolkit_wrapper)
        assert len(matches) == 1 # there should be a unique match, so one atom tuple is returned
        assert len(matches[0]) == 2 # it should have two tagged atoms
        assert set(matches[0]) == set([atom_C, atom_H])
        matches = molecule.chemical_environment_matches('[Cl:1]-[C:2]-[H:3]', toolkit_registry=toolkit_wrapper)
        assert len(matches) == 1 # there should be a unique match, so one atom tuple is returned
        assert len(matches[0]) == 3 # it should have three tagged atoms
        assert set(matches[0]) == set([atom_Cl, atom_C, atom_H])
        matches = molecule.chemical_environment_matches('[#6:1]~[*:2]', toolkit_registry=toolkit_wrapper)
        assert len(matches) == 4 # there should be four matches
        for match in matches:
            assert len(match) == 2 # each match should have two tagged atoms

    # TODO: Test forgive undef amide enol stereo
    # TODO: test forgive undef phospho linker stereo
    # TODO: test forgive undef C=NH stereo
    # TODO: test forgive undef phospho stereo
    # Potentially better OE stereo check: OEFlipper — Toolkits - - Python
    # https: // docs.eyesopen.com / toolkits / python / omegatk / OEConfGenFunctions / OEFlipper.html

    @requires_rdkit
    def test_chemical_environment_matches_RDKit(self):
        """Test chemical environment matches"""
        # Create chiral molecule
        from simtk.openmm.app import element
        toolkit_wrapper = RDKitToolkitWrapper()
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
        matches = molecule.chemical_environment_matches('[#6:1]', toolkit_registry=toolkit_wrapper)
        assert len(matches) == 1 # there should be a unique match, so one atom tuple is returned
        assert len(matches[0]) == 1 # it should have one tagged atom
        assert set(matches[0]) == set([atom_C])
        matches = molecule.chemical_environment_matches('[#6:1]~[#1:2]', toolkit_registry=toolkit_wrapper)
        assert len(matches) == 1 # there should be a unique match, so one atom tuple is returned
        assert len(matches[0]) == 2 # it should have two tagged atoms
        assert set(matches[0]) == set([atom_C, atom_H])
        matches = molecule.chemical_environment_matches('[Cl:1]-[C:2]-[H:3]', toolkit_registry=toolkit_wrapper)
        assert len(matches) == 1 # there should be a unique match, so one atom tuple is returned
        assert len(matches[0]) == 3 # it should have three tagged atoms
        assert set(matches[0]) == set([atom_Cl, atom_C, atom_H])
        matches = molecule.chemical_environment_matches('[#6:1]~[*:2]', toolkit_registry=toolkit_wrapper)
        assert len(matches) == 4 # there should be four matches
        for match in matches:
            assert len(match) == 2 # each match should have two tagged atoms

    @pytest.mark.slow
    def test_compute_partial_charges(self):
        """Test computation/retrieval of partial charges"""
        # TODO: Test only one molecule for speed?
        # TODO: Do we need to deepcopy each molecule, or is setUp called separately for each test method?
        from simtk import unit
        import numpy as np

        # Do not modify original molecules.
        molecules = copy.deepcopy(mini_drug_bank())

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
                for molecule in molecules[:1]: # Just test first molecule to save time
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

    @requires_openeye
    def test_assign_fractional_bond_orders(self):
        """Test assignment of fractional bond orders
        """
        # TODO: Test only one molecule for speed?
        # TODO: Do we need to deepcopy each molecule, or is setUp called separately for each test method?

        # Do not modify the original molecules.
        molecules = copy.deepcopy(mini_drug_bank())

        toolkit_to_bondorder_method = {OpenEyeToolkitWrapper:['am1','pm3']}
        for toolkit in list(toolkit_to_bondorder_method.keys()):
            toolkit_registry = ToolkitRegistry(toolkit_precedence=[toolkit])
            for charge_model in toolkit_to_bondorder_method[toolkit]:
                for molecule in molecules[:5]: # Just test first five molecules for speed
                    molecule.compute_wiberg_bond_orders(charge_model=charge_model, toolkit_registry=toolkit_registry)
                    fbo1 = [bond.fractional_bond_order for bond in molecule.bonds]
                    # Call should be faster the second time due to caching
                    molecule.compute_wiberg_bond_orders(charge_model=charge_model, toolkit_registry=toolkit_registry)
                    fbo2 = [bond.fractional_bond_order for bond in molecule.bonds]
                    assert fbo1 == fbo2
