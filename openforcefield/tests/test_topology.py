from functools import partial
from unittest import TestCase
import parmed
from openforcefield import utils
from openforcefield import topology
import pickle

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

class TestMolecule(TestCase):
    from openforcefield.topology import Molecule

    def setUp(self):
        self.molecules = pickle.load('zinc-subset-offmols.pkl')

    def test_serialize(self):
        serialized = pickle.dumps(self.molecules)
        molecules_copy = pickle.loads(serialized)

    # TODO: Only run tests if the RDKit is installed
    def test_rdkit_roundtrip(self):
        for molecule in self.molecules:
            rdmol = molecule.to_rdkit()
            molecule2 = Molecule.from_rdmol(rdmol)
            assert_molecule_is_equal(molecule, molecule2, "Molecule.to_rdmol()/from_rdmol() round trip failed")
            molecule3 = Molecule(rdmol)
            assert_molecule_is_equal(molecule, molecule3, "Molecule(rdmol) constructor failed")

    # TODO: Only run tests if the OpenEye toolkit is installed
    def test_oemol_roundtrip(self):
        """Test creation of Molecule object from OpenEye OEMol
        """
        for molecule in self.molecules:
            oemol = molecule.to_openeye()
            molecule2 = Molecule.from_openeye(oemol)
            assert_molecule_is_equal(molecule, molecule2, "Molecule.to_openeye()/from_openeye() round trip failed")
            molecule3 = Molecule(oemol)
            assert_molecule_is_equal(molecule, molecule3, "Molecule(oemol) constructor failed")

    def test_assign_partial_charges(self):
        """Test assignment of partial charges
        """
        # TODO: Only the OpenEye toolkit currently supports charge models
        # TODO: Does RDKit support other charge models?
        for molecule in self.molecules:
            for charge_model in topology.ALLOWED_CHARGE_MODELS:
                molecule.assign_partial_charges(method=charge_model)

    def test_assign_fractional_bond_orders(self):
        """Test assignment of fractional bond orders
        """
        for molecule in self.molecules:
            for charge_model in topology.ALLOWED_FRACTIONAL_BONDORDER_MODELS:
                molecule.assign_fractional_bond_orders(method=charge_model)

    def test_bonds(self):
        """Test iteration over bonds
        """
        for molecule in self.molecules:
            bonds = molecule.bonds()
            # TODO: Check known cases

    def test_angles(self):
        """Test iteration over angles
        """
        for molecule in self.molecules:
            angles = molecule.angles()
            # TODO: Check known cases

    def test_torsions(self):
        """Test iteration over torsions
        """
        for molecule in self.molecules:
            torsion = molecule.torsions()
            # TODO: Check known cases

    def test_propers(self):
        """Test iteration over proper torsions
        """
        for molecule in self.molecules:
            torsion = molecule.propers()
            # TODO: Check known cases

    def test_impropers(self):
        """Test iteration over improper torsions
        """
        for molecule in self.molecules:
            torsion = molecule.impropers()
            # TODO: Check known cases
