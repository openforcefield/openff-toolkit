import pytest

from openff.toolkit import Molecule
from openff.toolkit._tests.utils import requires_openeye


class TestProtomerEnumeration:
    @requires_openeye
    def test_enumerating_no_protomers(self):
        """Make sure the input molecule is returned when there is only one protomer."""

        mol = Molecule.from_smiles("CC")

        assert len(mol.enumerate_protomers()) == 1
        assert mol.enumerate_protomers()[0] == mol

    @requires_openeye
    def test_enumerating_protomers(self):
        """Test enumerating the formal charges."""

        # there should be three protomers, in addition to the input state
        mol = Molecule.from_smiles("Oc2ccc(c1ccncc1)cc2")

        protomers = mol.enumerate_protomers()

        assert mol in protomers
        assert len(protomers) == 4

        # make sure generating extra states produces the same result
        assert len(protomers) == len(mol.enumerate_protomers(max_states=10))

        # make sure each protomer is unique
        assert len(protomers) == len(set(protomers))

    @requires_openeye
    def test_tetracarboxylic_acid(self):
        acid = Molecule.from_smiles("C(C(=O)O)(C(=O)O)=C(C(=O)O)(C(=O)O)")

        protomers = acid.enumerate_protomers()

        assert len(protomers) == 16
        assert acid in protomers

    @pytest.mark.slow
    @requires_openeye
    def test_many_protomers(self):
        acid = Molecule.from_smiles(5 * "C(C(=O)O)(C(=O)O)")

        protomers = acid.enumerate_protomers()

        assert acid in protomers
        assert len(protomers) == 6400
