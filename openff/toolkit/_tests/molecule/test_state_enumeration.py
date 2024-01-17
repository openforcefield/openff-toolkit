import pytest

from openff.toolkit import Molecule
from openff.toolkit._tests.utils import requires_openeye


class TestProtomerEnumeration:
    @requires_openeye
    def test_enumerating_no_protomers(self):
        """Make sure the input molecule is returned when there is only one protomers."""

        mol = Molecule.from_smiles("CC")

        assert len(mol.enumerate_protomers()) == 1
        assert mol.enumerate_protomers()[0] == mol

    @requires_openeye
    def test_enumerating_protomers(self):
        """Test enumerating the formal charges."""

        mol = Molecule.from_smiles("Oc2ccc(c1ccncc1)cc2")

        # there should be three protomers for this molecule so restrict the output
        protomers = mol.enumerate_protomers(max_states=2)

        assert mol in protomers
        assert len(protomers) == 3

        # now make sure we can generate them all
        protomers = mol.enumerate_protomers(max_states=10)

        assert mol in protomers
        assert len(protomers) == 4

        # make sure each protomer is unique
        unique_protomers = set(protomers)
        assert len(protomers) == len(unique_protomers)

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
