import pytest

from openff.toolkit import Molecule
from openff.toolkit._tests.utils import requires_openeye
from openff.toolkit.utils.toolkits import OpenEyeToolkitWrapper, RDKitToolkitWrapper


class TestMoleculeStateEnumeration:
    tautomer_data = [
        {"molecule": "Oc1c(cccc3)c3nc2ccncc12", "tautomers": 2},
        {"molecule": "CN=c1nc[nH]cc1", "tautomers": 2},
        {"molecule": "c1[nH]c2c(=O)[nH]c(nc2n1)N", "tautomers": 14},
    ]

    @pytest.mark.parametrize(
        "toolkit_class", [OpenEyeToolkitWrapper, RDKitToolkitWrapper]
    )
    @pytest.mark.parametrize("molecule_data", tautomer_data)
    def test_enumerating_tautomers(self, molecule_data, toolkit_class):
        """Test the ability of each toolkit to produce tautomers of an input molecule."""

        if toolkit_class.is_available():
            toolkit = toolkit_class()
            mol = Molecule.from_smiles(
                molecule_data["molecule"],
                allow_undefined_stereo=True,
                toolkit_registry=toolkit,
            )

            tautomers = mol.enumerate_tautomers(toolkit_registry=toolkit)

            assert len(tautomers) == molecule_data["tautomers"]
            assert mol not in tautomers
            # check that the molecules are not isomorphic of the input
            for taut in tautomers:
                assert taut.n_conformers == 0
                assert mol.is_isomorphic_with(taut) is False

        else:
            pytest.skip("Required toolkit is unavailable")

    @pytest.mark.parametrize(
        "toolkit_class", [OpenEyeToolkitWrapper, RDKitToolkitWrapper]
    )
    def test_enumerating_tautomers_options(self, toolkit_class):
        """Test the enumeration options"""

        if toolkit_class.is_available():
            toolkit = toolkit_class()
            # test the max molecules option
            mol = Molecule.from_smiles(
                "c1[nH]c2c(=O)[nH]c(nc2n1)N",
                toolkit_registry=toolkit,
                allow_undefined_stereo=True,
            )

            tauts_no = 5
            tautomers = mol.enumerate_tautomers(
                max_states=tauts_no, toolkit_registry=toolkit
            )
            assert len(tautomers) <= tauts_no
            assert mol not in tautomers

    @pytest.mark.parametrize(
        "toolkit_class", [RDKitToolkitWrapper, OpenEyeToolkitWrapper]
    )
    def test_enumerating_no_tautomers(self, toolkit_class):
        """Test that the toolkits return an empty list if there are no tautomers to enumerate."""

        if toolkit_class.is_available():
            toolkit = toolkit_class()
            mol = Molecule.from_smiles("CC", toolkit_registry=toolkit)

            tautomers = mol.enumerate_tautomers(toolkit_registry=toolkit)
            assert tautomers == []

        else:
            pytest.skip("Required toolkit is unavailable")

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

    @pytest.mark.parametrize(
        "toolkit_class", [OpenEyeToolkitWrapper, RDKitToolkitWrapper]
    )
    def test_enumerating_stereobonds(self, toolkit_class):
        """Test the backend toolkits in enumerating the stereo bonds in a molecule."""

        if toolkit_class.is_available():
            toolkit = toolkit_class()
            mol = Molecule.from_smiles(
                "ClC=CCl", allow_undefined_stereo=True, toolkit_registry=toolkit
            )

            # use the default options
            isomers = mol.enumerate_stereoisomers()
            assert len(isomers) == 2

            assert mol not in isomers
            # make sure the input molecule is only different by bond stereo
            for ismol in isomers:
                assert (
                    Molecule.are_isomorphic(
                        mol,
                        ismol,
                        return_atom_map=False,
                        bond_stereochemistry_matching=False,
                    )[0]
                    is True
                )
                assert mol.is_isomorphic_with(ismol) is False

            # make sure the isomers are different
            assert isomers[0].is_isomorphic_with(isomers[1]) is False

        else:
            pytest.skip("Required toolkit is unavailable")

    @pytest.mark.parametrize(
        "toolkit_class", [OpenEyeToolkitWrapper, RDKitToolkitWrapper]
    )
    def test_enumerating_stereocenters(self, toolkit_class):
        """Test the backend toolkits in enumerating the stereo centers in a molecule."""

        if toolkit_class.is_available():
            toolkit = toolkit_class()
            mol = Molecule.from_smiles(
                "NC(Cl)(F)O", toolkit_registry=toolkit, allow_undefined_stereo=True
            )

            isomers = mol.enumerate_stereoisomers(toolkit_registry=toolkit)

            assert len(isomers) == 2
            # make sure the mol is not in the isomers and that they only differ by stereo chem
            assert mol not in isomers
            for ismol in isomers:
                assert ismol.n_conformers != 0
                assert (
                    Molecule.are_isomorphic(
                        mol,
                        ismol,
                        return_atom_map=False,
                        atom_stereochemistry_matching=False,
                    )[0]
                    is True
                )
                assert mol.is_isomorphic_with(ismol) is False

            # make sure the two isomers are different
            assert isomers[0].is_isomorphic_with(isomers[1]) is False

        else:
            pytest.skip("Required toolkit is unavailable")

    @pytest.mark.parametrize(
        "toolkit_class", [OpenEyeToolkitWrapper, RDKitToolkitWrapper]
    )
    def test_enumerating_stereo_options(self, toolkit_class):
        """Test the enumerating stereo chem options"""

        if toolkit_class.is_available():
            toolkit = toolkit_class()

            # test undefined only
            mol = Molecule.from_smiles(
                "ClC=CCl", toolkit_registry=toolkit, allow_undefined_stereo=True
            )
            isomers = mol.enumerate_stereoisomers(
                undefined_only=True, rationalise=False
            )

            assert len(isomers) == 2
            for isomer in isomers:
                assert isomer.n_conformers == 0

            mol = Molecule.from_smiles(
                r"Cl/C=C\Cl", toolkit_registry=toolkit, allow_undefined_stereo=True
            )
            isomers = mol.enumerate_stereoisomers(
                undefined_only=True, rationalise=False
            )

            assert isomers == []

            mol = Molecule.from_smiles(
                r"Cl/C=C\Cl", toolkit_registry=toolkit, allow_undefined_stereo=True
            )
            isomers = mol.enumerate_stereoisomers(
                undefined_only=False, rationalise=False
            )

            assert len(isomers) == 1

            # test max isomers
            mol = Molecule.from_smiles(
                "BrC=C[C@H]1OC(C2)(F)C2(Cl)C1",
                toolkit_registry=toolkit,
                allow_undefined_stereo=True,
            )
            isomers = mol.enumerate_stereoisomers(
                max_isomers=5,
                undefined_only=True,
                toolkit_registry=toolkit,
                rationalise=True,
            )

            assert len(isomers) <= 5
            for isomer in isomers:
                assert isomer.n_conformers == 1

        else:
            pytest.skip("Required toolkit is unavailable")

    @pytest.mark.parametrize(
        "toolkit_class", [OpenEyeToolkitWrapper, RDKitToolkitWrapper]
    )
    @pytest.mark.parametrize(
        "smiles, undefined_only, expected",
        [
            (
                "FC(Br)(Cl)[C@@H](Br)(Cl)",
                False,
                [
                    "F[C@](Br)(Cl)[C@@H](Br)(Cl)",
                    "F[C@](Br)(Cl)[C@H](Br)(Cl)",
                    "F[C@@](Br)(Cl)[C@@H](Br)(Cl)",
                    "F[C@@](Br)(Cl)[C@H](Br)(Cl)",
                ],
            ),
            (
                "FC(Br)(Cl)[C@@H](Br)(Cl)",
                True,
                ["F[C@](Br)(Cl)[C@@H](Br)(Cl)", "F[C@@](Br)(Cl)[C@@H](Br)(Cl)"],
            ),
            ("F[C@H](Cl)Br", False, ["F[C@@H](Cl)Br"]),
            ("F[C@H](Cl)Br", True, []),
        ],
    )
    def test_enumerating_stereo_partially_defined(
        self, toolkit_class, smiles, undefined_only, expected
    ):
        """Test the enumerating stereo of molecules with partially defined chirality"""

        if not toolkit_class.is_available():
            pytest.skip("Required toolkit is unavailable")

        toolkit = toolkit_class()

        # test undefined only
        mol = Molecule.from_smiles(
            smiles, toolkit_registry=toolkit, allow_undefined_stereo=True
        )
        stereoisomers = mol.enumerate_stereoisomers(
            undefined_only=undefined_only, rationalise=False
        )

        # Ensure that the results of the enumeration are what the test expects.
        # This roundtrips the expected output from SMILES --> OFFMol --> SMILES,
        # since the SMILES for stereoisomers generated in this test may change depending
        # on which cheminformatics toolkit is used.
        expected = {
            Molecule.from_smiles(stereoisomer, allow_undefined_stereo=True).to_smiles(
                explicit_hydrogens=True, isomeric=True, mapped=False
            )
            for stereoisomer in expected
        }
        actual = {
            stereoisomer.to_smiles(explicit_hydrogens=True, isomeric=True, mapped=False)
            for stereoisomer in stereoisomers
        }

        assert expected == actual
