from typing import List

import pytest

from openff.toolkit.typing.chemistry import (
    AngleChemicalEnvironment,
    AtomChemicalEnvironment,
    BondChemicalEnvironment,
    ChemicalEnvironment,
    ImproperChemicalEnvironment,
    TorsionChemicalEnvironment,
)
from openff.toolkit.utils.exceptions import SMIRKSMismatchError, SMIRKSParsingError
from openff.toolkit.utils.toolkits import OPENEYE_AVAILABLE

# TODO: Evaluate which tests in this file should be moved to test_toolkits
toolkits: List = []
if OPENEYE_AVAILABLE:
    from openff.toolkit.utils.toolkits import OpenEyeToolkitWrapper, RDKitToolkitWrapper

    toolkits.append("openeye")
    toolkits.append(OpenEyeToolkitWrapper())
else:
    from openff.toolkit.utils.toolkits import RDKitToolkitWrapper

    toolkits.append("rdkit")
    toolkits.append(RDKitToolkitWrapper())


class TestChemicalEnvironments:
    def test_deprecation_warning(self):
        from openff.toolkit.typing.chemistry.environment import (
            ChemicalEnvironmentDeprecationWarning,
        )

        with pytest.warns(ChemicalEnvironmentDeprecationWarning):
            ChemicalEnvironment("[*:1]")

    def test_createEnvironments(self):
        """
        Test all types of ChemicalEnvironment objects with defined atoms and bonds
        Each will be tetrahedral carbons connected by ring single bonds
        """
        AtomChemicalEnvironment("[#6X4:1]", "CT")
        BondChemicalEnvironment("[#6X4:1]-[#6X4:2]", "CT-CT")
        AngleChemicalEnvironment("[#6X4:1]-[#6X4:2]-[#6X4:3]", "CT-CT-CT")
        TorsionChemicalEnvironment("[#6X4:1]-[#6X4:2]-[#6X4:3]-[#6X4:4]", "CT-CT-CT-CT")
        ImproperChemicalEnvironment(
            "[#6X4:1]-[#6X4:2](-[#6X4:3])-[#6X4:4]", "CT-CT(-CT)-CT"
        )

    @pytest.mark.parametrize(
        ["smirks", "expected_valence", "expected_chemenv_class"],
        [
            ["[#6](-[#1])-[#8]", None, ChemicalEnvironment],
            ["[#6&X4&H0:1](-[#1])-[#6&X4]", "Atom", AtomChemicalEnvironment],
            ["[#6&X4&H0:1](-[#1])-[#6&X4:2]", "Bond", BondChemicalEnvironment],
            ["[*:1]-[*:2](-[#6&X4])-[*:3]", "Angle", AngleChemicalEnvironment],
            [
                "[#6&X4&H0:1](-[#1])-[#6&X4:2]-[#6&X4&H0:3](-[#1])-[#6&X4:4]",
                "ProperTorsion",
                TorsionChemicalEnvironment,
            ],
            [
                "[#1:1]-[#6&X4:2](-[#8:3])-[#1:4]",
                "ImproperTorsion",
                ImproperChemicalEnvironment,
            ],
            # Test that an improper smirks is also valid as a general ChemicalEnvironment
            [
                "[#1:1]-[#6&X4:2](-[#8:3])-[*:4](-[#6&H1])-[#8:5]",
                None,
                ChemicalEnvironment,
            ],
            ["[#6$(*~[#6]=[#8])$(*-,=[#7!-1,#8,#16,#7])]", None, ChemicalEnvironment],
            ["CCC", None, ChemicalEnvironment],
            ["[#6:1]1(-;!@[#1,#6])=;@[#6]-;@[#6]1", "Atom", ChemicalEnvironment],
            ["C(O-[#7,#8])CC=[*]", None, ChemicalEnvironment],
            [
                "[#6$([#6X4](~[#7!-1,#8!-1,#16!-1,#9,#17,#35,#53])(~[#8]~[#1])):1]-"
                "[#6X2H2;+0:2]-,=,:;!@;!#[#7!-1,#8,#16:3]-[#4:4]",
                "ProperTorsion",
                TorsionChemicalEnvironment,
            ],
            [
                "[#6$([#6X4](~[#7!-1,#8!-1,#16!-1,#9,#17,#35,#53])(~[#8]~[#1])):1]1=CCCC1",
                "Atom",
                AtomChemicalEnvironment,
            ],
            [
                "[*:1]-[#7X3:2](-[#6a$(*1ccc(-[#8-1X1])cc1):3])-[*:4]",
                "ImproperTorsion",
                ImproperChemicalEnvironment,
            ],
            ["[#6X4:1]1~[*:2]~[*$(*~[#1]):3]1", "Angle", AngleChemicalEnvironment],
            ["[$([#7]1~[#6]-CC1)]", None, ChemicalEnvironment],
            ["[$(c1ccccc1)]", None, ChemicalEnvironment],
            # The next two tests are for ring-closing bonds
            [
                "[H][C@:4]1(C(C([C:3]([N:2]1[C:1](=O)C([H])([H])[H])([H])[H])([H])[H])([H])[H])C=O",
                "ImproperTorsion",
                ChemicalEnvironment,
            ],
            ["[P:1]=1=[P]=[P]=[P]=[P:2]=1", "Bond", BondChemicalEnvironment],
        ],
    )
    @pytest.mark.parametrize("toolkit", toolkits)
    def test_parseSMIRKS(
        self, toolkit, smirks, expected_valence, expected_chemenv_class
    ):
        """
        Test creating environments with SMIRKS
        """
        env = expected_chemenv_class(smirks=smirks, toolkit_registry=toolkit)
        actual_type = env.get_type()
        assert (
            actual_type == expected_valence
        ), f"SMIRKS ({smirks}) classified as {actual_type} instead of {expected_valence} using {toolkit} toolkit"

    @pytest.mark.parametrize(
        ("smirks", "wrong_envs"),
        [
            (
                "[*]",
                [
                    AtomChemicalEnvironment,
                    BondChemicalEnvironment,
                    AngleChemicalEnvironment,
                    TorsionChemicalEnvironment,
                    ImproperChemicalEnvironment,
                ],
            ),
            (
                "[*:1]",
                [
                    BondChemicalEnvironment,
                    AngleChemicalEnvironment,
                    TorsionChemicalEnvironment,
                    ImproperChemicalEnvironment,
                ],
            ),
            (
                "[*:1]~[*:2]",
                [
                    AtomChemicalEnvironment,
                    AngleChemicalEnvironment,
                    TorsionChemicalEnvironment,
                    ImproperChemicalEnvironment,
                ],
            ),
            (
                "[*:3]~[*:2]~[*:1]",
                [
                    AtomChemicalEnvironment,
                    BondChemicalEnvironment,
                    TorsionChemicalEnvironment,
                    ImproperChemicalEnvironment,
                ],
            ),
            (
                "[*:1]~[*:2]~[*:3]~[*:4]",
                [
                    AtomChemicalEnvironment,
                    BondChemicalEnvironment,
                    AngleChemicalEnvironment,
                    ImproperChemicalEnvironment,
                ],
            ),
            (
                "[*:1]~[*:2](~[*:3])~[*:4]",
                [
                    AtomChemicalEnvironment,
                    BondChemicalEnvironment,
                    AngleChemicalEnvironment,
                    TorsionChemicalEnvironment,
                ],
            ),
            (
                "[*:1]~[*:2]~[*:3]~[*:4]~[*:5]",
                [
                    AtomChemicalEnvironment,
                    BondChemicalEnvironment,
                    AngleChemicalEnvironment,
                    TorsionChemicalEnvironment,
                    ImproperChemicalEnvironment,
                ],
            ),
        ],
    )
    def test_creating_wrong_environments(self, smirks, wrong_envs):
        """
        Test exceptions for making environments with the wrong smirks
        """
        for wrong_env in wrong_envs:
            with pytest.raises(SMIRKSMismatchError):
                wrong_env(smirks)

    @pytest.mark.parametrize("toolkit", toolkits)
    def test_wrong_smirks_error(self, toolkit):
        """
        Check that an imparseable SMIRKS raises errors
        """
        smirks = "[*;:1]"
        with pytest.raises(SMIRKSParsingError):
            ChemicalEnvironment(smirks, toolkit_registry=toolkit)

    def test_embedded_atoms_smirks(self):
        """
        Check embedded atom parsing works
        """
        smirks = "[#1$(*-[#6](-[#7,#8,#9,#16,#17,#35])-[#7,#8,#9,#16,#17,#35]):1]~[$([#1]~[#6])]"
        ChemicalEnvironment(smirks)
