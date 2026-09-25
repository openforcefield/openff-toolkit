import pytest
from rdkit import Chem

from openff.toolkit import Molecule
from openff.toolkit.utils.toolkits import OpenEyeToolkitWrapper

ATOM_MOLBLOCK = """
     RDKit          3D

 12 12  0  0  0  0  0  0  0  0999 V2000
    1.4501   -0.0070    0.2136 C   0  0  2  0  0  0  0  0  0  0  0  0
    0.4490   -1.1287    0.0101 C   0  0  2  0  0  0  0  0  0  0  0  0
   -0.2954    0.0069   -0.6415 C   0  0  2  0  0  0  0  0  0  0  0  0
    0.5002    1.1328   -0.0918 C   0  0  1  0  0  0  0  0  0  0  0  0
    1.9656    0.2213    1.1500 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.8561   -1.9011   -0.5958 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1142    0.1098   -1.8269 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.9522    1.9015   -0.8006 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.5907   -0.2572   -0.7380 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4653   -1.5490    1.1650 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.0448    1.6586    1.0771 F   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1067    0.0572   -0.4945 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  1  5  1  1
  1  9  1  0
  1  4  1  0
  1  2  1  0
  2  6  1  6
  2  3  1  0
  2 10  1  0
  3  7  1  6
  3  4  1  0
  3 12  1  0
  4  8  1  6
  4 11  1  0
M  END
"""

BOND_MOLBLOCK = """
     RDKit          3D

 13 13  0  0  0  0  0  0  0  0999 V2000
    1.3530    0.3119    0.6181 C   0  0  1  0  0  0  0  0  0  0  0  0
    0.7075    1.4537   -0.4274 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5398    1.5231   -0.3512 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1576    0.4603    0.4717 C   0  0  2  0  0  0  0  0  0  0  0  0
   -0.0095   -0.4814    0.2775 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6018    0.4601    1.6545 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.3024    2.2655   -1.0164 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0855    2.3208   -0.8148 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3336    0.7556    1.5274 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3591   -1.3768   -0.6865 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.2753   -0.7074   -0.0612 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3588   -0.4836    0.0369 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2645   -1.6240   -0.4028 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  6  1  1
  1 11  1  0
  1  5  1  0
  1  2  1  0
  2  7  1  0
  2  3  2  0
  3  8  1  0
  3  4  1  0
  4  9  1  1
  4  5  1  0
  4 12  1  0
  5 13  2  0
 10 13  1  0
M  CHG  3  11  -1  12  -1  13   1
M  END
"""


def molecule(molblock):
    rdmol = Chem.MolFromMolBlock(molblock, removeHs=False)
    return Molecule.from_rdkit(rdmol, allow_undefined_stereo=True)


def labels(mol):
    return (
        {a.molecule_atom_index: a.stereochemistry for a in mol.atoms if a.stereochemistry},
        {tuple(sorted((b.atom1_index, b.atom2_index))): b.stereochemistry for b in mol.bonds if b.stereochemistry},
    )


@pytest.mark.parametrize("molblock", [ATOM_MOLBLOCK, BOND_MOLBLOCK], ids=["atom_0237e4a4", "bond_b69936ff"])
class TestMinimalRecords:
    def test_input_has_stereochemistry(self, molblock):
        atoms, bonds = labels(molecule(molblock))
        assert atoms or bonds

    def test_to_openeye_does_not_raise_inconsistent_stereochemistry(self, molblock):
        OpenEyeToolkitWrapper().to_openeye(molecule(molblock))

    def test_stereochemistry_is_not_silently_lost(self, molblock):
        mol = molecule(molblock)
        back = OpenEyeToolkitWrapper().from_openeye(
            OpenEyeToolkitWrapper().to_openeye(mol), allow_undefined_stereo=True
        )
        assert labels(back) == labels(mol)
