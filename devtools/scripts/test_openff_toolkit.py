from packaging.version import Version

from openff.toolkit import __version__
from openff.toolkit.topology import Molecule
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.utils.toolkits import (
    GLOBAL_TOOLKIT_REGISTRY,
    NAGLToolkitWrapper,
    RDKitToolkitWrapper,
)


def test_versions_sane():
    # partially safeguard against a malformed version number;
    # 0.minor.patch and major.0.patch should be okay but not 0.0.patch
    found_version = Version(__version__)
    assert not (found_version.minor == 0 and found_version.major == 0), found_version

    print(f"found version {found_version}")

    print(f"{GLOBAL_TOOLKIT_REGISTRY.registered_toolkit_versions=}")

    assert RDKitToolkitWrapper().is_available()
    assert NAGLToolkitWrapper().is_available()


def test_basic_system_creation():
    offmol = Molecule.from_smiles("CCO")
    Molecule.from_rdkit(offmol.to_rdkit())
    ff = ForceField("openff-2.3.0.offxml")
    ff.create_openmm_system(offmol.to_topology())
