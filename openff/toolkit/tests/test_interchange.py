from openff.toolkit.topology import Molecule
from openff.toolkit.typing.engines.smirnoff.forcefield import ForceField


def test_basic_construction():
    top = Molecule.from_smiles("C").to_topology()
    parsley = ForceField("openff-1.0.0.offxml")
    out = parsley._create_openff_interchange(top)

    assert out.topology is not None
    assert out.topology.mdtop is not None
    assert out.positions is None
    for key in ["vdW", "Bonds", "ProperTorsions", "Electrostatics"]:
        assert key in out.handlers
