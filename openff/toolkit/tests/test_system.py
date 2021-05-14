from openff.system.stubs import ForceField

from openff.toolkit.topology import Molecule, Topology


def test_basic_construction():
    top = Molecule.from_smiles("C").to_topology()
    parsley = ForceField("openff-1.0.0.offxml")
    off_sys = parsley.create_openff_system(top)

    assert off_sys.topology is not None
    assert off_sys.topology.mdtop is not None
    assert off_sys.positions is None
    for key in ["vdW", "Bonds", "Electrostatics"]:
        assert key in off_sys.handlers
