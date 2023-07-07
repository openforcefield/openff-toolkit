from openff.interchange._tests import get_test_file_path
from openff.interchange.interop.gromacs._import._import import from_files
from openff.interchange.interop.gromacs.models.models import RyckaertBellemansDihedral


def test_load_rb_torsions(monkeypatch):
    monkeypatch.setenv("INTERCHANGE_EXPERIMENTAL", "1")

    system = from_files(
        get_test_file_path("ethanol_rb_torsions.top"),
        get_test_file_path("ethanol_rb_torsions.gro"),
    )

    assert len(system.molecule_types["Compound"].dihedrals) > 0

    for torsion in system.molecule_types["Compound"].dihedrals:
        assert isinstance(torsion, RyckaertBellemansDihedral)
