import pytest
from openff.units import unit

from openff.interchange import Interchange
from openff.interchange._tests import needs_gmx
from openff.interchange.drivers.gromacs import get_gromacs_energies
from openff.interchange.drivers.openmm import get_openmm_energies


@pytest.mark.slow()
@needs_gmx
def test_group_impropers(cb8_host, no_charges):
    out = Interchange.from_smirnoff(no_charges, [cb8_host], box=[4, 4, 4])

    openmm_torsions = get_openmm_energies(out)["Torsion"]
    gromacs_torsions = get_gromacs_energies(out)["Torsion"]

    assert abs(openmm_torsions - gromacs_torsions).m_as(unit.kilojoule_per_mole) < 1e-3
