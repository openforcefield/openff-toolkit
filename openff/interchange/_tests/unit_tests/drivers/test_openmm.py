"""\
def _process(
    raw_energies: Dict[int, openmm.Force],
    system,
    combine_nonbonded_forces: bool,
    detailed: bool,
) -> EnergyReport:
"""

import openmm
import pytest

from openff.interchange.constants import kj_mol
from openff.interchange.drivers.openmm import _process


class TestProcess:
    @pytest.fixture()
    def dummy_system(self):
        system = openmm.System()
        for force in [
            openmm.PeriodicTorsionForce,
            openmm.HarmonicAngleForce,
            openmm.HarmonicBondForce,
            openmm.RBTorsionForce,
            openmm.NonbondedForce,
        ]:
            system.addForce(force())

        return system

    @pytest.fixture()
    def dummy_system_split(self, dummy_system):
        dummy_system.addForce(openmm.CustomNonbondedForce("sigma*epsilon"))
        dummy_system.addForce(openmm.CustomBondForce("sigma*epsilon"))
        dummy_system.addForce(openmm.CustomBondForce("qq"))

        return dummy_system

    def test_simple(self, dummy_system):
        processed = _process(
            {
                0: 1.1 * kj_mol,
                1: 2.2 * kj_mol,
                2: 3.3 * kj_mol,
                3: 4.4 * kj_mol,
                4: 5.5 * kj_mol,
            },
            dummy_system,
            True,
            False,
        )

        assert processed["Bond"].m_as(kj_mol) == 3.3
        assert processed["Angle"].m_as(kj_mol) == 2.2
        assert processed["Torsion"].m_as(kj_mol) == 1.1
        assert processed["RBTorsion"].m_as(kj_mol) == 4.4
        assert processed["Nonbonded"].m_as(kj_mol) == 5.5

    def test_split_forces(self, dummy_system_split):
        processed = _process(
            {
                0: 1.1 * kj_mol,
                1: 2.2 * kj_mol,
                2: 3.3 * kj_mol,
                3: 4.4 * kj_mol,
                4: 0.5 * kj_mol,
                5: -1 * kj_mol,
                6: -2 * kj_mol,  # vdW
                7: -3 * kj_mol,  # Electrostatics
            },
            dummy_system_split,
            False,
            True,
        )

        assert processed["Electrostatics"].m_as(kj_mol) == 0.5
        assert processed["Electrostatics 1-4"].m_as(kj_mol) == -3
        assert processed["vdW"].m_as(kj_mol) == -1
        assert processed["vdW 1-4"].m_as(kj_mol) == -2
