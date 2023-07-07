import random

import openmm.unit
import pytest
from openff.toolkit.tests.create_molecules import create_ethanol
from openff.units import unit

from openff.interchange import Interchange
from openff.interchange._tests import _BaseTest
from openff.interchange.constants import kj_mol
from openff.interchange.drivers.openmm import get_openmm_energies
from openff.interchange.exceptions import UnsupportedImportError
from openff.interchange.interop.openmm._import import from_openmm
from openff.interchange.interop.openmm._import._import import _convert_nonbonded_force


class TestFromOpenMM(_BaseTest):
    def test_simple_roundtrip(self, monkeypatch, sage_unconstrained):
        monkeypatch.setenv("INTERCHANGE_EXPERIMENTAL", "1")

        molecule = create_ethanol()
        molecule.generate_conformers(n_conformers=1)

        interchange = Interchange.from_smirnoff(
            sage_unconstrained,
            [molecule],
            box=[4, 4, 4] * unit.nanometer,
        )

        system = interchange.to_openmm(combine_nonbonded_forces=True)

        converted = from_openmm(
            topology=interchange.topology.to_openmm(),
            system=system,
            positions=interchange.positions,
            box_vectors=interchange.box,
        )

        get_openmm_energies(interchange).compare(
            get_openmm_energies(converted),
            tolerances={
                "Bond": 1e-6 * kj_mol,
                "Angle": 1e-6 * kj_mol,
                "Torsion": 1e-6 * kj_mol,
                "Nonbonded": 1e-3 * kj_mol,
            },
        )


class TestConvertNonbondedForce:
    @pytest.mark.parametrize(
        "method",
        [
            openmm.NonbondedForce.NoCutoff,
            openmm.NonbondedForce.CutoffNonPeriodic,
            openmm.NonbondedForce.CutoffPeriodic,
        ],
    )
    def test_unsupported_method(self, method):
        force = openmm.NonbondedForce()
        force.setNonbondedMethod(method)

        with pytest.raises(UnsupportedImportError):
            _convert_nonbonded_force(force)

    def test_parse_switching_distance(self):
        force = openmm.NonbondedForce()
        force.setNonbondedMethod(openmm.NonbondedForce.PME)

        cutoff = 1 + random.random() * 0.1
        switch_width = random.random() * 0.1

        force.setCutoffDistance(cutoff)
        force.setUseSwitchingFunction(True)
        force.setSwitchingDistance(cutoff - switch_width)

        vdw, _ = _convert_nonbonded_force(force)

        assert vdw.cutoff.m_as(unit.nanometer) == pytest.approx(cutoff)
        assert vdw.switch_width.m_as(unit.nanometer) == pytest.approx(switch_width)

    def test_parse_switching_distance_unused(self):
        force = openmm.NonbondedForce()
        force.setNonbondedMethod(openmm.NonbondedForce.PME)

        cutoff = 1 + random.random() * 0.1

        force.setCutoffDistance(cutoff)

        vdw, _ = _convert_nonbonded_force(force)

        assert vdw.cutoff.m_as(unit.nanometer) == pytest.approx(cutoff)
        assert vdw.switch_width.m_as(unit.nanometer) == 0.0


class TestConvertConstraints(_BaseTest):
    def test_num_constraints(self, monkeypatch, sage, basic_top):
        """Test that the number of constraints is preserved when converting to and from OpenMM"""
        monkeypatch.setenv("INTERCHANGE_EXPERIMENTAL", "1")

        interchange = sage.create_interchange(basic_top)

        converted = from_openmm(
            topology=interchange.topology.to_openmm(),
            system=interchange.to_openmm(combine_nonbonded_forces=True),
        )

        assert "Constraints" in interchange.collections
        assert "Constraints" in converted.collections

        assert len(interchange["Constraints"].key_map) == len(
            converted["Constraints"].key_map,
        )
