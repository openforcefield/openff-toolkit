import pytest
from openff.units import unit
from openff.units.openmm import ensure_quantity

from openff.interchange._tests import _BaseTest
from openff.interchange.constants import kj_mol
from openff.interchange.drivers.report import EnergyReport
from openff.interchange.exceptions import (
    EnergyError,
    IncompatibleTolerancesError,
    InvalidEnergyError,
)


class TestEnergyReport(_BaseTest):
    @pytest.fixture()
    def report(self):
        return EnergyReport(
            energies={
                "Bond": 10 * kj_mol,
                "Angle": 10 * kj_mol,
                "Torsion": 2 * kj_mol,
                "vdW": 20 * kj_mol,
                "Electrostatics": -10 * kj_mol,
            },
        )

    def test_coerce_to_quantity(self):
        assert isinstance(
            EnergyReport(
                energies={
                    "Bond": ensure_quantity(10 * kj_mol, "openmm"),
                },
            )["Bond"],
            unit.Quantity,
        )

    def test_getitem(self, report):
        assert "Bond" in str(report)
        assert report["Bond"].units == kj_mol
        assert report["Nonbonded"] is None

        with pytest.raises(LookupError, match="type <class 'int'>"):
            report[0]

        assert report["Total"].m_as(kj_mol) == 32
        assert report.total_energy == report["Total"]

    def test_bad_constructor(self):
        with pytest.raises(InvalidEnergyError, match="foo not understood."):
            EnergyReport(
                energies={
                    "foo": 1 * kj_mol,
                },
            )

    def test_update(self, report):
        assert report["Bond"].m == 10.0
        report.update({"Bond": -200.0 * kj_mol})

        assert report["Bond"].m == -200.0

    @pytest.mark.parametrize(
        "unit_system",
        ["openff", "openmm"],
    )
    def test_update_unit_types(self, report, unit_system):
        report.update(
            {
                "Bond": ensure_quantity(55.55 * kj_mol, unit_system),
            },
        )

        assert isinstance(report["Bond"], unit.Quantity)
        assert report["Bond"].m_as(kj_mol) == pytest.approx(55.55)

    def test_bad_update(self, report):
        with pytest.raises(InvalidEnergyError, match="foo not understood."):
            report.update({"foo": 1 * kj_mol})

    def test_sub(self):
        a = EnergyReport(energies={"Torsion": 10 * kj_mol})
        b = EnergyReport(energies={"Torsion": 15 * kj_mol})
        c = EnergyReport(energies={"Torsion": 15 * kj_mol, "Nonbonded": 10 * kj_mol})

        assert (b - a)["Torsion"] == 5 * unit.kilojoule / unit.mol

        with pytest.warns(UserWarning, match="Did not find key Nonbonded"):
            c - a

    def test_diff_identical(self, report):
        differences = report.diff(report)

        assert differences == {
            "Bond": 0 * kj_mol,
            "Angle": 0 * kj_mol,
            "Torsion": 0 * kj_mol,
            "vdW": 0 * kj_mol,
            "Electrostatics": 0 * kj_mol,
        }

    def test_diff_combined_nonbonded(self, report):
        single_nonbonded = EnergyReport(
            energies={
                "Bond": 10 * kj_mol,
                "Angle": 10 * kj_mol,
                "Torsion": 2 * kj_mol,
                "Nonbonded": 8 * kj_mol,
            },
        )

        differences = report.diff(single_nonbonded)

        assert "vdW" not in differences
        assert "Electrostatics" not in differences

        assert abs(differences["Nonbonded"].m) == 2

    def test_diff_multiple(self, report):
        other = EnergyReport(
            energies={
                "Bond": -4 * kj_mol,
                "Angle": 10 * kj_mol,
                "Torsion": -2 * kj_mol,
                "vdW": 10 * kj_mol,
                "Electrostatics": -10 * kj_mol,
            },
        )

        diff = report.diff(other)

        assert len(diff) == 5
        assert len([val for val in diff.values() if val.m != 0]) == 3

    def test_compare_identical(self, report):
        report.compare(report)

    def test_compare_different(self, report):
        other = EnergyReport(
            energies={
                key: val + 0.2 * kj_mol
                for key, val in {
                    "Bond": 10 * kj_mol,
                    "Angle": 10 * kj_mol,
                    "Torsion": 2 * kj_mol,
                    "vdW": 20 * kj_mol,
                    "Electrostatics": -10 * kj_mol,
                }.items()
            },
        )

        with pytest.raises(EnergyError):
            report.compare(other)

    def test_compare_multiple(self, report):
        other = EnergyReport(
            energies={
                "Bond": -4 * kj_mol,
                "Angle": 10 * kj_mol,
                "Torsion": -2 * kj_mol,
                "vdW": 10 * kj_mol,
                "Electrostatics": -10 * kj_mol,
            },
        )

        with pytest.raises(
            EnergyError,
            match="Bond.*Torsion.*vdW",
        ):
            report.compare(other)

    def test_compare_incompatible_tolerances(self, report):
        with pytest.raises(
            IncompatibleTolerancesError,
            match="whether nonbonded",
        ):
            report.compare(report, {"Nonbonded": 0.0 * kj_mol})
