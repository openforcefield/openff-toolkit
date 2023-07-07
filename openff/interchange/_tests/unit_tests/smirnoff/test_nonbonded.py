import numpy
import pytest
from openff.toolkit import Topology
from openff.toolkit.typing.engines.smirnoff import (
    ChargeIncrementModelHandler,
    ElectrostaticsHandler,
    LibraryChargeHandler,
    ToolkitAM1BCCHandler,
)
from openff.units import unit

from openff.interchange import Interchange
from openff.interchange._tests import _BaseTest
from openff.interchange.smirnoff._nonbonded import SMIRNOFFElectrostaticsCollection


class TestNonbonded(_BaseTest):
    @pytest.mark.slow()
    def test_electrostatics_am1_handler(self, methane):
        methane.assign_partial_charges(partial_charge_method="am1bcc")

        # Explicitly store these, since results differ RDKit/AmberTools vs. OpenEye
        reference_charges = [c.m for c in methane.partial_charges]

        parameter_handlers = [
            ElectrostaticsHandler(version=0.3),
            ToolkitAM1BCCHandler(version=0.3),
        ]

        electrostatics_handler = SMIRNOFFElectrostaticsCollection.create(
            parameter_handlers,
            methane.to_topology(),
        )
        numpy.testing.assert_allclose(
            [charge.m_as(unit.e) for charge in electrostatics_handler.charges.values()],
            reference_charges,
        )

    def test_electrostatics_library_charges(self, methane):
        library_charge_handler = LibraryChargeHandler(version=0.3)
        library_charge_handler.add_parameter(
            {
                "smirks": "[#6X4:1]-[#1:2]",
                "charge1": -0.1 * unit.elementary_charge,
                "charge2": 0.025 * unit.elementary_charge,
            },
        )

        parameter_handlers = [
            ElectrostaticsHandler(version=0.3),
            library_charge_handler,
        ]

        electrostatics_handler = SMIRNOFFElectrostaticsCollection.create(
            parameter_handlers,
            methane.to_topology(),
        )

        numpy.testing.assert_allclose(
            [charge.m_as(unit.e) for charge in electrostatics_handler.charges.values()],
            [-0.1, 0.025, 0.025, 0.025, 0.025],
        )

    def test_electrostatics_charge_increments(self, hydrogen_chloride):
        hydrogen_chloride.assign_partial_charges(partial_charge_method="am1-mulliken")

        reference_charges = [c.m for c in hydrogen_chloride.partial_charges]
        reference_charges[0] += 0.1
        reference_charges[1] -= 0.1

        charge_increment_handler = ChargeIncrementModelHandler(version=0.3)
        charge_increment_handler.add_parameter(
            {
                "smirks": "[#17:1]-[#1:2]",
                "charge_increment1": 0.1 * unit.elementary_charge,
                "charge_increment2": -0.1 * unit.elementary_charge,
            },
        )

        parameter_handlers = [
            ElectrostaticsHandler(version=0.3),
            charge_increment_handler,
        ]

        electrostatics_handler = SMIRNOFFElectrostaticsCollection.create(
            parameter_handlers,
            hydrogen_chloride.to_topology(),
        )

        # AM1-Mulliken charges are [-0.168,  0.168], increments are [0.1, -0.1],
        # sum is [-0.068,  0.068]
        numpy.testing.assert_allclose(
            [charge.m_as(unit.e) for charge in electrostatics_handler.charges.values()],
            reference_charges,
        )

    @pytest.mark.slow()
    def test_toolkit_am1bcc_uses_elf10_if_oe_is_available(self, sage, hexane_diol):
        """
        Ensure that the ToolkitAM1BCCHandler assigns ELF10 charges if OpenEye is available.

        Taken from https://github.com/openforcefield/openff-toolkit/pull/1214,
        """
        try:
            hexane_diol.assign_partial_charges(partial_charge_method="am1bccelf10")
            uses_elf10 = True
        except ValueError:
            # This assumes that the ValueError stems from "am1bccelf10" and not other sources; the
            # toolkit should implement a failure mode that does not clash with other `ValueError`s
            hexane_diol.assign_partial_charges(partial_charge_method="am1bcc")
            uses_elf10 = False

        partial_charges = [c.m for c in hexane_diol.partial_charges]

        assigned_charges = [
            v.m
            for v in Interchange.from_smirnoff(sage, [hexane_diol])[
                "Electrostatics"
            ].charges.values()
        ]

        try:
            from openeye import oechem

            openeye_available = oechem.OEChemIsLicensed()
        except ImportError:
            openeye_available = False

        if openeye_available:
            assert uses_elf10
            numpy.testing.assert_allclose(partial_charges, assigned_charges)
        else:
            assert not uses_elf10
            numpy.testing.assert_allclose(partial_charges, assigned_charges)


class TestSMIRNOFFChargeIncrements(_BaseTest):
    @pytest.fixture()
    def hydrogen_cyanide_charge_increments(self):
        handler = ChargeIncrementModelHandler(
            version=0.4,
            partial_charge_method="formal_charge",
        )
        handler.add_parameter(
            {
                "smirks": "[H:1][C:2]",
                "charge_increment1": -0.111 * unit.elementary_charge,
                "charge_increment2": 0.111 * unit.elementary_charge,
            },
        )
        handler.add_parameter(
            {
                "smirks": "[C:1]#[N:2]",
                "charge_increment1": 0.5 * unit.elementary_charge,
                "charge_increment2": -0.5 * unit.elementary_charge,
            },
        )

        return handler

    def test_no_charge_increments_applied(self, sage, hexane_diol):
        gastiger_charges = [c.m for c in hexane_diol.partial_charges]
        sage.deregister_parameter_handler("ToolkitAM1BCC")

        no_increments = ChargeIncrementModelHandler(
            version=0.3,
            partial_charge_method="gasteiger",
        )
        sage.register_parameter_handler(no_increments)

        assert len(sage["ChargeIncrementModel"].parameters) == 0

        out = Interchange.from_smirnoff(sage, [hexane_diol])
        assert numpy.allclose(
            numpy.asarray([v.m for v in out["Electrostatics"].charges.values()]),
            gastiger_charges,
        )

    def test_overlapping_increments(self, sage, methane):
        """Test that separate charge increments can be properly applied to the same atom."""
        sage.deregister_parameter_handler("ToolkitAM1BCC")
        charge_handler = ChargeIncrementModelHandler(
            version=0.3,
            partial_charge_method="formal_charge",
        )
        charge_handler.add_parameter(
            {
                "smirks": "[C:1][H:2]",
                "charge_increment1": 0.111 * unit.elementary_charge,
                "charge_increment2": -0.111 * unit.elementary_charge,
            },
        )
        sage.register_parameter_handler(charge_handler)
        assert 0.0 == pytest.approx(
            sum(
                v.m
                for v in Interchange.from_smirnoff(sage, [methane])[
                    "Electrostatics"
                ].charges.values()
            ),
        )

    def test_charge_increment_forwawrd_reverse_molecule(
        self,
        sage,
        hydrogen_cyanide,
        hydrogen_cyanide_reversed,
        hydrogen_cyanide_charge_increments,
    ):
        sage.deregister_parameter_handler("ToolkitAM1BCC")
        sage.register_parameter_handler(hydrogen_cyanide_charge_increments)

        topology = Topology.from_molecules(
            [hydrogen_cyanide, hydrogen_cyanide_reversed],
        )

        out = Interchange.from_smirnoff(sage, topology)

        expected_charges = [-0.111, 0.611, -0.5, -0.5, 0.611, -0.111]

        # TODO: Fix get_charges to return the atoms in order
        found_charges = [0.0] * topology.n_atoms
        for key, val in out["Electrostatics"].charges.items():
            found_charges[key.atom_indices[0]] = val.m

        assert numpy.allclose(expected_charges, found_charges)
