from typing import Union

import pytest
from openff.toolkit import Topology
from openff.units import unit

from openff.interchange._tests import _BaseTest
from openff.interchange.components.mdconfig import (
    MDConfig,
    get_intermol_defaults,
    get_smirnoff_defaults,
)
from openff.interchange.constants import _PME


@pytest.fixture()
def system_no_constraints(sage_unconstrained, basic_top):
    return sage_unconstrained.create_interchange(basic_top)


@pytest.fixture()
def rigid_water_box(sage, water):
    topology = water.to_topology()
    topology.box_vectors = unit.Quantity([5, 5, 5], unit.nanometer)
    return sage.create_interchange(topology)


@pytest.fixture()
def constrained_ligand_rigid_water_box(sage, basic_top, water):
    topology = Topology(basic_top)
    topology.add_molecule(water)
    return sage.create_interchange(topology)


@pytest.fixture()
def unconstrained_ligand_rigid_water_box(sage_unconstrained, basic_top, water):
    topology = Topology(basic_top)
    topology.add_molecule(water)
    return sage_unconstrained.create_interchange(topology)


def parse_mdp(file: str) -> dict[str, str]:
    """Naively parse an MDP file into a dict structure."""
    options = dict()

    with open(file) as f:
        for line in f.readlines():
            split = [token.replace("-", "") for token in line.split()]

            if len(split) == 0:
                continue

            elif len(split) == 3:
                assert split[1] == "="

                options[split[0].lower()] = split[2].lower()

            else:
                raise Exception

    return options


def parse_sander(file: str) -> dict[str, Union[dict, str]]:
    """Naively parse (sections of) a sander input file into a dict structure."""
    options: dict[str, Union[dict, str]] = dict()
    current_level = options

    with open(file) as f:
        for number, line in enumerate(f.readlines()):
            if number == 0:
                continue

            line = line.strip().replace(" ", "").replace(",", "").lower()

            if line.startswith("/"):
                continue

            if line.startswith("&"):
                current_section = line[1:]
                current_level[current_section] = dict()
                current_level = current_level[current_section]  # type: ignore[assignment]
                continue

            if "=" in line:
                key, value = line.split("=")
                current_level[key.lower()] = value
                continue

            raise ValueError(f"Unexpected content {line} on line {number} of {file}.")

    return options


class TestMDConfigFromInterchange(_BaseTest):
    @pytest.mark.parametrize("switch", [True, False])
    def test_from_interchange(self, sage, basic_top, switch):
        from openff.units import unit

        from openff.interchange import Interchange

        if not switch:
            sage["vdW"].switch_width = 0.0 * unit.nanometer

        interchange = Interchange.from_smirnoff(sage, basic_top)
        config = MDConfig.from_interchange(interchange)

        if switch:
            assert config.switching_function
            assert config.switching_distance.m_as(unit.nanometer) == pytest.approx(0.8)
        else:
            assert not config.switching_function
            # No need to check the value of `switching_distance` ... right?


class TestSMIRNOFFDefaults(_BaseTest):
    @pytest.mark.parametrize("periodic", [True, False])
    def test_apply_smirnoff_defaults(self, sage, basic_top, periodic):
        from openff.units import unit

        from openff.interchange import Interchange

        interchange = Interchange.from_smirnoff(sage, basic_top)
        get_smirnoff_defaults(periodic=True).apply(interchange)

        for attr, value in zip(["cutoff", "switch_width"], [0.9, 0.1]):
            assert getattr(interchange["vdW"], attr).m_as(
                unit.nanometer,
            ) == pytest.approx(value)

        assert interchange["vdW"].method == "cutoff"

        if periodic:
            assert interchange["Electrostatics"].periodic_potential == _PME
        else:
            assert interchange["Electrostatics"].nonperiodic_potential == "Coulomb"


class TestIntermolDefaults(_BaseTest):
    @pytest.mark.parametrize("periodic", [True, False])
    def test_write_mdp(self, periodic):
        """
        https://github.com/shirtsgroup/InterMol/blob/master/intermol/tests/gromacs/grompp_vacuum.mdp
        https://github.com/shirtsgroup/InterMol/blob/master/intermol/tests/gromacs/grompp.mdp
        """

        get_intermol_defaults(
            periodic=periodic,
        ).write_mdp_file("tmp.mdp")

        options = parse_mdp("tmp.mdp")

        assert options["pbc"] == "xyz" if periodic else "no"

        assert options["coulombtype"] == "pme" if periodic else "cutoff"
        assert options["rcoulomb"] == "0.9" if periodic else "2.0"
        assert options["coulombmodifier"] == "none"

        assert options["vdwtype"] == "cutoff"
        assert options["rvdw"] == "0.9"
        assert options["vdwmodifier"] == "none"

        # This is set to "no" in the vacuum file, though it probably doesn't matter
        assert options["dispcorr"] == "ener"

        assert options["constraints"] == "none"

    def test_apply_intermol_defaults(self, sage, basic_top):
        from openff.units import unit

        from openff.interchange import Interchange

        interchange = Interchange.from_smirnoff(sage, basic_top)
        get_intermol_defaults(periodic=True).apply(interchange)

        assert interchange["vdW"].switch_width.m == 0.0

        for key in ["vdW", "Electrostatics"]:
            assert interchange[key].cutoff.m_as(unit.nanometer) == pytest.approx(0.9)


class TestWriteSanderInput(_BaseTest):
    def test_system_no_constraints(self, system_no_constraints):
        MDConfig.from_interchange(system_no_constraints).write_sander_input_file("1.in")

        options = parse_sander("1.in")

        assert options["cntrl"]["ntf"] == "1"
        assert options["cntrl"]["ntc"] == "1"

    def test_rigid_water_box(self, rigid_water_box):
        MDConfig.from_interchange(rigid_water_box).write_sander_input_file("2.in")

        options = parse_sander("2.in")

        # Unclear what's "supposed" to be the value here
        assert options["cntrl"]["ntf"] in ("2", "4")
        assert options["cntrl"]["ntc"] == "1"

    def test_constrained_ligand_rigid_water_box(
        self,
        constrained_ligand_rigid_water_box,
    ):
        MDConfig.from_interchange(
            constrained_ligand_rigid_water_box,
        ).write_sander_input_file("3.in")

        options = parse_sander("3.in")

        # Amber does not support this case (constrain all h-bonds and water angles, but no other angles),
        # this option seems to be the best one
        assert options["cntrl"]["ntf"] == "2"
        assert options["cntrl"]["ntc"] == "1"

    def test_unconstrained_ligand_rigid_water_box(
        self,
        unconstrained_ligand_rigid_water_box,
    ):
        MDConfig.from_interchange(
            unconstrained_ligand_rigid_water_box,
        ).write_sander_input_file("3.in")

        options = parse_sander("3.in")

        # Amber does not support this case (constrain wtaer bonds and water angles, but no other bonds or angles),
        # this option seems to be the best one
        assert options["cntrl"]["ntf"] == "2"
        assert options["cntrl"]["ntc"] == "1"

    def test_fswitch_negative_when_no_switching_function(
        self,
        unconstrained_ligand_rigid_water_box,
    ):
        """Reproduce issue 745."""
        MDConfig.from_interchange(
            unconstrained_ligand_rigid_water_box,
        ).write_sander_input_file("yes.in")

        # With OpenFF defaults, the switch starts at 8 A and the cutoff is at 9 A
        assert parse_sander("yes.in")["cntrl"]["fswitch"] == "8.0"

        unconstrained_ligand_rigid_water_box["vdW"].switch_width = 0.0 * unit.nanometer

        MDConfig.from_interchange(
            unconstrained_ligand_rigid_water_box,
        ).write_sander_input_file("no.in")

        assert parse_sander("no.in")["cntrl"]["fswitch"] == "-1.0"
