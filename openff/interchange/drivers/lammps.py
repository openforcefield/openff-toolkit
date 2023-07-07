"""Functions for running energy evluations with LAMMPS."""
import subprocess
from shutil import which
from typing import Optional

import numpy
from openff.units import unit

from openff.interchange import Interchange
from openff.interchange.components.mdconfig import MDConfig
from openff.interchange.drivers.report import EnergyReport
from openff.interchange.exceptions import LAMMPSNotFoundError, LAMMPSRunError


def _find_lammps_executable(raise_exception: bool = False) -> Optional[str]:
    """Attempt to locate a LAMMPS executable based on commonly-used names."""
    lammps_executable_names = ["lammps", "lmp_serial", "lmp_mpi"]

    for name in lammps_executable_names:
        if which(name):
            return name

    if raise_exception:
        raise LAMMPSNotFoundError
    else:
        return None


def get_lammps_energies(
    interchange: Interchange,
    round_positions: Optional[int] = None,
    detailed: bool = False,
) -> EnergyReport:
    """
    Given an OpenFF Interchange object, return single-point energies as computed by LAMMPS.

    .. warning :: This API is experimental and subject to change.

    .. todo :: Split out _running_ LAMMPS into a separate internal function

    Parameters
    ----------
    interchange : openff.interchange.Interchange
        An OpenFF Interchange object to compute the single-point energy of
    round_positions : int, optional
        The number of decimal places, in nanometers, to round positions. This can be useful when
        comparing to i.e. GROMACS energies, in which positions may be rounded.
    detailed : bool, optional
        If True, return a detailed energy report containing all energy components.

    Returns
    -------
    report : EnergyReport
        An `EnergyReport` object containing the single-point energies.

    """
    return _process(
        _get_lammps_energies(interchange, round_positions),
        detailed,
    )


def _get_lammps_energies(
    interchange: Interchange,
    round_positions: Optional[int] = None,
) -> dict[str, unit.Quantity]:
    lmp = _find_lammps_executable(raise_exception=True)

    if round_positions is not None:
        interchange.positions = numpy.round(interchange.positions, round_positions)

    interchange.to_lammps("out.lmp")
    mdconfig = MDConfig.from_interchange(interchange)
    mdconfig.write_lammps_input(
        input_file="tmp.in",
    )

    run_cmd = f"{lmp} -i tmp.in"

    proc = subprocess.Popen(
        run_cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )

    _, err = proc.communicate()

    if proc.returncode:
        raise LAMMPSRunError(err)

    # thermo_style custom ebond eangle edihed eimp epair evdwl ecoul elong etail pe
    parsed_energies = unit.kilocalorie_per_mole * _parse_lammps_log("log.lammps")

    return {
        "Bond": parsed_energies[0],
        "Angle": parsed_energies[1],
        "ProperTorsion": parsed_energies[2],
        "ImproperTorsion": parsed_energies[3],
        "vdW": parsed_energies[5],
        "DispersionCorrection": parsed_energies[8],
        "ElectrostaticsShort": parsed_energies[6],
        "ElectrostaticsLong": parsed_energies[7],
    }


def _process(
    energies: dict[str, unit.Quantity],
    detailed: bool = False,
) -> EnergyReport:
    if detailed:
        return EnergyReport(energies=energies)

    return EnergyReport(
        energies={
            "Bond": energies["Bond"],
            "Angle": energies["Angle"],
            "Torsion": energies["ProperTorsion"] + energies["ImproperTorsion"],
            "vdW": energies["vdW"] + energies["DispersionCorrection"],
            "Electrostatics": (
                energies["ElectrostaticsShort"] + energies["ElectrostaticsLong"]
            ),
        },
    )


def _parse_lammps_log(file_in: str) -> list[float]:
    """Parse a LAMMPS log file for energy components."""
    tag = False
    with open(file_in) as fi:
        for line in fi.readlines():
            if tag:
                data = [float(val) for val in line.split()]
                tag = False
            if line.strip().startswith("E_bond"):
                tag = True

    return data
