"""Functions for running energy evluations with all available engines."""
from collections.abc import Iterable

from openff.utilities.utilities import requires_package
from pandas import DataFrame

from openff.interchange import Interchange
from openff.interchange.drivers.amber import get_amber_energies
from openff.interchange.drivers.gromacs import get_gromacs_energies
from openff.interchange.drivers.lammps import get_lammps_energies
from openff.interchange.drivers.openmm import get_openmm_energies
from openff.interchange.drivers.report import EnergyReport
from openff.interchange.exceptions import (
    AmberError,
    GMXError,
    LAMMPSError,
    UnsupportedCutoffMethodError,
)


def get_all_energies(
    interchange: "Interchange",
    combine_nonbonded_forces: bool = False,
    _engines: Iterable[str] = ("OpenMM", "Amber", "GROMACS", "LAMMPS"),
) -> dict[str, EnergyReport]:
    """
    Given an Interchange object, return single-point energies as computed by all available engines.
    """
    # TODO: Return something nan-like if one driver fails, but still return others that succeed
    # TODO: Have each driver return the version of the engine that was used

    try:
        # TODO: Worth wiring this argument up to this function? kwargs complexity is not fun
        all_energies = {
            "OpenMM": get_openmm_energies(
                interchange,
                combine_nonbonded_forces=combine_nonbonded_forces,
            ),
        }
    except UnsupportedCutoffMethodError:
        all_energies = {
            "OpenMM": get_openmm_energies(interchange, combine_nonbonded_forces=False),
        }

    for engine_name, engine_driver, engine_exception in [
        ("Amber", get_amber_energies, AmberError),
        ("GROMACS", get_gromacs_energies, GMXError),
        ("LAMMPS", get_lammps_energies, LAMMPSError),
    ]:
        if engine_name not in _engines:
            continue
        try:
            all_energies[engine_name] = engine_driver(interchange)  # type: ignore[operator]
        except engine_exception:
            pass

    return all_energies


@requires_package("pandas")
def get_summary_data(
    interchange: "Interchange",
    combine_nonbonded_forces: bool = False,
    _engines: Iterable[str] = ("OpenMM", "Amber", "GROMACS", "LAMMPS"),
) -> "DataFrame":
    """Return a pandas DataFrame with summaries of energies from all available engines."""
    from openff.units import unit
    from pandas import DataFrame

    kj_mol = unit.kilojoule / unit.mol

    energies = get_all_energies(
        interchange,
        combine_nonbonded_forces=combine_nonbonded_forces,
        _engines=_engines,
    )

    for k, v in energies.items():
        for kk in v.energies:
            energies[k].energies[kk] = energies[k].energies[kk].m_as(kj_mol)  # type: ignore[union-attr]

    return DataFrame({k: v.energies for k, v in energies.items()}).T
