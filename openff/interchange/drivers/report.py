"""Storing and processing results of energy evaluations."""
import warnings
from typing import Optional

from openff.models.models import DefaultModel
from openff.models.types import FloatQuantity
from openff.units import unit
from pydantic import validator

from openff.interchange.constants import kj_mol
from openff.interchange.exceptions import (
    EnergyError,
    IncompatibleTolerancesError,
    InvalidEnergyError,
)

_KNOWN_ENERGY_TERMS: set[str] = {
    "Bond",
    "Angle",
    "Torsion",
    "RBTorsion",
    "Nonbonded",
    "vdW",
    "Electrostatics",
    "vdW 1-4",
    "Electrostatics 1-4",
}


class EnergyReport(DefaultModel):
    """A lightweight class containing single-point energies as computed by energy tests."""

    # TODO: Should the default be None or 0.0 kj_mol?
    energies: dict[str, Optional[FloatQuantity]] = {
        "Bond": None,
        "Angle": None,
        "Torsion": None,
        "vdW": None,
        "Electrostatics": None,
    }

    @validator("energies")
    def validate_energies(cls, v: dict) -> dict:
        """Validate the structure of a dict mapping keys to energies."""
        for key, val in v.items():
            if key not in _KNOWN_ENERGY_TERMS:
                raise InvalidEnergyError(f"Energy type {key} not understood.")
            if not isinstance(val, unit.Quantity):
                v[key] = FloatQuantity.validate_type(val)
        return v

    @property
    def total_energy(self):
        """Return the total energy."""
        return self["total"]

    def __getitem__(self, item: str) -> Optional[FloatQuantity]:
        if type(item) != str:
            raise LookupError(
                "Only str arguments can be currently be used for lookups.\n"
                f"Found item {item} of type {type(item)}",
            )
        if item in self.energies.keys():
            return self.energies[item]
        if item.lower() == "total":
            return sum(self.energies.values())  # type: ignore
        else:
            return None

    def update(self, new_energies: dict) -> None:
        """Update the energies in this report with new value(s)."""
        self.energies.update(self.validate_energies(new_energies))

    def compare(
        self,
        other: "EnergyReport",
        tolerances: Optional[dict[str, FloatQuantity]] = None,
    ):
        """
        Compare two energy reports.

        Parameters
        ----------
        other: EnergyReport
            The other `EnergyReport` to compare energies against

        tolerances: dict of str: `FloatQuantity`
            Per-key allowed differences in energies

        """
        default_tolerances = {
            "Bond": 1e-3 * kj_mol,
            "Angle": 1e-3 * kj_mol,
            "Torsion": 1e-3 * kj_mol,
            "vdW": 1e-3 * kj_mol,
            "Electrostatics": 1e-3 * kj_mol,
        }

        if tolerances:
            default_tolerances.update(tolerances)

        tolerances = default_tolerances

        energy_differences = self.diff(other)

        if ("Nonbonded" in tolerances) != ("Nonbonded" in energy_differences):
            raise IncompatibleTolerancesError(
                "Mismatch between energy reports and tolerances with respect to whether nonbonded "
                "interactions are collapsed into a single value.",
            )

        errors = dict()

        for key, diff in energy_differences.items():
            if abs(energy_differences[key]) > tolerances[key]:
                errors[key] = diff

        if errors:
            raise EnergyError(errors)

    def diff(
        self,
        other: "EnergyReport",
    ) -> dict[str, FloatQuantity]:
        """
        Return the per-key energy differences between these reports.

        Parameters
        ----------
        other: EnergyReport
            The other `EnergyReport` to compare energies against

        Returns
        -------
        energy_differences : dict of str: `FloatQuantity`
            Per-key energy differences

        """
        energy_differences: dict[str, FloatQuantity] = dict()

        nonbondeds_processed = False

        for key in self.energies:
            if key in ("Bond", "Angle", "Torsion"):
                energy_differences[key] = self[key] - other[key]  # type: ignore[operator]

                continue

            if key in ("Nonbonded", "vdW", "Electrostatics"):
                if nonbondeds_processed:
                    continue

                if (self["vdW"] and other["vdW"]) is not None and (
                    self["Electrostatics"] and other["Electrostatics"]
                ) is not None:
                    for key in ("vdW", "Electrostatics"):
                        energy_differences[key] = self[key] - other[key]  # type: ignore[operator]
                        energy_differences[key] = self[key] - other[key]  # type: ignore[operator]

                        nonbondeds_processed = True

                        continue

                else:
                    energy_differences["Nonbonded"] = (
                        self._get_nonbonded_energy() - other._get_nonbonded_energy()
                    )

                    nonbondeds_processed = True

                    continue

        return energy_differences

    def __sub__(self, other: "EnergyReport") -> dict[str, FloatQuantity]:
        diff = dict()
        for key in self.energies:
            if key not in other.energies:
                warnings.warn(f"Did not find key {key} in second report")
                continue
            diff[key]: FloatQuantity = self.energies[key] - other.energies[key]  # type: ignore

        return diff

    def __str__(self) -> str:
        return (
            "Energies:\n\n"
            f"Bond:          \t\t{self['Bond']}\n"
            f"Angle:         \t\t{self['Angle']}\n"
            f"Torsion:       \t\t{self['Torsion']}\n"
            f"RBTorsion:     \t\t{self['RBTorsion']}\n"
            f"Nonbonded:     \t\t{self['Nonbonded']}\n"
            f"vdW:           \t\t{self['vdW']}\n"
            f"Electrostatics:\t\t{self['Electrostatics']}\n"
        )

    def _get_nonbonded_energy(self) -> FloatQuantity:
        nonbonded_energy = 0.0 * kj_mol
        for key in ("Nonbonded", "vdW", "Electrostatics"):
            if key in self.energies is not None:
                nonbonded_energy += self.energies[key]

        return nonbonded_energy
