"""
Commonly-used constants.
"""
from openff.units import unit

_PME = "Ewald3D-ConductingBoundary"
kj_mol = unit.Unit("kilojoule / mol")
kcal_mol = unit.kilocalorie_per_mole

kcal_ang = kcal_mol / unit.angstrom**2
kcal_rad = kcal_mol / unit.radian**2

kj_nm = kj_mol / unit.nanometer**2
kj_rad = kj_mol / unit.radian**2

AMBER_COULOMBS_CONSTANT = 18.2223
kcal_mol_a2 = kcal_mol / unit.angstrom**2
kcal_mol_rad2 = kcal_mol / unit.radian**2
