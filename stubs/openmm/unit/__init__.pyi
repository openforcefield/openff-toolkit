from collections.abc import Generator

import numpy
from openmm.unit.baseunit import BaseUnit

class Quantity:
    def __init__(self, *args, **kwargs) -> None: ...
    @property
    def unit(self) -> Unit: ...
    def value_in_unit(self, unit: Unit | str) -> float | int | list | numpy.ndarray: ...

class Unit:
    def __init__(self, name: str) -> None: ...
    def iter_base_or_scaled_units(
        self,
    ) -> Generator[tuple[BaseUnit, int]]: ...
    def __mul__(self, other: float | int | list | numpy.ndarray) -> Quantity: ...
    def __rmul__(self, other: float | int | list | numpy.ndarray) -> Quantity: ...

dalton = Unit("dalton")
dimensionless = Unit("dimensionless")
atmosphere = Unit("atmosphere")
angstrom = Unit("angstrom")
angstroms = Unit("angstroms")
