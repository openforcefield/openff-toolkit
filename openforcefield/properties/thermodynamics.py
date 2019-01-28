# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Thermodynamics API.

Authors
-------
* John D. Chodera <john.chodera@choderalab.org>
* Levi N. Naden <levi.naden@choderalab.org>
* Simon Boothroyd <simon.boothroyd@choderalab.org>

"""
# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

from enum import Enum
from typing import Optional

from pydantic import BaseModel, validator
from simtk import unit

from openforcefield.utils.serialization import deserialize_quantity, serialize_quantity


# =============================================================================================
# THERMODYNAMIC STATE
# =============================================================================================

class Ensemble(Enum):
    """An enum describing the available thermodynamic ensembles.
    """
    NVT = "NVT"
    NPT = "NPT"


class ThermodynamicState(BaseModel):
    """
    Data specifying a physical thermodynamic state obeying Boltzmann statistics.

    Properties
    ----------
    temperature : simtk.unit.Quantity with units compatible with kelvin
        The external temperature
    pressure : simtk.unit.Quantity with units compatible with atmospheres
        The external pressure

    Examples
    --------
    Specify an NPT state at 298 K and 1 atm pressure.

    >>> state = ThermodynamicState(temperature=298.0*unit.kelvin, pressure=1.0*unit.atmospheres)

    Note that the pressure is only relevant for periodic systems.

    """

    temperature: Optional[unit.Quantity] = None
    pressure: Optional[unit.Quantity] = None

    class Config:

        arbitrary_types_allowed = True

        json_encoders = {
            unit.Quantity: lambda v: serialize_quantity(v),
        }

    @validator('temperature', pre=True, whole=True)
    def validate_temperature(cls, v):

        if isinstance(v, dict):
            v = deserialize_quantity(v)

        if isinstance(v, unit.Quantity):
            v = v.in_units_of(unit.kelvin)

        return v
    
    @validator('pressure', pre=True, whole=True)
    def validate_pressure(cls, v):

        if isinstance(v, dict):
            v = deserialize_quantity(v)

        if isinstance(v, unit.Quantity):
            v = v.in_units_of(unit.atmospheres)

        return v

    def __repr__(self):
        """
        Returns a string representation of a state.
        """

        return_value = "ThermodynamicState("

        if self.temperature is not None:
            return_value += "temperature={0:s}, ".format(repr(self.temperature))
        if self.pressure is not None:
            return_value += "pressure = {0:s}".format(repr(self.pressure))

        return_value += ")"

        return return_value

    def __str__(self):

        return_value = "<ThermodynamicState object"

        if self.temperature is not None:
            return_value += ", temperature = {0:s}".format(str(self.temperature))
        if self.pressure is not None:
            return_value += ", pressure = {0:s}".format(str(self.pressure))

        return_value += ">"

        return return_value

    def __hash__(self):
        return hash((str(self.temperature), str(self.pressure)))

    def __eq__(self, other):

        return (abs(self.temperature - other.temperature) < 0.0001 * unit.kelvin and
                abs(self.pressure - other.pressure) < 0.0001 * unit.atmosphere)

    def __ne__(self, other):
        return not (self == other)
