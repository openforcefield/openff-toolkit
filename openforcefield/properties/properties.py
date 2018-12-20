# !/usr/bin/env python

# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Properties base API.

Authors
-------
* Levi N. Naden <levi.naden@choderalab.org> (original layout)
* Simon Boothroyd <simon.boothroyd@choderalab.org>

"""
# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

from enum import Enum, IntFlag, unique

from pydantic import BaseModel, validator
from pydantic.validators import dict_validator
from typing import Optional

from simtk import unit

from openforcefield.utils.serialization import deserialize_quantity, serialize_quantity

from openforcefield.properties.thermodynamics import ThermodynamicState
from openforcefield.properties.substances import Substance


# =============================================================================================
# Property Registry
# =============================================================================================

registered_properties = {}


def register_property(cls):

    registered_properties[cls.__name__] = cls
    return cls


# =============================================================================================
# Property Enums
# =============================================================================================

@unique
class PropertyPhase(IntFlag):
    """An enum describing the phase a property was collected in.
    """

    Undefined = 0x00
    Solid = 0x01
    Liquid = 0x02
    Gas = 0x04

    def __str__(self):
        """
        Returns
        ---
        str
            A string representation of the PropertyPhase enum
        """
        phases = '|'.join([phase.name for phase in PropertyPhase if self & phase])
        return phases

    def __repr__(self):
        """
        Returns
        ---
        str
            A string representation of the PropertyPhase enum
        """
        return str(self)


@unique
class CalculationFidelity(IntFlag):
    """An enum describing the fidelity at which a property was measured.
    """
    SurrogateModel = 1
    Reweighting = 2
    DirectSimulation = 3


# =============================================================================================
# Property Sources
# =============================================================================================

class Source(BaseModel):
    """Container class for information about how a property was measured / calculated.
    """

    @classmethod
    def get_validators(cls):
        # yield dict_validator
        yield cls.validate

    @classmethod
    def validate(cls, value):
        # A dirty hack to ensure proper inheritance..

        if isinstance(value, cls):
            return value
        else:

            if 'doi' in value:
                return MeasurementSource(**value)
            elif 'fidelity' in value:
                return CalculationSource(**value)

            cls(**dict_validator(value))


class MeasurementSource(Source):
    """Contains any metadata about how a physical property was measured by experiment.

    This class contains either the DOI and/or the reference, but must contain at
    least one as the observable must have a source, even if it was measured in lab.

    Attributes
    ----------
    doi : str or None, default None
        The DOI for the source, preferred way to identify for source
    reference : str
        The long form description of the source if no DOI is available, or more
        information is needed or wanted.
    """

    doi: Optional[str] = None
    reference: Optional[str] = None

    @classmethod
    def get_validators(cls):
        # yield dict_validator
        yield cls.validate

    @classmethod
    def validate(cls, value):
        if isinstance(value, cls):
            return value
        else:
            return cls(**dict_validator(value))


class CalculationSource(Source):
    """Contains any metadata about how a physical property was calculated.

    This includes at which fidelity the property was calculated at (e.g Direct
    simulation, reweighting, ...) in addition to the parameters which were
    used as part of the calculations.

    Attributes
    ----------
    fidelity : CalculationFidelity
        The fidelity at which the property was calculated
    provenance : str
        A JSON string containing information about how the property was calculated.
    """

    fidelity: CalculationFidelity = CalculationFidelity.DirectSimulation
    provenance: str = None


# =============================================================================================
# Property Definitions
# =============================================================================================

class PhysicalProperty(BaseModel):
    """Represents the value of any physical property and it's uncertainty.

    It additionally stores the thermodynamic state at which the property
    was collected, the phase it was collected in, information about
    the composition of the observed system, and metadata about how the
    property was collected.
    """

    thermodynamic_state: ThermodynamicState = None
    phase: PropertyPhase = PropertyPhase.Undefined

    substance: Substance = None

    value: unit.Quantity = None
    uncertainty: unit.Quantity = None

    source: Source = None

    @validator('value', 'uncertainty', pre=True, whole=True)
    def validate_quantity(cls, v):

        if isinstance(v, dict):
            v = deserialize_quantity(v)

        return v

    class Config:

        # A dirty hack to allow simtk.unit.Quantities...
        # TODO: Should really invesitigate QCElemental as an alternative.
        arbitrary_types_allowed = True

        json_encoders = {
            unit.Quantity: lambda v: serialize_quantity(v),
        }

    @property
    def temperature(self):
        """simtk.unit.Quantity or None: The temperature at which the property was collected."""
        return None if self.thermodynamic_state is None else self.thermodynamic_state.temperature

    @property
    def pressure(self):
        """simtk.unit.Quantity or None: The pressure at which the property was collected."""
        return None if self.thermodynamic_state is None else self.thermodynamic_state.pressure

    def set_value(self, value, uncertainty):
        """Set the value and uncertainty of this property.

        Parameters
        ----------
        value : simtk.unit.Quantity
            The value of the property.
        uncertainty : simtk.unit.Quantity
            The uncertainty in the properties value.
        """
        self.value = value
        self.uncertainty = uncertainty

    @staticmethod
    def get_default_calculation_schema():
        """Returns the set of steps needed to calculate
        this property by direct simulation methods.

        Returns
        -------
        openforcefield.properties.estimator.CalculationSchema
            The calculation schema to follow.
        """
        raise NotImplementedError()
