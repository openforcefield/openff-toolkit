#!/usr/bin/env python

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

from simtk import unit


# =============================================================================================
# Property Enums
# =============================================================================================

@unique
class PropertyType(IntFlag):

    Undefined          = 0x00
    Density            = 0x01
    DielectricConstant = 0x02

    def __str__(self):

        phases = '|'.join([phase.name for phase in PropertyType if self & phase])
        return phases

    def __repr__(self):
        return str(self)


@unique
class PropertyPhase(IntFlag):

    Undefined = 0x00
    Solid     = 0x01
    Liquid    = 0x02
    Gas       = 0x04

    def __str__(self):

        phases = '|'.join([phase.name for phase in PropertyPhase if self & phase])
        return phases

    def __repr__(self):
        return str(self)


@unique
class CalculationFidelity(Enum):
    """
        An enum describing the fidelity at which a calculation was performed.
    """

    SurrogateModel =   1
    Reweighting =      2
    DirectSimulation = 3


# =============================================================================================
# Property Source
# =============================================================================================

class Source(object):
    """
    Container class for DOI and reference for a given observable

    This class contains either the DOI and/or the reference, but must contain at least one as the observable must
    have a source, even if it was measured in lab.

    Parameters
    ----------
    doi : str or None, default None
        The DOI for the source, preferred way to identify for source
    reference :
        The long form description of the source if no DOI is available, or more
        information is needed or wanted.

    """
    def __init__(self, doi=None, reference=None):

        if doi is None and reference is None:
            raise ValueError("Either a doi and / or a reference must be set")

        self.doi = doi
        self.reference = reference


# =============================================================================================
# Property Representations
# =============================================================================================

class MeasuredPhysicalProperty(object):
    """
    A Measured Physical Property is the implementation of property measured and reported in the ThermoML
    database. It is the combination of a :class:``MeasurementMethod`` method, a :class:``ThermodynamicState``,
    as :class:``Substance``, and then a measured value and uncertainty.


    Implement the :func:`measurement_method` to return the

    Parameters
    ----------
    substance : Substance
        Material/Substance/Chemical that the property was measured for.
    thermodynamic_state : ThermodynamicState
        Physical thermodynamic state with temperature and pressure the property was measured at.
    property_type : PropertyType
        The type of property (e.g density) that was measured.
    property_phase : PropertyPhase
        The phase in which the property was measured.
    value : simtk.unit.Quantity, float, or None
        Value the observable was measured at.
        If None, then this is an unobserved value and just a placeholder (e.g. if its not recorded in database)
        If a float, then the default units for the ThermoML database are assumed
    uncertainty : simtk.unit.Quantity, float, or None
        Uncertainty in the observed measurement
        If None, then it is assumed uncertainty was unreported for this observable
        If a float, then the default units for the ThermoML database are assumed
    source : Source or None, optional, default None
        Reference from which the measurement was taken,
    """
    def __init__(self, substance=None, thermodynamic_state=None,
                 property_type=PropertyType.Undefined, property_phase=PropertyPhase.Undefined,
                 value=None, uncertainty=None, source=None):

        self._thermodynamic_state = thermodynamic_state

        # self._measurement_method = measurement_method
        self.method_name = None

        self.type = property_type
        self.phase = property_phase

        self.substance = substance

        self.value = value
        self.uncertainty = uncertainty

        self._source = source

    @property
    def temperature(self):
        """The temperature which the property was measured at"""
        return None if self._thermodynamic_state is None else self._thermodynamic_state.temperature

    @property
    def pressure(self):
        """The pressure which the property was measured at"""
        return None if self._thermodynamic_state is None else self._thermodynamic_state.pressure

    @property
    def thermodynamic_state(self):
        """Get the thermodynamic state which the property was measured at"""
        return self._thermodynamic_state

    @thermodynamic_state.setter
    def thermodynamic_state(self, thermodynamic_state):
        """Set the thermodynamic state which the property was measured at"""
        self._thermodynamic_state = thermodynamic_state

    # @property
    # def measurement_method(self):
    #     """The method used to measure this physical property"""
    #     return self._measurement_method
    #
    # @measurement_method.setter
    # def measurement_method(self, value):
    #     raise ValueError("measurement_method is not a property which can be set, only implemented as subclass")

    @property
    def doi(self):
        """Get the doi of the article in which this property was published"""
        return None if self._source is None else self._source.doi

    @property
    def reference(self):
        """Get a reference to the location of where this property was published"""
        return None if self._source is None else self._source.reference

    @property
    def source(self):
        """Get the source of this published property"""
        return self._source

    @source.setter
    def source(self, source):
        """Set the source of this published property"""
        self._source = source

    def set_value(self, value, uncertainty):
        """Set the value and uncertainty of this property"""
        self.value = value
        self.uncertainty = uncertainty


class CalculatedPhysicalProperty(object):
    """
    A CalculatedPhysicalProperty represents a property estimated by either
    a surrogate model, reweighting existing data, or directly from a simulation.

    Parameters
    ----------
    substance : Substance
        Material/Substance/Chemical that the property was calculated for.
    thermodynamic_state : ThermodynamicState
        Physical thermodynamic state with temperature and pressure the property was calculated at.
    property_type : PropertyType
        The type of property (e.g density) that was calculated.
    property_phase : PropertyPhase
        The phase in which the property was calculated.
    fidelity : CalculationFidelity
        The fidelity (e.g direct simulation) at which the property was calculated.
    value : simtk.unit.Quantity, float, or None
        The value of the calculated property.
    uncertainty : simtk.unit.Quantity, float, or None
        The uncertainty of the calculated property.
    """
    def __init__(self, substance, thermodynamic_state, property_type,
                 property_phase, fidelity, value, uncertainty):

        self._substance = substance
        self._thermodynamic_state = thermodynamic_state

        self._type = property_type
        self._phase = property_phase

        self._fidelity = fidelity

        self._value = value
        self._uncertainty = uncertainty

    @property
    def temperature(self):
        """The temperature which the property was measured at"""
        return self._thermodynamic_state.temperature

    @property
    def pressure(self):
        """The pressure which the property was measured at"""
        return self._thermodynamic_state.pressure

    @property
    def thermodynamic_state(self):
        """The thermodynamic state which the property was measured at"""
        return self._thermodynamic_state

    @property
    def value(self):
        """The value of the calculated property"""
        return self._value

    @property
    def uncertainty(self):
        """The uncertainty in the calculated property"""
        return self._uncertainty

    @property
    def type(self):
        """The type of property calculated."""
        return self._type

    @property
    def phase(self):
        """The phase in which the property was calculated."""
        return self._phase

    def set_value(self, value, uncertainty):
        """Set the calculated value and uncertainty."""
        self._value = value
        self._uncertainty = uncertainty


class CalculatedPropertySet:
    """
    A collection of calculated physical properties, calculated using
    a specific parameter_set.

    Parameters
    ----------
    properties : dict(str, list(CalculatedPhysicalProperty))
        The list of calculated properties.
    parameter_set : openforcefield.typing.engines.smirnoff.ForceField
        The parameters used to calculate the properties.
    """

    def __init__(self, properties, parameter_set):

        self._properties = properties
        self._parameter_set = parameter_set

    @property
    def properties(self):
        """Returns the calculated properties"""
        return self._properties

    @property
    def parameter_set(self):
        """Returns the parameter set used to calculate the properties"""
        return self._parameter_set
