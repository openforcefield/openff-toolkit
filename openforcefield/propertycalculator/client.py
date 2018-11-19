#!/usr/bin/env python

# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Property calculator client side API.

Authors
-------
* Simon Boothroyd <simon.boothroyd@choderalab.org>

"""


# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

from enum import Enum, unique
from openforcefield.propertycalculator.runner import PropertyCalculationRunner


# =============================================================================================
# Helper Classes
# =============================================================================================

class CalculationFidelity(Enum):
    """
        An enum describing the fidelity at which a calculation was performed.
    """

    SurrogateModel:   0x01
    Reweighting:      0x02
    DirectSimulation: 0x03


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
    property_phase : PropertyType
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
    """

    def __init__(self, properties, parameter_set):

        self._calculated_properties = properties
        self._parameter_set = parameter_set

    @property
    def calculated_properties(self):
        """The calculated properties"""
        return self._calculated_properties

    @property
    def parameter_set(self):
        """The parameter set used to calculate the properties"""
        return self._parameter_set


# =============================================================================================
# Property Estimator
# =============================================================================================

class PropertyEstimator(object):
    """
    The object responsible for requesting a set of properties
    be calculated by the low-level property calculation backend,
    and for analysing the performance of the parameters.
    """

    @staticmethod
    def compute_properties(data_set, parameter_set):
        """
        Submit the property and parameter set for calculation.

        Parameters
        ----------
        data_set : PropertyDataSet
            The set of properties to attempt to compute.
        parameter_set : ParameterSet
            The OpenFF parameter set to use for the calculations.
        """

        if data_set is None or parameter_set is None:
            raise ValueError('Both a data set and parameter set must be '
                             'present to compute physical properties.')

        # In principle this method would simply push all properties and
        # to the backend which will decide what to do with them.

        # For now, just create the backend manually on the local device.
        calculation_runner = PropertyCalculationRunner()

        # In practice such a complicated runner will need to report back
        # detailed diagnostics of what ran and why, and what if anything
        # went wrong.
        calculated_properties = calculation_runner.run(data_set, parameter_set)

        return calculated_properties

    @staticmethod
    def produce_calculation_report(measured_data_set, calculated_data_set):
        """
        Produce a report detailing how well a measured and calculated data
        set match.

        Parameters
        ----------
        measured_data_set : PropertyDataSet
            The set of measured properties to compare against.
        calculated_data_set : ParameterSet
            The set of calculated properties to analyse.
        """

        return ''


