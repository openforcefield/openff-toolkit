#!/usr/bin/env python

# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Measurement Methods API.

Authors
-------
* Levi N. Naden <levi.naden@choderalab.org>

TODO
----

"""
# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

import re
import numpy as np
from openforcefield import substances, thermodynamics

from abc import ABC, abstractproperty, abstractmethod

from simtk import unit

from typing import Optional, Union


# =============================================================================================
# MODULE FUNCTIONS
# =============================================================================================

def cast_thermoml_string(thermoml_string: str) -> str:
    """
    Helper function that converts an ThemoML tag into the base PythonClass style.
    Obeys the following rules:

    1. Characters including and after the Regex match "[^A-Za-z0-9 -]" (non-word character) are dropped
    2. The word "method" at the end is dropped
    3. Remaining words are capitalized
    4. Remaining spaces are removed

    Examples
    --------
    >>> cast_thermoml_string("Capillary tube (Ostwald; Ubbelohde) method")
    "CapillaryTube"
    >>> cast_thermoml_string("Vibrating tube method")
    "VibratingTube"

    Parameters
    ----------
    thermoml_string : str
        Data associated with a ThermoML Database tag to be converted

    Returns
    -------
    converted_thermoml_string : str
    """
    stripped_string = re.sub("[^A-Za-z0-9 -].*", "", thermoml_string).lower()
    stripped_string = re.sub("method$", "", stripped_string)
    return stripped_string.title().replace(" ", "")


# =============================================================================================
# MEASURED PHYSICAL PROPERTY CLASS
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
    measurement_method : MeasurementMethod
        Observable and the method used to measure it. This object must be fully implemented
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
    def __init__(self, substance: substances.Substance,
                 thermodynamic_state: thermodynamics.ThermodynamicState,
                 measurement_method: MeasurementMethod,
                 value: Optional[Union[float, unit.Quantity]]=None,
                 uncertainty: Optional[Union[float, unit.Quantity]]=None,
                 source: Optional[Source]=None):
        self.substance = substance
        self._thermodynamic_state = thermodynamic_state
        self.value = value
        self.uncertainty = uncertainty
        self._measurement_method = measurement_method
        self._source = source

    @property
    def temperature(self) -> unit.Quantity:
        """Temperature which the property was measured at"""
        return self._thermodynamic_state.temperature

    @property
    def pressure(self) -> unit.Quantity:
        """Pressure the property was measured at"""
        return self._thermodynamic_state.pressure

    @property
    def measurement_method(self) -> MeasurementMethod:
        """The method used to measure this physical property"""
        return self._measurement_method

    @measurement_method.setter
    def measurement_method(self, value):
        raise ValueError("measurement_method is not a property which can be set, only implemented as subclass")

    @property
    def doi(self) -> Union[str, None]:
        try:
            self._source.doi
        except AttributeError:
            return None

    @property
    def reference(self) -> Union[str, None]:
        return self._fetch_source('reference')

    @property
    def source(self):
        return self._source

    def _fetch_source(self, source_attr):
        try:
            return getattr(self._source, source_attr)
        except AttributeError:
            return None


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
    def __init__(self, doi: Optional[str]=None, reference: Optional[str]=None):
        if doi is None and reference is None:
            raise ValueError("Must set either doi or reference")
        self.doi = doi
        self.reference = reference


# =============================================================================================
# ABSTRACT CLASSES
# =============================================================================================

class MeasurementMethod(ABC):
    """
    A Measurement Method is the abstract implementation of ThermoML Physical Measurement. Encoded within is the
    measurement method itself, and the physical observable the measurement is computing.

    This class knows the structure of the property in the ThermoML database, and how to compute the observable from
    simulation.

    Encoded in this class should be both the
    """
    def __init__(self):
        pass

    # ------------------
    # Abstract functions
    # ------------------

    @abstractproperty
    def method_name(self) -> str:
        """
        Observable calculation method as it appears in the ThermoML ``eMethodName`` tag.

        The ``eMethodName``'s are standardized in ThermoML, however, ``sMethodNames`` are not. The ``sMethodNames``
        are therefor not supported
        """
        pass

    @abstractproperty
    def prop_name(self) -> str:
        """Observable name as it appears in the ThermoML ``ePropName`` tag"""
        pass

    @abstractmethod
    def log_likelihood_data(self,
                            experimental_property,
                            calculated_property,
                            thermodynamic_state: thermodynamics.ThermodynamicState
                            ) -> float:
        """
        Likelihood of the observable is accurate given the experimental value of the observable, the calculated
        value, and the current thermodynamic state.

        Parameters
        ----------
        experimental_property
        calculated_property
        thermodynamic_state

        Returns
        -------
        likelihood : float
        """
        pass

    @abstractmethod
    def error_model(self, *args):
        """Function to compute the error model of the observable"""
        pass

    @abstractmethod
    def compute_property(self,
                         sim_data,
                         thermodynamic_state: thermodynamics.ThermodynamicState
                         ):
        """
        Compute the observable this class has bee

        Parameters
        ----------
        sim_data : Simulation Data
            Data to compute the property through simulation or other synthetic estimation technique
        thermodynamic_state : thermodynamics.ThermodynamicState
            Physical conditions to emulate

        Returns
        -------
        observable : float or unit.Quantity
            The observable being computed
        sim_data : simulation data
            Output simulation data to be analyzed for statistics after property estimation
        """
        pass

    # -----------------------
    # Sanity/Safety functions
    # -----------------------

    @method_name.setter
    def method_name(self, value):
        raise ValueError("method_name is not a property which can be set, only implemented as subclass")

    @prop_name.setter
    def prop_name(self, value):
        raise ValueError("s_method_name is not a property which can be set, only implemented as subclass")


# =============================================================================================
# IMPLEMENTED METHODS
# =============================================================================================

# =============================================================================================
# Implemented Measurement Methods
# =============================================================================================

class DummyMethods(MeasurementMethod):
    """
    Dummy class for testing functions
    """
    def error_model(self, *args):
        return 1

    def compute_property(self,
                         sim_data,
                         thermodynamic_state: thermodynamics.ThermodynamicState
                         ):
        dummy_data = np.zeros([10, 3]) * unit.nanometer
        dummy_observable = 1
        return dummy_observable, dummy_data

    def log_likelihood_data(self,
                            experimental_property,
                            calculated_property,
                            thermodynamic_state: thermodynamics.ThermodynamicState
                            ):
        dummy_likliehood = 1
        return dummy_likliehood


class DensityMeasurement(DummyMethods):

    @property
    def prop_name(self) -> str:
        return "Mass density, kg/m3"


class VibratingTube(DensityMeasurement):

    @property
    def method_name(self) -> str:
        return "Vibrating tube method"


class Pycnometer(DensityMeasurement):

    @property
    def method_name(self) -> str:
        return "Pycnometric method"


class ExcessEnthalpyMeasurements(DummyMethods):

    @property
    def prop_name(self) -> str:
        return "Excess molar enthalpy (molar enthalpy of mixing), kJ/mol"


class FlowCalorimiter(ExcessEnthalpyMeasurements):

    @property
    def method_name(self) -> str:
        return "Flow calorimetry"


class CalvetCalorimiter(ExcessEnthalpyMeasurements):

    @property
    def method_name(self) -> str:
        return "Calvet calorimetry"
