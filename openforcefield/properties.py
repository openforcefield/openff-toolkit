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


# =============================================================================================
# Property Enums
# =============================================================================================

@unique
class PropertyType(IntFlag):
    """An enum describing the type of property that was collected.
    """

    Undefined          = 0x00
    Density            = 0x01
    DielectricConstant = 0x02

    def __str__(self):
        """
        Returns
        ---
        str
            A string representation of the PropertyType enum
        """
        phases = '|'.join([phase.name for phase in PropertyType if self & phase])
        return phases

    def __repr__(self):
        """
        Returns
        ---
        str
            A string representation of the PropertyType enum
        """
        return str(self)


@unique
class PropertyPhase(IntFlag):
    """An enum describing the phase a property was collected in.
    """

    Undefined = 0x00
    Solid     = 0x01
    Liquid    = 0x02
    Gas       = 0x04

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
class CalculationFidelity(Enum):
    """An enum describing the fidelity at which a calculation was performed.
    """

    SurrogateModel =   1
    Reweighting =      2
    DirectSimulation = 3


# =============================================================================================
# Property Source
# =============================================================================================

class Source(object):
    """Container class for DOI and reference for a given observable

    This class contains either the DOI and/or the reference, but must contain at
    least one as the observable must have a source, even if it was measured in lab.

    Parameters
    ----------
    doi : str or None, default None
        The DOI for the source, preferred way to identify for source
    reference : str
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

class PhysicalProperty(object):
    """Represents the value of any physical property and it's uncertainty.

    It additionally stores the thermodynamic state at which the property
    was collected, the phase it was collected in, and information about
    the composition of the observed system.
    """
    def __init__(self):

        self._thermodynamic_state = None

        self.type = PropertyType.Undefined
        self.phase = PropertyPhase.Undefined

        self.substance = None

        self.value = None
        self.uncertainty = None

    @property
    def temperature(self):
        """simtk.unit.Quantity or None: The temperature at which the property was collected."""
        return None if self._thermodynamic_state is None else self._thermodynamic_state.temperature

    @property
    def pressure(self):
        """simtk.unit.Quantity or None: The pressure at which the property was collected."""
        return None if self._thermodynamic_state is None else self._thermodynamic_state.pressure

    @property
    def thermodynamic_state(self):
        """ThermodynamicState or None: The thermodynamic state at which the property was collected."""
        return self._thermodynamic_state

    @thermodynamic_state.setter
    def thermodynamic_state(self, thermodynamic_state):
        self._thermodynamic_state = thermodynamic_state

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


class MeasuredPhysicalProperty(PhysicalProperty):
    """A MeasuredPhysicalProperty is the implementation of property measured
    experimentally.

    It tracks the source from which the property was extracted (e.g ThermoML).
    """

    def __init__(self):

        super().__init__()

        self.method_name = None
        self._source = None

    @property
    def doi(self):
        """str or None: The doi of the article in which this property was published"""
        return None if self._source is None else self._source.doi

    @property
    def reference(self):
        """str or None: The reference to the location of where this property was published"""
        return None if self._source is None else self._source.reference

    @property
    def source(self):
        """Source or None: The source of this published property"""
        return self._source

    @source.setter
    def source(self, source):
        self._source = source


class CalculatedPhysicalProperty(PhysicalProperty):
    """
    A CalculatedPhysicalProperty represents a property estimated by either
    a surrogate model, reweighting existing data, or directly from a simulation.
    """

    def __init__(self):

        super().__init__()
        self.fidelity = CalculationFidelity.DirectSimulation
