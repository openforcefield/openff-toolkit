#!/usr/bin/env python

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

TODO
----

"""
# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

from simtk import unit
from typing import Optional


# =============================================================================================
# THERMODYNAMIC STATE
# =============================================================================================

class ThermodynamicState(object):
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

    def __init__(self,
                 temperature: Optional[unit.Quantity]=None,
                 pressure: Optional[unit.Quantity]=None):
        """
        Initialize the thermodynamic state.

        Parameters
        ----------
        temperature : simtk.unit.Quantity compatible with 'kelvin', optional, default=None
           The temperature for a system with constant temperature
        pressure : simtk.unit.Quantity compatible with 'atmospheres', optional, default=None
           The pressure for constant-pressure systems (default: None)

        """
        # Init values
        self._temperature = None
        self._pressure = None

        # Assign values using the setter wrappers
        # Check that units are compatible.
        # These throw a TypeException if this conversion cannot be done.
        self.temperature = temperature
        self.pressure = pressure

    @property
    def temperature(self):
        return self._temperature

    @temperature.setter
    def temperature(self, value):
        if value is not None:
            value = value.in_units_of(unit.kelvin)
        self._temperature = value

    @property
    def pressure(self):
        return self._pressure

    @pressure.setter
    def pressure(self, value):
        if value is not None:
            value = value.in_units_of(unit.atmospheres)
        self._pressure = value

    # TODO: This method seems to be in the wrong place?
    # def is_compatible_with(self, state):
    #     """
    #     Determine whether another state is in the same thermodynamic ensemble (e.g. NVT, NPT).
    #
    #     Parameters
    #     ----------
    #     state : ThermodynamicState
    #         Thermodynamic state whose compatibility is to be determined
    #
    #     Returns
    #     -------
    #     is_compatible : bool
    #         True if 'state' is of the same ensemble (e.g. both NVT, both NPT), False otherwise
    #
    #     Examples
    #     --------
    #
    #     Create NVT and NPT states.
    #
    #     >>> from simtk import unit
    #     >>> from openmmtools import testsystems
    #     >>> testsystem = testsystems.LennardJonesCluster()
    #     >>> [system, positions] = [testsystem.system, testsystem.positions]
    #     >>> nvt_state = ThermodynamicState(system=system, temperature=100.0*unit.kelvin)
    #     >>> npt_state = ThermodynamicState(system=system, temperature=100.0*unit.kelvin, pressure=1.0*unit.atmospheres)
    #
    #     Test compatibility.
    #
    #     >>> test1 = nvt_state.is_compatible_with(nvt_state)
    #     >>> test2 = nvt_state.is_compatible_with(npt_state)
    #     >>> test3 = npt_state.is_compatible_with(nvt_state)
    #     >>> test4 = npt_state.is_compatible_with(npt_state)
    #
    #     """
    #
    #     is_compatible = True
    #
    #     # Make sure systems have the same number of atoms.
    #     if ((self.system != None) and (state.system != None)):
    #         if (self.system.getNumParticles() != state.system.getNumParticles()):
    #             is_compatible = False
    #
    #     # Make sure other terms are defined for both states.
    #     # TODO: Use introspection to get list of parameters?
    #     for parameter in ['temperature', 'pressure']:
    #         if (parameter in dir(self)) is not (parameter in dir(state)):
    #             # parameter is not shared by both states
    #             is_compatible = False
    #
    #     return is_compatible

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
