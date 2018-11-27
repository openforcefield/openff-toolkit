#!/usr/bin/env python

# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Data sets API.

Authors
-------
* John D. Chodera <john.chodera@choderalab.org>
* Levi N. Naden <levi.naden@choderalab.org>
* Simon Boothroyd <simon.boothroyd@choderalab.org>

"""

# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================


# =============================================================================================
# DATASET
# =============================================================================================

class PhysicalPropertyDataSet(object):
    """
    A data set of physical property measurements / calculations.
    Contains functionality for merging multiple data sets and
    filtering existing ones.
    """

    def __init__(self):
        """
        Constructs a new PhysicalPropertyDataSet
        """
        self._properties = {}
        self._sources = []

    @property
    def properties(self):
        """
        dict(str, list(PhysicalProperty)): The property list which
        takes a substance as the key.
        """
        return self._properties

    @property
    def sources(self):
        """list of Source: The list of sources from which the properties were gathered"""
        return self._sources

    def merge(self, data_set):
        """Merge another data set into the current one.

        Parameters
        ----------
        data_set : PhysicalPropertyDataSet
            The secondary data set to merge into this one.
        """
        # TODO: Do we need to check whether merging the same data set here?
        for substance_hash in data_set.properties:

            if substance_hash not in self._properties:
                self._properties[substance_hash] = []

            self._properties[substance_hash].extend(
                data_set.properties[substance_hash])

        self._sources.extend(data_set.sources)

    def filter_by_function(self, filter_function):
        """Filter the data set using a given filter function.

        Parameters
        ----------
        filter_function : lambda
            The filter function.
        """

        # This works for now - if we wish to be able to undo a filter then
        # a 'filtered' list needs to be maintained separately to the main list.
        for substance_hash in self._properties:

            self._properties[substance_hash] = list(filter(
                filter_function, self._properties[substance_hash]))

    def filter_by_properties(self, property_types):
        """Filter the data set based on the type of property (e.g Density).

        Parameters
        ----------
        property_types : PropertyType
            The type of property which should be retained.
        """
        def filter_function(x): return x.type & property_types
        self.filter_by_function(filter_function)

    def filter_by_phases(self, phases):
        """Filter the data set based on the phase of the property (e.g Liquid).

        Parameters
        ----------
        phases : PropertyPhase
            The phase of property which should be retained.
        """
        def filter_function(x): return x.phase & phases
        self.filter_by_function(filter_function)

    def filter_by_temperature(self, min_temperature, max_temperature):
        """Filter the data set based on a minimum and maximum temperature.

        Parameters
        ----------
        min_temperature : float
            The minimum temperature.
        max_temperature : float
            The maximum temperature.
        """
        def filter_function(x):
            return min_temperature <= x.thermodynamic_state.temperature <= max_temperature

        self.filter_by_function(filter_function)


class CalculatedPropertySet(object):
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
        """
        dict(str, list(PhysicalProperty)): The calculated property list which
        takes a substance as the key.
        """
        return self._properties

    @property
    def parameter_set(self):
        """
        ForceField: The parameters used to calculate the property.
        """
        return self._parameter_set
