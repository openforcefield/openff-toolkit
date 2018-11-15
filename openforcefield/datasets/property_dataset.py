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
from openforcefield.properties import PropertyType


class PhysicalPropertyDataSet(object):
    """A data set of physical property measurements.

    Implements the container API.

    Properties
    ----------

    """

    def __init__(self):

        self._measured_properties = []
        self._sources = []

    @property
    def measured_properties(self):
        return self._measured_properties

    @property
    def sources(self):
        return self._sources

    def merge(self, data_set):

        # TODO: Do we need to check wether merging the same data set here?
        self._measured_properties.extend(data_set.measured_properties)
        self._sources.extend(data_set.sources)

    def filter_by_function(self, function):

        # This works for now - if we wish to be able to undo a filter then
        # a 'filtered' list needs to be maintained separately to the main list.
        self._measured_properties = list(filter(function, self._measured_properties))

    def filter_by_properties(self, property_types):

        def filter_function(x): return x.type & property_types
        self.filter_by_function(filter_function)

    def filter_by_phases(self, phases):

        def filter_function(x): return x.phase & phases
        self.filter_by_function(filter_function)
