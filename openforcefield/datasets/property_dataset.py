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

from enum import Enum, unique


# =============================================================================================
# DATASET
# =============================================================================================

class PhysicalPropertyDataSet(object):
    """A data set of physical property measurements.

    Implements the container API.

    Properties
    ----------

    """

    def __init__(self):

        self._measured_properties = []

    @property
    def measured_properties(self):
        return self._measured_properties

    def filter_by_property(self, property):

        pass

    def filter_by_phase(self, property):

        pass
