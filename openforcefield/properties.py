#!/usr/bin/env python

# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Properties base API.

Authors
-------
* Simon Boothroyd <simon.boothroyd@choderalab.org>

"""
# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

from enum import IntFlag, unique


# =============================================================================================
# PropertyType
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


# =============================================================================================
# PropertyPhase
# =============================================================================================

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
