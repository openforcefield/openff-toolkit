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

    Undefined   = 0x00
    MassDensity = 0x01

    def __str__(self):
        return self.name

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
        return self.name

    def __repr__(self):
        return str(self)
