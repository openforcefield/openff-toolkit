#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
forcefield_labeler.py

Class structure mirroring forcefield.py but for simply determining what parameter numbers would be assigned to specified SMIRKS.

AUTHORS

John D. Chodera <john.chodera@choderalab.org>
David L. Mobley <dmobley@mobleylab.org>

Baseed on simtk.openmm.app.forcefield written by Peter Eastman.

"""
#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import sys
import string

#import xml.etree.ElementTree as etree
import lxml.etree as etree

import os
import math
import copy
import re
import numpy
import random

import openeye.oechem
import openeye.oeomega
import openeye.oequacpac

from openeye import oechem

from simtk import openmm, unit

import time

import networkx

#=============================================================================================
# PRIVATE SUBROUTINES
#=============================================================================================

def _convertParameterToNumber(param):
    """
    Convert parameter to OpenMM units.
    """
    if unit.is_quantity(param):
        if param.unit.is_compatible(unit.bar):
            return param / unit.bar
        return param.value_in_unit_system(unit.md_unit_system)
    return float(param)

#=============================================================================================
# FORCEFIELD LABELER
#=============================================================================================

