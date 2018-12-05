#!/usr/bin/env python

# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Groups API.

Authors
-------
* Simon Boothroyd <simon.boothroyd@choderalab.org>
"""


# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

import copy
import logging

import mdtraj

import arch.bootstrap

import uuid

import numpy as np

from os import path

from enum import Enum

from pymbar import timeseries

from openeye import oechem, oeomega

from openforcefield.utils import packmol
from openforcefield.utils.exceptions import XmlNodeMissingException
from openforcefield.typing.engines import smirnoff

from simtk import openmm, unit
from simtk.openmm import app


# =============================================================================================
# Groups
# =============================================================================================


