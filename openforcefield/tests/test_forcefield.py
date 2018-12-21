#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Tests for forcefield class

"""


#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

from functools import partial
from unittest import TestCase
from openforcefield import utils
from simtk import unit
import numpy as np
from numpy.testing import assert_almost_equal

import pytest
from openforcefield.utils.toolkits import ToolkitWrapper, OpenEyeToolkitWrapper, RDKitToolkitWrapper, AmberToolsToolkitWrapper, ToolkitRegistry

from openforcefield.utils import get_data_filename

from openforcefield.topology.molecule import Molecule
from openforcefield.typing.engines.smirnoff import ForceField


#=============================================================================================
# TESTS
#=============================================================================================

class TestForceField(TestCase):
    """Test the ForceField class"""



    def test_create_forcefield(self):
        """Test empty constructor"""
        filename = get_data_filename('forcefield/smirnoff99Frosst.offxml')
        forcefield = ForceField(filename)
        raise Exception(forcefield._parameter_handlers['Angles']._parameters[10].__dict__)


# from_filename
# from xml_string
# from_xml_bytes
# from_url
# get_new_parameterhandler
# get_existing_parameterhandler
# get_parameter
# add_parameter
# add_parameter_fractional_bondorder
# create_force_fractional_bondorder
# store_cosmetic_attribs
# write_cosmetic_attribs
# store_cosmetic_elements
# write_cosmetic_elements
# add_handler_with_incompatible_kwargs (for example different scale14 vals)