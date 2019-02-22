#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Test classes and function in module openforcefield.typing.engines.smirnoff.io.

"""


#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import ast

import pytest
from simtk import unit

from openforcefield.typing.engines.smirnoff.io import (
    _ast_unit_eval,
    _extract_attached_units,
)


#=============================================================================================
# QUANTITY PARSING UTILITIES
#=============================================================================================

@pytest.mark.parametrize('unit_string,expected_unit',[
    ('kilocalories_per_mole', unit.kilocalories_per_mole),
    ('kilocalories_per_mole/angstrom**2', unit.kilocalories_per_mole/unit.angstrom**2),
    ('joule/(mole * nanometer**2)', unit.joule/(unit.mole * unit.nanometer**2)),
    ('picosecond**(-1)', unit.picosecond**(-1)),
    ('300.0 * kelvin', 300*unit.kelvin),
    ('1 * kilojoule + 500 * joule', 1.5*unit.kilojoule),
    ('1 / meter', 1.0 / unit.meter)
])
def test_ast_unit_eval(unit_string, expected_unit):
    """Test that _ast_unit_eval() correctly parses string quantities."""
    ast_root_node = ast.parse(unit_string, mode='eval').body
    parsed_units = _ast_unit_eval(ast_root_node)
    assert parsed_units == expected_unit


@pytest.mark.parametrize('attributes,expected',[
    ({'not_parsed': 'blabla'},
        {}),
    ({'not_parsed': 'blabla', 'attr_unit': 'angstrom/femtosecond'},
        {'attr': unit.angstrom/unit.femtosecond}),
    ({'not_parsed': 'blabla', 'attr1_unit': 'meter', 'attr2_unit': 'kilojoule_per_mole'},
        {'attr1': unit.meter, 'attr2': unit.kilojoule_per_mole}),
])
def test_extract_attached_units(attributes, expected):
    """Test that _extract_attached_units() correctly parses the correct."""
    assert _extract_attached_units(attributes) == expected


@pytest.mark.parametrize('attributes',[
    {'attr_unit': '300.0 * kelvin'},
    {'attr_unit': '1 / picosecond'}
])
def test_extract_attached_units_raises(attributes):
    """Test that _extract_attached_units() raises an error when a quantity is specified instead of a unit."""
    with pytest.raises(ValueError, match='associated to a quantity rather than only units'):
        _extract_attached_units(attributes)
