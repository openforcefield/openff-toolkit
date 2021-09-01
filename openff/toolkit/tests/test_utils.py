#!/usr/bin/env python

# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Tests for utility methods

"""

import ast
import os

import pytest
from simtk import unit

# =============================================================================================
# TESTS
# =============================================================================================


def test_requires_package():
    """Test the @requires_package decorator"""
    from openff.toolkit.utils.utils import MissingDependencyError, requires_package

    @requires_package("numpy")
    def fn_installed():
        pass

    fn_installed()

    @requires_package("foobar")
    def fn_missing():
        pass

    with pytest.raises(MissingDependencyError, match="foobar"):
        fn_missing()


def test_subclasses():
    """Test that all subclasses (and descendents) are correctly identified by all_subclasses()"""

    from openff.toolkit.utils.utils import all_subclasses

    class Foo:
        pass

    class FooSubclass1(Foo):
        pass

    class FooSubclass2(Foo):
        pass

    class FooSubSubclass(FooSubclass1):
        pass

    subclass_names = [cls.__name__ for cls in all_subclasses(Foo)]
    assert set(subclass_names) == set(
        ["FooSubclass1", "FooSubclass2", "FooSubSubclass"]
    )


def test_temporary_cd():
    """Test temporary_cd() context manager"""
    initial_dir = os.getcwd()
    temporary_dir = "/"
    from openff.toolkit.utils import temporary_cd

    with temporary_cd(temporary_dir):
        assert os.getcwd() == temporary_dir
    assert os.getcwd() == initial_dir


def test_get_data_file_path():
    """Test get_data_file_path()"""
    from openff.toolkit.utils import get_data_file_path

    filename = get_data_file_path("test_forcefields/tip3p.offxml")
    assert os.path.exists(filename)


@pytest.mark.parametrize(
    "unit_string,expected_unit",
    [
        ("kilocalories_per_mole", unit.kilocalories_per_mole),
        (
            "kilocalories_per_mole/angstrom**2",
            unit.kilocalories_per_mole / unit.angstrom ** 2,
        ),
        ("joule/(mole * nanometer**2)", unit.joule / (unit.mole * unit.nanometer ** 2)),
        ("picosecond**(-1)", unit.picosecond ** (-1)),
        ("300.0 * kelvin", 300 * unit.kelvin),
        ("1 * kilojoule + 500 * joule", 1.5 * unit.kilojoule),
        ("1 / meter", 1.0 / unit.meter),
    ],
)
def test_ast_eval(unit_string, expected_unit):
    """Test that _ast_eval() correctly parses string quantities."""
    from openff.toolkit.utils.utils import _ast_eval

    ast_root_node = ast.parse(unit_string, mode="eval").body
    parsed_units = _ast_eval(ast_root_node)
    assert parsed_units == expected_unit


def test_dimensionless_units():
    from openff.toolkit.utils.utils import string_to_unit, unit_to_string

    assert string_to_unit("dimensionless") == unit.dimensionless

    unit_string = unit_to_string(unit.dimensionless)
    unit_value = string_to_unit(unit_string)

    assert unit_value == unit.dimensionless


def test_sort_smirnoff_dict():
    from collections import OrderedDict

    from openff.toolkit.typing.engines.smirnoff import ForceField
    from openff.toolkit.utils.utils import sort_smirnoff_dict

    forcefield = ForceField("test_forcefields/test_forcefield.offxml")
    smirnoff_dict = forcefield._to_smirnoff_data()

    # Ensure data is not created or destroyed
    # dict.__eq__ does not check order
    assert smirnoff_dict == OrderedDict(
        sort_smirnoff_dict(forcefield._to_smirnoff_data())
    )


def test_import_message_exception_raises_warning(caplog):
    # TODO: Remove when removing MessageException
    msg = "DEPRECATED and will be removed in a future release of the OpenFF Toolkit"

    with pytest.warns(DeprecationWarning, match=msg):
        from openff.toolkit.utils.exceptions import MessageException

    with pytest.warns(None) as rec:
        from openff.toolkit.utils.exceptions import SMILESParseError

    assert len(rec) == 0
