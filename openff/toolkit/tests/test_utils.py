#!/usr/bin/env python

# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Tests for utility methods

"""

import os

import pytest
from openff.units import unit

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
