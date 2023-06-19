"""
Tests for utility methods

"""

import os

from openff.units import unit


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


def test_object_to_quantity_accepts_openmm():
    import openmm

    from openff.toolkit.utils.utils import object_to_quantity

    val = object_to_quantity(1.0 * openmm.unit.angstrom)
    assert val == 1.0 * unit.angstrom
    val = object_to_quantity(1.0 * openmm.unit.nanometer)
    assert val == 10.0 * unit.angstrom
    val = object_to_quantity(
        [
            2.0 * openmm.unit.angstrom,
            2.0 * openmm.unit.nanometer,
            2.0 * openmm.unit.dimensionless,
        ]
    )
    assert val == [2.0 * unit.angstrom, 20.0 * unit.angstrom, 2 * unit.dimensionless]


def test_sort_smirnoff_dict():
    from collections import OrderedDict

    from openff.toolkit.typing.engines.smirnoff import ForceField
    from openff.toolkit.utils.utils import get_data_file_path, sort_smirnoff_dict

    forcefield = ForceField(
        get_data_file_path("test_forcefields/test_forcefield.offxml")
    )
    smirnoff_dict = forcefield._to_smirnoff_data()

    # Ensure data is not created or destroyed
    # dict.__eq__ does not check order
    assert smirnoff_dict == OrderedDict(
        sort_smirnoff_dict(forcefield._to_smirnoff_data())
    )
