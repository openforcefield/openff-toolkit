#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Tests for utility methods

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import ast
import os

import pytest
from simtk import unit

from openforcefield import utils


#=============================================================================================
# TESTS
#=============================================================================================

def test_subclasses():
    """Test that all subclasses (and descendents) are correctly identified by all_subclasses()"""
    class Foo:
        pass
    class FooSubclass1(Foo):
        pass
    class FooSubclass2(Foo):
        pass
    class FooSubSubclass(FooSubclass1):
        pass

    subclass_names = [ cls.__name__ for cls in utils.all_subclasses(Foo) ]
    assert set(subclass_names) == set(['FooSubclass1', 'FooSubclass2', 'FooSubSubclass'])

def test_temporary_cd():
    """Test temporary_cd() context manager"""
    initial_dir = os.getcwd()
    temporary_dir = '/'
    from openforcefield.utils import temporary_cd
    with temporary_cd(temporary_dir):
        assert os.getcwd() == temporary_dir
    assert os.getcwd() == initial_dir

def test_temporary_directory():
    """Test temporary_directory() context manager"""
    from openforcefield.utils import temporary_directory
    initial_dir = os.getcwd()
    with temporary_directory() as tmp_dir:
        # Make sure the temporary directory is not the current directory
        assert tmp_dir != initial_dir

        # Make sure the temporary directory is writeable
        output_filename = os.path.join(tmp_dir, 'test.out')
        outfile = open(output_filename,'w')
        outfile.write('test')
        outfile.close()

        # Make sure the file we created exists
        assert os.path.exists(output_filename)

    # Make sure the directory has been cleaned up
    assert not os.path.exists(tmp_dir), "Temporary directory was not automatically cleaned up."
    assert not os.path.exists(output_filename), "Temporary directory was not automatically cleaned up."

def test_get_data_file_path():
    """Test get_data_file_path()"""
    from openforcefield.utils import get_data_file_path
    filename = get_data_file_path('test_forcefields/tip3p.offxml')
    assert os.path.exists(filename)


@pytest.mark.parametrize('unit_string,expected_unit',[
    ('kilocalories_per_mole', unit.kilocalories_per_mole),
    ('kilocalories_per_mole/angstrom**2', unit.kilocalories_per_mole/unit.angstrom**2),
    ('joule/(mole * nanometer**2)', unit.joule/(unit.mole * unit.nanometer**2)),
    ('picosecond**(-1)', unit.picosecond**(-1)),
    ('300.0 * kelvin', 300*unit.kelvin),
    ('1 * kilojoule + 500 * joule', 1.5*unit.kilojoule),
    ('1 / meter', 1.0 / unit.meter)
])
def test_ast_eval(unit_string, expected_unit):
    """Test that _ast_eval() correctly parses string quantities."""
    from openforcefield.utils.utils import _ast_eval
    ast_root_node = ast.parse(unit_string, mode='eval').body
    parsed_units = _ast_eval(ast_root_node)
    assert parsed_units == expected_unit
