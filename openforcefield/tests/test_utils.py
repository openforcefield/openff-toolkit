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

from functools import partial
from unittest import TestCase
import parmed
from openforcefield import utils
import os

#=============================================================================================
# TESTS
#=============================================================================================

class TestUtils(TestCase):
    def test_subclasses(self):
        """Test that all subclasses (and descendents) are correctly identified by all_subclasses()"""
        class Foo(object):
            pass
        class FooSubclass1(Foo):
            pass
        class FooSubclass2(Foo):
            pass
        class FooSubSubclass(FooSubclass1):
            pass

        subclass_names = [ cls.__name__ for cls in utils.all_subclasses(Foo) ]
        assert set(subclass_names) == set(['FooSubclass1', 'FooSubclass2', 'FooSubSubclass'])

    def test_temporary_cd(self):
        """Test temporary_cd() context manager"""
        initial_dir = os.getcwd()
        temporary_dir = '/'
        from openforcefield.utils import temporary_cd
        with temporary_cd(temporary_dir):
            assert os.getcwd() == temporary_dir
        assert os.getcwd() == initial_dir

    def test_temporary_direcftory(self):
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

    def test_get_data_filename(self):
        """Test get_data_filename()"""
        from openforcefield.utils import get_data_filename
        filename = get_data_filename('forcefield/tip3p.offxml')
        assert os.path.exists(filename)
