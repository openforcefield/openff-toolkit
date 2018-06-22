from functools import partial
from unittest import TestCase
import parmed
from openforcefield import utils

class TestUtils(TestCase):
    def test_toolkits(self):
        """Test that cheminformatics toolkits are properly detected.
        """
        assert utils.RDKIT_INSTALLED != utils.RDKIT_UNAVAILABLE
        assert utils.OPENEYE_INSTALLED != utils.OPENEYE_UNAVAILABLE
        assert utils.TOOLKIT_INSTALLED == (utils.RDKIT_INSTALLED or utils.OPENEYE_INSTALLED)

    def test_subclasses(self):
        """Test that all subclasses (and descendents) are correctly identified.
        """
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
