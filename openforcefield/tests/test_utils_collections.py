#!/usr/bin/env python

# =====================================================================
# MODULE DOCSTRING
# =====================================================================

"""
Tests for custom collections classes.

"""

# =====================================================================
# GLOBAL IMPORTS
# =====================================================================

import pytest

from openforcefield.utils.collections import ValidatedList


# =====================================================================
# Test Callbackable class
# =====================================================================

class TestValidatedList:
    """Test suite for the ValidatedList class."""

    def test_validators(self):
        """Validators of ValidatedList are called correctly."""
        def is_positive(value):
            if value <= 0:
                raise TypeError('value is not positive')

        # The constructor should check all elements.
        with pytest.raises(TypeError, match='value is not positive'):
            ValidatedList([1, -2], validator=is_positive)

        l = ValidatedList([1, 2, 3], validator=is_positive)

        # __setitem__()
        with pytest.raises(TypeError, match='value is not positive'):
            l[2] = -1
        with pytest.raises(TypeError, match='value is not positive'):
            l[0:2] = [2, 3, -1]

        # append()
        with pytest.raises(TypeError, match='value is not positive'):
            l.append(-4)

        # extend() and __iadd__()
        with pytest.raises(TypeError, match='value is not positive'):
            l.extend([6, -1])
        with pytest.raises(TypeError, match='value is not positive'):
            l += [6, -1]

        # insert()
        with pytest.raises(TypeError, match='value is not positive'):
            l.insert(1, -3)

    def test_converters(self):
        """Custom converters of ValidatedList are called correctly."""
        # All elements are converted on construction.
        l = ValidatedList([1, 2.0, '3'], converter=int)
        assert l == [1, 2, 3]

        # __setitem__()
        l[2] = '4'
        assert l[2] == 4
        l[0:3] = ['2', '3', 4]
        assert l == [2, 3, 4]

        # append()
        l.append('5')
        assert l[3] == 5

        # extend() and __iadd__()
        l.extend([6, '7'])
        assert l[5] == 7
        l += ['8', 9]
        assert l[6] == 8

        # insert()
        l.insert(5, '10')
        assert l[5] == 10

    def test_validators_and_converters(self):
        """Custom converters of ValidatedList are called correctly."""
        def is_positive(value):
            if value <= 0:
                raise TypeError('value is not positive')

        # Validators are run after converters.
        l = ValidatedList([1, 2, -3], converter=abs, validator=is_positive)
        assert l == [1, 2, 3]

        # __setitem__
        l[2] = -1
        assert l[2] == 1
        with pytest.raises(TypeError, match='value is not positive'):
            l[2] = 0
        l[0:3] = [2, 3, -1]
        assert l == [2, 3, 1]
        with pytest.raises(TypeError, match='value is not positive'):
            l[0:3] = [2, 3, 0]

        # append()
        l.append(-4)
        assert l[-1] == 4
        with pytest.raises(TypeError, match='value is not positive'):
            l.append(0)

        # extend() and __iadd__()
        l.extend([6, -1])
        assert l[-2:] == [6, 1]
        with pytest.raises(TypeError, match='value is not positive'):
            l.extend([6, 0])
        l += [6, -2]
        assert l[-2:] == [6, 2]
        with pytest.raises(TypeError, match='value is not positive'):
            l += [6, 0]

        # insert()
        l.insert(1, -3)
        assert l[1] == 3
        with pytest.raises(TypeError, match='value is not positive'):
            l.insert(1, 0)
