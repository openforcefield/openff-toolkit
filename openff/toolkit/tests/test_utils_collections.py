"""
Tests for custom collections classes.

"""

import copy

import pytest

from openff.toolkit.utils.collections import ValidatedDict, ValidatedList


class TestValidatedMixin:
    def test_pickle(self):
        """Test pickle roundtripping"""
        # TODO: implement this test
        pass


class TestValidatedList(TestValidatedMixin):
    """Test suite for the ValidatedList class."""

    def test_validators(self):
        """Validators of ValidatedList are called correctly."""

        def is_positive(value):
            if value <= 0:
                raise TypeError("value is not positive")

        # The constructor should check all elements.
        with pytest.raises(TypeError, match="value is not positive"):
            ValidatedList([1, -2], validator=is_positive)

        vlist = ValidatedList([1, 2, 3], validator=is_positive)

        # __setitem__()
        with pytest.raises(TypeError, match="value is not positive"):
            vlist[2] = -1
        with pytest.raises(TypeError, match="value is not positive"):
            vlist[0:2] = [2, 3, -1]

        # append()
        with pytest.raises(TypeError, match="value is not positive"):
            vlist.append(-4)

        # extend() and __iadd__()
        with pytest.raises(TypeError, match="value is not positive"):
            vlist.extend([6, -1])
        with pytest.raises(TypeError, match="value is not positive"):
            vlist += [6, -1]

        # insert()
        with pytest.raises(TypeError, match="value is not positive"):
            vlist.insert(1, -3)

    def test_converters(self):
        """Custom converters of ValidatedList are called correctly."""
        # All elements are converted on construction.
        vlist = ValidatedList([1, 2.0, "3"], converter=int)
        assert vlist == [1, 2, 3]

        # __setitem__()
        vlist[2] = "4"
        assert vlist[2] == 4
        vlist[0:3] = ["2", "3", 4]
        assert vlist == [2, 3, 4]

        # append()
        vlist.append("5")
        assert vlist[3] == 5

        # extend() and __iadd__()
        vlist.extend([6, "7"])
        assert vlist[5] == 7
        vlist += ["8", 9]
        assert vlist[6] == 8

        # insert()
        vlist.insert(5, "10")
        assert vlist[5] == 10

    def test_validators_and_converters(self):
        """Custom converters of ValidatedList are called correctly."""

        def is_positive(value):
            if value <= 0:
                raise TypeError("value is not positive")

        # Validators are run after converters.
        vlist = ValidatedList([1, 2, -3], converter=abs, validator=is_positive)
        assert vlist == [1, 2, 3]

        # __setitem__
        vlist[2] = -1
        assert vlist[2] == 1
        with pytest.raises(TypeError, match="value is not positive"):
            vlist[2] = 0
        vlist[0:3] = [2, 3, -1]
        assert vlist == [2, 3, 1]
        with pytest.raises(TypeError, match="value is not positive"):
            vlist[0:3] = [2, 3, 0]

        # append()
        vlist.append(-4)
        assert vlist[-1] == 4
        with pytest.raises(TypeError, match="value is not positive"):
            vlist.append(0)

        # extend() and __iadd__()
        vlist.extend([6, -1])
        assert vlist[-2:] == [6, 1]
        with pytest.raises(TypeError, match="value is not positive"):
            vlist.extend([6, 0])
        vlist += [6, -2]
        assert vlist[-2:] == [6, 2]
        with pytest.raises(TypeError, match="value is not positive"):
            vlist += [6, 0]

        # insert()
        vlist.insert(1, -3)
        assert vlist[1] == 3
        with pytest.raises(TypeError, match="value is not positive"):
            vlist.insert(1, 0)

    def test_multiple_converters(self):
        """Multiple converters of ValidatedList are called in order."""
        vlist = ValidatedList([1, 2, -3], converter=[abs, str])
        assert vlist == ["1", "2", "3"]

    def test_multiple_validators(self):
        """Multiple converters of ValidatedList are called in order."""

        def is_positive(value):
            if value <= 0:
                raise TypeError("value must be positive")

        def is_odd(value):
            if value % 2 == 0:
                raise TypeError("value must be odd")

        with pytest.raises(TypeError, match="value must be positive"):
            ValidatedList([-1, -3], validator=[is_positive, is_odd])

        with pytest.raises(TypeError, match="value must be odd"):
            ValidatedList([2, 4], validator=[is_positive, is_odd])

    def test_copy(self):
        """A copy of a ValidatedList returns another ValidatedList."""
        vlist = ValidatedList([1, 2, 3])
        assert isinstance(vlist.copy(), ValidatedList)
        assert isinstance(copy.copy(vlist), ValidatedList)
        assert isinstance(copy.deepcopy(vlist), ValidatedList)

    def test_slice(self):
        """A slice of a ValidatedList returns another ValidatedList."""
        vlist = ValidatedList([1, 2, 3])
        assert isinstance(vlist[:2], ValidatedList)


class TestValidatedDict(TestValidatedMixin):
    """Test suite for the ValidatedDict class."""

    def test_validators(self):
        """Validators of ValidatedDict are called correctly."""

        def is_positive(value):
            if value <= 0:
                raise TypeError("value is not positive")

        # The constructor should check all elements.
        with pytest.raises(TypeError, match="value is not positive"):
            ValidatedDict({"a": 1, "b": -2}, validator=is_positive)

        d = ValidatedDict({"x": 1, "y": 2, "z": 3}, validator=is_positive)

        # __setitem__()
        with pytest.raises(TypeError, match="value is not positive"):
            d["me!"] = -1

        # update()
        with pytest.raises(TypeError, match="value is not positive"):
            d.update((("moe", 2), ("larry", 3), ("curly", -4)))
        with pytest.raises(TypeError, match="value is not positive"):
            d.update({"moe": 2, "larry": 3, "curly": -4})

    def test_converters(self):
        """Custom converters of ValidatedDict are called correctly."""
        # All elements are converted on construction.
        d = ValidatedDict({"x": 1, "y": 2.0, "z": "3"}, converter=int)
        assert d == {"x": 1, "y": 2, "z": 3}

        # __setitem__()
        d["y"] = "4"
        assert d["y"] == 4
        d["w"] = "-20"
        assert d["w"] == -20

        # update()
        d.update({6: "7"})
        assert d[6] == 7

        d.update(((6, "8"), ("shemp", "-30")))
        assert d[6] == 8
        assert d["shemp"] == -30

    def test_validators_and_converters(self):
        """Custom converters of ValidatedDict are called correctly."""

        def is_positive(value):
            if value <= 0:
                raise TypeError("value is not positive")

        # Validators are run after converters.
        d = ValidatedDict(
            {"a": 1, "b": 2, "c": -3}, converter=abs, validator=is_positive
        )
        assert d == {"a": 1, "b": 2, "c": 3}

        # __setitem__
        d[2] = -1
        assert d[2] == 1
        with pytest.raises(TypeError, match="value is not positive"):
            d[2] = 0

        # update()
        d.update({"x": 6, "y": -1})
        assert d["y"] == 1
        with pytest.raises(TypeError, match="value is not positive"):
            d.update({6: 0})

        d.update((("x", 6), ("y", -1)))
        assert d["y"] == 1

    def test_multiple_converters(self):
        """Multiple converters of ValidatedDict are called in order."""
        d = ValidatedDict({"u": 1, "v": 2, "w": -3}, converter=[abs, str])
        assert d == {"u": "1", "v": "2", "w": "3"}

    def test_multiple_validators(self):
        """Multiple converters of ValidatedDict are called in order."""

        def is_positive(value):
            if value <= 0:
                raise TypeError("value must be positive")

        def is_odd(value):
            if value % 2 == 0:
                raise TypeError("value must be odd")

        with pytest.raises(TypeError, match="value must be positive"):
            ValidatedDict({"first": -1, "second": -3}, validator=[is_positive, is_odd])

        with pytest.raises(TypeError, match="value must be odd"):
            ValidatedDict({"first": 2, "second": 4}, validator=[is_positive, is_odd])

    def test_copy(self):
        """A copy of a ValidatedDict returns another ValidatedDict."""
        d = ValidatedDict({"a": 1, "b": 2, "c": 3})
        assert isinstance(d.copy(), ValidatedDict)
        assert isinstance(copy.copy(d), ValidatedDict)
        assert isinstance(copy.deepcopy(d), ValidatedDict)
