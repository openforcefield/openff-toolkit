#!/usr/bin/env python

"""
Custom collections classes.

"""

__all__ = [
    'ValidatedList',
]


# =====================================================================
# GLOBAL IMPORTS
# =====================================================================

from collections import abc


# =====================================================================
# VALIDATED LIST
# =====================================================================

class ValidatedList(list):
    """A list that runs custom converter and validators when new elements are added.

    Multiple converters and validators can be assigned to the list. These
    are executed in the given order with converters run before validators.

    Validators must take the new element as the first argument and raise
    an exception if validation fails.

        validator(new_element) -> None

    Converters must also take the new element as the first argument, but
    they have to return the converted value.

        converter(new_element) -> converted_value

    Examples
    --------
    We can define validator and converter functions that are run on each
    element of the list.

    >>> def is_positive_validator(value):
    ...     if value <= 0:
    ...         raise TypeError('value must be positive')
    ...
    >>> vl = ValidatedList([1, -1], validator=is_positive_validator)
    Traceback (most recent call last):
    ...
    TypeError: value must be positive

    Multiple converters that are run before the validators can be specified.

    >>> vl = ValidatedList([-1, '2', 3.0], converter=[float, abs],
    ...                    validator=is_positive_validator)
    >>> vl
    [1.0, 2.0, 3.0]

    """

    def __init__(self, seq=(), converter=None, validator=None):
        """
        Initialize the list.

        Parameters
        ----------
        seq : Iterable
            A sequence of elements.
        converter : callable or List[callable]
            Functions that will be used to convert each new element of
            the list.
        validator : callable or List[callable]
            Functions that will be used to convert each new element of
            the list.

        """
        # Make sure converter and validator are always iterables.
        if not (converter is None or isinstance(converter, abc.Iterable)):
            converter = [converter]
        if not (validator is None or isinstance(validator, abc.Iterable)):
            validator = [validator]
        self._converters = converter
        self._validators = validator

        # Validate and convert the whole sequence.
        seq = self._convert_and_validate(seq)
        super().__init__(seq)

    def extend(self, iterable):
        iterable = self._convert_and_validate(iterable)
        super().extend(iterable)

    def append(self, p_object):
        p_object = self._convert_and_validate([p_object])[0]
        super().append(p_object)

    def insert(self, index, p_object):
        p_object = self._convert_and_validate([p_object])[0]
        super().insert(index, p_object)

    def __iadd__(self, other):
        other = self._convert_and_validate(other)
        return super().__iadd__(other)

    def __setitem__(self, key, value):
        if isinstance(key, slice):
            value = self._convert_and_validate(value)
        else:
            value = self._convert_and_validate([value])[0]
        super().__setitem__(key, value)

    def copy(self):
        # Make sure a shallow copy still returns a ValidatedList.
        return self.__class__(self)

    def __getitem__(self, item):
        # Make sure a slice returns a ValidatedList.
        if isinstance(item, slice):
            return self.__class__(super().__getitem__(item))
        return super().__getitem__(item)

    # This is needed for pickling. See https://github.com/openforcefield/openforcefield/issues/411
    # for more details.
    # TODO: Is there a cleaner way (getstate/setstate perhaps?) to allow FFs to be
    #       pickled?
    def __reduce__(self):
        return (__class__, ( list( self),), self.__dict__)


    def _convert_and_validate(self, seq):
        """Run all converters and the validator on the given sequence."""
        # Run all element converters.
        if self._converters is not None:
            for converter in self._converters:
                seq = [converter(element) for element in seq]

        # Run all element validators.
        if self._validators is not None:
            for validator in self._validators:
                for element in seq:
                    validator(element)
        return seq

class ValidatedListMapping(ValidatedList):
    """A list of mappings that runs custom converter and validators when new
    elements are added.

    Multiple converters and validators can be assigned to the list of mappings.
    These are executed in the given order with converters run before
    validators.

    Validators must take the new element as the first argument and raise
    an exception if validation fails.

        validator(new_element) -> None

    Converters must also take the new element as the first argument, but
    they have to return the converted value.

        converter(new_element) -> converted_value

    Examples
    --------
    We can define validator and converter functions that are run on each
    element of the list.

    >>> def is_positive_validator(value):
    ...     if value <= 0:
    ...         raise TypeError('value must be positive')
    ...
    >>> vl = ValidatedList([1, -1], validator=is_positive_validator)
    Traceback (most recent call last):
    ...
    TypeError: value must be positive

    Multiple converters that are run before the validators can be specified.

    >>> vl = ValidatedList([-1, '2', 3.0], converter=[float, abs],
    ...                    validator=is_positive_validator)
    >>> vl
    [1.0, 2.0, 3.0]

    """

    def __init__(self, seq=(), converter=None, validator=None):
        """
        Initialize the list.

        Parameters
        ----------
        seq : Iterable
            A sequence of elements.
        converter : callable or List[callable]
            Functions that will be used to convert each new element of
            the list.
        validator : callable or List[callable]
            Functions that will be used to convert each new element of
            the list.

        """
        # Make sure converter and validator are always iterables.
        if not (converter is None or isinstance(converter, abc.Iterable)):
            converter = [converter]
        if not (validator is None or isinstance(validator, abc.Iterable)):
            validator = [validator]
        self._converters = converter
        self._validators = validator

        # Validate and convert the whole sequence.
        seq = self._convert_and_validate(seq)
        super().__init__(seq)

    def extend(self, iterable):
        iterable = self._convert_and_validate(iterable)
        super().extend(iterable)

    def append(self, p_object):
        p_object = self._convert_and_validate([p_object])[0]
        super().append(p_object)

    def insert(self, index, p_object):
        p_object = self._convert_and_validate([p_object])[0]
        super().insert(index, p_object)

    def __iadd__(self, other):
        other = self._convert_and_validate(other)
        return super().__iadd__(other)

    def __setitem__(self, key, value):
        if isinstance(key, slice):
            value = self._convert_and_validate(value)
        else:
            value = self._convert_and_validate([value])[0]
        super().__setitem__(key, value)

    def copy(self):
        # Make sure a shallow copy still returns a ValidatedList.
        return self.__class__(self)

    def __getitem__(self, item):
        # Make sure a slice returns a ValidatedList.
        if isinstance(item, slice):
            return self.__class__(super().__getitem__(item))
        return super().__getitem__(item)

    # This is needed for pickling. See https://github.com/openforcefield/openforcefield/issues/411
    # for more details.
    # TODO: Is there a cleaner way (getstate/setstate perhaps?) to allow FFs to be
    #       pickled?
    def __reduce__(self):
        return (__class__, ( list( self),), self.__dict__)


    def _convert_and_validate(self, seq):
        """Run all converters and the validator on the given sequence."""
        # Run all element converters.
        if self._converters is not None:
            for converter in self._converters:
                seq = [converter(element) for element in seq]

        # Run all element validators.
        if self._validators is not None:
            for validator in self._validators:
                for element in seq:
                    validator(element)
        return seq

if __name__ == '__main__':
    import doctest
    doctest.run_docstring_examples(ValidatedList, globals())
