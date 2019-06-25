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

    def __init__(self, seq=(), converter=None, validator=None):
        self._converter = converter
        self._validator = validator

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

    def _convert_and_validate(self, seq):
        if self._converter is not None:
            seq = [self._converter(element) for element in seq]
        if self._validator is not None:
            for element in seq:
                self._validator(element)
        return seq
