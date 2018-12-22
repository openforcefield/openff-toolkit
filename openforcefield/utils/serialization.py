#!/usr/bin/env python

# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
An API to store general exceptions that may be raised.

Authors
-------
* Simon Boothroyd <simon.boothroyd@choderalab.org>
"""


# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

import sys

from pydantic.validators import dict_validator
from simtk import unit

from pydantic import BaseModel, validator


# =============================================================================================
# Module Constants
# =============================================================================================

class TypedBaseModel(BaseModel):

    module_metadata: str = ''
    type_metadata: str = ''

    @classmethod
    def get_validators(cls):
        # yield dict_validator
        yield cls.validate

    @classmethod
    def validate(cls, value):
        if isinstance(value, cls):
            return value
        else:
            class_type = getattr(sys.modules[value['module_metadata']], value['type_metadata'])
            return class_type(**dict_validator(value))

    def __init__(self, **data):
        super().__init__(**data)

        self.module_metadata = type(self).__module__
        self.type_metadata = type(self).__name__


def serialize_quantity(quantity):
    """
    Serialized a simtk.unit.Quantity into a dict of {'unitless_value': X, 'unit': Y}

    Todo
    ----
    Currently duplicates Jeff Wagners implementation.

    Parameters
    ----------
    quantity : A simtk.unit.Quantity-wrapped value or iterator over values
        The object to serialize

    Returns
    -------
    serialzied : dict
        The serialized object
    """

    serialized = dict()

    # If it's None, just return None in all fields
    if quantity is None:
        serialized['unitless_value'] = None
        serialized['unit'] = None
        return serialized

    # If it's not None, make sure it's a simtk.unit.Quantity
    assert (hasattr(quantity, 'unit'))


    quantity_unit = list()
    for base_unit in quantity.unit.iter_all_base_units():
        quantity_unit.append((base_unit[0].name, base_unit[1]))

    unitless_value = quantity / quantity.unit
    serialized['unitless_value'] = unitless_value
    serialized['unit'] = quantity_unit
    return serialized


def deserialize_quantity(serialized):
    """
    Deserializes a simtk.unit.Quantity.

    Todo
    ----
    Currently duplicates Jeff Wagners implementation.

    Parameters
    ----------
    serialized : dict
        Serialized representation of a simtk.unit.Quantity. Must have keys ["unitless_value", "unit"]

    Returns
    -------
    simtk.unit.Quantity
    """
    if (serialized['unitless_value'] is None) and (serialized['unit'] is None):
        return None
    quantity_unit = None
    for unit_name, power in serialized['unit']:
        unit_name = unit_name.replace(' ','_') # Convert eg. 'elementary charge' to 'elementary_charge'
        if quantity_unit is None:
            quantity_unit = (getattr(unit, unit_name) ** power)
        else:
            quantity_unit *= (getattr(unit, unit_name) ** power)
    quantity = unit.Quantity(serialized['unitless_value'], quantity_unit)
    return quantity
