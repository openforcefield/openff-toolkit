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

import importlib
import json
import sys
from enum import Enum

from json import JSONEncoder, JSONDecoder, JSONDecodeError

from pydantic import BaseModel, ValidationError
from pydantic.validators import dict_validator
from simtk import unit


# =============================================================================================
# Module Constants
# =============================================================================================

class TypedBaseModel(BaseModel):

    module_metadata: str = ''
    type_metadata: str = ''

    @classmethod
    def __get_validators__(cls):
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


class PolymorphicEncoder(JSONEncoder):

    def default(self, obj):

        serializable_object = {}

        value_to_serialize = obj
        type_to_serialize = type(obj)

        from simtk import unit

        if isinstance(obj, PolymorphicDataType):

            value_to_serialize = obj.value
            type_to_serialize = obj.type

        if isinstance(value_to_serialize, BaseModel):

            serializable_object = value_to_serialize.dict()

        elif isinstance(value_to_serialize, Enum):

            serializable_object = value_to_serialize.value

        elif isinstance(value_to_serialize, unit.Quantity):

            serializable_object = serialize_quantity(value_to_serialize)

        else:

            try:

                json.dumps(value_to_serialize)  # Check if the value is natively serializable.
                serializable_object = value_to_serialize

            except TypeError as e:

                if isinstance(value_to_serialize, list):

                    serializable_object = []
                    list_object_type = None

                    for value_in_list in value_to_serialize:

                        if isinstance(value_in_list, PolymorphicDataType):
                            raise ValueError("Nested lists of PolymorphicDataType's are not serializable.")

                        if list_object_type is not None and type(value_in_list) != list_object_type:
                            raise ValueError("Nested lists of PolymorphicDataType's are not serializable.")

                        serializable_object.append(self.default(value_in_list))
                        list_object_type = type(value_in_list)

                else:

                    serializable_object = value_to_serialize.__getstate__()

        qualified_name = type_to_serialize.__qualname__
        qualified_name = qualified_name.replace('.', '->')

        return_value = {
            '@type': '{}.{}'.format(type_to_serialize.__module__,
                                    qualified_name),

            'value': serializable_object
        }

        return return_value


class PolymorphicDataType:
    """A helper object wrap values which have a type unknown
    ahead of time.
    """

    def __init__(self, value):
        """Creates a new PolymorphicDataType object.

        Parameters
        ----------
        value: Any
            The value to wrap.
        """

        self.value = value
        self.type = type(value)

    @classmethod
    def __get_validators__(cls):
        yield cls.validate

    @classmethod
    def validate(cls, value):
        """A pydantic helper method for deserializing the object.
        """
        if isinstance(value, PolymorphicDataType):
            return value

        return PolymorphicDataType.deserialize(value)

    @staticmethod
    def deserialize(json_dictionary):
        """A method to deserialize the polymorphic value from its
        JSON dictionary representation.

        Parameters
        ----------
        json_dictionary: Dict[str, Any]
            The JSON dictionary to deserialize.

        Returns
        -------
        Any
            The deserialized object.
        """

        if '@type' not in json_dictionary or 'value' not in json_dictionary:
            raise ValidationError('{} is not a valid PolymorphicDataType'.format(json_dictionary))

        type_string = json_dictionary['@type']
        last_period_index = type_string.rfind('.')

        if last_period_index < 0:
            raise ValidationError('{} is not a valid PolymorphicDataType'.format(json_dictionary))

        module_name = type_string[0:last_period_index]
        module = importlib.import_module(module_name)

        class_name = type_string[last_period_index + 1:]

        class_name_split = class_name.split('->')
        class_object = module

        while len(class_name_split) > 0:

            class_name_current = class_name_split.pop(0)
            class_object = getattr(class_object, class_name_current)

        value_object = json_dictionary['value']

        parsed_object = None

        from simtk import unit

        if issubclass(class_object, BaseModel):

            parsed_object = class_object.parse_obj(value_object)

        elif issubclass(class_object, Enum):

            parsed_object = class_object(value_object)

        elif issubclass(class_object, unit.Quantity):

            parsed_object = deserialize_quantity(value_object)

        elif issubclass(class_object, list):

            parsed_object = []

            for list_item in value_object:
                parsed_object.append(PolymorphicDataType.deserialize(list_item).value)

        else:

            parsed_object = value_object

            if hasattr(class_object, 'validate'):
                parsed_object = class_object.validate(parsed_object)

            elif not isinstance(parsed_object, class_object):

                created_object = class_object()
                created_object.__setstate__(parsed_object)

                parsed_object = created_object

        return PolymorphicDataType(parsed_object)

    @staticmethod
    def serialize(value_to_serialize):
        """A method to serialize a polymorphic value, along with its
        type in the form of a JSON dictionary.

        Parameters
        ----------
        value_to_serialize: PolymorphicDataType
            The value to serialize.

        Returns
        -------
        str
            The JSON serialized value.
        """

        # value_json = ''
        #
        # if isinstance(value_to_serialize.value, BaseModel):
        #
        #     value_json = value_to_serialize.value.json()
        #
        # elif isinstance(value_to_serialize.value, Enum):
        #
        #     value_json = value_to_serialize.value.value
        #
        # else:
        #     try:
        #         value_json = json.dumps(value_to_serialize.value)
        #     except TypeError as e:
        #         value_json = json.dumps(value_to_serialize.value.__getstate__())
        #
        # qualified_name = value_to_serialize.type.__qualname__
        # qualified_name = qualified_name.replace('.', '->')
        #
        # return_value = {
        #     '@type': '{}.{}'.format(value_to_serialize.type.__module__,
        #                             qualified_name),
        #
        #     'value': value_json
        # }

        encoder = PolymorphicEncoder()

        return_value = encoder.default(value_to_serialize)
        return return_value


def serialize_quantity(quantity):
    """
    Serialized a simtk.unit.Quantity into a dict of {'unitless_value': X, 'unit': Y}

    .. todo:: Currently duplicates Jeff Wagners implementation.

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

    conversion_factor = quantity.unit.get_conversion_factor_to_base_units()

    unitless_value = (quantity / quantity.unit) * conversion_factor
    serialized['unitless_value'] = unitless_value
    serialized['unit'] = quantity_unit
    return serialized


def deserialize_quantity(serialized):
    """
    Deserializes a simtk.unit.Quantity.

    .. todo:: Currently duplicates Jeff Wagners implementation.

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
