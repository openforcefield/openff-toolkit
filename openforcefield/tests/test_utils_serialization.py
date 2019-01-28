"""
Units tests for openforcefield.utils.serialization
"""

import logging
from enum import Enum, IntEnum
from typing import Dict

from pydantic import BaseModel

from openforcefield.utils.serialization import PolymorphicDataType


class Foo:

    def __init__(self):

        self.field1 = 'field1'
        self.field2 = 2

    def __getstate__(self):

        return {
            'field1': self.field1,
            'field2': self.field2
        }

    def __setstate__(self, state):

        self.field1 = state['field1']
        self.field2 = state['field2']


class Bar(BaseModel):

    field1: str = 'field1'
    field2: int = 2


class Baz(Enum):

    Option1 = "Option1"
    Option2 = "Option2"


class Qux(IntEnum):

    Option1 = 1
    Option2 = 2


class NestedParent:

    class NestedChild(Enum):

        Option1 = "Option1"
        Option2 = "Option2"


class PydanticTestClass(BaseModel):

    inputs: Dict[str, PolymorphicDataType] = None

    class Config:

        arbitrary_types_allowed = True

        json_encoders = {
            PolymorphicDataType: lambda value: PolymorphicDataType.serialize(value),
        }


def test_polymorphic_dictionary():
    """Test the polymorphic dictionary helper class."""

    test_dictionary = {
        "test_str": PolymorphicDataType(value='test1'),
        "test_int": PolymorphicDataType(value=1),
        "test_bool": PolymorphicDataType(value=True),
        "test_Foo": PolymorphicDataType(value=Foo()),
        "test_Bar": PolymorphicDataType(value=Bar()),
        "test_Baz": PolymorphicDataType(value=Baz.Option1),
        "test_Qux": PolymorphicDataType(value=Qux.Option1),
        "test_Nested": PolymorphicDataType(value=NestedParent.NestedChild.Option1)
    }

    pydantic_object = PydanticTestClass(inputs=test_dictionary)
    pydantic_json = pydantic_object.json()

    pydantic_recreated = PydanticTestClass.parse_raw(pydantic_json)
    pydantic_recreated_json = pydantic_recreated.json()

    assert pydantic_json == pydantic_recreated_json
