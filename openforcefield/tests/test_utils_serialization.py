"""
Units tests for openforcefield.utils.serialization
"""

import logging
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
        "test_Foo": PolymorphicDataType(value=Foo()),
        "test_Bar": PolymorphicDataType(value=Bar())
    }

    pydantic_object = PydanticTestClass(inputs=test_dictionary)
    pydantic_json = pydantic_object.json()

    pydantic_recreated = PydanticTestClass.parse_raw(pydantic_json)

    logging.info(pydantic_json)
