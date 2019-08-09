#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Tests for utility methods for serialization

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import pytest


#=============================================================================================
# TESTS
#=============================================================================================

from openforcefield.utils.serialization import Serializable

class Thing(Serializable):

    def __init__(self, description, mylist):
        self.description = description
        self.mylist = mylist

    def to_dict(self):
        return {
            'description' : self.description,
            'mylist' : self.mylist
        }

    @classmethod
    def from_dict(cls, d):
        return cls(d['description'], d['mylist'])

    def __eq__(self, other):
        """Comparator for asserting object field equality."""
        equality = self.description == other.description
        equality &= self.mylist == other.mylist
        return equality


# DEBUG
def write(filename, contents):
    if type(contents) == str:
        mode = 'w'
    elif type(contents) == bytes:
        mode = 'wb'
    else:
        raise Exception('Cannot handle contents of type {}'.format(type(contents)))
    with open(filename, mode) as outfile:
        outfile.write(contents)


class TestUtilsSerialization:
    """Test serialization and deserialization of a simple class."""

    @classmethod
    def setup_class(cls):
        cls.thing = Thing('blorb', [1, 2, 3])

    def test_json(self):
        """Test JSON serialization"""
        json_thing = self.thing.to_json()
        thing_from_json = self.thing.__class__.from_json(json_thing)
        assert self.thing == thing_from_json

    def test_yaml(self):
        """Test YAML serialization"""
        yaml_thing = self.thing.to_yaml()
        thing_from_yaml = self.thing.__class__.from_yaml(yaml_thing)
        assert self.thing == thing_from_yaml

    def test_bson(self):
        """Test BSON serialization"""
        bson_thing = self.thing.to_bson()
        thing_from_bson = self.thing.__class__.from_bson(bson_thing)
        assert self.thing == thing_from_bson

    @pytest.mark.wip(reason=("the current implementation of to_toml cannot handle dict "
                             "keys associated to None (e.g. the ToolkitAM1BCC tag in the "
                             "TestUtilsSMIRNOFFSerialization suite)."))
    def test_toml(self):
        """Test TOML serialization"""
        toml_thing = self.thing.to_toml()
        thing_from_toml = self.thing.__class__.from_toml(toml_thing)
        assert self.thing == thing_from_toml

    def test_messagepack(self):
        """Test MessagePack serialization"""
        messagepack_thing = self.thing.to_messagepack()
        thing_from_messagepack = self.thing.__class__.from_messagepack(messagepack_thing)
        assert self.thing == thing_from_messagepack

    @pytest.mark.wip(reason='from/to_xml is not implemented yet. This test fails '
                            'because the list of integers is saved in the XML as '
                            'a list of numeric strings.')
    def test_xml(self):
        """Test XML serialization"""
        xml_thing = self.thing.to_xml()
        thing_from_xml = self.thing.__class__.from_xml(xml_thing)
        assert self.thing == thing_from_xml

    def test_pickle(self):
        """Test pickle serialization"""
        pkl_thing = self.thing.to_pickle()
        thing_from_pkl = self.thing.__class__.from_pickle(pkl_thing)
        assert self.thing == thing_from_pkl


class DictionaryContainer(Serializable):
    def __init__(self, dictionary):
        import copy
        self.dictionary = copy.deepcopy(dictionary)

    def to_dict(self):
        return self.dictionary

    @staticmethod
    def from_dict(dictionary):
        return DictionaryContainer(dictionary)

    def __eq__(self, other):
        """Comparator for asserting object field equality."""
        return self.dictionary == other.dictionary


class TestUtilsSMIRNOFFSerialization(TestUtilsSerialization):
    """Test serialization and deserialization of smirnoff99Frosst."""

    @classmethod
    def setup_class(cls):
        cls.thing = Thing('blorb', [1, 2, 3])
        # Create an example object holding the SMIRNOFF xmltodict dictionary representation
        import xmltodict
        from openforcefield.utils import get_data_file_path
        filename = get_data_file_path('test_forcefields/smirnoff99Frosst.offxml')
        with open(filename, 'r') as f:
            xml = f.read()
            dictionary = xmltodict.parse(xml)
            cls.thing = DictionaryContainer(dictionary)
