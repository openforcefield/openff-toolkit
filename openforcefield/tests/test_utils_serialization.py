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

from unittest import TestCase

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

    @staticmethod
    def from_dict(d):
        return Thing(d['description'], d['mylist'])

    def __eq__(self, other):
        """Comparator for asserting object field equality."""
        equality = True
        equality == equality and (self.description == other.description)
        equality == equality and (self.mylist == other.mylist)
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

class TestUtilsSerialization(TestCase):
    """Test serialization and deserialization of a simple class."""

    def setUp(self):
        # Create an example object
        self.thing = Thing('blorb', [1,2,3])
        self.thing_class = Thing

    def test_json(self):
        """Test JSON serialization"""
        json_thing = self.thing.to_json()
        write('out.json', json_thing) # DEBUG
        thing_from_json = self.thing_class.from_json(json_thing)
        assert self.thing == thing_from_json

    def test_yaml(self):
        """Test YAML serialization"""
        yaml_thing = self.thing.to_yaml()
        write('out.yaml', yaml_thing) # DEBUG
        thing_from_yaml = self.thing_class.from_yaml(yaml_thing)
        assert self.thing == thing_from_yaml

    def test_bson(self):
        """Test BSON serialization"""
        bson_thing = self.thing.to_bson()
        write('out.bson', bson_thing) # DEBUG
        thing_from_bson = self.thing_class.from_bson(bson_thing)
        assert self.thing == thing_from_bson

    def test_toml(self):
        """Test TOML serialization"""
        toml_thing = self.thing.to_toml()
        write('out.toml', toml_thing) # DEBUG
        thing_from_toml = self.thing_class.from_toml(toml_thing)
        assert self.thing == thing_from_toml

    def test_messagepack(self):
        """Test MessagePack serialization"""
        messagepack_thing = self.thing.to_messagepack()
        write('out.messagepack', messagepack_thing) # DEBUG
        thing_from_messagepack = self.thing_class.from_messagepack(messagepack_thing)
        assert self.thing == thing_from_messagepack

    def test_xml(self):
        """Test XML serialization"""
        xml_thing = self.thing.to_xml()
        write('out.xml', xml_thing) # DEBUG
        thing_from_xml = self.thing_class.from_xml(xml_thing)
        assert self.thing == thing_from_xml

    def test_pickle(self):
        """Test pickle serialization"""
        pkl_thing = self.thing.to_pickle()
        write('out.pkl', pkl_thing) # DEBUG
        thing_from_pkl = self.thing_class.from_pickle(pkl_thing)
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

    def setUp(self):
        # Create an example object holding the SMIRNOFF xmltodict dictionary representation
        import xmltodict
        from openforcefield.utils import get_data_filename
        filename = get_data_filename('forcefield/smirnoff99Frosst.offxml')
        with open(filename, 'r') as f:
            xml = f.read()
            dictionary = xmltodict.parse(xml)
            self.thing = DictionaryContainer(dictionary)
            self.thing_class = DictionaryContainer
