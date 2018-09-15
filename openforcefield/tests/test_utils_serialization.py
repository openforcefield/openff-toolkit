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

class TestUtilsSerialization(TestCase):
    """Test serialization and deserialization of a simple class."""

    def setUp(self):
        # Create an example object
        self.thing = Thing('blorb', [1,2,3])

    def test_json(self):
        """Test JSON serialization"""
        json_thing = self.thing.to_json()
        print(json_thing)
        thing_from_json = Thing.from_json(json_thing)
        assert self.thing == thing_from_json

    def test_yaml(self):
        """Test YAML serialization"""
        yaml_thing = self.thing.to_yaml()
        print(yaml_thing)
        thing_from_yaml = Thing.from_yaml(yaml_thing)
        assert self.thing == thing_from_yaml

    def test_bson(self):
        """Test BSON serialization"""
        bson_thing = self.thing.to_bson()
        print(bson_thing)
        thing_from_bson = Thing.from_bson(bson_thing)
        assert self.thing == thing_from_bson

    def test_toml(self):
        """Test TOML serialization"""
        toml_thing = self.thing.to_toml()
        print(toml_thing)
        thing_from_toml = Thing.from_toml(toml_thing)
        assert self.thing == thing_from_toml

    def test_messagepack(self):
        """Test MessagePack serialization"""
        messagepack_thing = self.thing.to_messagepack()
        print(messagepack_thing)
        thing_from_messagepack = Thing.from_messagepack(messagepack_thing)
        assert self.thing == thing_from_messagepack

    def test_xml(self):
        """Test XML serialization"""
        xml_thing = self.thing.to_xml()
        print(xml_thing)
        thing_from_xml = Thing.from_xml(xml_thing)
        assert self.thing == thing_from_xml
