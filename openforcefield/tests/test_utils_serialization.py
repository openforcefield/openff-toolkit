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
    def __init__(self, description):
        self.description = description

    def to_dict(self):
        return { 'description' : self.description }

    @staticmethod
    def from_dict(d):
        return Thing(d['description'])

class TestUtilsSerialization(TestCase):


    def setUp(self):
        # Create an example object
        self.thing = Thing('blorb')

    def test_json(self):
        """Test JSON serialization"""
        json_thing = self.thing.to_json()
        thing_from_json = Thing.from_json(json_thing)

    def test_yaml(self):
        """Test YAML serialization"""
        yaml_thing = self.thing.to_yaml()
        thing_from_yaml = Thing.from_yaml(yaml_thing)

    def test_bson(self):
        """Test BSON serialization"""
        bson_thing = self.thing.to_bson()
        thing_from_bson = Thing.from_bson(bson_thing)

    def test_toml(self):
        """Test TOML serialization"""
        toml_thing = self.thing.to_toml()
        thing_from_toml = Thing.from_toml(toml_thing)

    def test_messagepack(self):
        """Test MessagePack serialization"""
        messagepack_thing = self.thing.to_messagepack()
        thing_from_messagepack = Thing.from_messagepack(messagepack_thing)

    def test_xml(self):
        """Test XML serialization"""
        xml_thing = self.thing.to_xml()
        thing_from_xml = Thing.from_xml(xml_thing)
