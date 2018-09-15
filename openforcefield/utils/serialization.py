#!/usr/bin/env python

"""
Serialization mix-in

.. todo ::

   Currently, the ``openforcefield`` toolkit package requires a number of dependencies to support all of these serialization protocols.
   Instead, should we not include these by default, and instead raise a helpful exception with installation instructions if one of the serialization schemes is called but the requisite library is not installed?

"""
#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

#=============================================================================================
# SERIALIZATION MIX-IN
#=============================================================================================

class Serializable(object):
    """Mix-in to add serialization and deserialization support via JSON, YAML, BSON, TOML, MessagePack, and XML.

    To use this mix-in, the class inheriting from this class must have implemented ``to_dict()`` and ``from_dict()`` methods
    that utilize dictionaries containing only serialiable Python objects.

    Examples
    --------

    Example class using :class:`Serializable` mix-in:

    >>> class Thing(Serializable):
    ...     def __init__(self, description):
    ...         self.description = description
    ...
    ...     def to_dict(self):
    ...         return { 'description' : self.description }
    ...
    ...     @staticmethod
    ...     def from_dict(d):
    ...         return Thing(d['description'])
    ...
    ... # Create an example object
    ... thing = Thing('blorb')

    Get JSON representation

    >>> json_thing = thing.to_json()
    >>> thing_from_json = Thing.from_json(json_thing)

    Get YAML representation

    >>> yaml_thing = thing.to_yaml()
    >>> thing_from_yaml = Thing.from_yaml(yaml_thing)

    Get BSON representation

    >>> bson_thing = thing.to_bson()
    >>> thing_from_bson = Thing.from_bson(bson_thing)

    Get TOML representation

    >>> toml_thing = thing.to_toml()
    >>> thing_from_toml = Thing.from_toml(toml_thing)

    Get MessagePack representation

    >>> messagepack_thing = thing.to_messagepack()
    >>> thing_from_messagepack = Thing.from_messagepack(messagepack_thing)

    Get XML representation

    >>> xml_thing = thing.to_xml()
    >>> thing_from_xml = Thing.from_xml(xml_thing)

    """
    def to_json(self):
        """
        Return a JSON serialized representation.

        Specification: https://www.json.org/

        Returns
        -------
        serialized : str
            A JSON serialized representation of the object

        """
        import json
        d = self.to_dict()
        return json.dumps(d)

    @classmethod
    def from_json(cls, serialized):
        """
        Instantiate an object from a JSON serialized representation.

        Specification: https://www.json.org/

        Parameters
        ----------
        serialized : str
            A JSON serialized representation of the object

        Returns
        -------
        instance : cls
            An instantiated object

        """
        import json
        d = json.loads(serialized)
        return cls.from_dict(d)

    def to_bson(self):
        """
        Return a BSON serialized representation.

        Specification: http://bsonspec.org/

        Returns
        -------
        serialized : bytes
            A BSON serialized representation of the objecft

        """
        import bson
        d = self.to_dict()
        return bson.dumps(d)

    @classmethod
    def from_bson(cls, serialized):
        """
        Instantiate an object from a BSON serialized representation.

        Specification: http://bsonspec.org/

        Parameters
        ----------
        serialized : bytes
            A BSON serialized representation of the object

        Returns
        -------
        instance : cls
            An instantiated object

        """
        import bson
        d = bson.loads(serialized)
        return cls.from_dict(d)

    def to_toml(self):
        """
        Return a TOML serialized representation.

        Specification: https://github.com/toml-lang/toml

        Returns
        -------
        serialized : str
            A TOML serialized representation of the object

        """
        import toml
        d = self.to_dict()
        return toml.dumps(d)

    @classmethod
    def from_toml(cls, serialized):
        """
        Instantiate an object from a TOML serialized representation.

        Specification: https://github.com/toml-lang/toml

        Parameters
        -------
        serlialized : str
            A TOML serialized representation of the object

        Returns
        -------
        instance : cls
            An instantiated object

        """
        import toml
        d = toml.loads(serialized)
        return cls.from_dict(d)

    def to_yaml(self):
        """
        Return a YAML serialized representation.

        Specification: http://yaml.org/

        Returns
        -------
        serialized : str
            A YAML serialized representation of the object

        """
        import yaml
        d = self.to_dict()
        return yaml.dump(d)

    @classmethod
    def from_yaml(cls, serialized):
        """
        Instantiate from a YAML serialized representation.

        Specification: http://yaml.org/

        Parameters
        -------
        serialized : str
            A YAML serialized representation of the object

        Returns
        -------
        instance : cls
            Instantiated object

        """
        import yaml
        d = yaml.load(serialized)
        return cls.from_dict(d)

    def to_messagepack(self):
        """
        Return a MessagePack representation.

        Specification: https://msgpack.org/index.html

        Returns
        -------
        serialized : bytes
            A MessagePack-encoded bytes serialized representation of the object

        """
        import msgpack
        d = self.to_dict()
        return msgpack.dumps(d, use_bin_type=True)

    @classmethod
    def from_messagepack(cls, serialized):
        """
        Instantiate an object from a MessagePack serialized representation.

        Specification: https://msgpack.org/index.html

        Parameters
        ----------
        serialized : bytes
            A MessagePack-encoded bytes serialized representation

        Returns
        -------
        instance : cls
            Instantiated object.

        """
        import msgpack
        d = msgpack.loads(serialized, raw=False)
        return cls.from_dict(d)

    def to_xml(self, pretty=True):
        """
        Return an XML representation.

        Specification: https://www.w3.org/XML/

        Parameters
        ----------
        pretty : bool, optional, default=True
            If True, will pretty-format the XML by inserting additional spaces

        Returns
        -------
        serialized : bytes
            A MessagePack-encoded bytes serialized representation.

        """
        import xmltodict
        d = self.to_dict()
        root_name = self.__class__.__name__
        return xmltodict.unparse({root_name : d}, pretty=pretty)

    @classmethod
    def from_xml(cls, serialized):
        """
        Instantiate an object from an XML serialized representation.

        Specification: https://www.w3.org/XML/

        Parameters
        ----------
        serialized : bytes
            An XML serialized representation

        Returns
        -------
        instance : cls
            Instantiated object.

        """
        import xmltodict
        d = xmltodict.parse(serialized)
        root_name = cls.__name__
        return cls.from_dict(d[root_name])
