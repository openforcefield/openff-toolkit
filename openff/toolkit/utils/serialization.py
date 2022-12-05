"""
Serialization mix-in

.. todo ::

   Currently, the ``openff-toolkit`` toolkit package requires a number
   of dependencies to support all of these serialization protocols.
   Instead, should we not include these by default, and instead raise
   a helpful exception with installation instructions if one of the
   serialization schemes is called but the requisite library is not
   installed?

"""
import abc
from typing import Dict

from openff.toolkit.utils.utils import requires_package


class Serializable(abc.ABC):
    """Mix-in to add serialization and deserialization support via JSON, YAML, BSON, TOML, MessagePack, and XML.

    For more information on these formats, see: `JSON <https://www.json.org/>`_,
    `BSON <http://bsonspec.org/>`_, `YAML <http://yaml.org/>`_, `TOML <https://github.com/toml-lang/toml>`_,
    `MessagePack <https://msgpack.org/index.html>`_, and `XML <https://www.w3.org/XML/>`_.

    To use this mix-in, the class inheriting from this class must have implemented ``to_dict()`` and ``from_dict()``
    methods that utilize dictionaries containing only serialiable Python objects.

    .. warning ::

       The serialization/deserialiation schemes used here place some strict constraints on what kinds of ``dict``
       objects can be serialized. No effort is made to add further protection to ensure serialization is possible.
       Use with caution.

    Examples
    --------

    Example class using :class:`Serializable` mix-in:

    >>> from openff.toolkit.utils.serialization import Serializable
    >>> class Thing(Serializable):
    ...     def __init__(self, description):
    ...         self.description = description
    ...
    ...     def to_dict(self):
    ...         return { 'description' : self.description }
    ...
    ...     @classmethod
    ...     def from_dict(cls, d):
    ...         return cls(d['description'])
    ...
    >>> # Create an example object
    >>> thing = Thing('blorb')

    Get `JSON <https://www.json.org/>`_ representation:

    >>> json_thing = thing.to_json()

    Reconstruct an object from its `JSON <https://www.json.org/>`_ representation:

    >>> thing_from_json = Thing.from_json(json_thing)

    Get `BSON <http://bsonspec.org/>`_ representation:

    >>> bson_thing = thing.to_bson()

    Reconstruct an object from its `BSON <http://bsonspec.org/>`_ representation:

    >>> thing_from_bson = Thing.from_bson(bson_thing)

    Get `YAML <http://yaml.org/>`_ representation:

    >>> yaml_thing = thing.to_yaml()

    Reconstruct an object from its `YAML <http://yaml.org/>`_ representation:

    >>> thing_from_yaml = Thing.from_yaml(yaml_thing)

    Get `MessagePack <https://msgpack.org/index.html>`_ representation:

    >>> messagepack_thing = thing.to_messagepack()

    Reconstruct an object from its `MessagePack <https://msgpack.org/index.html>`_ representation:

    >>> thing_from_messagepack = Thing.from_messagepack(messagepack_thing)

    Get `XML <https://www.w3.org/XML/>`_ representation:

    >>> xml_thing = thing.to_xml()

    """

    @abc.abstractmethod
    def to_dict(self):
        pass

    @classmethod
    @abc.abstractmethod
    def from_dict(cls, d):
        pass

    def to_json(self, indent=None) -> str:
        """
        Return a JSON serialized representation.

        Specification: https://www.json.org/

        Parameters
        ----------
        indent : int, optional, default=None
            If not None, will pretty-print with specified number of spaces for indentation

        Returns
        -------
        serialized : str
            A JSON serialized representation of the object

        """
        import json

        d = self.to_dict()

        if _contains_bytes(d):
            d = _prep_numpy_data_for_json(d)

        return json.dumps(d, indent=indent)

    @classmethod
    def from_json(cls, serialized: str):
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

    @requires_package("bson")
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
    @requires_package("bson")
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

    @requires_package("toml")
    def to_toml(self):
        """
        Return a TOML serialized representation.

        Specification: https://github.com/toml-lang/toml

        Returns
        -------
        serialized : str
            A TOML serialized representation of the object

        """
        raise NotImplementedError()
        # TODO: This implementation currently discards dict keys associated to the None value.
        #   See test_utils_serialization::TestUtilsSMIRNOFFSerialization::test_toml.
        # import toml
        # d = self.to_dict()
        # return toml.dumps(d)

    @classmethod
    @requires_package("toml")
    def from_toml(cls, serialized):
        """
        Instantiate an object from a TOML serialized representation.

        Specification: https://github.com/toml-lang/toml

        Parameters
        ----------
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

    @requires_package("yaml")
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
        return yaml.safe_dump(d, width=180)

    @classmethod
    @requires_package("yaml")
    def from_yaml(cls, serialized):
        """
        Instantiate from a YAML serialized representation.

        Specification: http://yaml.org/

        Parameters
        ----------
        serialized : str
            A YAML serialized representation of the object

        Returns
        -------
        instance : cls
            Instantiated object

        """
        import yaml

        d = yaml.safe_load(serialized)
        return cls.from_dict(d)

    @requires_package("msgpack")
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
    @requires_package("msgpack")
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

    def to_xml(self, indent=2):
        """
        Return an XML representation.

        Specification: https://www.w3.org/XML/

        Parameters
        ----------
        indent : int, optional, default=2
            If not None, will pretty-print with specified number of spaces for indentation

        Returns
        -------
        serialized : bytes
            A MessagePack-encoded bytes serialized representation.

        """
        import xmltodict

        # An XML document requires one and only one root node.
        root_name = self.__class__.__name__
        d = {root_name: self.to_dict()}
        # Configure indentation level.
        if indent is not None:
            pretty = True
            indent = " " * indent
        else:
            pretty = False
        # Convert data from dictionary to XML format.
        return xmltodict.unparse(d, pretty=pretty, indent=indent)

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
        raise NotImplementedError()
        # TODO: This implementation currently loads numbers as strings.
        #   See test_utils_serialization::TestUtilsSerialization::test_xml.
        # import xmltodict
        # d = xmltodict.parse(serialized)
        # root_name = cls.__name__
        # return cls.from_dict(d[root_name])

    def to_pickle(self):
        """
        Return a pickle serialized representation.

        .. warning ::

           This is not recommended for safe, stable storage since the pickle specification
           may change between Python versions.

        Returns
        -------
        serialized : str
            A pickled representation of the object

        """
        import pickle

        d = self.to_dict()
        return pickle.dumps(d)

    @classmethod
    def from_pickle(cls, serialized):
        """
        Instantiate an object from a pickle serialized representation.

        .. warning ::

           This is not recommended for safe, stable storage since the pickle specification
           may change between Python versions.

        Parameters
        ----------
        serialized : str
            A pickled representation of the object

        Returns
        -------
        instance : cls
            An instantiated object

        """
        import pickle

        d = pickle.loads(serialized)
        return cls.from_dict(d)


def _contains_bytes(val) -> bool:
    """Report if any values in list are bytes."""
    if val is None:
        return False
    elif isinstance(val, bytes):
        return True
    elif isinstance(val, (int, float, str, bool)):
        return False
    elif isinstance(val, (list, tuple)):
        return any([_contains_bytes(x) for x in val])
    elif isinstance(val, dict):
        return any([_contains_bytes(x) for x in val.values()])
    else:
        raise Exception(f"type {val}")


def _prep_numpy_data_for_json(data: Dict) -> Dict:
    """Recursively search through a dict and convert the bytes fields to lists"""
    import numpy as np

    big_endian_float = np.dtype("float").newbyteorder(">")

    for key, val in data.items():
        if isinstance(val, np.ndarray):
            data[key] = val.tolist()
        if isinstance(val, dict):
            data[key] = _prep_numpy_data_for_json(val)
        if isinstance(val, bytes):
            data[key] = np.frombuffer(val, dtype=big_endian_float).tolist()
        if isinstance(val, (list, tuple)):
            for i, element in enumerate(val):
                if isinstance(element, bytes):
                    # Handles case of List[np.array], like Molecule.conformers
                    data[key][i] = np.frombuffer(
                        element, dtype=big_endian_float
                    ).tolist()
                elif isinstance(element, dict):
                    # Handles case of List[Molecule], like Topology.molecules
                    data[key][i] = _prep_numpy_data_for_json(element)
    return data
