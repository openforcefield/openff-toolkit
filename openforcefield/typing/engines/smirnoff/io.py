#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================
"""
XML I/O parser for the SMIRNOFF (SMIRKS Native Open Force Field) format.

.. codeauthor:: John D. Chodera <john.chodera@choderalab.org>
.. codeauthor:: David L. Mobley <dmobley@mobleylab.org>
.. codeauthor:: Peter K. Eastman <peastman@stanford.edu>

"""

__all__ = [
    'ParameterIOHandler',
    'XMLParameterIOHandler',
]


#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import logging

import xmltodict
from simtk import unit


#=============================================================================================
# CONFIGURE LOGGER
#=============================================================================================

logger = logging.getLogger(__name__)


#=============================================================================================
# QUANTITY PARSING UTILITIES
#=============================================================================================

def _ast_unit_eval(node):
    """
    Performs a safe algebraic syntax tree evaluation of a unit.

    This will likely be replaced by the native implementation in Pint
    if/when we'll switch over to Pint units.
    """
    import ast
    import operator as op

    operators = {ast.Add: op.add, ast.Sub: op.sub, ast.Mult: op.mul,
        ast.Div: op.truediv, ast.Pow: op.pow, ast.BitXor: op.xor,
        ast.USub: op.neg}

    if isinstance(node, ast.Num):  # <number>
        return node.n
    elif isinstance(node, ast.BinOp):  # <left> <operator> <right>
        return operators[type(node.op)](_ast_unit_eval(node.left), _ast_unit_eval(node.right))
    elif isinstance(node, ast.UnaryOp):  # <operator>( <operand> ) e.g., -1
        return operators[type(node.op)](_ast_unit_eval(node.operand))
    elif isinstance(node, ast.Name):
        # Check if this is a simtk unit.
        u = getattr(unit, node.id)
        if not isinstance(u, unit.Unit):
            raise ValueError('No unit named {} found in simtk.unit.'.format(node.id))
        return u
    else:
        raise TypeError(node)


#=============================================================================================
# Base ParameterIOHandler
#=============================================================================================


class ParameterIOHandler:
    """
    Base class for handling serialization/deserialization of SMIRNOFF ForceField objects
    """
    _FORMAT = None

    def __init__(self):
        """
        Create a new ParameterIOHandler.

        """
        pass

    def parse_file(self, file_path):
        """

        Parameters
        ----------
        file_path

        Returns
        -------

        """
        pass

    def parse_string(self, data):
        """
        Parse a SMIRNOFF force field definition in a seriaized format

        Parameters
        ----------
        data

        Returns
        -------

        """
        pass

    def to_file(self, file_path, smirnoff_data):
        """
        Write the current forcefield parameter set to a file.

        Parameters
        ----------
        file_path : str
            The path to the file to write to.
        smirnoff_data : dict
            A dictionary structured in compliance with the SMIRNOFF spec

        Returns
        -------

        """
        pass

    def to_string(self, smirnoff_data):
        """
        Render the forcefield parameter set to a string

        Parameters
        ----------
        smirnoff_data : dict
            A dictionary structured in compliance with the SMIRNOFF spec

        Returns
        -------
        str
        """
        pass


#=============================================================================================
# XML I/O
#=============================================================================================

class XMLParameterIOHandler(ParameterIOHandler):
    """
    Handles serialization/deserialization of SMIRNOFF ForceField objects from OFFXML format.
    """
    # TODO: Come up with a better keyword for format
    _FORMAT = 'XML'

    def parse_file(self, source):
        """Parse a SMIRNOFF force field definition in XML format, read from a file.

        Parameters
        ----------
        source : str or io.RawIOBase
            File path of file-like object implementing a ``read()`` method
            specifying a SMIRNOFF force field definition in `the SMIRNOFF XML format <https://github.com/openforcefield/openforcefield/blob/master/The-SMIRNOFF-force-field-format.md>`_.

        Raises
        ------
        ParseError
            If the XML cannot be processed.
        FileNotFoundError
            If the file could not found.

        """
        # If this is a file-like object, we should be able to read it.
        try:
            raw_data = source.read()
        except AttributeError:
            # This raises FileNotFoundError if the file doesn't exist.
            raw_data = open(source).read()

        # Parse the data in string format.
        return self.parse_string(raw_data)

    def parse_string(self, data):
        """Parse a SMIRNOFF force field definition in XML format.

        A ``ParseError`` is raised if the XML cannot be processed.

        Parameters
        ----------
        data : str
            A SMIRNOFF force field definition in `the SMIRNOFF XML format <https://github.com/openforcefield/openforcefield/blob/master/The-SMIRNOFF-force-field-format.md>`_.

        """
        from pyexpat import ExpatError
        from openforcefield.typing.engines.smirnoff.forcefield import ParseError

        # Parse XML file
        try:
            smirnoff_data = xmltodict.parse(data, attr_prefix='')
            return smirnoff_data
        except ExpatError as e:
            raise ParseError(e)

    def to_file(self, file_path, smirnoff_data):
        """Write the current forcefield parameter set to a file.

        Parameters
        ----------
        file_path : str
            The path to the file to be written.
            The `.offxml` or `.xml` file extension must be present.
        smirnoff_data : dict
            A dict structured in compliance with the SMIRNOFF data spec.

        """
        xml_string = self.to_string(smirnoff_data)
        with open(file_path, 'w') as of:
            of.write(xml_string)

    def to_string(self, smirnoff_data):
        """
        Write the current forcefield parameter set to an XML string.

        Parameters
        ----------
        smirnoff_data : dict
            A dictionary structured in compliance with the SMIRNOFF spec

        Returns
        -------
        serialized_forcefield : str
            XML String representation of this forcefield.

        """
        def prepend_all_keys(d, char='@'):
            """
            Modify a dictionary in-place, prepending a specified string to each key
            that doesn't refer to a value that is list or dict.

            Parameters
            ----------
            d : dict
                Hierarchical dictionary to traverse and modify keys
            char : string, optional. Default='@'
                String to prepend onto each applicable dictionary key

            """
            if isinstance(d, dict):
                for key in list(d.keys()):
                    if isinstance(d[key], list) or isinstance(d[key], dict):
                        prepend_all_keys(d[key])
                    else:
                        new_key = char + key
                        d[new_key] = d[key]
                        del d[key]
                        prepend_all_keys(d[new_key])
            elif isinstance(d, list):
                for item in d:
                    prepend_all_keys(item)

        prepend_all_keys(smirnoff_data['SMIRNOFF'])
        return xmltodict.unparse(smirnoff_data, pretty=True)


