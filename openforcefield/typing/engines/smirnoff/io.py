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

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

from collections import OrderedDict
import logging
import os
import sys

import xmltodict

from simtk import openmm, unit
from openforcefield.utils.utils import get_data_filename


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


class ParameterIOHandler(object):
    """
    Handles serialization/deserialization of SMIRNOFF ForceField objects
    """
    _FORMAT = None

    def __init__(self):
        """
        Create a new ParameterIOHandler.

        """
        pass
        #self._forcefield = forcefield

    def parse_file(self, filename):
        """

        Parameters
        ----------
        filename

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

    def to_file(self, filename):
        """
        Write the current forcefield parameter set to a file.

        Parameters
        ----------
        filename

        Returns
        -------

        """
        pass

    def to_string(self):
        """
        Render the forcefield parameter set to a string

        Returns
        -------
        str
        """
        pass

    #def parameter_handler_tags_are_compatible(self, tags):
    #    """
    #
    #    Parameters
    #    ----------
    #    tags : dict of {'length_unit': 'angstroms', 'k_unit': 'kilocalories_per_mole/angstrom**2'}
    #
    #    Returns
    #    -------
    #
    #    """


#=============================================================================================
# XML I/O
#=============================================================================================



class XMLParameterIOHandler(ParameterIOHandler):
    # TODO: Come up with a better keyword for format
    _FORMAT = 'XML'



    # TODO: Fix this
    def parse_file(self, source):
        """Parse a SMIRNOFF force field definition in XML format, read from a file.

        A ``ParseError`` is raised if the XML cannot be processed.

        Parameters
        ----------
        source : str or file-like obj
            File path of file-like obj specifying a SMIRNOFF force field definition in `the SMIRNOFF XML format <https://github.com/openforcefield/openforcefield/blob/master/The-SMIRNOFF-force-field-format.md>`_.

        .. notes ::

           * New SMIRNOFF sections are handled independently, as if they were specified in the same file.
           * If a SMIRNOFF section that has already been read appears again, its definitions are appended to the end of the previously-read
             definitions if the sections are configured with compatible attributes; otherwise, an ``IncompatibleTagException`` is raised.
        """
        from pyexpat import ExpatError
        from openforcefield.typing.engines.smirnoff.forcefield import ParseError

        # this handles open file-like objects (and strings)
        try:
            smirnoff_dict = xmltodict.parse(source, attr_prefix='')
            print(smirnoff_dict)
            return smirnoff_dict
        except ExpatError:
            pass

        # This handles complete/local filenames
        try:
            # Check if the file exists in the data/forcefield directory
            data = open(source).read()
            smirnoff_data = xmltodict.parse(data, attr_prefix='')
            return smirnoff_data
        except ExpatError:
            pass
        except FileNotFoundError:
            pass

        # This handles nonlocal filenames
        try:
            # Check if the file exists in the data/forcefield directory
            temp_file = get_data_filename(os.path.join('forcefield', source))
            data = open(temp_file).read()
            smirnoff_data = xmltodict.parse(data, attr_prefix='')
            return smirnoff_data
        except Exception as e:
            # Fail with an error message about which file could not be read.
            # TODO: Also handle case where fallback to 'data' directory encounters problems,
            # but this is much less worrisome because we control those files.
            msg = str(e) + '\n'
            if hasattr(source, 'name'):
                filename = source.name
            else:
                filename = str(source)
            msg += "ForceField.loadFile() encountered an error reading file '%s'\n" % filename
            raise ParseError(msg)

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

    def to_file(self, filename, root):
        """Write the current forcefield parameter set to a file, autodetecting the type from the extension.

        Parameters
        ----------
        filename : str
            The name of the file to be written.
            The `.offxml` file extension is auto-detected.
        root : str, optional, default=None

        """
        (basename, extension) = os.path.splitext(filename)
        if extension == '.offxml':
            root.write(filename, xml_declaration=True, pretty_print=True)
        else:
            msg = "Cannot export forcefield parameters to file '{}'\n".format(
                filename)
            msg += 'Export to extension {} not implemented yet.\n'.format(
                extension)
            msg += "Supported choices are: ['.offxml']"
            raise NotImplementedError(msg)

    def to_string(self, smirnoff_data):
        """

        Parameters
        ----------
        smirnoff_data : dict
            A dictionary structured in compliance with the SMIRNOFF spec

        Returns
        -------

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
        print(smirnoff_data)
        print()
        return xmltodict.unparse(smirnoff_data, pretty=True)

    #@staticmethod
    #def from_xml(xml):
    #    """Construct a ForceField object from a SMIRNOFF XML file.
    #    """
    #    return ForceField(xml)


    # def from_lxml(self, root):
    #
    #
    #     cosmetic_tags = ['Date', 'Author']
    #
    #     try:
    #         exception_node = root  # node used for exception reporting
    #         if not (root.tag == 'SMIRNOFF'):
    #             raise ParseError("Root tag of tree is not 'SMIRNOFF'")
    #         self._parse_version(root)
    #         self._parse_aromaticity_model(root)
    #         # Process handlers
    #         #for section in root.iter(tag=etree.Element):
    #         for section in root:
    #             # Skip comment lines
    #             if not isinstance(section.tag, str):
    #                 continue
    #             if section.tag in cosmetic_tags:
    #                 continue
    #             exception_node = section  # node used for exception reporting
    #             # Extract parameter name from XML tag
    #             parameter_name = section.tag
    #
    #             # Split out attributes that refer to units
    #             handler_kwargs, attached_units = _extract_attached_units(
    #                 section.attrib)
    #             # TODO: Attach units to handler_kwargs
    #             # Make a copy of handler_kwargs that we can modify
    #             handler_kwargs_dict = dict(handler_kwargs)
    #             handler_kwargs_dict = _attach_units(handler_kwargs_dict,
    #                                                 attached_units)
    #
    #             # Retrieve or create parameter handler
    #             handler = self._forcefield.get_handler(parameter_name,
    #                                                    handler_kwargs_dict)
    #
    #             # Populate handler with parameter definitions
    #             for parameter in section:
    #                 # Skip comment lines
    #                 if not isinstance(parameter.tag, str):
    #                     continue
    #                 exception_node = parameter  # node used for exception reporting
    #
    #                 # parameter.attrib doesn't support assignment of Quantity type, so make a copy as dict
    #                 parameter_attrib_dict = dict(parameter.attrib)
    #
    #                 # Append units to parameters as needed
    #                 parameter_kwargs = _attach_units(parameter_attrib_dict,
    #                                                  attached_units)
    #                 # Add parameter definition
    #                 handler.add_parameter(parameter_kwargs)
    #
    #     except Exception as e:
    #         # Prepend the line and line text of XML file to aid debugging
    #         # TODO: Can we include the filename as well?
    #         print(self._get_sourceline(exception_node))
    #         #if not e.args:
    #         #    e.args = ('',)
    #         #e.args = e.args[0] + + e.args[1:]
    #         raise e

    def to_xml(self):
        """Render the forcefield parameter set to XML.

        Returns
        -------
        xml : str
            The SMIRNOFF parameter set rendered as XML.
        """
        return etree.tostring(self.to_lxml())
        # Test that this works

    # TODO: Do we need this? Should we use built-in dict-based serialization?
    def __getstate__(self):
        """Serialize to XML.
        """
        return self.to_xml()

    # TODO: Do we need this? Should we use built-in dict-based serialization?
    def __setstate__(self, state):
        """Deserialize from XML.
        """
        self._initialize()
        self.parse_xml(state)
