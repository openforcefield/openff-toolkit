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

import lxml.etree as etree

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


def _extract_attached_units(attrib):
    """Form a (potentially unit-bearing) quantity from the specified attribute name.

    The current implementation will likely be replaced by the native
    implementation in Pint if/when we'll switch over to Pint units.

    Parameters
    ----------
    attrib : dict
       Dictionary of XML node attributes.

    Returns
    -------
    attrib : dict
       XML node attributes with keys ending in ``_unit`` removed.
    attached_units : dict str : simtk.unit.Unit
       ``attached_units[parameter_name]`` is the simtk.unit.Unit combination
       that should be attached to corresponding parameter ``parameter_name``.

    """
    # TODO: Should this scheme also convert unit-bearing quantities such as '8*angstroms' to 8*unit.angstroms?
    # TODO: Should this scheme also convert "1" to int(1) and "8.0" to float(8.0)?

    import ast

    attached_units = OrderedDict()
    for key in attrib.keys():
        if key.endswith('_unit'):
            parameter_name = key[:-5]
            parameter_units_string = attrib[key]

            # Parse string expression and obtain units.
            try:
                ast_root_node = ast.parse(parameter_units_string, mode='eval').body
                parameter_units = _ast_unit_eval(ast_root_node)
            except Exception as e:
                # Re-raise parsing exception preserving the stack.
                err_msg = 'Could not parse units {}\n'.format(parameter_units_string) + str(e)
                raise type(e)(err_msg).with_traceback(sys.exc_info()[2])

            # Check that "attr_unit" was not associated to a Quantity that
            # could fail silently when multiplied with the value in "attr".
            if isinstance(parameter_units, unit.Quantity):
                err_msg = ('{key} was associated to a quantity rather than only units: '
                           '{expr}'.format(key=key, expr=parameter_units_string))
                raise ValueError(err_msg)
            attached_units[parameter_name] = parameter_units

    return attrib, attached_units


def _attach_units(attrib, attached_units):
    """Attach units to attributes for which units are specified.

    Parameters
    ----------
    attrib : dict
       Dictionary of XML node attributes.
    attached_units : dict str : simtk.unit.Unit
       ``attached_units[parameter_name]`` is the simtk.unit.Unit combination that should be attached to corresponding
       parameter ``parameter_name``

    Returns
    -------
    attrib : dict
       Updated XML node attributes with simtk.unit.Unit units attached to values for which units were specified for their keys

    """
    for parameter_name, units_to_attach in attached_units.items():
        if parameter_name in attrib.keys():
            parameter_attrib_string = attrib[parameter_name]
            try:
                attrib[parameter_name] = float(
                    parameter_attrib_string) * units_to_attach
            except ValueError as e:
                e.msg = (
                    "Expected numeric value for parameter '{}',"
                    "instead found '{}' when trying to attach units '{}'\n"
                ).format(parameter_name, parameter_attrib_string, units_to_attach)
                raise e

        # Now check for matches like "phase1", "phase2"
        c = 1
        while (parameter_name + str(c)) in attrib.keys():
            indexed_parameter_name = parameter_name + str(c)
            parameter_attrib_string = attrib[indexed_parameter_name]
            try:
                attrib[indexed_parameter_name] = float(
                    parameter_attrib_string) * units_to_attach
            except ValueError as e:
                e.msg = "Expected numeric value for parameter '{}', instead found '{}' when trying to attach units '{}'\n".format(
                    indexed_parameter_name, parameter_attrib_string,
                    units_to_attach)
                raise e
            c += 1
        #if parameter_name in attached_units:
        #units_to_attach = attached_units[parameter_name]
        # TODO: Do we have to worry about None or null values for parameters with attached units?
        #try:
        #    attrib[parameter_name_attrib] = float(parameter_value_string) * units_to_attach
        #except Exception as e:
        #    e.msg = "Expected numeric value for parameter '{}', instead found '{}' when trying to attach units '{}'\n".format(
        #        parameter_name, parameter_value_string, units_to_attach)
        #    raise e

    return attrib


class ParseError(Exception):
    """
    Exception for when a file is not parseable by a ParameterIOHandler
    """

    def __init__(self, msg):
        super(ParseError, self).__init__(self, msg)
        self.msg = msg


#=============================================================================================
# Base ParameterIOHandler
#=============================================================================================


class ParameterIOHandler(object):
    """
    Handles serialization/deserialization of SMIRNOFF ForceField objects
    """
    _FORMAT = None

    def __init__(self, forcefield):
        """
        Create a new ParameterIOHandler.

        Parameters
        ----------
        forcefield : openforcefield.typing.engines.smirnoff.ForceField
            The ForceField that this ParameterIOHandler belongs to. The ParameterIOHandler will read serialized
            ForceField representations, create ParameterType-derived objects, and add the ParameterType-derived objects
            into the ForceField's matching ParameterHandler.
        """
        self._forcefield = forcefield

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

# TODO: Instead of subclassing ForceField, what if we had an abstract base class for parameter IO?


class XMLParameterIOHandler(ParameterIOHandler):
    # TODO: Come up with a better keyword for format
    _FORMAT = 'XML'

    #def __init__(self, *args, **kwargs):
    #    super().__init__(*args, **kwargs)

    @staticmethod
    def _get_sourceline(node, filename=None):
        """Prepend the source line number to aid debugging if the XML parsing library supports itself.

        Parameters
        ----------
        node : lxml.etree.Node like, optional, default=None
            XML DOM node object, which may optionally have a ``sourceline`` field
        filename :

        """
        if filename:
            msg = "Error encountered parsing '{}'".format(filename)
        else:
            msg = "Error encountered parsing XML"
        if hasattr(node, 'sourceline'):
            msg += "line {}:\n{}\n".format(node.sourceline, str(node))
        else:
            msg += "line:\n{}\n".format(str(node))

        return msg

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

        # TODO: Search the right sequence of paths
        parser = etree.XMLParser(
            remove_blank_text=True)  # For pretty print on write
        try:
            # this handles either filenames or open file-like objects
            tree = etree.parse(source, parser)
        except IOError:
            # Check if the file exists in the data/forcefield directory
            temp_file = get_data_filename(os.path.join('forcefield', source))
            tree = etree.parse(temp_file, parser)
        #except Exception: # If it's a string
        #    string_data = source.read()
        #    tree = self.parse_string(source)
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
            raise Exception(msg)
        self.from_lxml(tree.getroot())

    def parse_string(self, data):
        """Parse a SMIRNOFF force field definition in XML format.

        A ``ParseError`` is raised if the XML cannot be processed.

        Parameters
        ----------
        data : str
            A SMIRNOFF force field definition in `the SMIRNOFF XML format <https://github.com/openforcefield/openforcefield/blob/master/The-SMIRNOFF-force-field-format.md>`_.

        """

        # Parse XML file
        root = etree.XML(data)
        #root = etree.fromstring(data)
        self.from_lxml(root)

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
        if (format == 'offxml') or (extension == '.offxml'):
            root.write(filename, xml_declaration=True, pretty_print=True)
        else:
            msg = "Cannot export forcefield parameters to file '{}'\n".format(
                filename)
            msg += 'Export to extension {} not implemented yet.\n'.format(
                extension)
            msg += "Supported choices are: ['.offxml']"
            raise NotImplementedError(msg)

    def to_string(self):
        """

        Returns
        -------

        """
        return self.to_xml()

    #@staticmethod
    #def from_xml(xml):
    #    """Construct a ForceField object from a SMIRNOFF XML file.
    #    """
    #    return ForceField(xml)

    def to_lxml(self):
        """Render the forcefield this ParameterIOHandler is registered to to an lxml.etree

        Returns
        -------
        root : lxml.etree
            Root node
        """
        root = etree.Element('SMIRNOFF', {'version': self._forcefield.version})
        for tagname, parameter_handler in self._forcefield._parameter_handlers.items(
        ):
            parameter_subtree = etree.SubElement(root, tagname,
                                                 parameter_handler.attribs)
            for parameter in parameter_handler.parameters:
                etree.SubElement(parameter_subtree, parameter.tagname,
                                 parameter.attribs)

        return root

    def from_lxml(self, root):
        if not (root.tag == 'SMIRNOFF'):
            raise Exception("Root tag of tree is not 'SMIRNOFF'")

        cosmetic_tags = ['Date', 'Author']
        #root = smirnoff_root[0]
        #raise Exception(root)
        try:
            exception_node = root  # node used for exception reporting

            # Process handlers
            #for section in root.iter(tag=etree.Element):
            for section in root:
                # Skip comment lines
                if not isinstance(section.tag, str):
                    continue
                if section.tag in cosmetic_tags:
                    continue
                exception_node = section  # node used for exception reporting
                # Extract parameter name from XML tag
                parameter_name = section.tag

                # Split out attributes that refer to units
                handler_kwargs, attached_units = _extract_attached_units(
                    section.attrib)
                # TODO: Attach units to handler_kwargs
                # Make a copy of handler_kwargs that we can modify
                handler_kwargs_dict = dict(handler_kwargs)
                handler_kwargs_dict = _attach_units(handler_kwargs_dict,
                                                    attached_units)

                # Retrieve or create parameter handler
                handler = self._forcefield.get_handler(parameter_name,
                                                       handler_kwargs_dict)

                # Populate handler with parameter definitions
                for parameter in section:
                    # Skip comment lines
                    if not isinstance(parameter.tag, str):
                        continue
                    exception_node = parameter  # node used for exception reporting

                    # parameter.attrib doesn't support assignment of Quantity type, so make a copy as dict
                    parameter_attrib_dict = dict(parameter.attrib)

                    # Append units to parameters as needed
                    parameter_kwargs = _attach_units(parameter_attrib_dict,
                                                     attached_units)
                    # Add parameter definition
                    handler.add_parameter(parameter_kwargs)

        except Exception as e:
            # Prepend the line and line text of XML file to aid debugging
            # TODO: Can we include the filename as well?
            print(self._get_sourceline(exception_node))
            #if not e.args:
            #    e.args = ('',)
            #e.args = e.args[0] + + e.args[1:]
            raise e

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
