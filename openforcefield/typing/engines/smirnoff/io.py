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

import os
import logging
from collections import OrderedDict

import lxml.etree as etree

from simtk import openmm, unit

#=============================================================================================
# CONFIGURE LOGGER
#=============================================================================================

logger = logging.getLogger(__name__)

#=============================================================================================
#
#=============================================================================================

@staticmethod
def _extract_attached_units(attrib):
    """Form a (potentially unit-bearing) quantity from the specified attribute name.

    Parameters
    ----------
    attrib : dict
       Dictionary of XML node attributes.

    Returns
    -------
    attrib : dict
       XML node attributes with keys ending in ``_unit`` removed.
    attached_units : dict str : simtk.unit.Unit
       ``attached_units[parameter_name]`` is the simtk.unit.Unit combination that should be attached to corresponding
       parameter ``parameter_name``

    """
    # TODO: Should this scheme also convert unit-bearing quantities such as '8*angstroms' to 8*unit.angstroms?
    # TODO: Should this scheme also convert "1" to int(1) and "8.0" to float(8.0)?

    attached_units = OrderedDict()
    for key in list(attrib.keys()):
        if key.endswith('_unit'):
            parameter_name = key[:-5]
            parameter_units_string = attrib[key]
            try:
                parameter_units = eval(parameter_units_string, unit.__dict__)
            except Exception as e:
                e.msg = "Could not parse units {}\n".format(parameter_units_string) + e.msg
                raise e
            attached_units[parameter_name] = parameter_units
            del attrib[parameter_name]

    return attrib, attached_units

@staticmethod
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
    try:
        for (parameter_name, parameter_value_string) in attrib.items():
            if parameter_name in attached_units:
                units_to_attach = attached_units[parameter_name]
                # TODO: Do we have to worry about None or null values for parameters with attached units?
                attrib[parameter_name] = float(parameter_value_string) * units_to_attach
    except Exception as e:
        e.msg = "Expected numeric value for parameter '{}', instead found '{}' when trying to attach units '{}'\n".format(parameter_name, parameter_value_string, units_to_attach)
        raise e

    return attrib

#=============================================================================================
# XML I/O
#=============================================================================================

# TODO: Instead of subclassing ForceField, what if we had an abstract base class for parameter IO?

class XMLForceField(ForceField):

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
    def parse_xml_filename(self, filename):
        """Parse a SMIRNOFF force field definition in XML format, read from a file.

        A ``ParseError`` is raised if the XML cannot be processed.

        Parameters
        ----------
        filename : str
            File path specifying a SMIRNOFF force field definition in `the SMIRNOFF XML format <https://github.com/openforcefield/openforcefield/blob/master/The-SMIRNOFF-force-field-format.md>`_.

        .. notes ::

           * New SMIRNOFF sections are handled independently, as if they were specified in the same file.
           * If a SMIRNOFF section that has already been read appears again, its definitions are appended to the end of the previously-read
             definitions if the sections are configured with compatible attributes; otherwise, an ``IncompatibleTagException`` is raised.
        """

        # TODO: Search the right sequence of paths

        parser = etree.XMLParser(remove_blank_text=True) # For pretty print on write
        try:
            # this handles either filenames or open file-like objects
            tree = etree.parse(file, parser)
        except IOError:
            # Check if the file exists in an installed directory
            temp_file = get_data_filename(file)
            tree = etree.parse(temp_file, parser)
        except Exception as e:
            # Fail with an error message about which file could not be read.
            # TODO: Also handle case where fallback to 'data' directory encounters problems,
            # but this is much less worrisome because we control those files.
            msg  = str(e) + '\n'
            if hasattr(file, 'name'):
                filename = file.name
            else:
                filename = str(file)
            msg += "ForceField.loadFile() encountered an error reading file '%s'\n" % filename
            raise Exception(msg)

        # TODO: Run through parsers

    def parse_xml_string(self, xml, filename=None):
        """Parse a SMIRNOFF force field definition in XML format.

        A ``ParseError`` is raised if the XML cannot be processed.

        Parameters
        ----------
        xml : str
            A SMIRNOFF force field definition in `the SMIRNOFF XML format <https://github.com/openforcefield/openforcefield/blob/master/The-SMIRNOFF-force-field-format.md>`_.
        filename : str, optional, default=None
            If specified, will include the filename in exceptions to aid debugging.

        """

        # Parse XML file
        root = etree.XML(xml)

        try:
            # Create new ForceField with provided attributes
            exception_node = root # node used for exception reporting
            forcefield = ForceField(**root.attrib)

            # Process handlers
            for section in root:
                exception_node = section # node used for exception reporting
                # Extract parameter name from XML tag
                parmeter_name = section.tag

                # Split out attributes that refer to units
                handler_kwargs, attached_units = _extract_attached_units(section.attrib)

                # Retrieve or create parameter handler
                handler = forcefield.get_handler(section.tag, handler_kwargs)

                # Populate handler with parameter definitions
                for parameter in section:
                    exception_node = parameter # node used for exception reporting
                    # Append units to parameters as needed
                    parameter_kwargs = _append_units(parameter.attrib, attached_units)
                    # Add parameter definition
                    handler.add_parameter(**parameter_kwargs)

        except Exception as e:
            # Prepend the line and line text of XML file to aid debugging
            # TODO: Can we include the filename as well?
            e.msg = self._get_sourceline(exception_node, filename=filename) + e.msg
            raise e

    def save(self, filename, format=None):
        """Write the current forcefield parameter set to a file, autodetecting the type from the extension.

        Parameters
        ----------
        filename : str
            The name of the file to be written.
            The `.offxml` file extension is auto-detected.
        format : str, optional, default=None
            If specified, will write in the specified format.
            Options: ['offxml']
        """
        (basename, extension) = os.path.splitext(filename)
        if (format == 'offxml') or (extension == '.offxml'):
            tree.write(filename, xml_declaration=True, pretty_print=True)
        else:
            msg = "Cannot export forcefield parameters to file '{}'\n".format(filename)
            msg += 'Export to extension {} not implemented yet.\n'.format(extension)
            msg += "Supported choices are: ['.offxml']"
            raise NotImplementedError(msg)

    @staticmethod
    def from_xml(xml):
        """Construct a ForceField object from a SMIRNOFF XML file.
        """
        return ForceField(xml)

    def to_lxml(self):
        """Render the forcefield parameter set to an lxml.etree

        Returns
        -------
        root : lxml.etree
            Root node
        """
        root = etree.Element('SMIRNOFF', self.attrib)
        for tagname, parameter_handler in self.parameter_handers.items():
            parameter_subtree = etree.SubElement(root, tagname, parameter_handler.attribs)
            for parameter in parameter_handler.parameters:
                etree.SubElement(parameter_subtree, parameter.tagname, parameter.attribs)

        return root

    def from_lxml(self):
        pass

    def to_xml(self):
        """Render the forcefield parameter set to XML.

        Returns
        -------
        xml : str
            The SMIRNOFF parameter set rendered as XML.
        """
        return etree.tostring(self.to_lxml())

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
