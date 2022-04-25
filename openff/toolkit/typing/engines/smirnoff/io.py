#!/usr/bin/env python

# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================
"""
XML I/O parser for the SMIRNOFF (SMIRKS Native Open Force Field) format.

.. codeauthor:: John D. Chodera <john.chodera@choderalab.org>
.. codeauthor:: David L. Mobley <dmobley@mobleylab.org>
.. codeauthor:: Peter K. Eastman <peastman@stanford.edu>

"""

__all__ = [
    "ParameterIOHandler",
    "XMLParameterIOHandler",
]

import logging
from typing import Optional

import xmltodict

# =============================================================================================
# CONFIGURE LOGGER
# =============================================================================================

logger = logging.getLogger(__name__)


# =============================================================================================
# Base ParameterIOHandler
# =============================================================================================


class ParameterIOHandler:
    """
    Base class for handling serialization/deserialization of SMIRNOFF ForceField objects
    """

    _FORMAT: Optional[str] = None

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
        Write the current force field parameter set to a file.

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
        Render the force field parameter set to a string

        Parameters
        ----------
        smirnoff_data : dict
            A dictionary structured in compliance with the SMIRNOFF spec

        Returns
        -------
        str
        """
        pass


# =============================================================================================
# XML I/O
# =============================================================================================


class XMLParameterIOHandler(ParameterIOHandler):
    """
    Handles serialization/deserialization of SMIRNOFF ForceField objects from OFFXML format.
    """

    # TODO: Come up with a better keyword for format
    _FORMAT = "XML"

    def parse_file(self, source):
        """Parse a SMIRNOFF force field definition in XML format, read from a file.

        Parameters
        ----------
        source : str or io.RawIOBase
            File path of file-like object implementing a ``read()`` method
            specifying a SMIRNOFF force field definition in `the SMIRNOFF XML format
            <https://openforcefield.github.io/standards/standards/smirnoff/#xml-representation>`_.

        Raises
        ------
        SMIRNOFFParseError
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

        A ``SMIRNOFFParseError`` is raised if the XML cannot be processed.

        Parameters
        ----------
        data : str
            A SMIRNOFF force field definition in `the SMIRNOFF XML format
            <https://openforcefield.github.io/standards/standards/smirnoff/#xml-representation>`_.

        """
        from pyexpat import ExpatError

        from openff.toolkit.utils.exceptions import SMIRNOFFParseError

        # Parse XML file
        try:
            smirnoff_data = xmltodict.parse(data, attr_prefix="")
            return smirnoff_data
        except ExpatError as e:
            raise SMIRNOFFParseError(str(e))

    def to_file(self, file_path, smirnoff_data):
        """Write the current force field parameter set to a file.

        Parameters
        ----------
        file_path : str
            The path to the file to be written.
            The `.offxml` or `.xml` file extension must be present.
        smirnoff_data : dict
            A dict structured in compliance with the SMIRNOFF data spec.

        """
        xml_string = self.to_string(smirnoff_data)
        with open(file_path, "w") as of:
            of.write(xml_string)

    def to_string(self, smirnoff_data):
        """
        Write the current force field parameter set to an XML string.

        Parameters
        ----------
        smirnoff_data : dict
            A dictionary structured in compliance with the SMIRNOFF spec

        Returns
        -------
        serialized_forcefield : str
            XML String representation of this force field.

        """

        def prepend_all_keys(d, char="@", ignore_keys=frozenset()):
            """
            Modify a dictionary in-place, prepending a specified string to each key
            that doesn't refer to a value that is list or dict.

            Parameters
            ----------
            d : dict
                Hierarchical dictionary to traverse and modify keys
            char : string, optional. Default='@'
                String to prepend onto each applicable dictionary key
            ignore_keys : iterable of str
                A set or list of strings, indicating keys not to prepend in the data structure

            """
            if isinstance(d, dict):
                for key in list(d.keys()):
                    if key in ignore_keys:
                        continue
                    if isinstance(d[key], list) or isinstance(d[key], dict):
                        prepend_all_keys(d[key], char=char, ignore_keys=ignore_keys)
                    else:
                        new_key = char + key
                        d[new_key] = d[key]
                        del d[key]
                        prepend_all_keys(d[new_key], char=char, ignore_keys=ignore_keys)
            elif isinstance(d, list):
                for item in d:
                    prepend_all_keys(item, char=char, ignore_keys=ignore_keys)

        # the "xmltodict" library defaults to print out all element attributes on separate lines
        # unless they're prepended by "@"
        prepend_all_keys(smirnoff_data["SMIRNOFF"], ignore_keys=["Author", "Date"])

        # Reorder parameter sections to put Author and Date at the top (this is the only
        # way to change the order of items in a dict, as far as I can tell)
        for key, value in list(smirnoff_data["SMIRNOFF"].items()):
            if key in ["Author", "Date"]:
                continue
            del smirnoff_data["SMIRNOFF"][key]
            smirnoff_data["SMIRNOFF"][key] = value

        return xmltodict.unparse(smirnoff_data, pretty=True, indent=" " * 4)
