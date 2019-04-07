#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Test classes and function in module openforcefield.typing.engines.smirnoff.io.

"""


#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import pytest


#=============================================================================================
# QUANTITY PARSING UTILITIES
#=============================================================================================

class TestXMLParameterIOHandler:

    def test_from_string(self):
        pass

    #
    # Tests for ForceField writing to XML files
    #

    # TODO: Remove ForceField from this whole file. All tests should be for converting between hierarchical SMIRNOFF
    #       dicts and XML
    @pytest.mark.skip(reason='Needs to be updated for 1.0.0 syntax')
    def test_save(self):
        """Test writing and reading of SMIRNOFF in XML format.
        """
        forcefield = ForceField(smirnoff99Frosst_offxml_filename)
        # Write XML to a file
        with TemporaryDirectory() as tmpdir:
            offxml_tmpfile = os.path.join(tmpdir, 'forcefield.offxml')
            forcefield.save(offxml_tmpfile)
            forcefield2 = ForceField(offxml_tmpfile)
            assert_forcefields_equal(cls.forcefield, forcefield2,
                                     "ForceField written to .offxml does not match original ForceField")

    @pytest.mark.skip(reason='Needs to be updated for 1.0.0 syntax')
    def test_to_xml(self):
        forcefield = ForceField(smirnoff99Frosst_offxml_filename)
        # Retrieve XML as a string
        xml = forcefield.to_xml()
        # Restore ForceField from XML
        forcefield2 = ForceField(xml)
        assert_forcefields_equal(cls.forcefield, forcefield2,
                                 "ForceField serialized to XML does not match original ForceField")


