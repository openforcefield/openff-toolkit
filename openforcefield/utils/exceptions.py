#!/usr/bin/env python

# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
An API to store general exceptions that may be raised.

Authors
-------
* Simon Boothroyd <simon.boothroyd@choderalab.org>
"""


# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================


# =============================================================================================
# Module Constants
# =============================================================================================

class XmlNodeMissingException(Exception):

    def __init__(self, node_name):

        message = 'The calculation template does not contain a <' + str(node_name) + '> node.'
        super().__init__(message)

