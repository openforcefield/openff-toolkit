#!/usr/bin/env python

"""
Utility subroutines for managing cheminformatics toolkits

.. todo::

   * Generalize this infrastructure to make it easier to support additional cheminformatics toolkits.

"""
#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import importlib
from functools import wraps

#=============================================================================================
# CHEMINFORMATICS TOOLKITS
#=============================================================================================

# TODO: Generalize this infrastructure to make it easier to support additional toolkits in future

# Control the precedence order in which cheminformatics toolkits are used
TOOLKIT_PRECEDENCE = ['openeye', 'rdkit']

# List of supported toolkits and messages indicating how they can be installed
SUPPORTED_TOOLKITS = {
    'rdkit' : 'A conda-installable version of the free and open source RDKit cheminformatics toolkit can be found at: https://anaconda.org/rdkit/rdkit',
    'openeye' : 'The OpenEye toolkit requires a (free for academics) license, and can be found at: https://docs.eyesopen.com/toolkits/python/quickstart-python/install.html'
}

class LicenseError(Exception):
    """This function requires a license that cannot be found."""
    pass

class MissingPackageError(Exception):
    """This function requires a package that is not installed."""
    pass

# TODO: Differentiate between is installed and is licensed.
def is_openeye_installed(oetools=('oechem', 'oequacpac', 'oeiupac', 'oeomega')):
    """
    Check if a given OpenEye tool is installed and Licensed

    If the OpenEye toolkit is not installed, returns False

    Parameters
    ----------
    oetools : str or iterable of strings, Optional, Default: ('oechem', 'oequacpac', 'oeiupac', 'oeomega')
        Set of tools to check by their string name. Defaults to the complete set that YANK *could* use, depending on
        feature requested.

        Only checks the subset of tools if passed. Also accepts a single tool to check as a string instead of an
        iterable of length 1.

    Returns
    -------
    all_installed : bool
        True if all tools in ``oetools`` are installed and licensed, False otherwise
    """
    # Complete list of module: License check
    tools_license = {
        'oechem': 'OEChemIsLicensed',
        'oequacpac': 'OEQuacPacIsLicensed',
        'oeiupac': 'OEIUPACIsLicensed',
        'oeomega': 'OEOmegaIsLicensed'
        }
    tool_keys = tools_license.keys()

    # Cast oetools to tuple if its a single string
    if type(oetools) is str:
        oetools = (oetools,)
    tool_set = set(oetools)
    valid_tool_set = set(tool_keys)
    if tool_set & valid_tool_set == set():
        # Check for empty set intersection
        raise ValueError("Expected OpenEye tools to have at least of the following {}, "
                         "but instead got {}".format(tool_keys, oetools))
    try:
        for tool in oetools:
            if tool in tool_keys:
                # Try loading the module
                try:
                    module = importlib.import_module('openeye', tool)
                except SystemError: # Python 3.4 relative import fix
                    module = importlib.import_module('openeye.' + tool)
                # Check that we have the license
                if not getattr(module, tools_license[tool])():
                    raise ImportError
    except ImportError:
        return False
    return True

# TODO: Add a method to raise the appropriate exception if OpenEye is not installed or licensed, like assert_openeye_available, assert_rdkit_available

def requires_openeye(*oetools):
    """
    Decorator to check that OpenEye licenses are found, raising LicenseError if valid license not found.

    """
    def decorator(func):
        @wraps(func)
        def wrapped_function(*args, **kwargs):
            if not is_openeye_installed(oetools=oetools):
                msg = 'This function requires the OpenEye toolkit with licenses for the following tools: {}'.format(oetools)
                raise LicenseError(msg)
        return wrapped_function
    return decorator

OPENEYE_INSTALLED = is_openeye_installed('oechem') and is_openeye_installed('oequacpac') and is_openeye_installed('oeiupac') and is_openeye_installed('oeomega')
OPENEYE_UNAVAILABLE = not OPENEYE_INSTALLED

def is_rdkit_installed():
    try:
        module = importlib.import_module('rdkit', 'Chem')
        return True
    except ImportError:
        return False

RDKIT_INSTALLED = is_rdkit_installed()
RDKIT_UNAVAILABLE = not RDKIT_INSTALLED

def requires_rdkit():
    """
    Decorator to check that rdkit is installed.

    """
    def decorator(func):
        @wraps(func)
        def wrapped_function(*args, **kwargs):
            if not is_rdkit_installed():
                msg = 'This function requires the RDKit toolkit'
                # TODO: Differentiate between is installed and is licensed with MissingPackageError and LicenseError
                raise MissingPackageError(msg)
        return wrapped_function
    return decorator

def is_toolkit_installed():
    if is_openeye_installed():
        return True
    if is_rdkit_installed():
        return True
    return False

TOOLKIT_INSTALLED = is_toolkit_installed()
TOOLKIT_UNAVAILABLE = not TOOLKIT_INSTALLED

def toolkit_is_available(toolkit_name):
    """Return True if the requested toolkit is available.
    """
    if toolkit_name == 'openeye':
        return is_openeye_installed()
    elif toolkit_name == 'rdkit':
        return is_rdkit_installed()
    else:
        raise Exception('Toolkit {} unknown; options are {}'.format(toolkit_name, TOOLKIT_PRECEDENCE))

# Filter toolkits by what is available
TOOLKIT_PRECEDENCE = [ toolkit for toolkit in TOOLKIT_PRECEDENCE if toolkit_is_available(toolkit) ]

# Warn if no toolkits are available
if TOOLKIT_UNAVAILABLE:
    msg = 'WARNING: No cheminfomatics toolkits are available.\n'
    for (toolkit_name, install_instructions) in SUPPORTED_TOOLKITS.items():
        msg += 'Please install one of the following toolkits:\n'
        msg += '{} : {}\n'.format(toolkit_name, install_instructions)
    print(msg)
