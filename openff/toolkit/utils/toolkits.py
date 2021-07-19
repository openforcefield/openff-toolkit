#!/usr/bin/env python
"""
Wrapper classes for providing a minimal consistent interface to cheminformatics toolkits

Currently supported toolkits:

* The `OpenEye Toolkit <https://docs.eyesopen.com/toolkits/python/quickstart-python/index.html>`_
* The `RDKit <http://www.rdkit.org/>`_
* `AmberTools <http://ambermd.org/AmberTools.php>`_

.. todo::

   * Add checks at the beginning of each toolkit method call to make sure toolkit is licened
   * Switch toolkit methods to object methods instead of static methods
   * Should this be under ``openff.toolkit.utils.toolkits`` or ``openff.toolkit.toolkits``?
   * Add singleton global toolkit registry that registers all available toolkits by default
        when this file is imported
   * Add description fields for each toolkit wrapper
   * Eliminate global variables in favor of a singleton pattern
   * Change global variables from _INSTALLED to _AVAILABLE

"""
__all__ = (
    # constants
    "DEFAULT_AROMATICITY_MODEL",
    "ALLOWED_AROMATICITY_MODELS",
    "DEFAULT_FRACTIONAL_BOND_ORDER_MODEL",
    "ALLOWED_FRACTIONAL_BOND_ORDER_MODELS",
    "DEFAULT_CHARGE_MODEL",
    "ALLOWED_CHARGE_MODELS",
    # exceptions and warnings
    "MessageException",
    "IncompatibleUnitError",
    "MissingDependencyError",
    "MissingPackageError",
    "ToolkitUnavailableException",
    "LicenseError",
    "InvalidToolkitError",
    "InvalidToolkitRegistryError",
    "UndefinedStereochemistryError",
    "GAFFAtomTypeWarning",
    "ChargeMethodUnavailableError",
    "IncorrectNumConformersError",
    "IncorrectNumConformersWarning",
    "ChargeCalculationError",
    "InvalidIUPACNameError",
    "AntechamberNotFoundError",
    "ParseError",
    # base_wrapper
    "ToolkitWrapper",
    # builtin_wrapper
    "BuiltInToolkitWrapper",
    # openeye_wrapper
    "OpenEyeToolkitWrapper",
    # rdkit_wrapper
    "RDKitToolkitWrapper",
    # ambertools_wrapper
    "AmberToolsToolkitWrapper",
    # toolkit_registry
    "ToolkitRegistry",
    # in this module
    "GLOBAL_TOOLKIT_REGISTRY",
    "OPENEYE_AVAILABLE",
    "RDKIT_AVAILABLE",
    "AMBERTOOLS_AVAILABLE",
    "BASIC_CHEMINFORMATICS_TOOLKITS",
)


# =============================================================================================
# IMPORTS
# =============================================================================================

import logging

from .ambertools_wrapper import AmberToolsToolkitWrapper
from .base_wrapper import ToolkitWrapper
from .builtin_wrapper import BuiltInToolkitWrapper
from .constants import (
    ALLOWED_AROMATICITY_MODELS,
    ALLOWED_CHARGE_MODELS,
    ALLOWED_FRACTIONAL_BOND_ORDER_MODELS,
    DEFAULT_AROMATICITY_MODEL,
    DEFAULT_CHARGE_MODEL,
    DEFAULT_FRACTIONAL_BOND_ORDER_MODEL,
)
from .exceptions import (
    AntechamberNotFoundError,
    ChargeCalculationError,
    ChargeMethodUnavailableError,
    GAFFAtomTypeWarning,
    IncompatibleUnitError,
    IncorrectNumConformersError,
    IncorrectNumConformersWarning,
    InvalidIUPACNameError,
    InvalidToolkitError,
    InvalidToolkitRegistryError,
    LicenseError,
    MessageException,
    MissingDependencyError,
    MissingPackageError,
    ParseError,
    ToolkitUnavailableException,
    UndefinedStereochemistryError,
)
from .openeye_wrapper import OpenEyeToolkitWrapper, requires_openeye_module
from .rdkit_wrapper import RDKitToolkitWrapper
from .toolkit_registry import ToolkitRegistry

# =============================================================================================
# CONFIGURE LOGGER
# =============================================================================================

logger = logging.getLogger(__name__)


# =============================================================================================
# GLOBAL TOOLKIT REGISTRY
# =============================================================================================

# Create global toolkit registry, where all available toolkits are registered
GLOBAL_TOOLKIT_REGISTRY = ToolkitRegistry(
    toolkit_precedence=[
        OpenEyeToolkitWrapper,
        RDKitToolkitWrapper,
        AmberToolsToolkitWrapper,
        BuiltInToolkitWrapper,
    ],
    exception_if_unavailable=False,
)

# =============================================================================================
# SET GLOBAL TOOLKIT-AVAIABLE VARIABLES
# =============================================================================================

OPENEYE_AVAILABLE = False
RDKIT_AVAILABLE = False
AMBERTOOLS_AVAILABLE = False

# Only available toolkits will have made it into the GLOBAL_TOOLKIT_REGISTRY
for toolkit in GLOBAL_TOOLKIT_REGISTRY.registered_toolkits:
    if type(toolkit) is OpenEyeToolkitWrapper:
        OPENEYE_AVAILABLE = True
    elif type(toolkit) is RDKitToolkitWrapper:
        RDKIT_AVAILABLE = True
    elif type(toolkit) is AmberToolsToolkitWrapper:
        AMBERTOOLS_AVAILABLE = True

# =============================================================================================
# WARN IF INSUFFICIENT TOOLKITS INSTALLED
# =============================================================================================

# Define basic toolkits that handle essential file I/O

BASIC_CHEMINFORMATICS_TOOLKITS = [RDKitToolkitWrapper, OpenEyeToolkitWrapper]

# Ensure we have at least one basic toolkit
if (
    sum(
        [
            tk.is_available()
            for tk in GLOBAL_TOOLKIT_REGISTRY.registered_toolkits
            if type(tk) in BASIC_CHEMINFORMATICS_TOOLKITS
        ]
    )
    == 0
):
    from .utils import all_subclasses

    msg = "WARNING: No basic cheminformatics toolkits are available.\n"
    msg += "At least one basic toolkit is required to handle SMARTS matching and file I/O. \n"
    msg += "Please install at least one of the following basic toolkits:\n"
    for wrapper in all_subclasses(ToolkitWrapper):
        if wrapper.toolkit_name is not None:
            msg += "{} : {}\n".format(
                wrapper._toolkit_name, wrapper._toolkit_installation_instructions
            )
    print(msg)
