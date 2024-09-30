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
    "DEFAULT_AROMATICITY_MODEL",
    "ALLOWED_AROMATICITY_MODELS",
    "DEFAULT_FRACTIONAL_BOND_ORDER_MODEL",
    "ALLOWED_FRACTIONAL_BOND_ORDER_MODELS",
    "DEFAULT_CHARGE_MODEL",
    "ALLOWED_CHARGE_MODELS",
    "IncompatibleUnitError",
    "MissingPackageError",
    "ToolkitUnavailableException",
    "LicenseError",
    "InvalidToolkitError",
    "InvalidToolkitRegistryError",
    "UndefinedStereochemistryError",
    "GAFFAtomTypeWarning",
    "ToolkitWrapper",
    "OpenEyeToolkitWrapper",
    "RDKitToolkitWrapper",
    "AmberToolsToolkitWrapper",
    "NAGLToolkitWrapper",
    "BuiltInToolkitWrapper",
    "ChargeMethodUnavailableError",
    "IncorrectNumConformersError",
    "IncorrectNumConformersWarning",
    "ChargeCalculationError",
    "InvalidIUPACNameError",
    "AntechamberNotFoundError",
    "SMILESParseError",
    "ToolkitRegistry",
    "GLOBAL_TOOLKIT_REGISTRY",
    "OPENEYE_AVAILABLE",
    "RDKIT_AVAILABLE",
    "AMBERTOOLS_AVAILABLE",
    "BASIC_CHEMINFORMATICS_TOOLKITS",
)


import logging

from openff.toolkit.utils.ambertools_wrapper import AmberToolsToolkitWrapper
from openff.toolkit.utils.base_wrapper import ToolkitWrapper
from openff.toolkit.utils.builtin_wrapper import BuiltInToolkitWrapper
from openff.toolkit.utils.constants import (
    ALLOWED_AROMATICITY_MODELS,
    ALLOWED_CHARGE_MODELS,
    ALLOWED_FRACTIONAL_BOND_ORDER_MODELS,
    DEFAULT_AROMATICITY_MODEL,
    DEFAULT_CHARGE_MODEL,
    DEFAULT_FRACTIONAL_BOND_ORDER_MODEL,
)
from openff.toolkit.utils.exceptions import (
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
    MissingPackageError,
    SMILESParseError,
    ToolkitUnavailableException,
    UndefinedStereochemistryError,
)
from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
from openff.toolkit.utils.openeye_wrapper import OpenEyeToolkitWrapper
from openff.toolkit.utils.rdkit_wrapper import RDKitToolkitWrapper
from openff.toolkit.utils.toolkit_registry import (
    ToolkitRegistry,
    toolkit_registry_manager,
)

logger = logging.getLogger(__name__)


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
"""
The toolkit registry used by default when no registry or wrapper is specified.

A toolkit registry is a list of toolkit wrappers used to provide functionality
not implemented by the OpenFF Toolkit itself. Any given functionality is
requested from each toolkit in order until one is found that provides it. For
more information, see :py:class:`ToolkitRegistry`.

To temporarily modify the ``GLOBAL_TOOLKIT_REGISTRY``,
see :py:func:`toolkit_registry_manager`.

Examples
========

Check the default toolkit precedence order:

>>> from openff.toolkit import GLOBAL_TOOLKIT_REGISTRY
>>> GLOBAL_TOOLKIT_REGISTRY.registered_toolkits() # doctest: +ELLIPSIS
[...]

Temporarily change the global toolkit registry:

>>> from openff.toolkit import ToolkitRegistry, RDKitToolkitWrapper
>>> from openff.toolkit.utils import toolkit_registry_manager
>>> with toolkit_registry_manager(ToolkitRegistry([RDKitToolkitWrapper]))):
>>>     GLOBAL_TOOLKIT_REGISTRY
<ToolkitRegistry containing The RDKit>
>>>     mol = Molecule.from_smiles("CCO")
>>> # when the context manager ends, the default registry is back to normal
>>> GLOBAL_TOOLKIT_REGISTRY
<ToolkitRegistry containing The RDKit, AmberTools, Built-in Toolkit>

To change the value of ``GLOBAL_TOOLKIT_REGISTRY`` for the remainder of the
Python process, individual toolkits can be added or removed with
the :py:meth:`ToolkitRegistry.deregister_toolkit`
and :py:meth:`ToolkitRegistry.register_toolkit` methods. This can be useful for
debugging and exploring subtly different behavior between toolkit wrappers:

>>> print(len(GLOBAL_TOOLKIT_REGISTRY.registered_toolkits))
4
>>> GLOBAL_TOOLKIT_REGISTRY.deregister_toolkit(RDKitToolkitWrapper)
>>> print(len(GLOBAL_TOOLKIT_REGISTRY.registered_toolkits))
3
>>> GLOBAL_TOOLKIT_REGISTRY.register_toolkit(RDKitToolkitWrapper)
>>> print(len(GLOBAL_TOOLKIT_REGISTRY.registered_toolkits))
4

Note that as with other global attributes in Python, assigning a new toolkit
registry to ``GLOBAL_TOOLKIT_REGISTRY`` is quite difficult to get right and very
likely to fail silently - we recommend modifying the existing value instead in
most cases.

See Also
========
toolkit_registry_manager, ToolkitRegistry, ToolkitWrapper
"""


OPENEYE_AVAILABLE = False
RDKIT_AVAILABLE = False
AMBERTOOLS_AVAILABLE = False

with toolkit_registry_manager(
    ToolkitRegistry(
        toolkit_precedence=[
            OpenEyeToolkitWrapper,
            RDKitToolkitWrapper,
            AmberToolsToolkitWrapper,
        ],
        exception_if_unavailable=False,
    )
):
    for toolkit in GLOBAL_TOOLKIT_REGISTRY.registered_toolkits:
        if type(toolkit) is OpenEyeToolkitWrapper:
            OPENEYE_AVAILABLE = True
        elif type(toolkit) is RDKitToolkitWrapper:
            RDKIT_AVAILABLE = True
        elif type(toolkit) is AmberToolsToolkitWrapper:
            AMBERTOOLS_AVAILABLE = True


# Define basic toolkits that handle essential file I/O

BASIC_CHEMINFORMATICS_TOOLKITS = [RDKitToolkitWrapper, OpenEyeToolkitWrapper]

# Ensure we have at least one basic toolkit
any_toolkits = False
for tk in GLOBAL_TOOLKIT_REGISTRY.registered_toolkits:
    if type(tk) in BASIC_CHEMINFORMATICS_TOOLKITS:
        if tk.is_available():
            any_toolkits = True
            break

if not any_toolkits:  # pragma: no cover
    from openff.toolkit.utils.utils import all_subclasses

    msg = "WARNING: No basic cheminformatics toolkits are available.\n"
    msg += "At least one basic toolkit is required to handle SMARTS matching and file I/O. \n"
    msg += "Please install at least one of the following basic toolkits:\n"
    for wrapper in all_subclasses(ToolkitWrapper):
        if wrapper.toolkit_name is not None:
            msg += f"{wrapper._toolkit_name} : {wrapper._toolkit_installation_instructions}\n"
    # TODO: Make this a warning!?
    print(msg)
