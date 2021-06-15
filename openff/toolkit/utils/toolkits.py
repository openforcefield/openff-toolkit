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
   * Add singleton global toolkit registry that registers all available toolkits by default when this file is imported
   * Add description fields for each toolkit wrapper
   * Eliminate global variables in favor of a singleton pattern
   * Change global variables from _INSTALLED to _AVAILABLE

"""

__all__ = [
    "DEFAULT_AROMATICITY_MODEL",
    "ALLOWED_AROMATICITY_MODELS",
    "DEFAULT_FRACTIONAL_BOND_ORDER_MODEL",
    "ALLOWED_FRACTIONAL_BOND_ORDER_MODELS",
    "DEFAULT_CHARGE_MODEL",
    "ALLOWED_CHARGE_MODELS",
    "LicenseError",
    "MissingPackageError",
    "ToolkitUnavailableException",
    "InvalidToolkitError",
    "InvalidToolkitRegistryError",
    "UndefinedStereochemistryError",
    "GAFFAtomTypeWarning",
    "ToolkitWrapper",
    "BuiltInToolkitWrapper",
    "OpenEyeToolkitWrapper",
    "RDKitToolkitWrapper",
    "AmberToolsToolkitWrapper",
    "BuiltInToolkitWrapper",
    "ToolkitRegistry",
    "GLOBAL_TOOLKIT_REGISTRY",
    "OPENEYE_AVAILABLE",
    "RDKIT_AVAILABLE",
    "AMBERTOOLS_AVAILABLE",
    "BASIC_CHEMINFORMATICS_TOOLKITS",
]


# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

import copy
import importlib
import inspect
import itertools
import logging
import re
import subprocess
import tempfile
from collections import defaultdict
from functools import wraps
from typing import TYPE_CHECKING, List, Optional, Tuple

import numpy as np
from simtk import unit

from openff.toolkit.utils.utils import (
    MessageException,
    all_subclasses,
    inherit_docstrings,
    temporary_cd,
)

if TYPE_CHECKING:
    from openforcefield.topology.molecule import Molecule

# =============================================================================================
# CONFIGURE LOGGER
# =============================================================================================

logger = logging.getLogger(__name__)

# =============================================================================================
# SUPPORTED MODELS
#
# TODO: We may no longer need these since we now require SMIRNOFF to specify these models explicitly.
# =============================================================================================

DEFAULT_AROMATICITY_MODEL = "OEAroModel_MDL"  # TODO: Is there a more specific name and reference for the aromaticity model?
ALLOWED_AROMATICITY_MODELS = ["OEAroModel_MDL"]

DEFAULT_FRACTIONAL_BOND_ORDER_MODEL = "Wiberg"  # TODO: Is there a more specific name and reference for the fractional bond order models?
ALLOWED_FRACTIONAL_BOND_ORDER_MODELS = ["Wiberg"]

DEFAULT_CHARGE_MODEL = "AM1-BCC"  # TODO: Should this be `AM1-BCC`, or should we encode BCCs explicitly via AM1-CM2 preprocessing?
ALLOWED_CHARGE_MODELS = ["AM1-BCC"]  # TODO: Which models do we want to support?


# =============================================================================================
# Exceptions
# =============================================================================================


class MissingPackageError(MessageException):
    """This function requires a package that is not installed."""


class ToolkitUnavailableException(MessageException):
    """The requested toolkit is unavailable."""

    # TODO: Allow toolkit to be specified and used in formatting/printing exception.


class LicenseError(ToolkitUnavailableException):
    """This function requires a license that cannot be found."""


class InvalidToolkitError(MessageException):
    """A non-toolkit object was received when a toolkit object was expected"""


class InvalidToolkitRegistryError(MessageException):
    """An object other than a ToolkitRegistry or toolkit wrapper was received"""


class UndefinedStereochemistryError(MessageException):
    """A molecule was attempted to be loaded with undefined stereochemistry"""


class GAFFAtomTypeWarning(RuntimeWarning):
    """A warning raised if a loaded mol2 file possibly uses GAFF atom types."""


class ChargeMethodUnavailableError(MessageException):
    """A toolkit does not support the requested partial_charge_method combination"""


class IncorrectNumConformersError(MessageException):
    """The requested partial_charge_method expects a different number of conformers than was provided"""


class IncorrectNumConformersWarning(Warning):
    """The requested partial_charge_method expects a different number of conformers than was provided"""


class ChargeCalculationError(MessageException):
    """An unhandled error occured in an external toolkit during charge calculation"""


class InvalidIUPACNameError(MessageException):
    """Failed to parse IUPAC name"""


class AntechamberNotFoundError(MessageException):
    """The antechamber executable was not found"""


# =============================================================================================
# TOOLKIT UTILITY DECORATORS
# =============================================================================================

# =============================================================================================
# UTILITY FUNCTIONS
# =============================================================================================

# =============================================================================================
# CHEMINFORMATICS TOOLKIT WRAPPERS
# =============================================================================================


class ToolkitWrapper:
    """
    Toolkit wrapper base class.

    .. warning :: This API is experimental and subject to change.
    """

    _is_available = None  # True if toolkit is available
    _toolkit_version = None
    _toolkit_name = None  # Name of the toolkit
    _toolkit_installation_instructions = (
        None  # Installation instructions for the toolkit
    )

    # @staticmethod
    # TODO: Right now, to access the class definition, I have to make this a classmethod
    # and thereby call it with () on the outermost decorator. Is this wasting time? Are we caching
    # the is_available results?
    @classmethod
    def requires_toolkit(cls):  # remember cls is a ToolkitWrapper subclass here
        def decorator(func):
            @wraps(func)
            def wrapped_function(*args, **kwargs):
                if not cls.is_available():
                    msg = "This function requires the {} toolkit".format(
                        cls._toolkit_name
                    )
                    raise ToolkitUnavailableException(msg)
                value = func(*args, **kwargs)
                return value

            return wrapped_function

        return decorator

    @property
    # @classmethod
    def toolkit_name(self):
        """
        Return the name of the toolkit wrapped by this class as a str

        .. warning :: This API is experimental and subject to change.

        Returns
        -------
        toolkit_name : str
            The name of the wrapped toolkit

        """
        return self.__class__._toolkit_name

    @property
    # @classmethod
    def toolkit_installation_instructions(self):
        """
        Instructions on how to install the wrapped toolkit.
        """
        return self._toolkit_installation_instructions

    # @classmethod
    @property
    def toolkit_file_read_formats(self):
        """
        List of file formats that this toolkit can read.
        """
        return self._toolkit_file_read_formats

    # @classmethod
    @property
    def toolkit_file_write_formats(self):
        """
        List of file formats that this toolkit can write.
        """
        return self._toolkit_file_write_formats

    @classmethod
    def is_available(cls):
        """
        Check whether the corresponding toolkit can be imported

        Returns
        -------
        is_installed : bool
            True if corresponding toolkit is installed, False otherwise.

        """
        return NotImplementedError

    @property
    def toolkit_version(self):
        """
        Return the version of the wrapped toolkit as a str

        .. warning :: This API is experimental and subject to change.

        Returns
        -------
        toolkit_version : str
            The version of the wrapped toolkit

        """
        return self._toolkit_version

    def from_file(self, file_path, file_format, allow_undefined_stereo=False):
        """
        Return an openff.toolkit.topology.Molecule from a file using this toolkit.

        Parameters
        ----------
        file_path : str
            The file to read the molecule from
        file_format : str
            Format specifier, usually file suffix (eg. 'MOL2', 'SMI')
            Note that not all toolkits support all formats. Check ToolkitWrapper.toolkit_file_read_formats for details.
        allow_undefined_stereo : bool, default=False
            If false, raises an exception if any molecules contain undefined stereochemistry.
        _cls : class
            Molecule constructor
        Returns
        -------
        molecules : Molecule or list of Molecules
            a list of Molecule objects is returned.

        """
        return NotImplementedError

    def from_file_obj(
        self, file_obj, file_format, allow_undefined_stereo=False, _cls=None
    ):
        """
        Return an openff.toolkit.topology.Molecule from a file-like object (an object with a ".read()" method using this
         toolkit.

        Parameters
        ----------
        file_obj : file-like object
            The file-like object to read the molecule from
        file_format : str
            Format specifier, usually file suffix (eg. 'MOL2', 'SMI')
            Note that not all toolkits support all formats. Check ToolkitWrapper.toolkit_file_read_formats for details.
        allow_undefined_stereo : bool, default=False
            If false, raises an exception if any molecules contain undefined stereochemistry. If false, the function
            skips loading the molecule.
        _cls : class
            Molecule constructor

        Returns
        -------
        molecules : Molecule or list of Molecules
            a list of Molecule objects is returned.
        """
        return NotImplementedError

    @staticmethod
    def _check_n_conformers(
        molecule,
        partial_charge_method,
        min_confs=None,
        max_confs=None,
        strict_n_conformers=False,
    ):
        """
        Private method for validating the number of conformers on a molecule prior to partial
        charge calculation

        Parameters
        ----------
        molecule : Molecule
            Molecule for which partial charges are to be computed
        partial_charge_method : str, optional, default=None
            The name of the charge method being used
        min_confs : int, optional, default=None
            The minimum number of conformers required to use this charge method
        max_confs : int, optional, default=None
            The maximum number of conformers required to use this charge method
        strict_n_conformers : bool, default=False
            Whether to raise an exception if an invalid number of conformers is provided.
            If this is False and an invalid number of conformers is found, a warning will be raised.

        Raises
        ------
        IncorrectNumConformersError
            If the wrong number of conformers is attached to the input molecule, and strict_n_conformers is True.
        """
        import warnings

        n_confs = molecule.n_conformers
        wrong_confs_msg = (
            f"Molecule '{molecule}' has {n_confs} conformers, "
            f"but charge method '{partial_charge_method}' expects"
        )
        exception_suffix = (
            "You can disable this error by setting `strict_n_conformers=False' "
            "when calling 'molecule.assign_partial_charges'."
        )
        # If there's no n_confs filter, then this molecule automatically passes
        if min_confs is None and max_confs is None:
            return
        # If there's constraints on both ends, check both limits
        elif min_confs is not None and max_confs is not None:
            if not (min_confs <= n_confs <= max_confs):
                if min_confs == max_confs:
                    wrong_confs_msg += f" exactly {min_confs}."
                else:
                    wrong_confs_msg += f" between {min_confs} and {max_confs}."

            else:
                return
        # If there's only a max constraint, check that
        elif min_confs is not None and max_confs is None:
            if not (min_confs <= n_confs):
                wrong_confs_msg += f" at least {min_confs}."
            else:
                return
        # If there's only a maximum constraint, check that
        elif min_confs is None and max_confs is not None:
            if not (n_confs <= max_confs):
                wrong_confs_msg += f" at most {max_confs}."
            else:
                return
        # If we've made it this far, the molecule has the wrong number of conformers
        if strict_n_conformers:
            wrong_confs_msg += exception_suffix
            raise IncorrectNumConformersError(wrong_confs_msg)
        else:
            warnings.warn(wrong_confs_msg, IncorrectNumConformersWarning)

    def __repr__(self):
        return (
            f"ToolkitWrapper around {self.toolkit_name} version {self.toolkit_version}"
        )


@inherit_docstrings
class BuiltInToolkitWrapper(ToolkitWrapper):
    """
    Built-in ToolkitWrapper for very basic functionality. This is intended for use in testing and not much more.

    .. warning :: This API is experimental and subject to change.
    """

    _toolkit_name = "Built-in Toolkit"
    _toolkit_installation_instructions = (
        "This toolkit is installed with the Open Force Field Toolkit and does "
        "not require additional dependencies."
    )

    def __init__(self):
        super().__init__()

        self._toolkit_file_read_formats = []
        self._toolkit_file_write_formats = []

    def assign_partial_charges(
        self,
        molecule,
        partial_charge_method=None,
        use_conformers=None,
        strict_n_conformers=False,
        _cls=None,
    ):
        """
        Compute partial charges with the built-in toolkit using simple arithmetic operations, and assign
        the new values to the partial_charges attribute.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        molecule : openff.toolkit.topology.Molecule
            Molecule for which partial charges are to be computed
        partial_charge_method: str, optional, default=None
            The charge model to use. One of ['zeros', 'formal_charge']. If None, 'formal_charge' will be used.
        use_conformers : iterable of simtk.unit.Quantity-wrapped numpy arrays, each with shape (n_atoms, 3) and dimension of distance. Optional, default = None
            Coordinates to use for partial charge calculation. If None, an appropriate number of conformers
            will be generated.
        strict_n_conformers : bool, default=False
            Whether to raise an exception if an invalid number of conformers is provided for the given charge method.
            If this is False and an invalid number of conformers is found, a warning will be raised
            instead of an Exception.
        _cls : class
            Molecule constructor

        Raises
        ------
        ChargeMethodUnavailableError if the requested charge method can not be handled by this toolkit

        IncorrectNumConformersError if strict_n_conformers is True and use_conformers is provided and specifies an
        invalid number of conformers for the requested method

        ChargeCalculationError if the charge calculation is supported by this toolkit, but fails
        """

        PARTIAL_CHARGE_METHODS = {
            "zeros": {"rec_confs": 0, "min_confs": 0, "max_confs": 0},
            "formal_charge": {"rec_confs": 0, "min_confs": 0, "max_confs": 0},
        }

        if partial_charge_method is None:
            partial_charge_method = "formal_charge"

        if _cls is None:
            from openff.toolkit.topology.molecule import Molecule

            _cls = Molecule

        # Make a temporary copy of the molecule, since we'll be messing with its conformers
        mol_copy = _cls(molecule)

        partial_charge_method = partial_charge_method.lower()
        if partial_charge_method not in PARTIAL_CHARGE_METHODS:
            raise ChargeMethodUnavailableError(
                f'Partial charge method "{partial_charge_method}"" is not supported by '
                f"the Built-in toolkit. Available charge methods are "
                f"{list(PARTIAL_CHARGE_METHODS.keys())}"
            )

        if use_conformers is None:
            # Note that this refers back to the GLOBAL_TOOLKIT_REGISTRY by default, since
            # BuiltInToolkitWrapper can't generate conformers
            mol_copy.generate_conformers(
                n_conformers=PARTIAL_CHARGE_METHODS[partial_charge_method]["rec_confs"]
            )
        else:
            mol_copy._conformers = None
            for conformer in use_conformers:
                mol_copy._add_conformer(conformer)
            self._check_n_conformers(
                mol_copy,
                partial_charge_method=partial_charge_method,
                min_confs=0,
                max_confs=0,
                strict_n_conformers=strict_n_conformers,
            )

        partial_charges = unit.Quantity(
            np.zeros((molecule.n_particles)), unit.elementary_charge
        )
        if partial_charge_method == "zeroes":
            pass
        elif partial_charge_method == "formal_charge":
            for part_idx, particle in enumerate(molecule.particles):
                partial_charges[part_idx] = particle.formal_charge

        molecule.partial_charges = partial_charges


def publish(inherit_docstring=False, replacements={}):
    def publish_decorator(f):
        f._openff_info = {"published": True}
        if inherit_docstring:
            name = f.__name__
            base_func = getattr(ToolkitWrapper, name)
            docstring = inspect.getdoc(base_func)
            for from_s, to_s in replacements.items():
                docstring = docstring.replace(from_s, to_s)
            f.__doc__ = docstring
        return f
    return publish_decorator
    
class ModuleToolkitWrapper(ToolkitWrapper):
    def __init__(self):
        self._load_module()
        super().__init__()
        
    def _load_module(self):
        mod = importlib.import_module(self.module_name)
        if not mod.is_available():
            msg = (
                f"The required toolkit {mod.toolkit_name} is not "
                f"available. {mod.toolkit_installation_instructions}"
                )
            if not mod.is_installed():
                raise ToolkitUnavailableException(msg)
            if not mod.is_licensed():
                raise LicenseError(msg)
                
        self._copy_module_properties(mod)

    def _copy_module_properties(self, mod):
        #self._toolkit_name = mod.toolkit_name
        self._toolkit_version = mod.get_toolkit_version()
        self._toolkit_installation_instructions = mod.toolkit_installation_instructions
        self._toolkit_file_read_formats = mod.get_file_read_formats()
        self._toolkit_file_write_formats = mod.get_file_write_formats()
        
        for name, obj in mod.__dict__.items():
            if hasattr(obj, "_openff_info") and obj._openff_info["published"]:
                setattr(self, name, obj)

class OpenEyeToolkitWrapper(ModuleToolkitWrapper):
    """
    OpenEye toolkit wrapper

    .. warning :: This API is experimental and subject to change.
    """
    _toolkit_name = "OpenEye Toolkit"
    module_name = "openff.toolkit.utils.openeye_wrapper"
            
    

def requires_openeye_module(module_name):
    def inner_decorator(function):
        @wraps(function)
        def wrapper(*args, **kwargs):
            try:
                module = importlib.import_module("openeye." + module_name)
            except (ImportError, ModuleNotFoundError):
                # TODO: Custom exception
                raise Exception("openeye." + module_name)
            try:
                license_func = OpenEyeToolkitWrapper._license_functions[module_name]
            except KeyError:
                # TODO: Custom exception
                raise Exception(f"we do not currently use {module_name}")

            # TODO: Custom exception
            assert getattr(module, license_func)()

            return function(*args, **kwargs)

        return wrapper

    return inner_decorator



class RDKitToolkitWrapper(ModuleToolkitWrapper):
    """
    RDKit toolkit wrapper

    .. warning :: This API is experimental and subject to change.
    """
    _toolkit_name = "The RDKit"
    module_name = "openff.toolkit.utils.rdkit_wrapper"

class AmberToolsToolkitWrapper(ModuleToolkitWrapper):
    """
    AmberTools toolkit wrapper

    .. warning :: This API is experimental and subject to change.
    """
    _toolkit_name = "AmberTools"
    module_name = "openff.toolkit.utils.ambertools_wrapper"



# =============================================================================================
# Toolkit registry
# =============================================================================================


class ToolkitRegistry:
    """
    Registry for ToolkitWrapper objects

    Examples
    --------

    Register toolkits in a specified order, skipping if unavailable

    >>> from openff.toolkit.utils.toolkits import ToolkitRegistry
    >>> toolkit_precedence = [OpenEyeToolkitWrapper, RDKitToolkitWrapper, AmberToolsToolkitWrapper]
    >>> toolkit_registry = ToolkitRegistry(toolkit_precedence)
    >>> toolkit_registry
    ToolkitRegistry containing OpenEye Toolkit, The RDKit, AmberTools

    Register all available toolkits (in the order OpenEye, RDKit, AmberTools, built-in)

    >>> toolkits = [OpenEyeToolkitWrapper, RDKitToolkitWrapper, AmberToolsToolkitWrapper, BuiltInToolkitWrapper]
    >>> toolkit_registry = ToolkitRegistry(toolkit_precedence=toolkits)
    >>> toolkit_registry
    ToolkitRegistry containing OpenEye Toolkit, The RDKit, AmberTools, Built-in Toolkit

    Retrieve the global singleton toolkit registry, which is created when this module is imported from all available
    toolkits:

    >>> from openff.toolkit.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY as toolkit_registry
    >>> toolkit_registry
    ToolkitRegistry containing OpenEye Toolkit, The RDKit, AmberTools, Built-in Toolkit

    Note that this will contain different ToolkitWrapper objects based on what toolkits
    are currently installed.

    .. warning :: This API is experimental and subject to change.
    """

    def __init__(
        self,
        toolkit_precedence=[],
        exception_if_unavailable=True,
        _register_imported_toolkit_wrappers=False,
    ):
        """
        Create an empty toolkit registry.

        Parameters
        ----------
        toolkit_precedence : list, default=[]
            List of toolkit wrapper classes, in order of desired precedence when performing molecule operations. If
            None, no toolkits will be registered.

        exception_if_unavailable : bool, optional, default=True
            If True, an exception will be raised if the toolkit is unavailable

        _register_imported_toolkit_wrappers : bool, optional, default=False
            If True, will attempt to register all imported ToolkitWrapper subclasses that can be
            found in the order of toolkit_precedence, if specified. If toolkit_precedence is not
            specified, the default order is [OpenEyeToolkitWrapper, RDKitToolkitWrapper,
            AmberToolsToolkitWrapper, BuiltInToolkitWrapper].

        """
        self._toolkits = list()

        toolkits_to_register = list()

        if _register_imported_toolkit_wrappers:
            if toolkit_precedence is None:
                toolkit_precedence = [
                    OpenEyeToolkitWrapper,
                    RDKitToolkitWrapper,
                    AmberToolsToolkitWrapper,
                    BuiltInToolkitWrapper,
                ]
            all_importable_toolkit_wrappers = all_subclasses(ToolkitWrapper)
            for toolkit in toolkit_precedence:
                if toolkit in all_importable_toolkit_wrappers:
                    toolkits_to_register.append(toolkit)
        else:
            if toolkit_precedence:
                toolkits_to_register = toolkit_precedence

        if toolkits_to_register:
            for toolkit in toolkits_to_register:
                self.register_toolkit(
                    toolkit, exception_if_unavailable=exception_if_unavailable
                )

    @property
    def registered_toolkits(self):
        """
        List registered toolkits.

        .. warning :: This API is experimental and subject to change.

        .. todo :: Should this return a generator? Deep copies? Classes? Toolkit names?

        Returns
        -------
        toolkits : iterable of toolkit objects
        """
        return list(self._toolkits)

    @property
    def registered_toolkit_versions(self):
        """
        Return a dict containing the version of each registered toolkit.

        .. warning :: This API is experimental and subject to change.

        Returns
        -------
        toolkit_versions : dict[str, str]
            A dictionary mapping names and versions of wrapped toolkits

        """
        return dict(
            (tk.toolkit_name, tk.toolkit_version) for tk in self.registered_toolkits
        )

    def register_toolkit(self, toolkit_wrapper, exception_if_unavailable=True):
        """
        Register the provided toolkit wrapper class, instantiating an object of it.

        .. warning :: This API is experimental and subject to change.

        .. todo ::

           This method should raise an exception if the toolkit is unavailable, unless an optional argument
           is specified that silently avoids registration of toolkits that are unavailable.

        Parameters
        ----------
        toolkit_wrapper : instance or subclass of ToolkitWrapper
            The toolkit wrapper to register or its class.
        exception_if_unavailable : bool, optional, default=True
            If True, an exception will be raised if the toolkit is unavailable

        """
        # Instantiate class if class, or just add if already instantiated.
        if isinstance(toolkit_wrapper, type):
            try:
                toolkit_wrapper = toolkit_wrapper()
            except ToolkitUnavailableException:
                msg = "Unable to load toolkit '{}'. ".format(
                    toolkit_wrapper._toolkit_name
                )
                if exception_if_unavailable:
                    raise ToolkitUnavailableException(msg)
                else:
                    if "OpenEye" in msg:
                        msg += (
                            "The Open Force Field Toolkit does not require the OpenEye Toolkits, and can "
                            "use RDKit/AmberTools instead. However, if you have a valid license for the "
                            "OpenEye Toolkits, consider installing them for faster performance and additional "
                            "file format support: "
                            "https://docs.eyesopen.com/toolkits/python/quickstart-python/linuxosx.html "
                            "OpenEye offers free Toolkit licenses for academics: "
                            "https://www.eyesopen.com/academic-licensing"
                        )
                    logger.warning(f"Warning: {msg}")
                return

        # Add toolkit to the registry.
        self._toolkits.append(toolkit_wrapper)

    def deregister_toolkit(self, toolkit_wrapper):
        """
        Remove a ToolkitWrapper from the list of toolkits in this ToolkitRegistry

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        toolkit_wrapper : instance or subclass of ToolkitWrapper
            The toolkit wrapper to remove from the registry

        Raises
        ------
        InvalidToolkitError
            If toolkit_wrapper is not a ToolkitWrapper or subclass
        ToolkitUnavailableException
            If toolkit_wrapper is not found in the registry
        """
        # If passed a class, instantiate it
        if inspect.isclass(toolkit_wrapper):
            toolkit_wrapper = toolkit_wrapper()

        if not isinstance(toolkit_wrapper, ToolkitWrapper):
            msg = (
                f"Argument {toolkit_wrapper} must an ToolkitWrapper "
                f"or subclass of it. Found type {type(toolkit_wrapper)}."
            )
            raise InvalidToolkitError(msg)

        toolkits_to_remove = []

        for toolkit in self._toolkits:
            if type(toolkit) == type(toolkit_wrapper):
                toolkits_to_remove.append(toolkit)

        if not toolkits_to_remove:
            msg = (
                f"Did not find {toolkit_wrapper} in registry. "
                f"Currently registered toolkits are {self._toolkits}"
            )
            raise ToolkitUnavailableException(msg)

        for toolkit_to_remove in toolkits_to_remove:
            self._toolkits.remove(toolkit_to_remove)

    def add_toolkit(self, toolkit_wrapper):
        """
        Append a ToolkitWrapper onto the list of toolkits in this ToolkitRegistry

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        toolkit_wrapper : openff.toolkit.utils.ToolkitWrapper
            The ToolkitWrapper object to add to the list of registered toolkits

        Raises
        ------
        InvalidToolkitError
            If toolkit_wrapper is not a ToolkitWrapper or subclass
        """
        if not isinstance(toolkit_wrapper, ToolkitWrapper):
            msg = "Something other than a ToolkitWrapper object was passed to ToolkitRegistry.add_toolkit()\n"
            msg += "Given object {} of type {}".format(
                toolkit_wrapper, type(toolkit_wrapper)
            )
            raise InvalidToolkitError(msg)
        self._toolkits.append(toolkit_wrapper)

    # TODO: Can we automatically resolve calls to methods that are not explicitly defined using some Python magic?

    def resolve(self, method_name):
        """
        Resolve the requested method name by checking all registered toolkits in
        order of precedence for one that provides the requested method.

        Parameters
        ----------
        method_name : str
            The name of the method to resolve

        Returns
        -------
        method
            The method of the first registered toolkit that provides the requested method name

        Raises
        ------
        NotImplementedError if the requested method cannot be found among the registered toolkits

        Examples
        --------

        Create a molecule, and call the toolkit ``to_smiles()`` method directly

        >>> from openff.toolkit.topology import Molecule
        >>> molecule = Molecule.from_smiles('Cc1ccccc1')
        >>> toolkit_registry = ToolkitRegistry([OpenEyeToolkitWrapper, RDKitToolkitWrapper, AmberToolsToolkitWrapper])
        >>> method = toolkit_registry.resolve('to_smiles')
        >>> smiles = method(molecule)

        .. todo :: Is there a better way to figure out which toolkits implement given methods by introspection?

        """
        for toolkit in self._toolkits:
            if hasattr(toolkit, method_name):
                method = getattr(toolkit, method_name)
                return method

        # No toolkit was found to provide the requested capability
        # TODO: Can we help developers by providing a check for typos in expected method names?
        msg = 'No registered toolkits can provide the capability "{}".\n'.format(
            method_name
        )
        msg += "Available toolkits are: {}\n".format(self.registered_toolkits)
        raise NotImplementedError(msg)

    # TODO: Can we instead register available methods directly with `ToolkitRegistry`, so we can just use `ToolkitRegistry.method()`?
    def call(self, method_name, *args, raise_exception_types=None, **kwargs):
        """
        Execute the requested method by attempting to use all registered toolkits in order of precedence.

        ``*args`` and ``**kwargs`` are passed to the desired method, and return values of the method are returned

        This is a convenient shorthand for ``toolkit_registry.resolve_method(method_name)(*args, **kwargs)``

        Parameters
        ----------
        method_name : str
            The name of the method to execute
        raise_exception_types : list of Exception subclasses, default=None
            A list of exception-derived types to catch and raise immediately. If None, this will be set to [Exception],
            which will raise an error immediately if the first ToolkitWrapper in the registry fails. To try each
            ToolkitWrapper that provides a suitably-named method, set this to the empty list ([]). If all
            ToolkitWrappers run without raising any exceptions in this list, a single ValueError will be raised
            containing the each ToolkitWrapper that was tried and the exception it raised.

        Raises
        ------
        NotImplementedError if the requested method cannot be found among the registered toolkits

        ValueError if no exceptions in the raise_exception_types list were raised by ToolkitWrappers, and
        all ToolkitWrappers in the ToolkitRegistry were tried.

        Other forms of exceptions are possible if raise_exception_types is specified.
        These are defined by the ToolkitWrapper method being called.

        Examples
        --------

        Create a molecule, and call the toolkit ``to_smiles()`` method directly

        >>> from openff.toolkit.topology import Molecule
        >>> molecule = Molecule.from_smiles('Cc1ccccc1')
        >>> toolkit_registry = ToolkitRegistry([OpenEyeToolkitWrapper, RDKitToolkitWrapper])
        >>> smiles = toolkit_registry.call('to_smiles', molecule)

        """
        if raise_exception_types is None:
            raise_exception_types = [Exception]

        errors = list()
        for toolkit in self._toolkits:
            if hasattr(toolkit, method_name):
                method = getattr(toolkit, method_name)
                try:
                    return method(*args, **kwargs)
                except Exception as e:
                    for exception_type in raise_exception_types:
                        if isinstance(e, exception_type):
                            raise e
                    errors.append((toolkit, e))

        # No toolkit was found to provide the requested capability
        # TODO: Can we help developers by providing a check for typos in expected method names?
        msg = (
            f'No registered toolkits can provide the capability "{method_name}" '
            f'for args "{args}" and kwargs "{kwargs}"\n'
        )

        msg += "Available toolkits are: {}\n".format(self.registered_toolkits)
        # Append information about toolkits that implemented the method, but could not handle the provided parameters
        for toolkit, error in errors:
            msg += " {} {} : {}\n".format(toolkit, type(error), error)
        raise ValueError(msg)

    def __repr__(self):
        return f"ToolkitRegistry containing " + ", ".join(
            [tk.toolkit_name for tk in self._toolkits]
        )


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
    msg = "WARNING: No basic cheminformatics toolkits are available.\n"
    msg += "At least one basic toolkit is required to handle SMARTS matching and file I/O. \n"
    msg += "Please install at least one of the following basic toolkits:\n"
    for wrapper in all_subclasses(ToolkitWrapper):
        if wrapper.toolkit_name is not None:
            msg += "{} : {}\n".format(
                wrapper._toolkit_name, wrapper._toolkit_installation_instructions
            )
    print(msg)
