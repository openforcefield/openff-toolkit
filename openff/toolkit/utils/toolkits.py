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
from distutils.spawn import find_executable
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
        self._toolkit_name = mod.toolkit_name
        self._toolkit_version = mod.toolkit_version
        self._toolkit_installation_instructions = mod.toolkit_installation_instructions
        self._toolkit_file_read_formats = mod.toolkit_file_read_formats
        self._toolkit_file_write_formats = mod.toolkit_file_write_formats
        
        for name, obj in mod.__dict__.items():
            if hasattr(obj, "_openff_info") and obj._openff_info["published"]:
                assert not name.startswith("_")
                setattr(self, name, obj)

class OpenEyeToolkitWrapper(ModuleToolkitWrapper):
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
    module_name = "openff.toolkit.utils.rdkit_wrapper"


class AmberToolsToolkitWrapper(ToolkitWrapper):
    """
    AmberTools toolkit wrapper

    .. warning :: This API is experimental and subject to change.
    """

    _toolkit_name = "AmberTools"
    _toolkit_installation_instructions = (
        "The AmberTools toolkit (free and open source) can be found at "
        "https://anaconda.org/conda-forge/ambertools"
    )

    def __init__(self):
        super().__init__()

        self._toolkit_file_read_formats = []
        self._toolkit_file_write_formats = []

        if not self.is_available():
            raise ToolkitUnavailableException(
                f"The required toolkit {self._toolkit_name} is not "
                f"available. {self._toolkit_installation_instructions}"
            )

        # TODO: More reliable way to extract AmberTools version
        out = subprocess.check_output(["antechamber", "-L"])
        ambertools_version = out.decode("utf-8").split("\n")[1].split()[3].strip(":")
        self._toolkit_version = ambertools_version

        # TODO: Find AMBERHOME or executable home, checking miniconda if needed
        # Store an instance of an RDKitToolkitWrapper for file I/O
        self._rdkit_toolkit_wrapper = RDKitToolkitWrapper()

    @staticmethod
    def is_available():
        """
        Check whether the AmberTools toolkit is installed

        Returns
        -------
        is_installed : bool
            True if AmberTools is installed, False otherwise.

        """
        # TODO: Check all tools needed
        # TODO: How should we implement find_executable?
        ANTECHAMBER_PATH = find_executable("antechamber")
        if ANTECHAMBER_PATH is None:
            return False
        # AmberToolsToolkitWrapper needs RDKit to do basically anything, since its interface requires SDF I/O
        if not (RDKitToolkitWrapper.is_available()):
            return False
        return True

    def assign_partial_charges(
        self,
        molecule,
        partial_charge_method=None,
        use_conformers=None,
        strict_n_conformers=False,
        _cls=None,
    ):
        """
        Compute partial charges with AmberTools using antechamber/sqm, and assign
        the new values to the partial_charges attribute.

        .. warning :: This API experimental and subject to change.

        .. todo ::

           * Do we want to also allow ESP/RESP charges?

        Parameters
        ----------
        molecule : openff.toolkit.topology.Molecule
            Molecule for which partial charges are to be computed
        partial_charge_method : str, optional, default=None
            The charge model to use. One of ['gasteiger', 'am1bcc', 'am1-mulliken']. If None, 'am1-mulliken' will be used.
        use_conformers : iterable of simtk.unit.Quantity-wrapped numpy arrays, each with shape (n_atoms, 3) and dimension of distance. Optional, default = None
            List of (n_atoms x 3) simtk.unit.Quantities to use for partial charge calculation.
            If None, an appropriate number of conformers will be generated.
        strict_n_conformers : bool, default=False
            Whether to raise an exception if an invalid number of conformers is provided for the given charge method.
            If this is False and an invalid number of conformers is found, a warning will be raised.
        _cls : class
            Molecule constructor

        Raises
        ------
        ChargeMethodUnavailableError if the requested charge method can not be handled by this toolkit

        ChargeCalculationError if the charge method is supported by this toolkit, but fails
        """

        import os
        import subprocess

        from openff.toolkit.topology import Molecule

        if partial_charge_method is None:
            partial_charge_method = "am1-mulliken"
        else:
            # Standardize method name for string comparisons
            partial_charge_method = partial_charge_method.lower()

        SUPPORTED_CHARGE_METHODS = {
            "am1bcc": {
                "antechamber_keyword": "bcc",
                "min_confs": 1,
                "max_confs": 1,
                "rec_confs": 1,
            },
            "am1-mulliken": {
                "antechamber_keyword": "mul",
                "min_confs": 1,
                "max_confs": 1,
                "rec_confs": 1,
            },
            "gasteiger": {
                "antechamber_keyword": "gas",
                "min_confs": 0,
                "max_confs": 0,
                "rec_confs": 0,
            },
        }

        if partial_charge_method not in SUPPORTED_CHARGE_METHODS:
            raise ChargeMethodUnavailableError(
                f"partial_charge_method '{partial_charge_method}' is not available from AmberToolsToolkitWrapper. "
                f"Available charge methods are {list(SUPPORTED_CHARGE_METHODS.keys())} "
            )

        charge_method = SUPPORTED_CHARGE_METHODS[partial_charge_method]

        if _cls is None:
            from openff.toolkit.topology.molecule import Molecule

            _cls = Molecule

        # Make a temporary copy of the molecule, since we'll be messing with its conformers
        mol_copy = _cls(molecule)

        if use_conformers is None:
            if charge_method["rec_confs"] == 0:
                mol_copy._conformers = None
            else:
                mol_copy.generate_conformers(
                    n_conformers=charge_method["rec_confs"],
                    rms_cutoff=0.25 * unit.angstrom,
                    toolkit_registry=RDKitToolkitWrapper(),
                )
            # TODO: What's a "best practice" RMS cutoff to use here?
        else:
            mol_copy._conformers = None
            for conformer in use_conformers:
                mol_copy._add_conformer(conformer)
            self._check_n_conformers(
                mol_copy,
                partial_charge_method=partial_charge_method,
                min_confs=charge_method["min_confs"],
                max_confs=charge_method["max_confs"],
                strict_n_conformers=strict_n_conformers,
            )

        # Find the path to antechamber
        # TODO: How should we implement find_executable?
        ANTECHAMBER_PATH = find_executable("antechamber")
        if ANTECHAMBER_PATH is None:
            raise AntechamberNotFoundError(
                "Antechamber not found, cannot run charge_mol()"
            )

        # Compute charges
        with tempfile.TemporaryDirectory() as tmpdir:
            with temporary_cd(tmpdir):
                net_charge = mol_copy.total_charge / unit.elementary_charge
                # Write out molecule in SDF format
                ## TODO: How should we handle multiple conformers?
                self._rdkit_toolkit_wrapper.to_file(
                    mol_copy, "molecule.sdf", file_format="sdf"
                )
                # Compute desired charges
                # TODO: Add error handling if antechamber chokes
                short_charge_method = charge_method["antechamber_keyword"]
                subprocess.check_output(
                    [
                        "antechamber",
                        "-i",
                        "molecule.sdf",
                        "-fi",
                        "sdf",
                        "-o",
                        "charged.mol2",
                        "-fo",
                        "mol2",
                        "-pf",
                        "yes",
                        "-dr",
                        "n",
                        "-c",
                        short_charge_method,
                        "-nc",
                        str(net_charge),
                    ]
                )
                # Write out just charges
                subprocess.check_output(
                    [
                        "antechamber",
                        "-dr",
                        "n",
                        "-i",
                        "charged.mol2",
                        "-fi",
                        "mol2",
                        "-o",
                        "charges2.mol2",
                        "-fo",
                        "mol2",
                        "-c",
                        "wc",
                        "-cf",
                        "charges.txt",
                        "-pf",
                        "yes",
                    ]
                )
                # Check to ensure charges were actually produced
                if not os.path.exists("charges.txt"):
                    # TODO: copy files into local directory to aid debugging?
                    raise ChargeCalculationError(
                        "Antechamber/sqm partial charge calculation failed on "
                        "molecule {} (SMILES {})".format(
                            molecule.name, molecule.to_smiles()
                        )
                    )
                # Read the charges
                with open("charges.txt", "r") as infile:
                    contents = infile.read()
                text_charges = contents.split()
                charges = np.zeros([molecule.n_atoms], np.float64)
                for index, token in enumerate(text_charges):
                    charges[index] = float(token)
                # TODO: Ensure that the atoms in charged.mol2 are in the same order as in molecule.sdf
        charges = unit.Quantity(charges, unit.elementary_charge)
        molecule.partial_charges = charges

    def compute_partial_charges_am1bcc(
        self, molecule, use_conformers=None, strict_n_conformers=False
    ):
        """
        Compute partial charges with AmberTools using antechamber/sqm. This will calculate AM1-BCC charges on the first
        conformer only.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        molecule : Molecule
            Molecule for which partial charges are to be computed
        use_conformers : iterable of simtk.unit.Quantity-wrapped numpy arrays, each with shape (n_atoms, 3) and dimension of distance. Optional, default = None
            Coordinates to use for partial charge calculation. If None, an appropriate number of conformers
            will be generated.
        strict_n_conformers : bool, default=False
            Whether to raise an exception if an invalid number of conformers is provided.
            If this is False and an invalid number of conformers is found, a warning will be raised
            instead of an Exception.

        Returns
        -------
        charges : numpy.array of shape (natoms) of type float
            The partial charges
        """

        import warnings

        warnings.warn(
            "compute_partial_charges_am1bcc will be deprecated in an upcoming release. "
            "Use assign_partial_charges(partial_charge_method='am1bcc') instead.",
            DeprecationWarning,
        )

        self.assign_partial_charges(
            molecule,
            partial_charge_method="AM1BCC",
            use_conformers=use_conformers,
            strict_n_conformers=strict_n_conformers,
        )
        return molecule.partial_charges

    def _modify_sqm_in_to_request_bond_orders(self, file_path):
        """
        Modify a sqm.in file produced by antechamber to include the "printbondorders=1" directive
        in the header. This method will overwrite the original file.

        Parameters
        ----------
        file_path : str
            The path to sqm.in
        """

        data = open(file_path).read()

        # Original sqm.in file headerlooks like:

        # Run semi-empirical minimization
        #  &qmmm
        #    qm_theory='AM1', grms_tol=0.0005,
        #    scfconv=1.d-10, ndiis_attempts=700,   qmcharge=0,
        #  /
        # ... (atom coordinates in something like XYZ format) ...

        # To get WBOs, we need to add "printbondorders=1" to the list of keywords

        # First, split the sqm.in text at the "/" mark at the end of the header
        datasp = data.split("/")
        # Insert the "printbondorders" directive in a new line and re-add the "/"
        datasp.insert(1, "printbondorders=1, \n /")
        # Reassemble the file text
        new_data = "".join(datasp)
        # Write the new file contents, overwriting the original file.
        with open(file_path, "w") as of:
            of.write(new_data)

    def _get_fractional_bond_orders_from_sqm_out(
        self, file_path, validate_elements=None
    ):
        """
        Process a SQM output file containing bond orders, and return a dict of the form
        dict[atom_1_index, atom_2_index] = fractional_bond_order

        Parameters
        ----------
        file_path : str
            File path for sqm output file
        validate_elements : iterable of str
            The element symbols expected in molecule index order. A ValueError will be raised
            if the elements are not found in this order.

        Returns
        -------
        bond_orders : dict[(int, int)]: float
            A dictionary where the keys are tuples of two atom indices and the values are
            floating-point bond orders. The keys are sorted in ascending order, such that
            the lower atom index is key[0] and the higher is key[1].
        """

        # Example sqm.out section with WBOs:
        #  Bond Orders
        #
        #   QMMM:    NUM1 ELEM1 NUM2 ELEM2      BOND_ORDER
        #   QMMM:       2   C      1   C        1.41107532
        #   QMMM:       3   C      1   C        1.41047804
        # ...
        #   QMMM:      15   H     13   H        0.00000954
        #   QMMM:      15   H     14   H        0.00000813
        #
        #            --------- Calculation Completed ----------

        data = open(file_path).read()

        begin_sep = """ Bond Orders
 
  QMMM:    NUM1 ELEM1 NUM2 ELEM2      BOND_ORDER
"""
        end_sep = """

           --------- Calculation Completed ----------
"""
        # Extract the chunk of text between begin_sep and end_sep, and split it by newline
        fbo_lines = data.split(begin_sep)[1].split(end_sep)[0].split("\n")

        # Iterate over the lines and populate the dict to return
        bond_orders = dict()
        for line in fbo_lines:
            linesp = line.split()
            atom_index_1 = int(linesp[1])
            atom_element_1 = linesp[2]
            atom_index_2 = int(linesp[3])
            atom_element_2 = linesp[4]
            bond_order = float(linesp[5])

            # If validate_elements was provided, ensure that the ordering of element symbols is what we expected
            if validate_elements is not None:
                if (atom_element_1 != validate_elements[atom_index_1 - 1]) or (
                    atom_element_2 != validate_elements[atom_index_2 - 1]
                ):
                    # raise ValueError('\n'.join(fbo_lines))
                    raise ValueError(
                        f"Elements or indexing in sqm output differ from expectation. "
                        f"Expected {validate_elements[atom_index_1]} with index {atom_index_1} and "
                        f"{validate_elements[atom_index_2]} with index {atom_index_2}, "
                        f"but SQM output has {atom_element_1} and {atom_element_2} for the same atoms."
                    )

            # To make lookup easier, we identify bonds as integer tuples with the lowest atom index
            # first and the highest second.
            index_tuple = tuple(sorted([atom_index_1, atom_index_2]))
            bond_orders[index_tuple] = bond_order
        return bond_orders

    def assign_fractional_bond_orders(
        self, molecule, bond_order_model=None, use_conformers=None, _cls=None
    ):
        """
        Update and store list of bond orders this molecule. Bond orders are stored on each
        bond, in the `bond.fractional_bond_order` attribute.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        molecule : openff.toolkit.topology.molecule Molecule
            The molecule to assign wiberg bond orders to
        bond_order_model : str, optional, default=None
            The charge model to use. Only allowed value is 'am1-wiberg'. If None, 'am1-wiberg' will be used.
        use_conformers : iterable of simtk.unit.Quantity(np.array) with shape (n_atoms, 3) and dimension of distance, optional, default=None
            The conformers to use for fractional bond order calculation. If None, an appropriate number
            of conformers will be generated by an available ToolkitWrapper.
        _cls : class
            Molecule constructor
        """
        from openff.toolkit.topology import Molecule

        # Find the path to antechamber
        # TODO: How should we implement find_executable?
        ANTECHAMBER_PATH = find_executable("antechamber")
        if ANTECHAMBER_PATH is None:
            raise AntechamberNotFoundError(
                "Antechamber not found, cannot run "
                "AmberToolsToolkitWrapper.assign_fractional_bond_orders()"
            )

        if _cls is None:
            from openff.toolkit.topology.molecule import Molecule

            _cls = Molecule

        # Make a copy since we'll be messing with this molecule's conformers
        temp_mol = _cls(molecule)

        if use_conformers is None:
            temp_mol.generate_conformers(
                n_conformers=1,
                toolkit_registry=self._rdkit_toolkit_wrapper,
            )
        else:
            temp_mol._conformers = None
            for conformer in use_conformers:
                temp_mol._add_conformer(conformer)

        if len(temp_mol.conformers) == 0:
            raise ValueError(
                "No conformers present in molecule submitted for fractional bond order calculation. Consider "
                "loading the molecule from a file with geometry already present or running "
                "molecule.generate_conformers() before calling molecule.assign_fractional_bond_orders"
            )

        # Compute bond orders
        bond_order_model_to_antechamber_keyword = {"am1-wiberg": "mul"}
        supported_bond_order_models = list(
            bond_order_model_to_antechamber_keyword.keys()
        )
        if bond_order_model is None:
            bond_order_model = "am1-wiberg"

        bond_order_model = bond_order_model.lower()

        if bond_order_model not in supported_bond_order_models:
            raise ValueError(
                f"Bond order model '{bond_order_model}' is not supported by AmberToolsToolkitWrapper. "
                f"Supported models are {supported_bond_order_models}"
            )
        ac_charge_keyword = bond_order_model_to_antechamber_keyword[bond_order_model]

        bond_orders = defaultdict(list)

        for conformer in [*temp_mol.conformers]:

            with tempfile.TemporaryDirectory() as tmpdir:

                with temporary_cd(tmpdir):
                    net_charge = temp_mol.total_charge
                    # Write out molecule in SDF format
                    temp_mol._conformers = [conformer]
                    self._rdkit_toolkit_wrapper.to_file(
                        temp_mol, "molecule.sdf", file_format="sdf"
                    )
                    # Prepare sqm.in file as if we were going to run charge calc
                    # TODO: Add error handling if antechamber chokes
                    subprocess.check_output(
                        [
                            "antechamber",
                            "-i",
                            "molecule.sdf",
                            "-fi",
                            "sdf",
                            "-o",
                            "sqm.in",
                            "-fo",
                            "sqmcrt",
                            "-pf",
                            "yes",
                            "-c",
                            ac_charge_keyword,
                            "-nc",
                            str(net_charge),
                        ]
                    )
                    # Modify sqm.in to request bond order calculation
                    self._modify_sqm_in_to_request_bond_orders("sqm.in")
                    # Run sqm to get bond orders
                    subprocess.check_output(
                        ["sqm", "-i", "sqm.in", "-o", "sqm.out", "-O"]
                    )
                    # Ensure that antechamber/sqm did not change the indexing by checking against
                    # an ordered list of element symbols for this molecule
                    expected_elements = [at.element.symbol for at in molecule.atoms]
                    conformer_bond_orders = (
                        self._get_fractional_bond_orders_from_sqm_out(
                            "sqm.out", validate_elements=expected_elements
                        )
                    )

                    for bond_indices, value in conformer_bond_orders.items():
                        bond_orders[bond_indices].append(value)

        # Note that sqm calculate WBOs for ALL PAIRS of atoms, not just those that have
        # bonds defined in the original molecule. So here we iterate over the bonds in
        # the original molecule and only nab the WBOs for those.
        for bond in molecule.bonds:
            # The atom index tuples that act as bond indices are ordered from lowest to highest by
            # _get_fractional_bond_orders_from_sqm_out, so here we make sure that we look them up in
            # sorted order as well
            sorted_atom_indices = sorted(
                tuple([bond.atom1_index + 1, bond.atom2_index + 1])
            )
            bond.fractional_bond_order = np.mean(
                bond_orders[tuple(sorted_atom_indices)]
            )


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
