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
   * Should this be under ``openforcefield.utils.toolkits`` or ``openforcefield.toolkits``?
   * Add singleton global toolkit registry that registers all available toolkits by default when this file is imported
   * Add description fields for each toolkit wrapper
   * Eliminate global variables in favor of a singleton pattern
   * Change global variables from _INSTALLED to _AVAILABLE

"""

__all__ = [
    'DEFAULT_AROMATICITY_MODEL',
    'ALLOWED_AROMATICITY_MODELS',
    'DEFAULT_FRACTIONAL_BOND_ORDER_MODEL',
    'ALLOWED_FRACTIONAL_BOND_ORDER_MODELS',
    'DEFAULT_CHARGE_MODEL',
    'ALLOWED_CHARGE_MODELS',
    'LicenseError',
    'MissingPackageError',
    'ToolkitUnavailableException',
    'InvalidToolkitError',
    'UndefinedStereochemistryError',
    'GAFFAtomTypeWarning',
    'ToolkitWrapper',
    'OpenEyeToolkitWrapper',
    'RDKitToolkitWrapper',
    'AmberToolsToolkitWrapper',
    'ToolkitRegistry',
    'GLOBAL_TOOLKIT_REGISTRY',
    'OPENEYE_AVAILABLE',
    'RDKIT_AVAILABLE',
    'AMBERTOOLS_AVAILABLE',
    'BASIC_CHEMINFORMATICS_TOOLKITS'
]


#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import copy
from distutils.spawn import find_executable
from functools import wraps
import importlib
import logging
import subprocess

from simtk import unit
import numpy as np

from openforcefield.utils import all_subclasses, MessageException, inherit_docstrings


#=============================================================================================
# CONFIGURE LOGGER
#=============================================================================================

logger = logging.getLogger(__name__)

#=============================================================================================
# SUPPORTED MODELS
#
# TODO: We may no longer need these since we now require SMIRNOFF to specify these models explicitly.
#=============================================================================================

DEFAULT_AROMATICITY_MODEL = 'OEAroModel_MDL'  # TODO: Is there a more specific name and reference for the aromaticity model?
ALLOWED_AROMATICITY_MODELS = ['OEAroModel_MDL']

DEFAULT_FRACTIONAL_BOND_ORDER_MODEL = 'Wiberg'  # TODO: Is there a more specific name and reference for the fractional bond order models?
ALLOWED_FRACTIONAL_BOND_ORDER_MODELS = ['Wiberg']

DEFAULT_CHARGE_MODEL = 'AM1-BCC'  # TODO: Should this be `AM1-BCC`, or should we encode BCCs explicitly via AM1-CM2 preprocessing?
ALLOWED_CHARGE_MODELS = ['AM1-BCC'
                         ]  # TODO: Which models do we want to support?


#=============================================================================================
# Exceptions
#=============================================================================================

class LicenseError(Exception):
    """This function requires a license that cannot be found."""
    pass


class MissingPackageError(MessageException):
    """This function requires a package that is not installed."""
    pass


class ToolkitUnavailableException(MessageException):
    """The requested toolkit is unavailable."""
    # TODO: Allow toolkit to be specified and used in formatting/printing exception.
    pass


class InvalidToolkitError(MessageException):
    """A non-toolkit object was received when a toolkit object was expected"""


class UndefinedStereochemistryError(MessageException):
    """A molecule was attempted to be loaded with undefined stereochemistry"""
    pass


class GAFFAtomTypeWarning(RuntimeWarning):
    """A warning raised if a loaded mol2 file possibly uses GAFF atom types."""
    pass

class ChargeMethodUnavailableError(MessageException):
    """A toolkit does not support the requested partial_charge_method combination"""
    pass

class IncorrectNumConformersError(MessageException):
    """The requested partial_charge_method expects a different number of conformers than was provided"""
    pass

class IncorrectNumConformersWarning(Warning):
    """The requested partial_charge_method expects a different number of conformers than was provided"""
    pass

class ChargeCalculationError(MessageException):
    """An unhandled error occured in an external toolkit during charge calculation"""
    pass

#=============================================================================================
# TOOLKIT UTILITY DECORATORS
#=============================================================================================

#=============================================================================================
# UTILITY FUNCTIONS
#=============================================================================================

#=============================================================================================
# CHEMINFORMATICS TOOLKIT WRAPPERS
#=============================================================================================

class ToolkitWrapper:
    """
    Toolkit wrapper base class.

    .. warning :: This API is experimental and subject to change.
    """
    _is_available = None  # True if toolkit is available
    _toolkit_name = None  # Name of the toolkit
    _toolkit_installation_instructions = None  # Installation instructions for the toolkit
    _toolkit_file_read_formats = None  # The file types that this toolkit can read
    _toolkit_file_write_formats = None  # The file types that this toolkit can write

    #@staticmethod
    # TODO: Right now, to access the class definition, I have to make this a classmethod
    # and thereby call it with () on the outermost decorator. Is this wasting time? Are we caching
    # the is_available results?
    @classmethod
    def requires_toolkit(cls):  #remember cls is a ToolkitWrapper subclass here
        def decorator(func):
            @wraps(func)
            def wrapped_function(*args, **kwargs):
                if not cls.is_available():
                    msg = 'This function requires the {} toolkit'.format(
                        cls._toolkit_name)
                    raise LicenseError(msg)
                value = func(*args, **kwargs)
                return value

            return wrapped_function

        return decorator

    @property
    #@classmethod
    def toolkit_name(self):
        """
        The name of the toolkit wrapped by this class.
        """
        return self.__class__._toolkit_name

    @property
    @classmethod
    def toolkit_installation_instructions(cls):
        """
        Instructions on how to install the wrapped toolkit.
        """
        return cls._toolkit_installation_instructions

    #@classmethod
    @property
    def toolkit_file_read_formats(self):
        """
        List of file formats that this toolkit can read.
        """
        return self._toolkit_file_read_formats

    #@classmethod
    @property
    def toolkit_file_write_formats(self):
        """
        List of file formats that this toolkit can write.
        """
        return self._toolkit_file_write_formats

    @staticmethod
    def is_available():
        """
        Check whether the corresponding toolkit can be imported

        Returns
        -------
        is_installed : bool
            True if corresponding toolkit is installed, False otherwise.

        """
        return NotImplementedError

    def from_file(self,
                  file_path,
                  file_format,
                  allow_undefined_stereo=False):
        """
        Return an openforcefield.topology.Molecule from a file using this toolkit.

        Parameters
        ----------
        file_path : str
            The file to read the molecule from
        file_format : str
            Format specifier, usually file suffix (eg. 'MOL2', 'SMI')
            Note that not all toolkits support all formats. Check ToolkitWrapper.toolkit_file_read_formats for details.
        allow_undefined_stereo : bool, default=False
            If false, raises an exception if any molecules contain undefined stereochemistry.
        Returns
        -------
        molecules : Molecule or list of Molecules
            a list of Molecule objects is returned.

        """
        return NotImplementedError

    def from_file_obj(self,
                      file_obj,
                      file_format,
                      allow_undefined_stereo=False):
        """
        Return an openforcefield.topology.Molecule from a file-like object (an object with a ".read()" method using this
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

        Returns
        -------
        molecules : Molecule or list of Molecules
            a list of Molecule objects is returned.
        """
        return NotImplementedError


    def _check_n_conformers(self,
                            molecule,
                            partial_charge_method,
                            min_confs=None,
                            max_confs=None,
                            strict_n_conformers=False):
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
        ValueError
            If the wrong number of conformers is attached to the input molecule, and strict_n_conformers is True.
        """
        import warnings
        n_confs = len(molecule.conformers)
        wrong_confs_msg = f"Molecule '{molecule}' has {n_confs} conformers, " \
            f"but charge method '{partial_charge_method}' expects"
        exception_suffix = "You can disable this error by setting `strict_n_conformers=False' " \
                           "when calling 'molecule.assign_partial_charges'."
        # If there's no n_confs filter, then this molecule automatically passes
        if min_confs is None and max_confs is None:
            return
        # If there's constraints on both ends, check both limits
        elif min_confs is not None and max_confs is not None:
            if not(min_confs <= n_confs <= max_confs):
                if min_confs == max_confs:
                    wrong_confs_msg += f" exactly {min_confs}."
                else:
                    wrong_confs_msg += f" between {min_confs} and {max_confs}."

            else:
                return
        # If there's only a max constraint, check that
        elif min_confs is not None and max_confs is None:
            if not(min_confs <= n_confs):
                wrong_confs_msg += f" at least {min_confs}."
            else:
                return
        # If there's only a maximum constraint, check that
        elif min_confs is None and max_confs is not None:
            if not(n_confs <= max_confs):
                wrong_confs_msg += f" at most {max_confs}."
            else:
                return
        # If we've made it this far, the molecule has the wrong number of conformers
        if strict_n_conformers:
            wrong_confs_msg += exception_suffix
            raise IncorrectNumConformersError(wrong_confs_msg)
        else:
            warnings.warn(wrong_confs_msg, IncorrectNumConformersWarning)





@inherit_docstrings
class BuiltInToolkitWrapper(ToolkitWrapper):
    """
    Built-in ToolkitWrapper for very basic functionality. This is intended for use in testing and not much more.

    .. warning :: This API is experimental and subject to change.
    """
    _toolkit_name = 'Built-in Toolkit'
    _toolkit_installation_instructions = 'This toolkit is installed with the Open Force Field Toolkit and does' \
                                         'not require additional dependencies.'
    _toolkit_file_read_formats = []
    _toolkit_file_write_formats = []


    def assign_partial_charges(self,
                               molecule,
                               partial_charge_method=None,
                               use_conformers=None,
                               strict_n_conformers=False):
        """
        Compute partial charges with the built-in toolkit using simple arithmetic operations, and assign
        the new values to the partial_charges attribute.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        molecule : openforcefield.topology.Molecule
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

        Raises
        ------
        ChargeMethodUnavailableError if the requested charge method can not be handled by this toolkit

        IncorrectNumConformersError if strict_n_conformers is True and use_conformers is provided and specifies an
        invalid number of conformers for the requested method

        ChargeCalculationError if the charge calculation is supported by this toolkit, but fails
        """
        from openforcefield.topology import Molecule
        PARTIAL_CHARGE_METHODS = {'zeros':         {'rec_confs':0,
                                                    'min_confs':0,
                                                    'max_confs':0},
                                  'formal_charge': {'rec_confs':0,
                                                    'min_confs':0,
                                                    'max_confs':0}
                                  }

        if partial_charge_method is None:
            partial_charge_method = 'formal_charge'

        # Make a temporary copy of the molecule, since we'll be messing with its conformers
        mol_copy = Molecule(molecule)

        partial_charge_method = partial_charge_method.lower()
        if partial_charge_method not in PARTIAL_CHARGE_METHODS:
            raise ChargeMethodUnavailableError(f'Partial charge method "{partial_charge_method}"" is not supported by '
                                               f'the Built-in toolkit. Available charge methods are '
                                               f'{list(PARTIAL_CHARGE_METHODS.keys())}')

        if use_conformers is None:
            # Note that this refers back to the GLOBAL_TOOLKIT_REGISTRY by default, since
            # BuiltInToolkitWrapper can't generate conformers
            mol_copy.generate_conformers(n_conformers=PARTIAL_CHARGE_METHODS[partial_charge_method]['rec_confs'])
        else:
            mol_copy._conformers = None
            for conformer in use_conformers:
                mol_copy.add_conformer(conformer)
            self._check_n_conformers(mol_copy, partial_charge_method=partial_charge_method, min_confs=0,
                                     max_confs=0,strict_n_conformers=strict_n_conformers)

        partial_charges = unit.Quantity(np.zeros((molecule.n_particles)), unit.elementary_charge)
        if partial_charge_method == 'zeroes':
            pass
        elif partial_charge_method == 'formal_charge':
            for part_idx, particle in enumerate(molecule.particles):
                partial_charges[part_idx] = particle.formal_charge

        molecule.partial_charges = partial_charges



@inherit_docstrings
class OpenEyeToolkitWrapper(ToolkitWrapper):
    """
    OpenEye toolkit wrapper

    .. warning :: This API is experimental and subject to change.
    """
    _toolkit_name = 'OpenEye Toolkit'
    _toolkit_installation_instructions = 'The OpenEye toolkit requires a (free for academics) license, and can be ' \
                                         'found at: ' \
                                         'https://docs.eyesopen.com/toolkits/python/quickstart-python/install.html'

    def __init__(self):

        self._toolkit_file_read_formats = [
            'CAN', 'CDX', 'CSV', 'FASTA', 'INCHI', 'INCHIKEY', 'ISM', 'MDL', 'MF',
            'MMOD', 'MOL2', 'MOL2H', 'MOPAC', 'OEB', 'PDB', 'RDF', 'SDF', 'SKC',
            'SLN', 'SMI', 'USM', 'XYC'
        ]
        self._toolkit_file_write_formats = [
            'CAN', 'CDX', 'CSV', 'FASTA', 'INCHI', 'INCHIKEY', 'ISM', 'MDL', 'MF',
            'MMOD', 'MOL2', 'MOL2H', 'MOPAC', 'OEB', 'PDB', 'RDF', 'SDF', 'SKC',
            'SLN', 'SMI', 'USM', 'XYC'
        ]

        # check if the toolkit can be loaded
        if not self.is_available():
            raise ToolkitUnavailableException(f'The required toolkit {self._toolkit_name} is not '
                                              f'available. {self._toolkit_installation_instructions}')

    @staticmethod
    def is_available(
            oetools=('oechem', 'oequacpac', 'oeiupac', 'oeomega')):
        """
        Check if the given OpenEye toolkit components are available.

        If the OpenEye toolkit is not installed or no license is found
        for at least one the given toolkits , ``False`` is returned.

        Parameters
        ----------
        oetools : str or iterable of strings, optional, default=('oechem', 'oequacpac', 'oeiupac', 'oeomega')
            Set of tools to check by their Python module name. Defaults
            to the complete set of tools supported by this function.
            Also accepts a single tool to check as a string instead of
            an iterable of length 1.

        Returns
        -------
        all_installed : bool
            ``True`` if all tools in ``oetools`` are installed and licensed,
            ``False`` otherwise

        """
        # Complete list of module -> license function to check.
        license_function_names = {
            'oechem': 'OEChemIsLicensed',
            'oequacpac': 'OEQuacPacIsLicensed',
            'oeiupac': 'OEIUPACIsLicensed',
            'oeomega': 'OEOmegaIsLicensed'
        }
        supported_tools = set(license_function_names.keys())

        # Make sure oetools is a set.
        if isinstance(oetools, str):
            oetools = {oetools}
        else:
            oetools = set(oetools)

        # Check for unkown tools.
        unknown_tools = oetools.difference(supported_tools)
        if len(unknown_tools) > 0:
            raise ValueError("Found unkown OpenEye tools: {}. Supported values are: {}".format(
                sorted(unknown_tools), sorted(supported_tools)))

        # Check license of all tools.
        all_licensed = True
        for tool in oetools:
            try:
                module = importlib.import_module('openeye.' + tool)
            except (ImportError, ModuleNotFoundError):
                return False
            else:
                all_licensed &= getattr(module, license_function_names[tool])()
        return all_licensed

    def from_object(self, object, allow_undefined_stereo=False):
        """
        If given an OEMol (or OEMol-derived object), this function will load it into an openforcefield.topology.molecule

        Parameters
        ----------
        object : A molecule-like object
            An object to by type-checked.
        allow_undefined_stereo : bool, default=False
            Whether to accept molecules with undefined stereocenters. If False,
            an exception will be raised if a molecule with undefined stereochemistry
            is passed into this function.

        Returns
        -------
        Molecule
            An openforcefield.topology.molecule Molecule.

        Raises
        ------
        ValueError
            If the object could not be converted into a Molecule.
        """
        # TODO: Add tests for the from_object functions
        from openeye import oechem
        if isinstance(object, oechem.OEMolBase):
            return self.from_openeye(object, allow_undefined_stereo=allow_undefined_stereo)
        raise TypeError('Cannot create Molecule from {} object'.format(type(object)))

    def from_file(self,
                  file_path,
                  file_format,
                  allow_undefined_stereo=False):
        """
        Return an openforcefield.topology.Molecule from a file using this toolkit.

        Parameters
        ----------
        file_path : str
            The file to read the molecule from
        file_format : str
            Format specifier, usually file suffix (eg. 'MOL2', 'SMI')
            Note that not all toolkits support all formats. Check ToolkitWrapper.toolkit_file_read_formats for details.
        allow_undefined_stereo : bool, default=False
            If false, raises an exception if oemol contains undefined stereochemistry.

        Returns
        -------
        molecules : List[Molecule]
            The list of ``Molecule`` objects in the file.

        Raises
        ------
        GAFFAtomTypeWarning
            If the loaded mol2 file possibly uses GAFF atom types, which
            are not supported.

        Examples
        --------

        Load a mol2 file into an OpenFF ``Molecule`` object.

        >>> from openforcefield.utils import get_data_file_path
        >>> mol2_file_path = get_data_file_path('molecules/cyclohexane.mol2')
        >>> toolkit = OpenEyeToolkitWrapper()
        >>> molecule = toolkit.from_file(mol2_file_path, file_format='mol2')

        """
        from openeye import oechem
        ifs = oechem.oemolistream(file_path)
        return self._read_oemolistream_molecules(ifs, allow_undefined_stereo, file_path=file_path)

    def from_file_obj(self,
                      file_obj,
                      file_format,
                      allow_undefined_stereo=False):
        """
        Return an openforcefield.topology.Molecule from a file-like object (an object with a ".read()" method using
        this toolkit.

        Parameters
        ----------
        file_obj : file-like object
            The file-like object to read the molecule from
        file_format : str
            Format specifier, usually file suffix (eg. 'MOL2', 'SMI')
            Note that not all toolkits support all formats. Check ToolkitWrapper.toolkit_file_read_formats for details.
        allow_undefined_stereo : bool, default=False
            If false, raises an exception if oemol contains undefined stereochemistry.

        Returns
        -------
        molecules : List[Molecule]
            The list of Molecule objects in the file object.

        Raises
        ------
        GAFFAtomTypeWarning
            If the loaded mol2 file possibly uses GAFF atom types, which
            are not supported.

        """
        from openeye import oechem

        # Configure input molecule stream.
        ifs = oechem.oemolistream()
        ifs.openstring(file_obj.read())
        oeformat = getattr(oechem, 'OEFormat_' + file_format)
        ifs.SetFormat(oeformat)

        return self._read_oemolistream_molecules(ifs, allow_undefined_stereo)

    def to_file_obj(self, molecule, file_obj, file_format):
        """
        Writes an OpenFF Molecule to a file-like object

        Parameters
        ----------
        molecule : an OpenFF Molecule
            The molecule to write
        file_obj
            The file-like object to write to
        file_format
            The format for writing the molecule data

        """
        from openeye import oechem
        from openforcefield.utils import temporary_directory, temporary_cd

        oemol = self.to_openeye(molecule)

        # TODO: This is inefficiently implemented. Is there any way to attach a file-like object to an oemolstream?
        with temporary_directory() as tmpdir:
            with temporary_cd(tmpdir):
                outfile = 'temp_molecule.' + file_format
                ofs = oechem.oemolostream(outfile)
                openeye_format = getattr(oechem, 'OEFormat_' + file_format)
                ofs.SetFormat(openeye_format)
                oechem.OEWriteMolecule(ofs, oemol)
                ofs.close()
                file_data = open(outfile).read()
        file_obj.write(file_data)

    def to_file(self, molecule, file_path, file_format):
        """
        Writes an OpenFF Molecule to a file-like object

        Parameters
        ----------
        molecule : an OpenFF Molecule
            The molecule to write
        file_path
            The file path to write to.
        file_format
            The format for writing the molecule data

        """
        from openeye import oechem
        oemol = self.to_openeye(molecule)
        ofs = oechem.oemolostream(file_path)
        openeye_format = getattr(oechem, 'OEFormat_' + file_format)
        ofs.SetFormat(openeye_format)
        oechem.OEWriteMolecule(ofs, oemol)
        ofs.close()

    @classmethod
    def _read_oemolistream_molecules(cls, oemolistream, allow_undefined_stereo, file_path=None):
        """
        Reads and return the Molecules in a OEMol input stream.

        Parameters
        ----------
        oemolistream : oechem.oemolistream
            The OEMol input stream to read from.
        allow_undefined_stereo : bool
            If false, raises an exception if oemol contains undefined stereochemistry.
        file_path : str, optional
            The path to the mol2 file. This is used exclusively to make
            the error message more meaningful when the mol2 files doesn't
            use Tripos atom types.

        Returns
        -------
        molecules : List[Molecule]
            The list of Molecule objects in the stream.

        """
        from openforcefield.topology import Molecule
        from openeye import oechem

        mols = list()
        oemol = oechem.OEMol()
        while oechem.OEReadMolecule(oemolistream, oemol):
            oechem.OEPerceiveChiral(oemol)
            oechem.OEAssignAromaticFlags(oemol, oechem.OEAroModel_MDL)
            oechem.OE3DToInternalStereo(oemol)
            mol = cls.from_openeye(
                oemol,
                allow_undefined_stereo=allow_undefined_stereo)
            mols.append(mol)

            # Check if this file may be using GAFF atom types.
            if oemolistream.GetFormat() == oechem.OEFormat_MOL2:
                cls._check_mol2_gaff_atom_type(mol, file_path)

        return mols

    @staticmethod
    def _check_mol2_gaff_atom_type(molecule, file_path=None):
        """Attempts to detect the presence of GAFF atom types in a molecule loaded from a mol2 file.

        For now, this raises a ``GAFFAtomTypeWarning`` if the molecule
        include Osmium and Holmium atoms, which have GAFF types OS and
        HO respectively.

        Parameters
        ----------
        molecule : openforcefield.topology.molecule.Molecule
            The loaded molecule.
        file_path : str, optional
            The path to the mol2 file. This is used exclusively to make
            the error message more meaningful.

        """
        # Handle default.
        if file_path is None:
            file_path = ''
        else:
            # Append a ':' character that will separate the file
            # path from the molecule string representation.
            file_path = file_path + ':'
        # atomic_number: (GAFF_type, element_name)
        warning_atomic_numbers = {
            76: ('OS', 'Osmium'),
            67: ('HO', 'Holmium')
        }

        for atom in molecule.atoms:
            try:
                atom_type, element_name = warning_atomic_numbers[atom.atomic_number]
            except KeyError:
                pass
            else:
                import warnings
                warn_msg = (f'OpenEye interpreted the type "{atom_type}" in {file_path}{molecule.name}'
                            f' as {element_name}. Does your mol2 file uses Tripos SYBYL atom types?'
                            ' Other atom types such as GAFF are not supported.')
                warnings.warn(warn_msg, GAFFAtomTypeWarning)

    @staticmethod
    def _openeye_cip_atom_stereochemistry(oemol, oeatom):
        """
        Determine CIP stereochemistry (R/S) for the specified atom

        Parameters
        ----------
        oemol : openeye.oechem.OEMolBase
            The molecule of interest
        oeatom : openeye.oechem.OEAtomBase
            The atom whose stereochemistry is to be computed

        Returns
        -------
        stereochemistry : str
            'R', 'S', or None if no stereochemistry is specified or the atom is not a stereocenter
        """
        from openeye import oechem

        if not oeatom.HasStereoSpecified():
            # No stereochemical information has been stored, so this could be unknown stereochemistry
            # TODO: Should we raise an exception?
            return None

        cip = oechem.OEPerceiveCIPStereo(oemol, oeatom)

        if cip == oechem.OECIPAtomStereo_S:
            return 'S'
        elif cip == oechem.OECIPAtomStereo_R:
            return 'R'
        elif cip == oechem.OECIPAtomStereo_NotStereo:
            # Not a stereocenter
            # TODO: Should this be a different case from ``None``?
            return None

    @staticmethod
    def _openeye_cip_bond_stereochemistry(oemol, oebond):
        """
        Determine CIP stereochemistry (E/Z) for the specified bond

        Parameters
        ----------
        oemol : openeye.oechem.OEMolBase
            The molecule of interest
        oebond : openeye.oechem.OEBondBase
            The bond whose stereochemistry is to be computed

        Returns
        -------
        stereochemistry : str
            'E', 'Z', or None if stereochemistry is unspecified or the bond is not a stereo bond

        """
        from openeye import oechem

        if not oebond.HasStereoSpecified():
            # No stereochemical information has been stored, so this could be unknown stereochemistry
            # TODO: Should we raise an exception?
            return None

        cip = oechem.OEPerceiveCIPStereo(oemol, oebond)

        if cip == oechem.OECIPBondStereo_E:
            return 'E'
        elif cip == oechem.OECIPBondStereo_Z:
            return 'Z'
        elif cip == oechem.OECIPBondStereo_NotStereo:
            return None

    @staticmethod
    def from_openeye(oemol, allow_undefined_stereo=False):
        """
        Create a Molecule from an OpenEye molecule. If the OpenEye molecule has
        implicit hydrogens, this function will make them explicit.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        oemol : openeye.oechem.OEMol
            An OpenEye molecule
        allow_undefined_stereo : bool, default=False
            If false, raises an exception if oemol contains undefined stereochemistry.

        Returns
        -------
        molecule : openforcefield.topology.Molecule
            An openforcefield molecule

        Examples
        --------

        Create a Molecule from an OpenEye OEMol

        >>> from openeye import oechem
        >>> from openforcefield.tests.utils import get_data_file_path
        >>> ifs = oechem.oemolistream(get_data_file_path('systems/monomers/ethanol.mol2'))
        >>> oemols = list(ifs.GetOEGraphMols())

        >>> toolkit_wrapper = OpenEyeToolkitWrapper()
        >>> molecule = toolkit_wrapper.from_openeye(oemols[0])

        """
        from openeye import oechem
        from openforcefield.topology.molecule import Molecule


        # Add explicit hydrogens if they're implicit
        if oechem.OEHasImplicitHydrogens(oemol):
            oechem.OEAddExplicitHydrogens(oemol)

        # TODO: Is there any risk to perceiving aromaticity here instead of later?
        oechem.OEAssignAromaticFlags(oemol, oechem.OEAroModel_MDL)

        oechem.OEPerceiveChiral(oemol)

        # Check that all stereo is specified
        # Potentially better OE stereo check: OEFlipper — Toolkits - - Python
        # https: // docs.eyesopen.com / toolkits / python / omegatk / OEConfGenFunctions / OEFlipper.html

        unspec_chiral = False
        unspec_db = False
        problematic_atoms = list()
        problematic_bonds = list()

        for oeatom in oemol.GetAtoms():
            if oeatom.IsChiral():
                if not (oeatom.HasStereoSpecified()):
                    unspec_chiral = True
                    problematic_atoms.append(oeatom)
        for oebond in oemol.GetBonds():
            if oebond.IsChiral():
                if not (oebond.HasStereoSpecified()):
                    unspec_db = True
                    problematic_bonds.append(oebond)
        if (unspec_chiral or unspec_db):

            def oeatom_to_str(oeatom):
                return 'atomic num: {}, name: {}, idx: {}, aromatic: {}, chiral: {}'.format(
                    oeatom.GetAtomicNum(), oeatom.GetName(), oeatom.GetIdx(),
                    oeatom.IsAromatic(), oeatom.IsChiral())

            def oebond_to_str(oebond):
                return "order: {}, chiral: {}".format(oebond.GetOrder(),
                                                      oebond.IsChiral())

            def describe_oeatom(oeatom):
                description = "Atom {} with bonds:".format(
                    oeatom_to_str(oeatom))
                for oebond in oeatom.GetBonds():
                    description += "\nbond {} to atom {}".format(
                        oebond_to_str(oebond),
                        oeatom_to_str(oebond.GetNbr(oeatom)))
                return description

            msg = "OEMol has unspecified stereochemistry. " \
                  "oemol.GetTitle(): {}\n".format(oemol.GetTitle())
            if len(problematic_atoms) != 0:
                msg += "Problematic atoms are:\n"
                for problematic_atom in problematic_atoms:
                    msg += describe_oeatom(problematic_atom) + '\n'
            if len(problematic_bonds) != 0:
                msg += "Problematic bonds are: {}\n".format(problematic_bonds)
            if allow_undefined_stereo:
                msg = 'Warning (not error because allow_undefined_stereo=True): ' + msg
                logger.warning(msg)
            else:
                msg = 'Unable to make OFFMol from OEMol: ' + msg
                raise UndefinedStereochemistryError(msg)

        # TODO: What other information should we preserve besides name?
        # TODO: How should we preserve the name?

        molecule = Molecule()
        molecule.name = oemol.GetTitle()

        # Copy any attached SD tag information
        # TODO: Should we use an API for this?
        molecule._properties = dict()
        for dp in oechem.OEGetSDDataPairs(oemol):
            molecule._properties[dp.GetTag()] = dp.GetValue()

        map_atoms = dict()  # {oemol_idx: molecule_idx}
        atom_mapping = {}
        for oeatom in oemol.GetAtoms():
            oe_idx = oeatom.GetIdx()
            map_id = oeatom.GetMapIdx()
            atomic_number = oeatom.GetAtomicNum()
            formal_charge = oeatom.GetFormalCharge() * unit.elementary_charge
            is_aromatic = oeatom.IsAromatic()
            stereochemistry = OpenEyeToolkitWrapper._openeye_cip_atom_stereochemistry(
                oemol, oeatom)
            #stereochemistry = self._openeye_cip_atom_stereochemistry(oemol, oeatom)
            name = ''
            if oeatom.HasData('name'):
                name = oeatom.GetData('name')
            atom_index = molecule.add_atom(
                atomic_number,
                formal_charge,
                is_aromatic,
                stereochemistry=stereochemistry,
                name=name)
            map_atoms[
                oe_idx] = atom_index  # store for mapping oeatom to molecule atom indices below
            atom_mapping[atom_index] = map_id

        # if we have a full atom map add it to the molecule, 0 indicates a missing mapping or no mapping
        if 0 not in atom_mapping.values():
            molecule._properties['atom_map'] = atom_mapping

        for oebond in oemol.GetBonds():
            atom1_index = map_atoms[oebond.GetBgnIdx()]
            atom2_index = map_atoms[oebond.GetEndIdx()]
            bond_order = oebond.GetOrder()
            is_aromatic = oebond.IsAromatic()
            stereochemistry = OpenEyeToolkitWrapper._openeye_cip_bond_stereochemistry(
                oemol, oebond)
            if oebond.HasData('fractional_bond_order'):
                fractional_bond_order = oebond.GetData('fractional_bond_order')
            else:
                fractional_bond_order = None

            molecule.add_bond(
                atom1_index,
                atom2_index,
                bond_order,
                is_aromatic=is_aromatic,
                stereochemistry=stereochemistry,
                fractional_bond_order=fractional_bond_order)

        # TODO: Copy conformations, if present
        # TODO: Come up with some scheme to know when to import coordinates
        # From SMILES: no
        # From MOL2: maybe
        # From other: maybe
        if hasattr(oemol, 'GetConfs'):
            for conf in oemol.GetConfs():
                n_atoms = molecule.n_atoms
                positions = unit.Quantity(
                    np.zeros([n_atoms, 3], np.float), unit.angstrom)
                for oe_id in conf.GetCoords().keys():
                    off_atom_coords = unit.Quantity(conf.GetCoords()[oe_id],
                                                    unit.angstrom)
                    off_atom_index = map_atoms[oe_id]
                    positions[off_atom_index, :] = off_atom_coords
                if (positions == 0 * unit.angstrom).all() and n_atoms > 1:
                    continue
                molecule.add_conformer(positions)

        # Copy partial charges, if present
        partial_charges = unit.Quantity(
            np.zeros(molecule.n_atoms, dtype=np.float),
            unit=unit.elementary_charge)
        for oe_atom in oemol.GetAtoms():
            oe_idx = oe_atom.GetIdx()
            off_idx = map_atoms[oe_idx]
            charge = oe_atom.GetPartialCharge() * unit.elementary_charge
            partial_charges[off_idx] = charge

        molecule.partial_charges = partial_charges

        return molecule

    @staticmethod
    def to_openeye(molecule, aromaticity_model=DEFAULT_AROMATICITY_MODEL):
        """
        Create an OpenEye molecule using the specified aromaticity model

        .. todo ::

           * Should the aromaticity model be specified in some other way?

       .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        molecule : openforcefield.topology.molecule.Molecule object
            The molecule to convert to an OEMol
        aromaticity_model : str, optional, default=DEFAULT_AROMATICITY_MODEL
            The aromaticity model to use

        Returns
        -------
        oemol : openeye.oechem.OEMol
            An OpenEye molecule

        Examples
        --------

        Create an OpenEye molecule from a Molecule

        >>> from openforcefield.topology import Molecule
        >>> toolkit_wrapper = OpenEyeToolkitWrapper()
        >>> molecule = Molecule.from_smiles('CC')
        >>> oemol = toolkit_wrapper.to_openeye(molecule)

        """
        from openeye import oechem

        oemol = oechem.OEMol()
        #if not(molecule.name is None):
        oemol.SetTitle(molecule.name)
        map_atoms = {}  # {off_idx : oe_idx}
        # Add atoms
        oemol_atoms = list()  # list of corresponding oemol atoms
        for atom in molecule.atoms:
            oeatom = oemol.NewAtom(atom.atomic_number)
            oeatom.SetFormalCharge(atom.formal_charge.value_in_unit(unit.elementary_charge)) # simtk.unit.Quantity(1, unit.elementary_charge)
            oeatom.SetAromatic(atom.is_aromatic)
            oeatom.SetData('name', atom.name)
            oemol_atoms.append(oeatom)
            map_atoms[atom.molecule_atom_index] = oeatom.GetIdx()

        # Add bonds
        oemol_bonds = list()  # list of corresponding oemol bonds
        for bond in molecule.bonds:
            #atom1_index = molecule.atoms.index(bond.atom1)
            #atom2_index = molecule.atoms.index(bond.atom2)
            atom1_index = bond.atom1_index
            atom2_index = bond.atom2_index
            oebond = oemol.NewBond(oemol_atoms[atom1_index],
                                   oemol_atoms[atom2_index])
            oebond.SetOrder(bond.bond_order)
            oebond.SetAromatic(bond.is_aromatic)
            if not (bond.fractional_bond_order is None):
                oebond.SetData('fractional_bond_order',
                               bond.fractional_bond_order)
            oemol_bonds.append(oebond)

        # Set atom stereochemistry now that all connectivity is in place
        for atom, oeatom in zip(molecule.atoms, oemol_atoms):
            if not atom.stereochemistry:
                continue

            # Set arbitrary initial stereochemistry
            neighs = [n for n in oeatom.GetAtoms()]
            oeatom.SetStereo(neighs, oechem.OEAtomStereo_Tetra,
                             oechem.OEAtomStereo_Right)

            # Flip chirality if stereochemistry isincorrect
            oeatom_stereochemistry = OpenEyeToolkitWrapper._openeye_cip_atom_stereochemistry(
                oemol, oeatom)
            if oeatom_stereochemistry != atom.stereochemistry:
                # Flip the stereochemistry
                oeatom.SetStereo(neighs, oechem.OEAtomStereo_Tetra,
                                 oechem.OEAtomStereo_Left)
                # Verify it matches now as a sanity check
                oeatom_stereochemistry = OpenEyeToolkitWrapper._openeye_cip_atom_stereochemistry(
                    oemol, oeatom)
                if oeatom_stereochemistry != atom.stereochemistry:
                    raise Exception(
                        'Programming error: OpenEye atom stereochemistry assumptions failed.'
                    )

        # Set bond stereochemistry
        for bond, oebond in zip(molecule.bonds, oemol_bonds):
            if not bond.stereochemistry:
                continue

            atom1_index = bond.molecule.atoms.index(bond.atom1)
            atom2_index = bond.molecule.atoms.index(bond.atom2)
            # Set arbitrary initial stereochemistry
            oeatom1, oeatom2 = oemol_atoms[atom1_index], oemol_atoms[
                atom2_index]
            oeatom1_neighbor = [
                n for n in oeatom1.GetAtoms() if not n == oeatom2
            ][0]
            oeatom2_neighbor = [
                n for n in oeatom2.GetAtoms() if not n == oeatom1
            ][0]
            #oebond.SetStereo([oeatom1, oeatom2], oechem.OEBondStereo_CisTrans, oechem.OEBondStereo_Cis)
            oebond.SetStereo([oeatom1_neighbor, oeatom2_neighbor],
                             oechem.OEBondStereo_CisTrans,
                             oechem.OEBondStereo_Cis)

            # Flip stereochemistry if incorrect
            oebond_stereochemistry = OpenEyeToolkitWrapper._openeye_cip_bond_stereochemistry(
                oemol, oebond)
            #print('AAAA', oebond_stereochemistry, bond.stereochemistry)
            if oebond_stereochemistry != bond.stereochemistry:
                # Flip the stereochemistry
                #oebond.SetStereo([oeatom1, oeatom2], oechem.OEBondStereo_CisTrans, oechem.OEBondStereo_Trans)
                oebond.SetStereo([oeatom1_neighbor, oeatom2_neighbor],
                                 oechem.OEBondStereo_CisTrans,
                                 oechem.OEBondStereo_Trans)
                # Verify it matches now as a sanity check
                oebond_stereochemistry = OpenEyeToolkitWrapper._openeye_cip_bond_stereochemistry(
                    oemol, oebond)
                if oebond_stereochemistry != bond.stereochemistry:
                    raise Exception(
                        'Programming error: OpenEye bond stereochemistry assumptions failed.'
                    )

        # Retain conformations, if present
        if molecule.n_conformers != 0:
            oemol.DeleteConfs()
            for conf in molecule._conformers:
                # OE needs a 1 x (3*n_Atoms) double array as input
                flat_coords = np.zeros((oemol.NumAtoms() * 3),
                                       dtype=np.float32)
                for index, oe_idx in map_atoms.items():
                    (x, y, z) = conf[index, :] / unit.angstrom
                    flat_coords[(3 * oe_idx)] = x
                    flat_coords[(3 * oe_idx) + 1] = y
                    flat_coords[(3 * oe_idx) + 2] = z

                # TODO: Do we need to do these internal unit checks?
                # TODO: Is there any risk that the atom indexing systems will change?
                #flat_coords = (conf.in_units_of(unit.angstrom) / unit.angstrom).flatten()
                oecoords = oechem.OEFloatArray(flat_coords)
                oemol.NewConf(oecoords)

        # Retain charges, if present
        if not (molecule._partial_charges is None):
            # for off_atom, oe_atom in zip(molecule.atoms, oemol_atoms):
            #    charge_unitless = molecule._partial_charges

            oe_indexed_charges = np.zeros((molecule.n_atoms), dtype=np.float)
            for off_idx, charge in enumerate(molecule._partial_charges):
                oe_idx = map_atoms[off_idx]
                charge_unitless = charge / unit.elementary_charge
                oe_indexed_charges[oe_idx] = charge_unitless
            for oe_idx, oe_atom in enumerate(oemol.GetAtoms()):
                oe_atom.SetPartialCharge(oe_indexed_charges[oe_idx])

        # TODO: Retain properties, if present
        # Clean Up phase
        # The only feature of a molecule that wasn't perceived above seemed to be ring connectivity, better to run it
        # here then for someone to inquire about ring sizes and get 0 when it shouldn't be
        oechem.OEFindRingAtomsAndBonds(oemol)

        return oemol

    @staticmethod
    def to_smiles(molecule):
        """
        Uses the OpenEye toolkit to convert a Molecule into a SMILES string.

        Parameters
        ----------
        molecule : An openforcefield.topology.Molecule
            The molecule to convert into a SMILES.

        Returns
        -------
        smiles : str
            The SMILES of the input molecule.
        """
        from openeye import oechem
        oemol = OpenEyeToolkitWrapper.to_openeye(molecule)
        smiles = oechem.OECreateSmiString(
            oemol, oechem.OESMILESFlag_DEFAULT | oechem.OESMILESFlag_Hydrogens
            | oechem.OESMILESFlag_Isotopes | oechem.OESMILESFlag_BondStereo
            | oechem.OESMILESFlag_AtomStereo)
        return smiles

    def to_inchi(self, molecule, fixed_hydrogens=False):
        """
        Create an InChI string for the molecule using the RDKit Toolkit.
        InChI is a standardised representation that does not capture tautomers unless specified using the fixed hydrogen
        layer.

        For information on InChi see here https://iupac.org/who-we-are/divisions/division-details/inchi/

        Parameters
        ----------
        molecule : An openforcefield.topology.Molecule
            The molecule to convert into a SMILES.

        fixed_hydrogens: bool, default=False
            If a fixed hydrogen layer should be added to the InChI, if `True` this will produce a non standard specific
            InChI string of the molecule.

        Returns
        --------
        inchi: str
            The InChI string of the molecule.
        """

        from openeye import oechem
        oemol = self.to_openeye(molecule)

        if fixed_hydrogens:
            opts = oechem.OEInChIOptions()
            opts.SetFixedHLayer(True)
            inchi = oechem.OEMolToInChI(oemol)

        else:
            inchi = oechem.OEMolToSTDInChI(oemol)

        return inchi

    def to_inchikey(self, molecule, fixed_hydrogens=False):
        """
        Create an InChIKey for the molecule using the RDKit Toolkit.
        InChIKey is a standardised representation that does not capture tautomers unless specified using the fixed hydrogen
        layer.

        For information on InChi see here https://iupac.org/who-we-are/divisions/division-details/inchi/

        Parameters
        ----------
        molecule : An openforcefield.topology.Molecule
            The molecule to convert into a SMILES.

        fixed_hydrogens: bool, default=False
            If a fixed hydrogen layer should be added to the InChI, if `True` this will produce a non standard specific
            InChI string of the molecule.

        Returns
        --------
        inchi_key: str
            The InChIKey representation of the molecule.
        """

        from openeye import oechem
        oemol = self.to_openeye(molecule)

        if fixed_hydrogens:
            opts = oechem.OEInChIOptions()
            opts.SetFixedHLayer(True)
            inchi_key = oechem.OEMolToInChIKey(oemol)

        else:
            inchi_key = oechem.OEMolToSTDInChIKey(oemol)

        return inchi_key

    def canonical_order_atoms(self, molecule):
        """
        Canonical order the atoms in the molecule using the OpenEye toolkit.

        Parameters
        ----------
        molecule: openforcefield.topology.Molecule
            The input molecule

         Returns
        -------
        molecule : openforcefield.topology.Molecule
            The input molecule, with canonically-indexed atoms and bonds.
        """

        from openeye import oechem
        oemol = self.to_openeye(molecule)

        oechem.OECanonicalOrderAtoms(oemol)
        oechem.OECanonicalOrderBonds(oemol)

        # reorder the iterator
        vatm = []
        for atom in oemol.GetAtoms():
            if atom.GetAtomicNum() != oechem.OEElemNo_H:
                vatm.append(atom)
        oemol.OrderAtoms(vatm)

        vbnd = []
        for bond in oemol.GetBonds():
            if bond.GetBgn().GetAtomicNum() != oechem.OEElemNo_H and bond.GetEnd().GetAtomicNum() != oechem.OEElemNo_H:
                vbnd.append(bond)
        oemol.OrderBonds(vbnd)

        oemol.Sweep()

        for bond in oemol.GetBonds():
            if bond.GetBgnIdx() > bond.GetEndIdx():
                bond.SwapEnds()

        return self.from_openeye(oemol)

    def from_smiles(self, smiles, hydrogens_are_explicit=False, allow_undefined_stereo=False):
        """
        Create a Molecule from a SMILES string using the OpenEye toolkit.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        smiles : str
            The SMILES string to turn into a molecule
        hydrogens_are_explicit : bool, default = False
            If False, OE will perform hydrogen addition using OEAddExplicitHydrogens
        allow_undefined_stereo : bool, default=False
            Whether to accept SMILES with undefined stereochemistry. If False,
            an exception will be raised if a SMILES with undefined stereochemistry
            is passed into this function.

        Returns
        -------
        molecule : openforcefield.topology.Molecule
            An openforcefield-style molecule.
        """

        from openeye import oechem
        oemol = oechem.OEGraphMol()
        oechem.OESmilesToMol(oemol, smiles)
        if not (hydrogens_are_explicit):
            result = oechem.OEAddExplicitHydrogens(oemol)
            if result == False:
                raise ValueError(
                    "Addition of explicit hydrogens failed in from_openeye")
        elif hydrogens_are_explicit and oechem.OEHasImplicitHydrogens(oemol):
            raise ValueError(
                f"'hydrogens_are_explicit' was specified as True, but OpenEye Toolkit interpreted "
                f"SMILES '{smiles}' as having implicit hydrogen. If this SMILES is intended to "
                f"express all explicit hydrogens in the molecule, then you should construct the "
                f"desired molecule as an OEMol (where oechem.OEHasImplicitHydrogens(oemol) returns "
                f"False), and then use Molecule.from_openeye() to create the desired OFFMol.")
        molecule = self.from_openeye(oemol,
                                     allow_undefined_stereo=allow_undefined_stereo)
        return molecule

    def from_inchi(self, inchi, allow_undefined_stereo=False):
        """
        Construct a Molecule from a InChI representation

        Parameters
        ----------
        inchi : str
            The InChI representation of the molecule.

        allow_undefined_stereo : bool, default=False
            Whether to accept InChI with undefined stereochemistry. If False,
            an exception will be raised if a InChI with undefined stereochemistry
            is passed into this function.

        Returns
        -------
        molecule : openforcefield.topology.Molecule
        """

        from openeye import oechem
        # This calls the same functions as OESmilesToMol
        oemol = oechem.OEGraphMol()
        oechem.OEInChIToMol(oemol, inchi)

        # try and catch InChI parsing fails
        # if there are no atoms don't build the molecule
        if oemol.NumAtoms() == 0:
            raise RuntimeError('There was an issue parsing the InChI string, please check and try again.')

        molecule = self.from_openeye(oemol, allow_undefined_stereo=allow_undefined_stereo)

        return molecule

    def generate_conformers(self, molecule, n_conformers=1, rms_cutoff=None, clear_existing=True):
        """
        Generate molecule conformers using OpenEye Omega.

        .. warning :: This API is experimental and subject to change.

        .. todo ::

           * which parameters should we expose? (or can we implement a general system with **kwargs?)
           * will the coordinates be returned in the OpenFF Molecule's own indexing system? Or is there a chance that
           they'll get reindexed when we convert the input into an OEmol?

        Parameters
        ---------
        molecule : a :class:`Molecule`
            The molecule to generate conformers for.
        n_conformers : int, default=1
            The maximum number of conformers to generate.
        rms_cutoff : simtk.Quantity-wrapped float, in units of distance, optional, default=None
            The minimum RMS value at which two conformers are considered redundant and one is deleted.
            If None, the cutoff is set to 1 Angstrom
        clear_existing : bool, default=True
            Whether to overwrite existing conformers for the molecule

        """
        from openeye import oeomega
        oemol = self.to_openeye(molecule)
        omega = oeomega.OEOmega()
        omega.SetMaxConfs(n_conformers)
        omega.SetCanonOrder(False)
        omega.SetSampleHydrogens(True)
        omega.SetEnergyWindow(15.0)  #unit?
        if rms_cutoff is None:
            omega.SetRMSThreshold(1.0)
        else:
            omega.SetRMSThreshold(rms_cutoff.value_in_unit(unit.angstrom))
        # Don't generate random stereoisomer if not specified
        omega.SetStrictStereo(True)
        status = omega(oemol)

        if status is False:
            omega.SetStrictStereo(False)
            new_status = omega(oemol)
            if new_status is False:
                raise Exception("OpenEye Omega conformer generation failed")


        molecule2 = self.from_openeye(oemol, allow_undefined_stereo=True)

        if clear_existing:
            molecule._conformers = list()

        for conformer in molecule2._conformers:
            molecule._add_conformer(conformer)


    def assign_partial_charges(self,
                               molecule,
                               partial_charge_method=None,
                               use_conformers=False,
                               strict_n_conformers=False):
        """
        Compute partial charges with OpenEye quacpac, and assign
        the new values to the partial_charges attribute.

        .. warning :: This API is experimental and subject to change.

        .. todo ::

           * Should the default be ELF?
           * Can we expose more charge models?


        Parameters
        ----------
        molecule : Molecule
            Molecule for which partial charges are to be computed
        partial_charge_method : str, optional, default=None
            The charge model to use. One of ['amberff94', 'mmff', 'mmff94', `am1-mulliken`, 'am1bcc',
            'am1bccnosymspt', 'am1bccelf10']
            If None, 'am1-mulliken' will be used.
        use_conformers : iterable of simtk.unit.Quantity-wrapped numpy arrays, each with shape (n_atoms, 3) and dimension of distance. Optional, default = None
            Coordinates to use for partial charge calculation. If None, an appropriate number of conformers will be generated.
        strict_n_conformers : bool, default=False
            Whether to raise an exception if an invalid number of conformers is provided for the given charge method.
            If this is False and an invalid number of conformers is found, a warning will be raised.

        Raises
        ------
        ChargeMethodUnavailableError if the requested charge method can not be handled by this toolkit

        ChargeCalculationError if the charge calculation is supported by this toolkit, but fails
        """

        from openeye import oequacpac, oechem
        import numpy as np
        from openforcefield.topology import Molecule

        SUPPORTED_CHARGE_METHODS = {'am1bcc': {'oe_charge_method': oequacpac.OEAM1BCCCharges,
                                               'min_confs': 1,
                                               'max_confs': 1,
                                               'rec_confs': 1,
                                               },
                                    'am1-mulliken': {'oe_charge_method': oequacpac.OEAM1Charges,
                                                     'min_confs': 1,
                                                     'max_confs': 1,
                                                     'rec_confs': 1,
                                                     },
                                    'gasteiger': {'oe_charge_method': oequacpac.OEGasteigerCharges,
                                                  'min_confs': 0,
                                                  'max_confs': 0,
                                                  'rec_confs': 0,
                                                  },
                                    'mmff94': {'oe_charge_method': oequacpac.OEMMFF94Charges,
                                                  'min_confs': 0,
                                                  'max_confs': 0,
                                                  'rec_confs': 0,
                                                  },
                                    'am1bccnosymspt': {'oe_charge_method': oequacpac.OEAM1BCCCharges,
                                                       'min_confs': 1,
                                                       'max_confs': 1,
                                                       'rec_confs': 1,
                                                       },
                                    'am1bccelf10': {'oe_charge_method': oequacpac.OEAM1BCCELF10Charges,
                                                  'min_confs': 1,
                                                  'max_confs': None,
                                                  'rec_confs': 500,
                                                  },
                                    }

        if partial_charge_method is None:
            partial_charge_method = "am1-mulliken"

        partial_charge_method = partial_charge_method.lower()

        if partial_charge_method not in SUPPORTED_CHARGE_METHODS:
            raise ChargeMethodUnavailableError(
                f"partial_charge_method '{partial_charge_method}' is not available from OpenEyeToolkitWrapper. "
                f"Available charge methods are {list(SUPPORTED_CHARGE_METHODS.keys())} "
            )

        charge_method = SUPPORTED_CHARGE_METHODS[partial_charge_method]
        # AVAILABLE_CHARGE_METHODS = ['am1-mulliken', 'gasteiger', 'mmff94', 'am1bcc', 'am1bccnosymspt', 'am1bccelf10']

        # Make a temporary copy of the molecule, since we'll be messing with its conformers
        mol_copy = Molecule(molecule)

        if use_conformers is None:
            if charge_method['rec_confs'] == 0:
                mol_copy._conformers = None
            else:
                self.generate_conformers(mol_copy,
                                         n_conformers=charge_method['rec_confs'],
                                         rms_cutoff=0.25*unit.angstrom)
                # TODO: What's a "best practice" RMS cutoff to use here?
        else:
            mol_copy._conformers = None
            for conformer in use_conformers:
                mol_copy.add_conformer(conformer)
            self._check_n_conformers(mol_copy,
                                     partial_charge_method=partial_charge_method,
                                     min_confs=charge_method['min_confs'],
                                     max_confs=charge_method['max_confs'],
                                     strict_n_conformers=strict_n_conformers)

        oemol = mol_copy.to_openeye()

        errfs = oechem.oeosstream()
        oechem.OEThrow.SetOutputStream(errfs)
        oechem.OEThrow.Clear()

        if partial_charge_method == "am1bccnosymspt":
            optimize = False
            symmetrize = False
            quacpac_status = oequacpac.OEAssignCharges(oemol, charge_method['oe_charge_method'](optimize, symmetrize))
        else:
            quacpac_status = oequacpac.OEAssignCharges(oemol, charge_method['oe_charge_method']())

        oechem.OEThrow.SetOutputStream(oechem.oeerr)  # restoring to original state
        # This logic handles errors encountered in #34, which can occur when using ELF10 conformer selection
        if not quacpac_status:
            if "SelectElfPop: issue with removing trans COOH conformers" in (errfs.str().decode("UTF-8")):
                logger.warning("Warning: OEAM1BCCELF10 charge assignment failed due to a known bug (toolkit issue "
                               "#346). Downgrading to OEAM1BCC charge assignment for this molecule. More information"
                               "is available at https://github.com/openforcefield/openforcefield/issues/346")
                quacpac_status = oequacpac.OEAssignCharges(oemol, oequacpac.OEAM1BCCCharges())

        if quacpac_status is False:
            raise ChargeCalculationError(f'Unable to assign charges: {errfs.str().decode("UTF-8")}')

        # Extract and return charges
        ## TODO: Behavior when given multiple conformations?
        ## TODO: Make sure atom mapping remains constant

        charges = unit.Quantity(
            np.zeros([oemol.NumAtoms()], np.float64), unit.elementary_charge)
        for index, atom in enumerate(oemol.GetAtoms()):
            charge = atom.GetPartialCharge()
            charge = charge * unit.elementary_charge
            charges[index] = charge

        if ((charges / unit.elementary_charge) == 0.).all():
            raise ChargeCalculationError(
                "Partial charge calculation failed. Charges from assign_partial_charges() are all 0."
            )
        molecule.partial_charges = charges


    def compute_partial_charges_am1bcc(self, molecule, use_conformers=None, strict_n_conformers=False):
        """
        Compute AM1BCC partial charges with OpenEye quacpac. This function will attempt to use
        the OEAM1BCCELF10 charge generation method, but may print a warning and fall back to
        normal OEAM1BCC if an error is encountered. This error is known to occur with some
        carboxylic acids, and is under investigation by OpenEye.


        .. warning :: This API is experimental and subject to change.

        .. todo ::

           * Should the default be ELF?
           * Can we expose more charge models?


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


        """
        from openeye import oequacpac, oechem
        import numpy as np

        import warnings
        warnings.warn("compute_partial_charges_am1bcc will be deprecated in an upcoming release. "
                      "Use assign_partial_charges(partial_charge_method='am1bcc') instead.",
                      DeprecationWarning)

        self.assign_partial_charges(molecule,
                                    partial_charge_method='am1bcc',
                                    use_conformers=use_conformers,
                                    strict_n_conformers=strict_n_conformers)

    def assign_fractional_bond_orders(self, molecule, bond_order_model=None, use_conformers=None):
        """
        Update and store list of bond orders this molecule. Bond orders are stored on each
        bond, in the `bond.fractional_bond_order` attribute.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        molecule : openforcefield.topology.molecule Molecule
            The molecule to assign wiberg bond orders to
        bond_order_model : str, optional, default=None
            The charge model to use. One of ['am1-wiberg', 'pm3-wiberg']. If None, 'am1-wiberg' will be used.
        use_conformers : iterable of simtk.unit.Quantity(np.array) with shape (n_atoms, 3) and dimension of distance, optional, default=None
            The conformers to use for fractional bond order calculation. If None, an appropriate number
            of conformers will be generated by an available ToolkitWrapper.

         """
        from openeye import oequacpac
        from openforcefield.topology import Molecule
        # Make a copy since we'll be messing with this molecule's conformers
        temp_mol = Molecule(molecule)
        if use_conformers is None:
            temp_mol.generate_conformers(n_conformers=1)
        else:
            temp_mol._conformers = None
            for conformer in use_conformers:
                temp_mol.add_conformer(conformer)


        if temp_mol.n_conformers == 0:
            raise Exception(
                "No conformers present in molecule submitted for fractional bond order calculation. Consider "
                "loading the molecule from a file with geometry already present or running "
                "molecule.generate_conformers() before calling molecule.compute_wiberg_bond_orders()"
            )

        if bond_order_model is None:
            bond_order_model = 'am1-wiberg'

        # Based on example at https://docs.eyesopen.com/toolkits/python/quacpactk/examples_summary_wibergbondorders.html
        oemol = self.to_openeye(temp_mol)
        am1 = oequacpac.OEAM1()
        am1results = oequacpac.OEAM1Results()
        am1options = am1.GetOptions()
        if bond_order_model == "am1-wiberg":
            am1options.SetSemiMethod(oequacpac.OEMethodType_AM1)
        elif bond_order_model == "pm3-wiberg":
            # TODO: Make sure that modifying am1options actually works
            am1options.SetSemiMethod(oequacpac.OEMethodType_PM3)
        else:
            raise ValueError(f"Bond order model '{bond_order_model}' is not supported by OpenEyeToolkitWrapper. "
                             f"Supported models are ['am1-wiberg', 'pm3-wiberg']")

        #for conf in oemol.GetConfs():
        #TODO: How to handle multiple confs here?
        status = am1.CalcAM1(am1results, oemol)

        if status is False:
            raise Exception(
                'Unable to assign charges (in the process of calculating fractional bond orders)'
            )

        # TODO: Will bonds always map back to the same index? Consider doing a topology mapping.
        # Loop over bonds
        for idx, bond in enumerate(oemol.GetBonds()):
            # Get bond order
            order = am1results.GetBondOrder(bond.GetBgnIdx(), bond.GetEndIdx())
            mol_bond = molecule._bonds[idx]
            mol_bond.fractional_bond_order = order

    @staticmethod
    def _find_smarts_matches(oemol, smarts, aromaticity_model=None):
        """Find all sets of atoms in the provided OpenEye molecule that match the provided SMARTS string.

        Parameters
        ----------
        oemol : openeye.oechem.OEMol or similar
            oemol to process with the SMIRKS in order to find matches
        smarts : str
            SMARTS string with any number of sequentially tagged atoms.
            If there are N tagged atoms numbered 1..N, the resulting matches will be N-tuples of atoms that match the corresponding tagged atoms.
        aromaticity_model : str, optional, default=None
            OpenEye aromaticity model designation as a string, such as ``OEAroModel_MDL``.
            If ``None``, molecule is processed exactly as provided; otherwise it is prepared with this aromaticity model prior to querying.

        Returns
        -------
        matches : list of tuples of atoms indices within the ``oemol``
            matches[index] is an N-tuple of atom numbers from the ``oemol``
            Matches are returned in no guaranteed order.
            # TODO: What is returned if no matches are found? An empty list, or None?
            # TODO: Ensure that SMARTS numbers 1, 2, 3... are rendered into order of returnd matches indexed by 0, 1, 2...

        .. notes ::

           * Raises ``LicenseError`` if valid OpenEye tools license is not found, rather than causing program to terminate
           * Raises ``ValueError`` if ``smarts`` query is malformed

        """
        from openeye import oechem
        from openeye.oechem import OESubSearch
        # Make a copy of molecule so we don't influence original (probably safer than deepcopy per C Bayly)
        mol = oechem.OEMol(oemol)

        # Set up query
        qmol = oechem.OEQMol()
        if not oechem.OEParseSmarts(qmol, smarts):
            raise ValueError(f"Error parsing SMARTS '{smarts}'")

        # Determine aromaticity model
        if aromaticity_model:
            if type(aromaticity_model) == str:
                # Check if the user has provided a manually-specified aromaticity_model
                if hasattr(oechem, aromaticity_model):
                    oearomodel = getattr(oechem,
                                         'OEAroModel_' + aromaticity_model)
                else:
                    raise ValueError(
                        "Error: provided aromaticity model not recognized by oechem."
                    )
            else:
                raise ValueError(
                    "Error: provided aromaticity model must be a string.")

            # If aromaticity model was provided, prepare molecule
            oechem.OEClearAromaticFlags(mol)
            oechem.OEAssignAromaticFlags(mol, oearomodel)
            # Avoid running OEPrepareSearch or we lose desired aromaticity, so instead:
            oechem.OEAssignHybridization(mol)

        # Build list of matches
        # TODO: The MoleculeImage mapping should preserve ordering of template molecule for equivalent atoms
        #       and speed matching for larger molecules.
        unique = False  # We require all matches, not just one of each kind
        substructure_search = OESubSearch(qmol)
        substructure_search.SetMaxMatches(0)
        matches = list()
        for match in substructure_search.Match(mol, unique):
            # Compile list of atom indices that match the pattern tags
            atom_indices = dict()
            for matched_atom in match.GetAtoms():
                if matched_atom.pattern.GetMapIdx() != 0:
                    atom_indices[matched_atom.pattern.GetMapIdx() -
                                 1] = matched_atom.target.GetIdx()
            # Compress into list
            atom_indices = [
                atom_indices[index] for index in range(len(atom_indices))
            ]
            # Convert to tuple
            matches.append(tuple(atom_indices))
        return matches

    def find_smarts_matches(self,
                            molecule,
                            smarts,
                            aromaticity_model='OEAroModel_MDL'):
        """
        Find all SMARTS matches for the specified molecule, using the specified aromaticity model.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        molecule : openforcefield.topology.Molecule
            The molecule for which all specified SMARTS matches are to be located
        smarts : str
            SMARTS string with optional SMIRKS-style atom tagging
        aromaticity_model : str, optional, default='OEAroModel_MDL'
            Aromaticity model to use during matching

        .. note :: Currently, the only supported ``aromaticity_model`` is ``OEAroModel_MDL``

        """
        oemol = self.to_openeye(molecule)
        return self._find_smarts_matches(oemol, smarts)


class RDKitToolkitWrapper(ToolkitWrapper):
    """
    RDKit toolkit wrapper

    .. warning :: This API is experimental and subject to change.
    """

    _toolkit_name = 'The RDKit'
    _toolkit_installation_instructions = 'A conda-installable version of the free and open source RDKit cheminformatics ' \
                                         'toolkit can be found at: https://anaconda.org/rdkit/rdkit'

    def __init__(self):
        super().__init__()

        self._toolkit_file_read_formats = ['SDF', 'MOL', 'SMI']  # TODO: Add TDT support

        if not self.is_available():
            raise ToolkitUnavailableException(f'The required toolkit {self._toolkit_name} is not '
                                              f'available. {self._toolkit_installation_instructions}')
        else:
            from rdkit import Chem
            # we have to make sure the toolkit can be loaded before formatting this dict
            # Note any new file write formats should be added here only
            self._toolkit_file_write_formats = {'SDF': Chem.SDWriter, 'MOL': Chem.SDWriter, 'SMI': Chem.SmilesWriter,
                                                'PDB': Chem.PDBWriter, 'TDT': Chem.TDTWriter}

    @property
    def toolkit_file_write_formats(self):
        """
        List of file formats that this toolkit can write.
        """
        return list(self._toolkit_file_write_formats.keys())

    @staticmethod
    def is_available():
        """
        Check whether the RDKit toolkit can be imported

        Returns
        -------
        is_installed : bool
            True if RDKit is installed, False otherwise.

        """
        try:
            importlib.import_module('rdkit', 'Chem')
            return True
        except ImportError:
            return False

    def from_object(self, object, allow_undefined_stereo=False):
        """
        If given an rdchem.Mol (or rdchem.Mol-derived object), this function will load it into an
        openforcefield.topology.molecule. Otherwise, it will return False.

        Parameters
        ----------
        object : A rdchem.Mol-derived object
            An object to be type-checked and converted into a Molecule, if possible.
        allow_undefined_stereo : bool, default=False
            Whether to accept molecules with undefined stereocenters. If False,
            an exception will be raised if a molecule with undefined stereochemistry
            is passed into this function.

        Returns
        -------
        Molecule or False
            An openforcefield.topology.molecule Molecule.

        Raises
        ------
        NotImplementedError
            If the object could not be converted into a Molecule.
        """
        # TODO: Add tests for the from_object functions
        from rdkit import Chem
        if isinstance(object, Chem.rdchem.Mol):
            return self.from_rdkit(object,
                                   allow_undefined_stereo=allow_undefined_stereo)
        raise TypeError('Cannot create Molecule from {} object'.format(type(object)))

    def from_pdb_and_smiles(self, file_path, smiles, allow_undefined_stereo=False):
        """
        Create a Molecule from a pdb file and a SMILES string using RDKit.

        Requires RDKit to be installed.

        The molecule is created and sanitised based on the SMILES string, we then find a mapping
        between this molecule and one from the PDB based only on atomic number and connections.
        The SMILES molecule is then reindex to match the PDB, the conformer is attached and the
        molecule returned.

        Parameters
        ----------
        file_path: str
            PDB file path
        smiles : str
            a valid smiles string for the pdb, used for seterochemistry and bond order

        allow_undefined_stereo : bool, default=False
            If false, raises an exception if oemol contains undefined stereochemistry.

        Returns
        --------
        molecule : openforcefield.Molecule
            An OFFMol instance with ordering the same as used in the PDB file.

        Raises
        ------
        InvalidConformerError : if the SMILES and PDB molecules are not isomorphic.
        """

        from rdkit import Chem
        from openforcefield.topology.molecule import Molecule, InvalidConformerError

        # Make the molecule from smiles
        offmol = self.from_smiles(smiles,
                                  allow_undefined_stereo=allow_undefined_stereo)

        # Make another molecule from the PDB, allow stero errors here they are expected
        pdbmol = self.from_rdkit(Chem.MolFromPDBFile(file_path, removeHs=False), allow_undefined_stereo=True)

        # check isomorphic and get the mapping if true the mapping will be
        # Dict[pdb_index: offmol_index] sorted by pdb_index
        isomorphic, mapping = Molecule.are_isomorphic(pdbmol, offmol, return_atom_map=True,
                                                      aromatic_matching=False,
                                                      formal_charge_matching=False,
                                                      bond_order_matching=False,
                                                      atom_stereochemistry_matching=False,
                                                      bond_stereochemistry_matching=False)

        if mapping is not None:
            new_mol = offmol.remap(mapping)

            # the pdb conformer is in the correct order so just attach it here
            new_mol.add_conformer(pdbmol.conformers[0])

            return new_mol

        else:
            raise InvalidConformerError('The PDB and SMILES structures do not match.')

    def from_file(self,
                  file_path,
                  file_format,
                  allow_undefined_stereo=False):
        """
        Create an openforcefield.topology.Molecule from a file using this toolkit.



        Parameters
        ----------
        file_path : str
            The file to read the molecule from
        file_format : str
            Format specifier, usually file suffix (eg. 'MOL2', 'SMI')
            Note that not all toolkits support all formats. Check ToolkitWrapper.toolkit_file_read_formats for details.
        allow_undefined_stereo : bool, default=False
            If false, raises an exception if oemol contains undefined stereochemistry.

        Returns
        -------
        molecules : iterable of Molecules
            a list of Molecule objects is returned.

        """
        from rdkit import Chem
        mols = list()
        if (file_format == 'MOL') or (file_format == 'SDF'):
            for rdmol in Chem.SupplierFromFilename(file_path, removeHs=False, sanitize=False, strictParsing=True):
                if rdmol is None:
                    continue

                # Sanitize the molecules (fails on nitro groups)
                try:
                    Chem.SanitizeMol(rdmol, Chem.SANITIZE_ALL ^ Chem.SANITIZE_SETAROMATICITY ^ Chem.SANITIZE_ADJUSTHS)
                    Chem.AssignStereochemistryFrom3D(rdmol)
                except ValueError as e:
                    print(rdmol.GetProp('_Name'), e)
                    continue
                Chem.SetAromaticity(rdmol, Chem.AromaticityModel.AROMATICITY_MDL)
                mol = self.from_rdkit(rdmol, allow_undefined_stereo=allow_undefined_stereo)
                mols.append(mol)

        elif (file_format == 'SMI'):
            # TODO: We have to do some special stuff when we import SMILES (currently
            # just adding H's, but could get fancier in the future). It might be
            # worthwhile to parse the SMILES file ourselves and pass each SMILES
            # through the from_smiles function instead
            for rdmol in Chem.SmilesMolSupplier(file_path, titleLine=False):
                rdmol = Chem.AddHs(rdmol)
                mol = self.from_rdkit(rdmol, allow_undefined_stereo=allow_undefined_stereo)
                mols.append(mol)

        elif (file_format == 'PDB'):
            raise Exception(
                "RDKit can not safely read PDBs on their own. Information about bond order and aromaticity "
                "is likely to be lost. To read a PDB using RDKit use Molecule.from_pdb_and_smiles()")
            # TODO: See if we can implement PDB+mol/smi combinations to get complete bond information.
            #  testing to see if we can make a molecule from smiles and then use the PDB conformer as the geometry
            #  and just reorder the molecule
            # https://github.com/openforcefield/openforcefield/issues/121
            # rdmol = Chem.MolFromPDBFile(file_path, removeHs=False)
            # mol = Molecule.from_rdkit(rdmol)
            # mols.append(mol)
            # TODO: Add SMI, TDT(?) support

        return mols

    def from_file_obj(self,
                      file_obj,
                      file_format,
                      allow_undefined_stereo=False):
        """
        Return an openforcefield.topology.Molecule from a file-like object (an object with a ".read()" method using
        this toolkit.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        file_obj : file-like object
            The file-like object to read the molecule from
        file_format : str
            Format specifier, usually file suffix (eg. 'MOL2', 'SMI')
            Note that not all toolkits support all formats. Check ToolkitWrapper.toolkit_file_read_formats for details.
        allow_undefined_stereo : bool, default=False
            If false, raises an exception if oemol contains undefined stereochemistry.

        Returns
        -------
        molecules : Molecule or list of Molecules
            a list of Molecule objects is returned.

        """
        from rdkit import Chem

        mols = []

        if (file_format == "MOL") or (file_format == "SDF"):
            # TODO: Iterate over all mols in file_data
            for rdmol in Chem.ForwardSDMolSupplier(file_obj):
                mol = self.from_rdkit(rdmol)
                mols.append(mol)

        if (file_format == 'SMI'):
            # TODO: Find a cleaner way to parse SMILES lines
            file_data = file_obj.read()
            lines = [line.strip() for line in file_data.split('\n')]
            # remove blank lines
            lines.remove('')
            for line in lines:
                mol = self.from_smiles(line)
                mols.append(mol)

        elif file_format == 'PDB':
            raise Exception(
                "RDKit can not safely read PDBs on their own. Information about bond order and aromaticity "
                "is likely to be lost. To read a PDB using RDKit use Molecule.from_pdb_and_smiles()")
            # TODO: See if we can implement PDB+mol/smi combinations to get complete bond information.
            # https://github.com/openforcefield/openforcefield/issues/121
            # file_data = file_obj.read()
            # rdmol = Chem.MolFromPDBBlock(file_data)
            # mol = Molecule.from_rdkit(rdmol)
            # mols.append(mol)
        # TODO: TDT file support
        return mols

    def to_file_obj(self, molecule, file_obj, file_format):
        """
        Writes an OpenFF Molecule to a file-like object

        Parameters
        ----------
        molecule : an OpenFF Molecule
            The molecule to write
        file_obj
            The file-like object to write to
        file_format
            The format for writing the molecule data

        Returns
        -------

        """

        file_format = file_format.upper()
        rdmol = self.to_rdkit(molecule)
        try:
            writer = self._toolkit_file_write_formats[file_format](file_obj)
            writer.write(rdmol)
            writer.close()
        # if we can not write to that file type catch the error here
        except KeyError:
            raise ValueError(f'The requested file type ({file_format}) is not supported to be written using '
                             f'RDKitToolkitWrapper.')

    def to_file(self, molecule, file_path, file_format):
        """
        Writes an OpenFF Molecule to a file-like object

        Parameters
        ----------
        molecule : an OpenFF Molecule
            The molecule to write
        file_path
            The file path to write to
        file_format
            The format for writing the molecule data

        Returns
        ------

        """

        # open a file object and pass to the object writer
        with open(file_path, 'w') as file_obj:
            self.to_file_obj(molecule=molecule, file_obj=file_obj, file_format=file_format)

    def canonical_order_atoms(self, molecule):
        """
        Canonical order the atoms in the molecule using the RDKit.

        Parameters
        ----------
        molecule: openforcefield.topology.Molecule
            The input molecule

         Returns
        -------
        molecule : openforcefield.topology.Molecule
            The input molecule, with canonically-indexed atoms and bonds.
        """

        from rdkit import Chem
        rdmol = self.to_rdkit(molecule)

        # get the canonical ordering with hydrogens first
        # this is the default behaviour of RDKit
        atom_order = list(Chem.CanonicalRankAtoms(rdmol, breakTies=True))

        heavy_atoms = rdmol.GetNumHeavyAtoms()
        hydrogens = rdmol.GetNumAtoms() - heavy_atoms

        # now go through and change the rankings to get the heavy atoms first if hydrogens are present
        if hydrogens != 0:
            for i in range(len(atom_order)):
                if rdmol.GetAtomWithIdx(i).GetAtomicNum() != 1:
                    atom_order[i] -= hydrogens
                else:
                    atom_order[i] += heavy_atoms

        # make an atom mapping from the atom_order and remap the molecule
        atom_mapping = dict((i, rank) for i, rank in enumerate(atom_order))

        return molecule.remap(atom_mapping, current_to_new=True)

    @classmethod
    def to_smiles(cls, molecule):
        """
        Uses the RDKit toolkit to convert a Molecule into a SMILES string.

        Parameters
        ----------
        molecule : An openforcefield.topology.Molecule
            The molecule to convert into a SMILES.

        Returns
        -------
        smiles : str
            The SMILES of the input molecule.
        """
        from rdkit import Chem
        rdmol = cls.to_rdkit(molecule)
        return Chem.MolToSmiles(rdmol, isomericSmiles=True, allHsExplicit=True)

    def from_smiles(self, smiles, hydrogens_are_explicit=False, allow_undefined_stereo=False):
        """
        Create a Molecule from a SMILES string using the RDKit toolkit.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        smiles : str
            The SMILES string to turn into a molecule
        hydrogens_are_explicit : bool, default=False
            If False, RDKit will perform hydrogen addition using Chem.AddHs
        allow_undefined_stereo : bool, default=False
            Whether to accept SMILES with undefined stereochemistry. If False,
            an exception will be raised if a SMILES with undefined stereochemistry
            is passed into this function.

        Returns
        -------
        molecule : openforcefield.topology.Molecule
            An openforcefield-style molecule.
        """
        from rdkit import Chem

        rdmol = Chem.MolFromSmiles(smiles, sanitize=False)
        # strip the atom map from the molecule if it has one
        # so we don't affect the sterochemistry tags
        for atom in rdmol.GetAtoms():
            if atom.GetAtomMapNum() != 0:
                # set the map back to zero but hide the index in the atom prop data
                atom.SetProp('_map_idx', str(atom.GetAtomMapNum()))
                # set it back to zero
                atom.SetAtomMapNum(0)

        # TODO: I think UpdatePropertyCache(strict=True) is called anyway in Chem.SanitizeMol().
        rdmol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(rdmol, Chem.SANITIZE_ALL ^ Chem.SANITIZE_ADJUSTHS ^ Chem.SANITIZE_SETAROMATICITY)
        Chem.SetAromaticity(rdmol, Chem.AromaticityModel.AROMATICITY_MDL)

        # Chem.MolFromSmiles adds bond directions (i.e. ENDDOWNRIGHT/ENDUPRIGHT), but
        # doesn't set bond.GetStereo(). We need to call AssignStereochemistry for that.
        Chem.AssignStereochemistry(rdmol)

        # Throw an exception/warning if there is unspecified stereochemistry.
        if allow_undefined_stereo==False:
            self._detect_undefined_stereo(rdmol, err_msg_prefix='Unable to make OFFMol from SMILES: ')

        # Add explicit hydrogens if they aren't there already
        if not hydrogens_are_explicit:
            rdmol = Chem.AddHs(rdmol)
        elif hydrogens_are_explicit:
            for atom_idx in range(rdmol.GetNumAtoms()):
                atom = rdmol.GetAtomWithIdx(atom_idx)
                if atom.GetNumImplicitHs() != 0:
                    raise ValueError(
                        f"'hydrogens_are_explicit' was specified as True, but RDKit toolkit interpreted "
                        f"SMILES '{smiles}' as having implicit hydrogen. If this SMILES is intended to "
                        f"express all explicit hydrogens in the molecule, then you should construct the "
                        f"desired molecule as an RDMol with no implicit hydrogens, and then use "
                        f"Molecule.from_rdkit() to create the desired OFFMol.")

        molecule = self.from_rdkit(rdmol,
                                   allow_undefined_stereo=allow_undefined_stereo)

        return molecule

    def from_inchi(self, inchi, allow_undefined_stereo=False):
        """
        Construct a Molecule from a InChI representation

        Parameters
        ----------
        inchi : str
            The InChI representation of the molecule.

        allow_undefined_stereo : bool, default=False
            Whether to accept InChI with undefined stereochemistry. If False,
            an exception will be raised if a InChI with undefined stereochemistry
            is passed into this function.

        Returns
        -------
        molecule : openforcefield.topology.Molecule
        """

        from rdkit import Chem
        # this seems to always remove the hydrogens
        rdmol = Chem.MolFromInchi(inchi, sanitize=False, removeHs=False)

        # try and catch an InChI parsing error
        if rdmol is None:
            raise RuntimeError('There was an issue parsing the InChI string, please check and try again.')

        # process the molecule
        # TODO do we need this with inchi?
        rdmol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(rdmol, Chem.SANITIZE_ALL ^ Chem.SANITIZE_ADJUSTHS ^ Chem.SANITIZE_SETAROMATICITY)
        Chem.SetAromaticity(rdmol, Chem.AromaticityModel.AROMATICITY_MDL)

        # add hydrogens back here
        rdmol = Chem.AddHs(rdmol)

        molecule = self.from_rdkit(rdmol, allow_undefined_stereo=allow_undefined_stereo)

        return molecule

    def generate_conformers(self, molecule, n_conformers=1, rms_cutoff=None, clear_existing=True):
        """
        Generate molecule conformers using RDKit.

        .. warning :: This API is experimental and subject to change.

        .. todo ::

           * which parameters should we expose? (or can we implement a general system with **kwargs?)
           * will the coordinates be returned in the OpenFF Molecule's own indexing system? Or is there a chance that they'll get reindexed when we convert the input into an RDMol?

        Parameters
        ---------
        molecule : a :class:`Molecule`
            The molecule to generate conformers for.
        n_conformers : int, default=1
            Maximum number of conformers to generate.
        rms_cutoff : simtk.Quantity-wrapped float, in units of distance, optional, default=None
            The minimum RMS value at which two conformers are considered redundant and one is deleted. Precise
            implementation of this cutoff may be toolkit-dependent.
        clear_existing : bool, default=True
            Whether to overwrite existing conformers for the molecule.


        """
        from rdkit.Chem import AllChem
        if rms_cutoff is None:
            rms_cutoff = 1. * unit.angstrom
        rdmol = self.to_rdkit(molecule)
        # TODO: This generates way more conformations than omega, given the same nConfs and RMS threshold. Is there some way to set an energy cutoff as well?
        AllChem.EmbedMultipleConfs(
            rdmol,
            numConfs=n_conformers,
            pruneRmsThresh=rms_cutoff / unit.angstrom,
            randomSeed=1,
            #params=AllChem.ETKDG()
        )
        molecule2 = self.from_rdkit(rdmol, allow_undefined_stereo=True)

        if clear_existing:
            molecule._conformers = list()

        for conformer in molecule2._conformers:
            molecule._add_conformer(conformer)

    def from_rdkit(self, rdmol, allow_undefined_stereo=False):
        """
        Create a Molecule from an RDKit molecule.

        Requires the RDKit to be installed.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        rdmol : rkit.RDMol
            An RDKit molecule
        allow_undefined_stereo : bool, default=False
            If false, raises an exception if rdmol contains undefined stereochemistry.

        Returns
        -------
        molecule : openforcefield.Molecule
            An openforcefield molecule

        Examples
        --------

        Create a molecule from an RDKit molecule

        >>> from rdkit import Chem
        >>> from openforcefield.tests.utils import get_data_file_path
        >>> rdmol = Chem.MolFromMolFile(get_data_file_path('systems/monomers/ethanol.sdf'))

        >>> toolkit_wrapper = RDKitToolkitWrapper()
        >>> molecule = toolkit_wrapper.from_rdkit(rdmol)

        """
        from rdkit import Chem
        from openforcefield.topology.molecule import Molecule

        # Make a copy of the RDKit Mol as we'll need to change it (e.g. assign stereo).
        rdmol = Chem.Mol(rdmol)

        # Sanitizing the molecule. We handle aromaticity and chirality manually.
        # This SanitizeMol(...) calls cleanUp, updatePropertyCache, symmetrizeSSSR,
        # assignRadicals, setConjugation, and setHybridization.
        Chem.SanitizeMol(rdmol, (Chem.SANITIZE_ALL ^ Chem.SANITIZE_SETAROMATICITY ^
                                 Chem.SANITIZE_ADJUSTHS ^ Chem.SANITIZE_CLEANUPCHIRALITY ^
                                 Chem.SANITIZE_KEKULIZE))
        Chem.SetAromaticity(rdmol, Chem.AromaticityModel.AROMATICITY_MDL)
        # SetAromaticity set aromatic bonds to 1.5, but Molecule.bond_order is an
        # integer (contrarily to fractional_bond_order) so we need the Kekule order.
        Chem.Kekulize(rdmol)

        # Make sure the bond stereo tags are set before checking for
        # undefined stereo. RDKit can figure out bond stereo from other
        # information in the Mol object like bond direction properties.
        # Do not overwrite eventual chiral tags provided by the user.
        Chem.AssignStereochemistry(rdmol, cleanIt=False)

        # Check for undefined stereochemistry.
        self._detect_undefined_stereo(rdmol, raise_warning=allow_undefined_stereo,
                                      err_msg_prefix="Unable to make OFFMol from RDMol: ")

        # Create a new openforcefield Molecule
        offmol = Molecule()

        # If RDMol has a title save it
        if rdmol.HasProp("_Name"):
            #raise Exception('{}'.format(rdmol.GetProp('name')))
            offmol.name = rdmol.GetProp("_Name")
        else:
            offmol.name = ""

        # Store all properties
        # TODO: Should there be an API point for storing properties?
        properties = rdmol.GetPropsAsDict()
        offmol._properties = properties

        # setting chirality in openeye requires using neighbor atoms
        # therefore we can't do it until after the atoms and bonds are all added
        map_atoms = {}
        map_bonds = {}
        # if we are loading from a mapped smiles extract the mapping
        atom_mapping = {}
        for rda in rdmol.GetAtoms():
            rd_idx = rda.GetIdx()
            # if the molecule was made from a mapped smiles this has been hidden
            # so that it does not affect the sterochemistry tags
            try:
                map_id = int(rda.GetProp('_map_idx'))
            except KeyError:
                map_id = rda.GetAtomMapNum()

            # create a new atom
            #atomic_number = oemol.NewAtom(rda.GetAtomicNum())
            atomic_number = rda.GetAtomicNum()
            formal_charge = rda.GetFormalCharge() * unit.elementary_charge
            is_aromatic = rda.GetIsAromatic()
            if rda.HasProp('_Name'):
                name = rda.GetProp('_Name')
            else:
                # check for PDB names
                try:
                    name = rda.GetMonomerInfo().GetName().strip()
                except AttributeError:
                    name = ''

            # If chiral, store the chirality to be set later
            stereochemistry = None
            #tag = rda.GetChiralTag()
            if rda.HasProp('_CIPCode'):
                stereo_code = rda.GetProp('_CIPCode')
                #if tag == Chem.CHI_TETRAHEDRAL_CCW:
                if stereo_code == 'R':
                    stereochemistry = 'R'
                #if tag == Chem.CHI_TETRAHEDRAL_CW:
                elif stereo_code == 'S':
                    stereochemistry = 'S'
                else:
                    raise UndefinedStereochemistryError("In from_rdkit: Expected atom stereochemistry of R or S. "
                                                        "Got {} instead.".format(stereo_code))

            atom_index = offmol.add_atom(
                atomic_number,
                formal_charge,
                is_aromatic,
                name=name,
                stereochemistry=stereochemistry)
            map_atoms[rd_idx] = atom_index
            atom_mapping[atom_index] = map_id

        # if we have a full atom map add it to the molecule, 0 indicates a missing mapping or no mapping
        if 0 not in atom_mapping.values():
            offmol._properties['atom_map'] = atom_mapping

        # Similar to chirality, stereochemistry of bonds in OE is set relative to their neighbors
        for rdb in rdmol.GetBonds():
            rdb_idx = rdb.GetIdx()
            a1 = rdb.GetBeginAtomIdx()
            a2 = rdb.GetEndAtomIdx()

            # Determine bond aromaticity and Kekulized bond order
            is_aromatic = rdb.GetIsAromatic()
            order = rdb.GetBondTypeAsDouble()
            # Convert floating-point bond order to integral bond order
            order = int(order)

            # create a new bond
            bond_index = offmol.add_bond(map_atoms[a1], map_atoms[a2], order,
                                         is_aromatic)
            map_bonds[rdb_idx] = bond_index

        # Now fill in the cached (structure-dependent) properties. We have to have the 2D structure of the molecule
        # in place first, because each call to add_atom and add_bond invalidates all cached properties
        for rdb in rdmol.GetBonds():
            rdb_idx = rdb.GetIdx()
            offb_idx = map_bonds[rdb_idx]
            offb = offmol.bonds[offb_idx]
            # determine if stereochemistry is needed
            stereochemistry = None
            tag = rdb.GetStereo()
            if tag == Chem.BondStereo.STEREOZ:
                stereochemistry = 'Z'
            elif tag == Chem.BondStereo.STEREOE:
                stereochemistry = 'E'
            elif tag == Chem.BondStereo.STEREOTRANS or tag == Chem.BondStereo.STEREOCIS:
                raise ValueError(
                    "Expected RDKit bond stereochemistry of E or Z, got {} instead"
                    .format(tag))
            offb._stereochemistry = stereochemistry
            fractional_bond_order = None
            if rdb.HasProp("fractional_bond_order"):
                fractional_bond_order = rdb.GetDoubleProp(
                    'fractional_bond_order')
            offb.fractional_bond_order = fractional_bond_order

        # TODO: Save conformer(s), if present
        # If the rdmol has a conformer, store its coordinates
        if len(rdmol.GetConformers()) != 0:
            for conf in rdmol.GetConformers():
                n_atoms = offmol.n_atoms
                # TODO: Will this always be angstrom when loading from RDKit?
                positions = unit.Quantity(
                    np.zeros((n_atoms, 3)), unit.angstrom)
                for rd_idx, off_idx in map_atoms.items():
                    atom_coords = conf.GetPositions()[rd_idx, :] * unit.angstrom
                    positions[off_idx, :] = atom_coords
                offmol.add_conformer(positions)

        partial_charges = unit.Quantity(
            np.zeros(offmol.n_atoms, dtype=np.float), unit=unit.elementary_charge)

        any_atom_has_partial_charge = False
        for rd_idx, rd_atom in enumerate(rdmol.GetAtoms()):
            off_idx = map_atoms[rd_idx]
            if rd_atom.HasProp("partial_charge"):
                charge = rd_atom.GetDoubleProp(
                    "partial_charge") * unit.elementary_charge
                partial_charges[off_idx] = charge
                any_atom_has_partial_charge = True
            else:
                # If some other atoms had partial charges but this one doesn't, raise an Exception
                if any_atom_has_partial_charge:
                    raise ValueError(
                        "Some atoms in rdmol have partial charges, but others do not."
                    )

            offmol.partial_charges = partial_charges
        return offmol

    @classmethod
    def to_rdkit(cls, molecule, aromaticity_model=DEFAULT_AROMATICITY_MODEL):
        """
        Create an RDKit molecule

        Requires the RDKit to be installed.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        aromaticity_model : str, optional, default=DEFAULT_AROMATICITY_MODEL
            The aromaticity model to use

        Returns
        -------
        rdmol : rkit.RDMol
            An RDKit molecule

        Examples
        --------

        Convert a molecule to RDKit

        >>> from openforcefield.topology import Molecule
        >>> ethanol = Molecule.from_smiles('CCO')
        >>> rdmol = ethanol.to_rdkit()

        """
        from rdkit import Chem, Geometry

        # Create an editable RDKit molecule
        rdmol = Chem.RWMol()

        # Set name
        # TODO: What is the best practice for how this should be named?
        if not (molecule.name is None):
            rdmol.SetProp('_Name', molecule.name)

        # TODO: Set other properties
        for name, value in molecule.properties.items():
            if type(value) == str:
                rdmol.SetProp(name, value)
            elif type(value) == int:
                rdmol.SetIntProp(name, value)
            elif type(value) == float:
                rdmol.SetDoubleProp(name, value)
            elif type(value) == bool:
                rdmol.SetBoolProp(name, value)
            else:
                # Shove everything else into a string
                rdmol.SetProp(name, str(value))

        _bondtypes = {
            1: Chem.BondType.SINGLE,
            1.5: Chem.BondType.AROMATIC,
            2: Chem.BondType.DOUBLE,
            3: Chem.BondType.TRIPLE,
            4: Chem.BondType.QUADRUPLE,
            5: Chem.BondType.QUINTUPLE,
            6: Chem.BondType.HEXTUPLE,
            7: Chem.BondType.ONEANDAHALF,
        }

        for index, atom in enumerate(molecule.atoms):
            rdatom = Chem.Atom(atom.atomic_number)
            rdatom.SetFormalCharge(atom.formal_charge.value_in_unit(unit.elementary_charge))
            rdatom.SetIsAromatic(atom.is_aromatic)
            rdatom.SetProp('_Name', atom.name)

            ## Stereo handling code moved to after bonds are added
            if atom.stereochemistry == 'S':
                rdatom.SetChiralTag(Chem.CHI_TETRAHEDRAL_CW)
            elif atom.stereochemistry == 'R':
                rdatom.SetChiralTag(Chem.CHI_TETRAHEDRAL_CCW)

            rd_index = rdmol.AddAtom(rdatom)

            # Let's make sure al the atom indices in the two molecules
            # are the same, otherwise we need to create an atom map.
            assert index == atom.molecule_atom_index
            assert index == rd_index

        for bond in molecule.bonds:
            atom_indices = (bond.atom1.molecule_atom_index, bond.atom2.molecule_atom_index)
            rdmol.AddBond(*atom_indices)
            rdbond = rdmol.GetBondBetweenAtoms(*atom_indices)
            if not (bond.fractional_bond_order is None):
                rdbond.SetDoubleProp("fractional_bond_order",
                                     bond.fractional_bond_order)
            # Assign bond type, which is based on order unless it is aromatic
            if bond.is_aromatic:
                rdbond.SetBondType(_bondtypes[1.5])
                rdbond.SetIsAromatic(True)
            else:
                rdbond.SetBondType(_bondtypes[bond.bond_order])
                rdbond.SetIsAromatic(False)

        Chem.SanitizeMol(rdmol, Chem.SANITIZE_ALL ^ Chem.SANITIZE_ADJUSTHS ^ Chem.SANITIZE_SETAROMATICITY)

        # Fix for aromaticity being lost
        if aromaticity_model == 'OEAroModel_MDL':
            Chem.SetAromaticity(rdmol, Chem.AromaticityModel.AROMATICITY_MDL)
        else:
            raise ValueError(f"Aromaticity model {aromaticity_model} not recognized")

        # Assign atom stereochemsitry and collect atoms for which RDKit
        # can't figure out chirality. The _CIPCode property of these atoms
        # will be forcefully set to the stereo we want (see #196).
        undefined_stereo_atoms = {}
        for index, atom in enumerate(molecule.atoms):
            rdatom = rdmol.GetAtomWithIdx(index)

            # Skip non-chiral atoms.
            if atom.stereochemistry is None:
                continue

            # Let's randomly assign this atom's (local) stereo to CW
            # and check if this causes the (global) stereo to be set
            # to the desired one (S or R).
            rdatom.SetChiralTag(Chem.CHI_TETRAHEDRAL_CW)
            # We need to do force and cleanIt to recalculate CIP stereo.
            Chem.AssignStereochemistry(rdmol, force=True, cleanIt=True)
            # If our random initial assignment worked, then we're set.
            if rdatom.HasProp('_CIPCode') and rdatom.GetProp("_CIPCode") == atom.stereochemistry:
                continue

            # Otherwise, set it to CCW.
            rdatom.SetChiralTag(Chem.CHI_TETRAHEDRAL_CCW)
            # We need to do force and cleanIt to recalculate CIP stereo.
            Chem.AssignStereochemistry(rdmol, force=True, cleanIt=True)
            # Hopefully this worked, otherwise something's wrong
            if rdatom.HasProp('_CIPCode') and rdatom.GetProp("_CIPCode") == atom.stereochemistry:
                continue

            # Keep track of undefined stereo atoms. We'll force stereochemistry
            # at the end to avoid the next AssignStereochemistry to overwrite.
            if not rdatom.HasProp('_CIPCode'):
                undefined_stereo_atoms[rdatom] = atom.stereochemistry
                continue

            # Something is wrong.
            err_msg = ("Unknown atom stereochemistry encountered in to_rdkit. "
                       "Desired stereochemistry: {}. Set stereochemistry {}".format(
                atom.stereochemistry, rdatom.GetProp("_CIPCode")))
            raise RuntimeError(err_msg)

        # Copy bond stereo info from molecule to rdmol.
        cls._assign_rdmol_bonds_stereo(molecule, rdmol)

        # Set coordinates if we have them
        if molecule._conformers:
            for conformer in molecule._conformers:
                rdmol_conformer = Chem.Conformer()
                for atom_idx in range(molecule.n_atoms):
                    x, y, z = conformer[atom_idx, :].value_in_unit(unit.angstrom)
                    rdmol_conformer.SetAtomPosition(atom_idx,
                                                    Geometry.Point3D(x, y, z))
                rdmol.AddConformer(rdmol_conformer)

        # Retain charges, if present
        if not (molecule._partial_charges is None):

            rdk_indexed_charges = np.zeros((molecule.n_atoms), dtype=np.float)
            for atom_idx, charge in enumerate(molecule._partial_charges):
                charge_unitless = charge.value_in_unit(unit.elementary_charge)
                rdk_indexed_charges[atom_idx] = charge_unitless
            for atom_idx, rdk_atom in enumerate(rdmol.GetAtoms()):
                rdk_atom.SetDoubleProp('partial_charge',
                                       rdk_indexed_charges[atom_idx])

        # Cleanup the rdmol
        rdmol.UpdatePropertyCache(strict=False)
        Chem.GetSSSR(rdmol)

        # Forcefully assign stereo information on the atoms that RDKit
        # can't figure out. This must be done last as calling AssignStereochemistry
        # again will delete these properties (see #196).
        for rdatom, stereochemistry in undefined_stereo_atoms.items():
            rdatom.SetProp('_CIPCode', stereochemistry)

        # Return non-editable version
        return Chem.Mol(rdmol)

    def to_inchi(self, molecule, fixed_hydrogens=False):
        """
        Create an InChI string for the molecule using the RDKit Toolkit.
        InChI is a standardised representation that does not capture tautomers unless specified using the fixed hydrogen
        layer.

        For information on InChi see here https://iupac.org/who-we-are/divisions/division-details/inchi/

        Parameters
        ----------
        molecule : An openforcefield.topology.Molecule
            The molecule to convert into a SMILES.

        fixed_hydrogens: bool, default=False
            If a fixed hydrogen layer should be added to the InChI, if `True` this will produce a non standard specific
            InChI string of the molecule.

        Returns
        --------
        inchi: str
            The InChI string of the molecule.
        """

        from rdkit import Chem
        rdmol = self.to_rdkit(molecule)
        if fixed_hydrogens:
            inchi = Chem.MolToInchi(rdmol, options='-FixedH')
        else:
            inchi = Chem.MolToInchi(rdmol)
        return inchi

    def to_inchikey(self, molecule, fixed_hydrogens=False):
        """
        Create an InChIKey for the molecule using the RDKit Toolkit.
        InChIKey is a standardised representation that does not capture tautomers unless specified using the fixed hydrogen
        layer.

        For information on InChi see here https://iupac.org/who-we-are/divisions/division-details/inchi/

        Parameters
        ----------
        molecule : An openforcefield.topology.Molecule
            The molecule to convert into a SMILES.

        fixed_hydrogens: bool, default=False
            If a fixed hydrogen layer should be added to the InChI, if `True` this will produce a non standard specific
            InChI string of the molecule.

        Returns
        --------
        inchi_key: str
            The InChIKey representation of the molecule.
        """

        from rdkit import Chem
        rdmol = self.to_rdkit(molecule)
        if fixed_hydrogens:
            inchi_key = Chem.MolToInchiKey(rdmol, options='-FixedH')
        else:
            inchi_key = Chem.MolToInchiKey(rdmol)
        return inchi_key

    @staticmethod
    def _find_smarts_matches(rdmol, smirks,
                             aromaticity_model='OEAroModel_MDL'):
        """Find all sets of atoms in the provided RDKit molecule that match the provided SMARTS string.

        Parameters
        ----------
        rdmol : rdkit.Chem.Mol
            rdmol to process with the SMIRKS in order to find matches
        smarts : str
            SMARTS string with any number of sequentially tagged atoms.
            If there are N tagged atoms numbered 1..N, the resulting matches will be N-tuples of atoms that match the corresponding tagged atoms.
        aromaticity_model : str, optional, default='OEAroModel_MDL'
            OpenEye aromaticity model designation as a string, such as ``OEAroModel_MDL``.
            If ``None``, molecule is processed exactly as provided; otherwise it is prepared with this aromaticity model prior to querying.

        Returns
        -------
        matches : list of tuples of atoms indices within the ``rdmol``
            matches[index] is an N-tuple of atom numbers from the ``rdmol``
            Matches are returned in no guaranteed order.
            # TODO: What is returned if no matches are found? An empty list, or None?
            # TODO: Ensure that SMARTS numbers 1, 2, 3... are rendered into order of returnd matches indexed by 0, 1, 2...

        .. notes ::

           * Raises ``ValueError`` if ``smarts`` query is malformed

        """
        from rdkit import Chem

        # Make a copy of the molecule
        rdmol = Chem.Mol(rdmol)
        # Use designated aromaticity model
        if aromaticity_model == 'OEAroModel_MDL':
            Chem.SanitizeMol(rdmol,
                             Chem.SANITIZE_ALL ^ Chem.SANITIZE_SETAROMATICITY)
            Chem.SetAromaticity(rdmol, Chem.AromaticityModel.AROMATICITY_MDL)
        else:
            # Only the OEAroModel_MDL is supported for now
            raise ValueError(
                'Unknown aromaticity model: {}'.aromaticity_models)

        # Set up query.
        qmol = Chem.MolFromSmarts(smirks)  #cannot catch the error
        if qmol is None:
            raise ValueError(
                'RDKit could not parse the SMIRKS string "{}"'.format(smirks))

        # Create atom mapping for query molecule
        idx_map = dict()
        for atom in qmol.GetAtoms():
            smirks_index = atom.GetAtomMapNum()
            if smirks_index != 0:
                idx_map[smirks_index - 1] = atom.GetIdx()
        map_list = [idx_map[x] for x in sorted(idx_map)]

        # Perform matching
        matches = list()

        # choose the largest unsigned int without overflow
        # since the C++ signature is a uint
        max_matches = np.iinfo(np.uintc).max
        for match in rdmol.GetSubstructMatches(qmol, uniquify=False,
                maxMatches=max_matches):
            mas = [match[x] for x in map_list]
            matches.append(tuple(mas))

        return matches

    def find_smarts_matches(self,
                            molecule,
                            smarts,
                            aromaticity_model='OEAroModel_MDL'):
        """
        Find all SMARTS matches for the specified molecule, using the specified aromaticity model.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        molecule : openforcefield.topology.Molecule
            The molecule for which all specified SMARTS matches are to be located
        smarts : str
            SMARTS string with optional SMIRKS-style atom tagging
        aromaticity_model : str, optional, default='OEAroModel_MDL'
            Aromaticity model to use during matching

        .. note :: Currently, the only supported ``aromaticity_model`` is ``OEAroModel_MDL``

        """
        rdmol = self.to_rdkit(molecule, aromaticity_model=aromaticity_model)
        return self._find_smarts_matches(
            rdmol, smarts, aromaticity_model='OEAroModel_MDL')

    # --------------------------------
    # Stereochemistry RDKit utilities.
    # --------------------------------

    @staticmethod
    def _find_undefined_stereo_atoms(rdmol, assign_stereo=False):
        """Find the chiral atoms with undefined stereochemsitry in the RDMol.

        Parameters
        ----------
        rdmol : rdkit.RDMol
            The RDKit molecule.
        assign_stereo : bool, optional, default=False
            As a side effect, this function calls ``Chem.AssignStereochemistry()``
            so by default we work on a molecule copy. Set this to ``True`` to avoid
            making a copy and assigning the stereochemistry to the Mol object.

        Returns
        -------
        undefined_atom_indices : List[int]
            A list of atom indices that are chiral centers with undefined
            stereochemistry.

        See Also
        --------
        rdkit.Chem.FindMolChiralCenters

        """
        from rdkit import Chem

        if not assign_stereo:
            # Avoid modifying the original molecule.
            rdmol = copy.deepcopy(rdmol)

        # Flag possible chiral centers with the "_ChiralityPossible".
        Chem.AssignStereochemistry(rdmol, force=True, flagPossibleStereoCenters=True)

        # Find all atoms with undefined stereo.
        undefined_atom_indices = []
        for atom_idx, atom in enumerate(rdmol.GetAtoms()):
            if (atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED and
                    atom.HasProp('_ChiralityPossible')):
                undefined_atom_indices.append(atom_idx)
        return undefined_atom_indices

    @staticmethod
    def _find_undefined_stereo_bonds(rdmol):
        """Find the chiral atoms with undefined stereochemsitry in the RDMol.

        Parameters
        ----------
        rdmol : rdkit.RDMol
            The RDKit molecule.

        Returns
        -------
        undefined_bond_indices : List[int]
            A list of bond indices with undefined stereochemistry.

        See Also
        --------
        Chem.EnumerateStereoisomers._getFlippers

        Links
        -----
        https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/Chirality.cpp#L1509-L1515
            This comment in FindPotentialStereoBonds mention that the method
            ignores ring bonds.
        https://github.com/DrrDom/rdk/blob/master/gen_stereo_rdkit3.py
            The function get_unspec_double_bonds() in this module looks like
            may solve the problem with the rings.

        """
        from rdkit import Chem

        # Copy the molecule to avoid side effects. Chem.FindPotentialStereoBonds
        # assign Bond.STEREOANY to unspecific bond, which make subsequent calls
        # of Chem.AssignStereochemistry ignore the bond even if there are
        # ENDDOWNRIGHT/ENDUPRIGHT bond direction indications.
        rdmol = copy.deepcopy(rdmol)

        # This function assigns Bond.GetStereo() == Bond.STEREOANY to bonds with
        # undefined stereochemistry.
        Chem.FindPotentialStereoBonds(rdmol)

        undefined_bond_indices = []
        for bond_idx, bond in enumerate(rdmol.GetBonds()):
            if bond.GetStereo() == Chem.BondStereo.STEREOANY:
                undefined_bond_indices.append(bond_idx)
        return undefined_bond_indices

    @classmethod
    def _detect_undefined_stereo(cls, rdmol, err_msg_prefix='', raise_warning=False):
        """Raise UndefinedStereochemistryError if the RDMol has undefined stereochemistry.

        Parameters
        ----------
        rdmol : rdkit.Chem.Mol
            The RDKit molecule.
        err_msg_prefix : str, optional
            A string to prepend to the error/warning message.
        raise_warning : bool, optional, default=False
            If True, a warning is issued instead of an exception.

        Raises
        ------
        UndefinedStereochemistryError
            If the RDMol has undefined atom or bond stereochemistry.

        """
        # Find undefined atom/bond stereochemistry.
        undefined_atom_indices = cls._find_undefined_stereo_atoms(rdmol)
        undefined_bond_indices = cls._find_undefined_stereo_bonds(rdmol)

        # Build error message.
        if len(undefined_atom_indices) == 0 and len(undefined_bond_indices) == 0:
            msg = None
        else:
            msg = err_msg_prefix + "RDMol has unspecified stereochemistry. "
            # The "_Name" property is not always assigned.
            if rdmol.HasProp("_Name"):
                msg += "RDMol name: " + rdmol.GetProp("_Name")

        # Details about undefined atoms.
        if len(undefined_atom_indices) > 0:
            msg += "Undefined chiral centers are:\n"
            for undefined_atom_idx in undefined_atom_indices:
                msg += ' - Atom {symbol} (index {index})\n'.format(
                    symbol=rdmol.GetAtomWithIdx(undefined_atom_idx).GetSymbol(),
                    index=undefined_atom_idx)

        # Details about undefined bond.
        if len(undefined_bond_indices) > 0:
            msg += "Bonds with undefined stereochemistry are:\n"
            for undefined_bond_idx in undefined_bond_indices:
                bond = rdmol.GetBondWithIdx(undefined_bond_idx)
                atom1, atom2 = bond.GetBeginAtom(), bond.GetEndAtom()
                msg += ' - Bond {bindex} (atoms {aindex1}-{aindex2} of element ({symbol1}-{symbol2})\n'.format(
                    bindex=undefined_bond_idx,
                    aindex1=atom1.GetIdx(), aindex2=atom2.GetIdx(),
                    symbol1=atom1.GetSymbol(), symbol2=atom2.GetSymbol())

        if msg is not None:
            if raise_warning:
                msg = 'Warning (not error because allow_undefined_stereo=True): '
                logger.warning(msg)
            else:
                msg = 'Unable to make OFFMol from RDMol: ' + msg
                raise UndefinedStereochemistryError(msg)

    @staticmethod
    def _flip_rdbond_direction(rdbond, paired_rdbonds):
        """Flip the rdbond and all those paired to it.

        Parameters
        ----------
        rdbond : rdkit.Chem.Bond
            The Bond whose direction needs to be flipped.
        paired_rdbonds : Dict[Tuple[int], List[rdkit.Chem.Bond]]
            Maps bond atom indices that are assigned a bond direction to
            the bonds on the other side of the double bond.
        """
        from rdkit import Chem

        # The function assumes that all bonds are either up or down.
        supported_directions = {Chem.BondDir.ENDUPRIGHT, Chem.BondDir.ENDDOWNRIGHT}

        def _flip(b, paired, flipped, ignored):
            # The function assumes that all bonds are either up or down.
            assert b.GetBondDir() in supported_directions
            bond_atom_indices = (b.GetBeginAtomIdx(), b.GetEndAtomIdx())

            # Check that we haven't flipped this bond already.
            if bond_atom_indices in flipped:
                # This should never happen.
                raise RuntimeError('Cannot flip the bond direction consistently.')

            # Flip the bond.
            if b.GetBondDir() == Chem.BondDir.ENDUPRIGHT:
                b.SetBondDir(Chem.BondDir.ENDDOWNRIGHT)
            else:
                b.SetBondDir(Chem.BondDir.ENDUPRIGHT)
            flipped.add(bond_atom_indices)

            # Flip all the paired bonds as well (if there are any).
            if bond_atom_indices in paired:
                for paired_rdbond in paired[bond_atom_indices]:
                    # Don't flip the bond that was flipped in the upper-level recursion.
                    if (paired_rdbond.GetBeginAtomIdx(), paired_rdbond.GetEndAtomIdx()) != ignored:
                        # Don't flip this bond in the next recursion.
                        _flip(paired_rdbond, paired, flipped, ignored=bond_atom_indices)

        _flip(rdbond, paired_rdbonds, flipped=set(), ignored=None)

    @classmethod
    def _assign_rdmol_bonds_stereo(cls, offmol, rdmol):
        """Copy the info about bonds stereochemistry from the OFF Molecule to RDKit Mol."""
        from rdkit import Chem

        # Map the bonds indices that are assigned bond direction
        # to the bond on the other side of the double bond.
        # (atom_index1, atom_index2) -> List[rdkit.Chem.Bond]
        paired_bonds = {}

        for bond in offmol.bonds:
            # No need to do anything with bonds without stereochemistry.
            if not bond.stereochemistry:
                continue

            # Isolate stereo RDKit bond object.
            rdbond_atom_indices = (bond.atom1.molecule_atom_index,
                                   bond.atom2.molecule_atom_index)
            stereo_rdbond = rdmol.GetBondBetweenAtoms(*rdbond_atom_indices)

            # Collect all neighboring rdbonds of atom1 and atom2.
            neighbor_rdbonds1 = [rdmol.GetBondBetweenAtoms(n.molecule_atom_index,
                                                           bond.atom1.molecule_atom_index)
                                 for n in bond.atom1.bonded_atoms if n != bond.atom2]
            neighbor_rdbonds2 = [rdmol.GetBondBetweenAtoms(bond.atom2.molecule_atom_index,
                                                           n.molecule_atom_index)
                                 for n in bond.atom2.bonded_atoms if n != bond.atom1]

            # Select only 1 neighbor bond per atom out of the two.
            neighbor_rdbonds = []
            for i, rdbonds in enumerate([neighbor_rdbonds1, neighbor_rdbonds2]):
                # If there are no neighbors for which we have already
                # assigned the bond direction, just pick the first one.
                neighbor_rdbonds.append(rdbonds[0])
                # Otherwise, pick neighbor that was already assigned to
                # avoid inconsistencies and keep the tree non-cyclic.
                for rdb in rdbonds:
                    if (rdb.GetBeginAtomIdx(), rdb.GetBeginAtomIdx()) in paired_bonds:
                        neighbor_rdbonds[i] = rdb
                        break

            # Assign a random direction to the bonds that were not already assigned
            # keeping track of which bond would be best to flip later (i.e. does that
            # are not already determining the stereochemistry of another double bond).
            flipped_rdbond = neighbor_rdbonds[0]
            for rdb in neighbor_rdbonds:
                if (rdb.GetBeginAtomIdx(), rdb.GetEndAtomIdx()) not in paired_bonds:
                    rdb.SetBondDir(Chem.BondDir.ENDUPRIGHT)
                    # Set this bond as a possible bond to flip.
                    flipped_rdbond = rdb

            Chem.AssignStereochemistry(rdmol, cleanIt=True, force=True)

            # Verify that the current directions give us the desired stereochemistries.
            assert bond.stereochemistry in {'E', 'Z'}
            if bond.stereochemistry == 'E':
                desired_rdk_stereo_code = Chem.rdchem.BondStereo.STEREOE
            else:
                desired_rdk_stereo_code = Chem.rdchem.BondStereo.STEREOZ

            # If that doesn't work, flip the direction of one bond preferring
            # those that are not already determining the stereo of another bond.
            if stereo_rdbond.GetStereo() != desired_rdk_stereo_code:
                cls._flip_rdbond_direction(flipped_rdbond, paired_bonds)
                Chem.AssignStereochemistry(rdmol, cleanIt=True, force=True)

                # The stereo should be set correctly here.
                assert stereo_rdbond.GetStereo() == desired_rdk_stereo_code

            # Update paired bonds map.
            neighbor_bond_indices = [(rdb.GetBeginAtomIdx(), rdb.GetEndAtomIdx()) for rdb in neighbor_rdbonds]
            for i, bond_indices in enumerate(neighbor_bond_indices):
                try:
                    paired_bonds[bond_indices].append(neighbor_rdbonds[1-i])
                except KeyError:
                    paired_bonds[bond_indices] = [neighbor_rdbonds[1-i]]


class AmberToolsToolkitWrapper(ToolkitWrapper):
    """
    AmberTools toolkit wrapper

    .. warning :: This API is experimental and subject to change.
    """
    _toolkit_name = 'AmberTools'
    _toolkit_installation_instructions = 'The AmberTools toolkit (free and open source) can be found at ' \
                                         'https://anaconda.org/omnia/ambertools'

    def __init__(self):
        super().__init__()

        self._toolkit_file_read_formats = []
        self._toolkit_file_write_formats = []

        if not self.is_available():
            raise ToolkitUnavailableException(f'The required toolkit {self._toolkit_name} is not '
                                              f'available. {self._toolkit_installation_instructions}')

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
        else:
            return True

    def __init__(self):
        # TODO: Find AMBERHOME or executable home, checking miniconda if needed
        # Store an instance of an RDKitToolkitWrapper for file I/O
        self._rdkit_toolkit_wrapper = RDKitToolkitWrapper()

    def assign_partial_charges(self,
                               molecule,
                               partial_charge_method=None,
                               use_conformers=None,
                               strict_n_conformers=False):
        """
        Compute partial charges with AmberTools using antechamber/sqm, and assign
        the new values to the partial_charges attribute.

        .. warning :: This API experimental and subject to change.

        .. todo ::

           * Do we want to also allow ESP/RESP charges?

        Parameters
        ----------
        molecule : Molecule
            Molecule for which partial charges are to be computed
        partial_charge_method : str, optional, default=None
            The charge model to use. One of ['gasteiger', 'am1-mulliken']. If None, 'am1-mulliken' will be used.
        use_conformers : iterable of simtk.unit.Quantity-wrapped numpy arrays, each with shape (n_atoms, 3) and dimension of distance. Optional, default = None
            Coordinates to use for partial charge calculation. If None, an appropriate number of conformers will be generated.
        strict_n_conformers : bool, default=False
            Whether to raise an exception if an invalid number of conformers is provided for the given charge method.
            If this is False and an invalid number of conformers is found, a warning will be raised.

        Raises
        ------
        ChargeMethodUnavailableError if the requested charge method can not be handled by this toolkit

        ChargeCalculationError if the charge calculation is supported by this toolkit, but fails
    """

        import os
        import subprocess
        from simtk import unit
        from openforcefield.topology import Molecule

        if partial_charge_method is None:
            partial_charge_method = 'am1-mulliken'
        else:
            # Standardize method name for string comparisons
            partial_charge_method = partial_charge_method.lower()

        SUPPORTED_CHARGE_METHODS = {'am1bcc': {'antechamber_keyword':'bcc',
                                               'min_confs': 1,
                                               'max_confs': 1,
                                               'rec_confs': 1,
                                               },
                                    'am1-mulliken': {'antechamber_keyword': 'mul',
                                                     'min_confs': 1,
                                                     'max_confs': 1,
                                                     'rec_confs': 1,
                                                     },
                                    'gasteiger': {'antechamber_keyword': 'gas',
                                                  'min_confs': 0,
                                                  'max_confs': 0,
                                                  'rec_confs': 0,
                                                  }
                                    }


        if partial_charge_method not in SUPPORTED_CHARGE_METHODS:
            raise ChargeMethodUnavailableError(
                f"partial_charge_method '{partial_charge_method}' is not available from AmberToolsToolkitWrapper. "
                f"Available charge methods are {list(SUPPORTED_CHARGE_METHODS.keys())} "
            )

        charge_method = SUPPORTED_CHARGE_METHODS[partial_charge_method]

        # Make a temporary copy of the molecule, since we'll be messing with its conformers
        mol_copy = Molecule(molecule)

        if use_conformers is None:
            if charge_method['rec_confs'] == 0:
                mol_copy._conformers = None
            else:
                mol_copy.generate_conformers(n_conformers=charge_method['rec_confs'],
                                             rms_cutoff=0.25*unit.angstrom,
                                             toolkit_registry=RDKitToolkitWrapper())
            # TODO: What's a "best practice" RMS cutoff to use here?
        else:
            mol_copy._conformers = None
            for conformer in use_conformers:
                mol_copy.add_conformer(conformer)
            self._check_n_conformers(mol_copy,
                                     partial_charge_method=partial_charge_method,
                                     min_confs=charge_method['min_confs'],
                                     max_confs=charge_method['max_confs'],
                                     strict_n_conformers=strict_n_conformers)

        # Find the path to antechamber
        # TODO: How should we implement find_executable?
        ANTECHAMBER_PATH = find_executable("antechamber")
        if ANTECHAMBER_PATH is None:
            raise IOError("Antechamber not found, cannot run charge_mol()")

        # Compute charges
        from openforcefield.utils import temporary_directory, temporary_cd
        with temporary_directory() as tmpdir:
            with temporary_cd(tmpdir):
                net_charge = mol_copy.total_charge
                # Write out molecule in SDF format
                ## TODO: How should we handle multiple conformers?
                self._rdkit_toolkit_wrapper.to_file(
                    mol_copy, 'molecule.sdf', file_format='sdf')
                #os.system('ls')
                #os.system('cat molecule.sdf')
                # Compute desired charges
                # TODO: Add error handling if antechamber chokes
                # TODO: Add something cleaner than os.system
                short_charge_method = charge_method['antechamber_keyword']
                subprocess.check_output([
                    "antechamber", "-i", "molecule.sdf", "-fi", "sdf", "-o", "charged.mol2", "-fo",
                    "mol2", "-pf", "yes", "-dr", "n", "-c", short_charge_method, "-nc", str(net_charge)]

                )
                #os.system('cat charged.mol2')

                # Write out just charges
                subprocess.check_output([
                    "antechamber", "-dr", "n", "-i", "charged.mol2", "-fi", "mol2", "-o", "charges2.mol2", "-fo",
                    "mol2", "-c", "wc",  "-cf", "charges.txt", "-pf", "yes"])
                #os.system('cat charges.txt')
                # Check to ensure charges were actually produced
                if not os.path.exists('charges.txt'):
                    # TODO: copy files into local directory to aid debugging?
                    raise ChargeCalculationError(
                        "Antechamber/sqm partial charge calculation failed on "
                        "molecule {} (SMILES {})".format(
                            molecule.name, molecule.to_smiles()))
                # Read the charges
                with open('charges.txt', 'r') as infile:
                    contents = infile.read()
                text_charges = contents.split()
                charges = np.zeros([molecule.n_atoms], np.float64)
                for index, token in enumerate(text_charges):
                    charges[index] = float(token)
                # TODO: Ensure that the atoms in charged.mol2 are in the same order as in molecule.sdf

        charges = unit.Quantity(charges, unit.elementary_charge)
        molecule.partial_charges = charges


    def compute_partial_charges_am1bcc(self, molecule, use_conformers=None, strict_n_conformers=False):
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

        Raises
        ------
        ValueError if the requested charge method could not be handled

        """

        import warnings
        warnings.warn("compute_partial_charges_am1bcc will be deprecated in an upcoming release. "
                      "Use assign_partial_charges(partial_charge_method='am1bcc') instead.",
                      DeprecationWarning)

        self.assign_partial_charges(molecule,
                                    partial_charge_method="AM1BCC",
                                    use_conformers=use_conformers,
                                    strict_n_conformers=strict_n_conformers)

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
        datasp.insert(1, 'printbondorders=1, \n /')
        # Reassemble the file text
        new_data = ''.join(datasp)
        # Write the new file contents, overwriting the original file.
        with open(file_path, 'w') as of:
            of.write(new_data)

    def _get_fractional_bond_orders_from_sqm_out(self, file_path, validate_elements=None):
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

        begin_sep = ''' Bond Orders
 
  QMMM:    NUM1 ELEM1 NUM2 ELEM2      BOND_ORDER
'''
        end_sep = '''

           --------- Calculation Completed ----------
'''
        # Extract the chunk of text between begin_sep and end_sep, and split it by newline
        fbo_lines = data.split(begin_sep)[1].split(end_sep)[0].split('\n')

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
                if ((atom_element_1 != validate_elements[atom_index_1-1]) or
                    (atom_element_2 != validate_elements[atom_index_2-1])):
                    #raise ValueError('\n'.join(fbo_lines))
                    raise ValueError(f'Elements or indexing in sqm output differ from expectation. '
                                     f'Expected {validate_elements[atom_index_1]} with index {atom_index_1} and '
                                     f'{validate_elements[atom_index_2]} with index {atom_index_2}, '
                                     f'but SQM output has {atom_element_1} and {atom_element_2} for the same atoms.')


            # To make lookup easier, we identify bonds as integer tuples with the lowest atom index
            # first and the highest second.
            index_tuple = tuple(sorted([atom_index_1, atom_index_2]))
            bond_orders[index_tuple] = bond_order
        return bond_orders

    def assign_fractional_bond_orders(self, molecule, bond_order_model=None, use_conformers=None):
        """
        Update and store list of bond orders this molecule. Bond orders are stored on each
        bond, in the `bond.fractional_bond_order` attribute.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        molecule : openforcefield.topology.molecule Molecule
            The molecule to assign wiberg bond orders to
        bond_order_model : str, optional, default=None
            The charge model to use. Only allowed value is 'am1-wiberg'. If None, 'am1-wiberg' will be used.
        use_conformers : iterable of simtk.unit.Quantity(np.array) with shape (n_atoms, 3) and dimension of distance, optional, default=None
            The conformers to use for fractional bond order calculation. If None, an appropriate number
            of conformers will be generated by an available ToolkitWrapper.
         """
        from openforcefield.topology import Molecule
        # Find the path to antechamber
        # TODO: How should we implement find_executable?
        ANTECHAMBER_PATH = find_executable("antechamber")
        if ANTECHAMBER_PATH is None:
            raise (IOError("Antechamber not found, cannot run "
                           "AmberToolsToolkitWrapper.assign_fractional_bond_orders()"))

        # Make a copy since we'll be messing with this molecule's conformers
        temp_mol = Molecule(molecule)

        if use_conformers is None:
            temp_mol.generate_conformers(n_conformers=1)
        else:
            temp_mol._conformers = None
            for conformer in use_conformers:
                temp_mol.add_conformer(conformer)

        if len(temp_mol.conformers) == 0:
            raise ValueError(
                "No conformers present in molecule submitted for fractional bond order calculation. Consider "
                "loading the molecule from a file with geometry already present or running "
                "molecule.generate_conformers() before calling molecule.assign_fractional_bond_orders"
            )
        if len(temp_mol.conformers) > 1:
            logger.warning(f"Warning: In AmberToolsToolkitWrapper.assign_fractional_bond_orders: "
                           f"Molecule '{molecule.name}' has more than one conformer, but this function "
                           f"will only generate fractional bond orders for the first one.")


        # Compute bond orders
        bond_order_model_to_antechamber_keyword = {'am1-wiberg': 'mul'}
        supported_bond_order_models = list(bond_order_model_to_antechamber_keyword.keys())
        bond_order_model = bond_order_model.lower()
        if bond_order_model is None:
            bond_order_model = 'am1-wiberg'
        if bond_order_model not in supported_bond_order_models:
            raise ValueError(f"Bond order model '{bond_order_model}' is not supported by AmberToolsToolkitWrapper. "
                             f"Supported models are {supported_bond_order_models}")
        ac_charge_keyword = bond_order_model_to_antechamber_keyword[bond_order_model]

        from openforcefield.utils import temporary_directory, temporary_cd
        with temporary_directory() as tmpdir:
            with temporary_cd(tmpdir):
                net_charge = temp_mol.total_charge
                # Write out molecule in SDF format
                self._rdkit_toolkit_wrapper.to_file(
                    temp_mol, 'molecule.sdf', file_format='sdf')
                # Prepare sqm.in file as if we were going to run charge calc
                # TODO: Add error handling if antechamber chokes
                subprocess.check_output([
                    "antechamber", "-i", "molecule.sdf", "-fi", "sdf", "-o",
                    "sqm.in", "-fo", "sqmcrt", "-pf", "yes", "-c", ac_charge_keyword,
                    "-nc", str(net_charge)
                ])
                # Modify sqm.in to request bond order calculation
                self._modify_sqm_in_to_request_bond_orders("sqm.in")
                # Run sqm to get bond orders
                subprocess.check_output(["sqm", "-i", "sqm.in", "-o", "sqm.out", "-O"])
                # Ensure that antechamber/sqm did not change the indexing by checking against
                # an ordered list of element symbols for this molecule
                expected_elements = [at.element.symbol for at in molecule.atoms]
                bond_orders = self._get_fractional_bond_orders_from_sqm_out('sqm.out',
                                                                            validate_elements=expected_elements)

        # Note that sqm calculate WBOs for ALL PAIRS of atoms, not just those that have
        # bonds defined in the original molecule. So here we iterate over the bonds in
        # the original molecule and only nab the WBOs for those.
        for bond in molecule.bonds:
            # The atom index tuples that act as bond indices are ordered from lowest to highest by
            # _get_fractional_bond_orders_from_sqm_out, so here we make sure that we look them up in
            # sorted order as well
            sorted_atom_indices = sorted(tuple([bond.atom1_index+1, bond.atom2_index+1]))
            bond.fractional_bond_order = bond_orders[tuple(sorted_atom_indices)]



#=============================================================================================
# Toolkit registry
#=============================================================================================


class ToolkitRegistry:
    """
    Registry for ToolkitWrapper objects

    Examples
    --------

    Register toolkits in a specified order, skipping if unavailable

    >>> from openforcefield.utils.toolkits import ToolkitRegistry
    >>> toolkit_registry = ToolkitRegistry()
    >>> toolkit_precedence = [OpenEyeToolkitWrapper, RDKitToolkitWrapper, AmberToolsToolkitWrapper]
    >>> for toolkit in toolkit_precedence:
    ...     if toolkit.is_available():
    ...         toolkit_registry.register_toolkit(toolkit)

    Register specified toolkits, raising an exception if one is unavailable

    >>> toolkit_registry = ToolkitRegistry()
    >>> toolkits = [OpenEyeToolkitWrapper, AmberToolsToolkitWrapper]
    >>> for toolkit in toolkits:
    ...     toolkit_registry.register_toolkit(toolkit)

    Register all available toolkits in arbitrary order

    >>> from openforcefield.utils import all_subclasses
    >>> toolkits = all_subclasses(ToolkitWrapper)
    >>> for toolkit in toolkit_precedence:
    ...     if toolkit.is_available():
    ...         toolkit_registry.register_toolkit(toolkit)

    Retrieve the global singleton toolkit registry, which is created when this module is imported from all available
    toolkits:

    >>> from openforcefield.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY as toolkit_registry
    >>> available_toolkits = toolkit_registry.registered_toolkits

    .. warning :: This API is experimental and subject to change.
    """

    def __init__(self,
                 register_imported_toolkit_wrappers=False,
                 toolkit_precedence=None,
                 exception_if_unavailable=True):
        """
        Create an empty toolkit registry.

        Parameters
        ----------
        register_imported_toolkit_wrappers : bool, optional, default=False
            If True, will attempt to register all imported ToolkitWrapper subclasses that can be found, in no particular
             order.
        toolkit_precedence : list, optional, default=None
            List of toolkit wrapper classes, in order of desired precedence when performing molecule operations. If
            None, defaults to [OpenEyeToolkitWrapper, RDKitToolkitWrapper, AmberToolsToolkitWrapper].
        exception_if_unavailable : bool, optional, default=True
            If True, an exception will be raised if the toolkit is unavailable

        """

        self._toolkits = list()

        if toolkit_precedence is None:
            toolkit_precedence = [
                OpenEyeToolkitWrapper, RDKitToolkitWrapper,
                AmberToolsToolkitWrapper
            ]

        if register_imported_toolkit_wrappers:
            # TODO: The precedence ordering of any non-specified remaining wrappers will be arbitrary.
            # How do we fix this?
            # Note: The precedence of non-specifid wrappers may be determined by the order in which
            # they were defined
            all_importable_toolkit_wrappers = all_subclasses(ToolkitWrapper)
            for toolkit in all_importable_toolkit_wrappers:
                if toolkit in toolkit_precedence:
                    continue
                toolkit_precedence.append(toolkit)

        for toolkit in toolkit_precedence:
            self.register_toolkit(toolkit, exception_if_unavailable=exception_if_unavailable)

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

    def register_toolkit(self,
                         toolkit_wrapper,
                         exception_if_unavailable=True):
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
                msg = "Unable to load toolkit '{}'. ".format(toolkit_wrapper._toolkit_name)
                if exception_if_unavailable:
                    raise ToolkitUnavailableException(msg)
                else:
                    if 'OpenEye' in msg:
                        msg += "The Open Force Field Toolkit does not require the OpenEye Toolkits, and can " \
                               "use RDKit/AmberTools instead. However, if you have a valid license for the " \
                               "OpenEye Toolkits, consider installing them for faster performance and additional " \
                               "file format support: " \
                               "https://docs.eyesopen.com/toolkits/python/quickstart-python/linuxosx.html " \
                               "OpenEye offers free Toolkit licenses for academics: " \
                               "https://www.eyesopen.com/academic-licensing"
                    logger.warning(f"Warning: {msg}")
                return

        # Add toolkit to the registry.
        self._toolkits.append(toolkit_wrapper)

    def add_toolkit(self, toolkit_wrapper):
        """
        Append a ToolkitWrapper onto the list of toolkits in this ToolkitRegistry

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        toolkit_wrapper : openforcefield.utils.ToolkitWrapper
            The ToolkitWrapper object to add to the list of registered toolkits

        """
        if not isinstance(toolkit_wrapper, ToolkitWrapper):
            msg = "Something other than a ToolkitWrapper object was passed to ToolkitRegistry.add_toolkit()\n"
            msg += "Given object {} of type {}".format(toolkit_wrapper,
                                                       type(toolkit_wrapper))
            raise Exception(msg)
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

        >>> from openforcefield.topology import Molecule
        >>> molecule = Molecule.from_smiles('Cc1ccccc1')
        >>> toolkit_registry = ToolkitRegistry(register_imported_toolkit_wrappers=True)
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
            method_name)
        msg += 'Available toolkits are: {}\n'.format(self.registered_toolkits)
        raise NotImplementedError(msg)

    # TODO: Can we instead register available methods directly with `ToolkitRegistry`, so we can just use `ToolkitRegistry.method()`?
    def call(self, method_name, *args, **kwargs):
        """
        Execute the requested method by attempting to use all registered toolkits in order of precedence.

        ``*args`` and ``**kwargs`` are passed to the desired method, and return values of the method are returned

        This is a convenient shorthand for ``toolkit_registry.resolve_method(method_name)(*args, **kwargs)``

        Parameters
        ----------
        method_name : str
            The name of the method to execute

        Raises
        ------
        NotImplementedError if the requested method cannot be found among the registered toolkits

        Examples
        --------

        Create a molecule, and call the toolkit ``to_smiles()`` method directly

        >>> from openforcefield.topology import Molecule
        >>> molecule = Molecule.from_smiles('Cc1ccccc1')
        >>> toolkit_registry = ToolkitRegistry(register_imported_toolkit_wrappers=True)
        >>> smiles = toolkit_registry.call('to_smiles', molecule)

        """
        # TODO: catch ValueError and compile list of methods that exist but rejected the specific parameters because they did not implement the requested methods

        errors = list()
        for toolkit in self._toolkits:
            if hasattr(toolkit, method_name):
                method = getattr(toolkit, method_name)
                try:
                    return method(*args, **kwargs)
                except NotImplementedError as not_implemented_error:
                    errors.append((toolkit, not_implemented_error))
                except TypeError as type_error:
                    errors.append((toolkit, type_error))
                except ValueError as value_error:
                    errors.append((toolkit, value_error))

        # No toolkit was found to provide the requested capability
        # TODO: Can we help developers by providing a check for typos in expected method names?
        msg = f'No registered toolkits can provide the capability "{method_name}" ' \
              f'for args "{args}" and kwargs "{kwargs}"\n'

        msg += 'Available toolkits are: {}\n'.format(self.registered_toolkits)
        # Append information about toolkits that implemented the method, but could not handle the provided parameters
        for toolkit, error in errors:
            msg += ' {} : {}\n'.format(toolkit, error)
        # for toolkit, value_error in value_errors:
        #     msg += ' {} : {}\n'.format(toolkit, value_error)
        raise ValueError(msg)


#=============================================================================================
# GLOBAL TOOLKIT REGISTRY
#=============================================================================================

# Create global toolkit registry, where all available toolkits are registered
# TODO: Should this be all lowercase since it's not a constant?
GLOBAL_TOOLKIT_REGISTRY = ToolkitRegistry(
    register_imported_toolkit_wrappers=True,
    exception_if_unavailable=False)

#=============================================================================================
# SET GLOBAL TOOLKIT-AVAIABLE VARIABLES
#=============================================================================================

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

#=============================================================================================
# WARN IF INSUFFICIENT TOOLKITS INSTALLED
#=============================================================================================

# Define basic toolkits that handle essential file I/O

BASIC_CHEMINFORMATICS_TOOLKITS = [RDKitToolkitWrapper, OpenEyeToolkitWrapper]

# Ensure we have at least one basic toolkit
if sum([
        tk.is_available()
        for tk in GLOBAL_TOOLKIT_REGISTRY.registered_toolkits
        if type(tk) in BASIC_CHEMINFORMATICS_TOOLKITS
]) == 0:
    msg = 'WARNING: No basic cheminformatics toolkits are available.\n'
    msg += 'At least one basic toolkit is required to handle SMARTS matching and file I/O. \n'
    msg += 'Please install at least one of the following basic toolkits:\n'
    for wrapper in all_subclasses(ToolkitWrapper):
        if wrapper.toolkit_name is not None:
            msg += '{} : {}\n'.format(
                wrapper._toolkit_name,
                wrapper._toolkit_installation_instructions)
    print(msg)
