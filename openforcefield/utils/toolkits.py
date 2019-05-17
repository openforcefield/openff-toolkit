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
    _toolkit_file_read_formats = [
        'CAN', 'CDX', 'CSV', 'FASTA', 'INCHI', 'INCHIKEY', 'ISM', 'MDL', 'MF',
        'MMOD', 'MOL2', 'MOL2H', 'MOPAC', 'OEB', 'PDB', 'RDF', 'SDF', 'SKC',
        'SLN', 'SMI', 'USM', 'XYC'
    ]
    _toolkit_file_write_formats = [
        'CAN', 'CDX', 'CSV', 'FASTA', 'INCHI', 'INCHIKEY', 'ISM', 'MDL', 'MF',
        'MMOD', 'MOL2', 'MOL2H', 'MOPAC', 'OEB', 'PDB', 'RDF', 'SDF', 'SKC',
        'SLN', 'SMI', 'USM', 'XYC'
    ]

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

    def from_object(self, object):
        """
        If given an OEMol (or OEMol-derived object), this function will load it into an openforcefield.topology.molecule
        Otherwise, it will return False.

        Parameters
        ----------
        object : A molecule-like object
            An object to by type-checked.

        Returns
        -------
        Molecule
            An openforcefield.topology.molecule Molecule.

        Raises
        ------
        NotImplementedError
            If the object could not be converted into a Molecule.
        """
        # TODO: Add tests for the from_object functions
        from openeye import oechem
        if isinstance(object, oechem.OEMolBase):
            return self.from_openeye(object)
        raise NotImplementedError('Cannot create Molecule from {} object'.format(type(object)))

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
            mol = Molecule.from_openeye(
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
        Create a Molecule from an OpenEye molecule.

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
        for oeatom in oemol.GetAtoms():
            oe_idx = oeatom.GetIdx()
            atomic_number = oeatom.GetAtomicNum()
            formal_charge = oeatom.GetFormalCharge()
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
        for oe_idx, oe_atom in enumerate(oemol.GetAtoms()):
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
            oeatom.SetFormalCharge(atom.formal_charge)
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

    def from_smiles(self, smiles, hydrogens_are_explicit=False):
        """
        Create a Molecule from a SMILES string using the OpenEye toolkit.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        smiles : str
            The SMILES string to turn into a molecule
        hydrogens_are_explicit : bool, default = False
            If False, OE will perform hydrogen addition using OEAddExplicitHydrogens

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
                raise Exception(
                    "Addition of explicit hydrogens failed in from_openeye")
        # TODO: Add allow_undefined_stereo to this function, and pass to from_openeye?
        molecule = self.from_openeye(oemol)
        return molecule

    def generate_conformers(self, molecule, n_conformers=1, clear_existing=True):
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
        omega.SetRMSThreshold(1.0)
        #Don't generate random stereoisomer if not specified
        omega.SetStrictStereo(True)
        status = omega(oemol)

        if status is False:
            raise Exception("OpenEye Omega conformer generation failed")

        molecule2 = self.from_openeye(oemol)

        if clear_existing:
            molecule._conformers = list()

        for conformer in molecule2._conformers:
            molecule._add_conformer(conformer)

    def compute_partial_charges(self, molecule, quantum_chemical_method="AM1-BCC", partial_charge_method='None'):
        #charge_model="am1bcc"):
        """
        Compute partial charges with OpenEye quacpac

        .. warning :: This API is experimental and subject to change.

        .. todo ::

           * Should the default be ELF?
           * Can we expose more charge models?


        Parameters
        ----------
        molecule : Molecule
            Molecule for which partial charges are to be computed
        charge_model : str, optional, default=None
            The charge model to use. One of ['noop', 'mmff', 'mmff94', 'am1bcc', 'am1bccnosymspt', 'amber',
            'amberff94', 'am1bccelf10']
            If None, 'am1bcc' will be used.

        Returns
        -------
        charges : numpy.array of shape (natoms) of type float
            The partial charges

        """
        raise NotImplementedError
        # TODO: Implement this in a way that's compliant with SMIRNOFF's <ChargeIncrementModel> tag when the spec gets finalized

        # from openeye import oequacpac
        # import numpy as np
        #
        # if molecule.n_conformers == 0:
        #     raise Exception(
        #         "No conformers present in molecule submitted for partial charge calculation. Consider "
        #         "loading the molecule from a file with geometry already present or running "
        #         "molecule.generate_conformers() before calling molecule.compute_partial_charges"
        #     )
        # oemol = molecule.to_openeye()
        #
        # ## This seems like a big decision. Implemented a simple solution here. Not to be considered final.
        # ## Some discussion at https://github.com/openforcefield/openforcefield/pull/86#issuecomment-350111236
        #
        # if charge_model is None:
        #     charge_model = "am1bcc"
        #
        # if charge_model == "noop":
        #     result = oequacpac.OEAssignCharges(oemol,
        #                                        oequacpac.OEChargeEngineNoOp())
        # elif charge_model == "mmff" or charge_model == "mmff94":
        #     result = oequacpac.OEAssignCharges(oemol,
        #                                        oequacpac.OEMMFF94Charges())
        # elif charge_model == "am1bcc":
        #     result = oequacpac.OEAssignCharges(oemol,
        #                                        oequacpac.OEAM1BCCCharges())
        # elif charge_model == "am1bccnosymspt":
        #     optimize = True
        #     symmetrize = True
        #     result = oequacpac.OEAssignCharges(
        #         oemol, oequacpac.OEAM1BCCCharges(not optimize, not symmetrize))
        # elif charge_model == "amber" or charge_model == "amberff94":
        #     result = oequacpac.OEAssignCharges(oemol,
        #                                        oequacpac.OEAmberFF94Charges())
        # elif charge_model == "am1bccelf10":
        #     result = oequacpac.OEAssignCharges(
        #         oemol, oequacpac.OEAM1BCCELF10Charges())
        # else:
        #     raise ValueError('charge_model {} unknown'.format(charge_model))
        #
        # if result is False:
        #     raise Exception('Unable to assign charges')
        #
        # # Extract and return charges
        # ## TODO: Behavior when given multiple conformations?
        # ## TODO: Make sure atom mapping remains constant
        #
        # charges = unit.Quantity(
        #     np.zeros([oemol.NumAtoms()], np.float64), unit.elementary_charge)
        # for index, atom in enumerate(oemol.GetAtoms()):
        #     charge = atom.GetPartialCharge()
        #     charge = charge * unit.elementary_charge
        #     charges[index] = charge
        #
        # if ((charges / unit.elementary_charge) == 0.
        #     ).all() and not (charge_model == 'noop'):
        #     # TODO: These will be 0 if the charging failed. What behavior do we want in that case?
        #     raise Exception(
        #         "Partial charge calculation failed. Charges from compute_partial_charges() are all 0."
        #     )
        # molecule.set_partial_charges(charges)



    def compute_partial_charges_am1bcc(self, molecule):
        """
        Compute AM1BCC partial charges with OpenEye quacpac

        .. warning :: This API is experimental and subject to change.

        .. todo ::

           * Should the default be ELF?
           * Can we expose more charge models?


        Parameters
        ----------
        molecule : Molecule
            Molecule for which partial charges are to be computed

        Returns
        -------
        charges : numpy.array of shape (natoms) of type float
            The partial charges

        """
        from openeye import oequacpac
        import numpy as np

        if molecule.n_conformers == 0:
            raise Exception(
                "No conformers present in molecule submitted for partial charge calculation. Consider "
                "loading the molecule from a file with geometry already present or running "
                "molecule.generate_conformers() before calling molecule.compute_partial_charges"
            )
        oemol = molecule.to_openeye()

        result = oequacpac.OEAssignCharges(oemol, oequacpac.OEAM1BCCELF10Charges())

        if result is False:
            raise Exception('Unable to assign charges')

        # Extract and return charges
        ## TODO: Make sure atom mapping remains constant

        charges = unit.Quantity(
            np.zeros([oemol.NumAtoms()], np.float64), unit.elementary_charge)
        for index, atom in enumerate(oemol.GetAtoms()):
            charge = atom.GetPartialCharge()
            charge = charge * unit.elementary_charge
            charges[index] = charge

        if ((charges / unit.elementary_charge) == 0.).all():
            # TODO: These will be 0 if the charging failed. What behavior do we want in that case?
            raise Exception(
                "Partial charge calculation failed. Charges from compute_partial_charges() are all 0."
            )
        return charges

    def compute_wiberg_bond_orders(self, molecule, charge_model=None):
        """
        Update and store list of bond orders this molecule. Can be used for initialization of bondorders list, or
        for updating bond orders in the list.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        molecule : openforcefield.topology.molecule Molecule
            The molecule to assign wiberg bond orders to
        charge_model : str, optional, default=None
            The charge model to use. One of ['am1', 'pm3']. If None, 'am1' will be used.

         """
        # TODO: Cache charged molecule so we don't have to redo the computation (Can we do this given the different
        # AM1 interfaces?)

        from openeye import oequacpac

        oemol = self.to_openeye(molecule)
        if molecule.n_conformers == 0:
            raise Exception(
                "No conformers present in molecule submitted for wiberg bond order calculation. Consider "
                "loading the molecule from a file with geometry already present or running "
                "molecule.generate_conformers() before calling molecule.compute_wiberg_bond_orders()"
            )

        if charge_model is None:
            charge_model = 'am1'

        # Based on example at https://docs.eyesopen.com/toolkits/python/quacpactk/examples_summary_wibergbondorders.html
        am1 = oequacpac.OEAM1()
        am1results = oequacpac.OEAM1Results()
        am1options = am1.GetOptions()
        if charge_model == "am1":
            am1options.SetSemiMethod(oequacpac.OEMethodType_AM1)
        elif charge_model == "pm3":
            # TODO: Make sure that modifying am1options actually works
            am1options.SetSemiMethod(oequacpac.OEMethodType_PM3)
        else:
            raise ValueError('charge_model {} unknown'.format(charge_model))

        #for conf in oemol.GetConfs():
        #TODO: How to handle multiple confs here?
        status = am1.CalcAM1(am1results, oemol)

        if status is False:
            raise Exception(
                'Unable to assign charges (in the process of calculating Wiberg bond orders)'
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
        # Make a copy of molecule so we don't influence original (probably safer than deepcopy per C Bayly)
        mol = oechem.OEMol(oemol)

        # Set up query
        qmol = oechem.OEQMol()
        if not oechem.OEParseSmarts(qmol, smarts):
            raise ValueError("Error parsing SMARTS '%s'" % smarts)

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
        substructure_search = oechem.OESubSearch(qmol)
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
        oemol = molecule.to_openeye()
        return self._find_smarts_matches(oemol, smarts)


class RDKitToolkitWrapper(ToolkitWrapper):
    """
    RDKit toolkit wrapper

    .. warning :: This API is experimental and subject to change.
    """
    _toolkit_name = 'The RDKit'
    _toolkit_installation_instructions = 'A conda-installable version of the free and open source RDKit cheminformatics ' \
                                         'toolkit can be found at: https://anaconda.org/rdkit/rdkit'
    _toolkit_file_read_formats = ['SDF', 'MOL', 'SMI']  #TODO: Add TDT support
    _toolkit_file_write_formats = ['SDF', 'MOL', 'SMI', 'PDB']

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

    def from_object(self, object):
        """
        If given an rdchem.Mol (or rdchem.Mol-derived object), this function will load it into an
        openforcefield.topology.molecule. Otherwise, it will return False.

        Parameters
        ----------
        object : A rdchem.Mol-derived object
            An object to be type-checked and converted into a Molecule, if possible.

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
            return self.from_rdkit(object)
        raise NotImplementedError('Cannot create Molecule from {} object'.format(type(object)))

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
        from openforcefield.topology import Molecule
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
                mol = Molecule.from_rdkit(
                    rdmol,
                    allow_undefined_stereo=allow_undefined_stereo
                )

                mols.append(mol)
        elif (file_format == 'SMI'):
            # TODO: We have to do some special stuff when we import SMILES (currently
            # just adding H's, but could get fancier in the future). It might be
            # worthwhile to parse the SMILES file ourselves and pass each SMILES
            # through the from_smiles function instead
            for rdmol in Chem.SmilesMolSupplier(file_path):
                rdmol = Chem.AddHs(rdmol)
                mol = Molecule.from_rdkit(rdmol)
                mols.append(mol)

        elif (file_format == 'PDB'):
            raise Exception(
                "RDKit can not safely read PDBs on their own. Information about bond order and aromaticity "
                "is likely to be lost.")
            # TODO: See if we can implement PDB+mol/smi combinations to get complete bond information.
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
        from openforcefield.topology import Molecule
        from rdkit import Chem

        mols = []

        if (file_format == "MOL") or (file_format == "SDF"):
            # TODO: Iterate over all mols in file_data
            for rdmol in Chem.ForwardSDMolSupplier(file_obj):
                mol = Molecule.from_rdkit(rdmol)
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
                "is likely to be lost.")
            # TODO: See if we can implement PDB+mol/smi combinations to get complete bond information.
            # https://github.com/openforcefield/openforcefield/issues/121
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
        from rdkit import Chem
        file_format = file_format.upper()
        rdmol = self.to_rdkit(molecule)
        rdkit_writers = {
            'SDF': Chem.SDWriter,
            'PDB': Chem.PDBWriter,
            'SMI': Chem.SmilesWriter,
            'TDT': Chem.TDTWriter
        }
        writer = rdkit_writers[file_format](file_obj)
        writer.write(rdmol)
        writer.close()

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
        from rdkit import Chem
        file_format = file_format.upper()
        with open(file_path, 'w') as file_obj:
            rdmol = self.to_rdkit(molecule)
            rdkit_writers = {
                'SDF': Chem.SDWriter,
                'PDB': Chem.PDBWriter,
                'SMI': Chem.SmilesWriter,
                'TDT': Chem.TDTWriter
            }
            writer = rdkit_writers[file_format](file_obj)
            writer.write(rdmol)
            writer.close()

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

    def from_smiles(self, smiles, hydrogens_are_explicit=False):
        """
        Create a Molecule from a SMILES string using the RDKit toolkit.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        smiles : str
            The SMILES string to turn into a molecule
        hydrogens_are_explicit : bool, default=False
            If False, RDKit will perform hydrogen addition using Chem.AddHs

        Returns
        -------
        molecule : openforcefield.topology.Molecule
            An openforcefield-style molecule.
        """
        from openforcefield.topology.molecule import Molecule
        from rdkit import Chem

        rdmol = Chem.MolFromSmiles(smiles, sanitize=False)
        # TODO: I think UpdatePropertyCache(strict=True) is called anyway in Chem.SanitizeMol().
        rdmol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(rdmol, Chem.SANITIZE_ALL ^ Chem.SANITIZE_ADJUSTHS ^ Chem.SANITIZE_SETAROMATICITY)
        Chem.SetAromaticity(rdmol, Chem.AromaticityModel.AROMATICITY_MDL)

        # Chem.MolFromSmiles adds bond directions (i.e. ENDDOWNRIGHT/ENDUPRIGHT), but
        # doesn't set bond.GetStereo(). We need to call AssignStereochemistry for that.
        Chem.AssignStereochemistry(rdmol)

        # Throw an exception/warning if there is unspecified stereochemistry.
        self._detect_undefined_stereo(rdmol, err_msg_prefix='Unable to make OFFMol from SMILES: ')

        # Add explicit hydrogens if they aren't there already
        if not hydrogens_are_explicit:
            rdmol = Chem.AddHs(rdmol)

        # TODO: Add allow_undefined_stereo to this function, and pass to from_rdkit?
        molecule = Molecule.from_rdkit(rdmol)

        return molecule

    def generate_conformers(self, molecule, n_conformers=1, clear_existing=True):
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
        clear_existing : bool, default=True
            Whether to overwrite existing conformers for the molecule.
        
        
        """
        from rdkit.Chem import AllChem
        rdmol = self.to_rdkit(molecule)
        # TODO: This generates way more conformations than omega, given the same nConfs and RMS threshold. Is there some way to set an energy cutoff as well?
        AllChem.EmbedMultipleConfs(
            rdmol,
            numConfs=n_conformers,
            pruneRmsThresh=1.0,
            randomSeed=1,
            #params=AllChem.ETKDG()
        )
        molecule2 = self.from_rdkit(rdmol)

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
        for rda in rdmol.GetAtoms():
            rd_idx = rda.GetIdx()

            # create a new atom
            #atomic_number = oemol.NewAtom(rda.GetAtomicNum())
            atomic_number = rda.GetAtomicNum()
            formal_charge = rda.GetFormalCharge()
            is_aromatic = rda.GetIsAromatic()
            if rda.HasProp('_Name'):
                name = rda.GetProp('_Name')
            else:
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
                    raise Exception(
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
            rdatom.SetFormalCharge(atom.formal_charge)
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
        for match in rdmol.GetSubstructMatches(qmol, uniquify=False):
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
    _toolkit_file_read_formats = []
    _toolkit_file_write_formats = []

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

    def compute_partial_charges(self, molecule, charge_model=None):
        """
        Compute partial charges with AmberTools using antechamber/sqm

        .. warning :: This API experimental and subject to change.

        .. todo ::

           * Do we want to also allow ESP/RESP charges?

        Parameters
        ----------
        molecule : Molecule
            Molecule for which partial charges are to be computed
        charge_model : str, optional, default=None
            The charge model to use. One of ['gas', 'mul', 'bcc']. If None, 'bcc' will be used.


        Raises
        ------
        ValueError if the requested charge method could not be handled

        Notes
        -----
        Currently only sdf file supported as input and mol2 as output
        https://github.com/choderalab/openmoltools/blob/master/openmoltools/packmol.py

        """
        raise NotImplementedError
        # TODO: Implement this in a way that's compliant with SMIRNOFF's <ChargeIncrementModel> tag when the spec gets finalized

        # import os
        # from simtk import unit
        #
        # if charge_model is None:
        #     charge_model = 'bcc'
        #
        # # Check that the requested charge method is supported
        # # Needs to be fixed: 'cm1', 'cm2',
        # SUPPORTED_ANTECHAMBER_CHARGE_MODELS = ['gas', 'mul', 'bcc']
        # if charge_model not in SUPPORTED_ANTECHAMBER_CHARGE_MODELS:
        #     raise ValueError(
        #         'Requested charge method {} not among supported charge '
        #         'methods {}'.format(charge_model,
        #                             SUPPORTED_ANTECHAMBER_CHARGE_MODELS))
        #
        # # Find the path to antechamber
        # # TODO: How should we implement find_executable?
        # ANTECHAMBER_PATH = find_executable("antechamber")
        # if ANTECHAMBER_PATH is None:
        #     raise (IOError("Antechamber not found, cannot run charge_mol()"))
        #
        # if len(molecule._conformers) == 0:
        #     raise Exception(
        #         "No conformers present in molecule submitted for partial charge calculation. Consider "
        #         "loading the molecule from a file with geometry already present or running "
        #         "molecule.generate_conformers() before calling molecule.compute_partial_charges"
        #     )
        #
        #
        # # Compute charges
        # from openforcefield.utils import temporary_directory, temporary_cd
        # with temporary_directory() as tmpdir:
        #     with temporary_cd(tmpdir):
        #         net_charge = molecule.total_charge
        #         # Write out molecule in SDF format
        #         ## TODO: How should we handle multiple conformers?
        #         self._rdkit_toolkit_wrapper.to_file(
        #             molecule, 'molecule.sdf', file_format='sdf')
        #         #os.system('ls')
        #         #os.system('cat molecule.sdf')
        #         # Compute desired charges
        #         # TODO: Add error handling if antechamber chokes
        #         # TODO: Add something cleaner than os.system
        #         os.system(
        #             "antechamber -i molecule.sdf -fi sdf -o charged.mol2 -fo mol2 -pf "
        #             "yes -c {} -nc {}".format(charge_model, net_charge))
        #         #os.system('cat charged.mol2')
        #
        #         # Write out just charges
        #         os.system(
        #             "antechamber -i charged.mol2 -fi mol2 -o charges2.mol2 -fo mol2 -c wc "
        #             "-cf charges.txt -pf yes")
        #         #os.system('cat charges.txt')
        #         # Check to ensure charges were actually produced
        #         if not os.path.exists('charges.txt'):
        #             # TODO: copy files into local directory to aid debugging?
        #             raise Exception(
        #                 "Antechamber/sqm partial charge calculation failed on "
        #                 "molecule {} (SMILES {})".format(
        #                     molecule.name, molecule.to_smiles()))
        #         # Read the charges
        #         with open('charges.txt', 'r') as infile:
        #             contents = infile.read()
        #         text_charges = contents.split()
        #         charges = np.zeros([molecule.n_atoms], np.float64)
        #         for index, token in enumerate(text_charges):
        #             charges[index] = float(token)
        #         # TODO: Ensure that the atoms in charged.mol2 are in the same order as in molecule.sdf
        #
        # charges = unit.Quantity(charges, unit.elementary_charge)
        #
        # molecule.set_partial_charges(charges)

    def compute_partial_charges_am1bcc(self, molecule):
        """
        Compute partial charges with AmberTools using antechamber/sqm. This will calculate AM1-BCC charges on the first
        conformer only.

        .. warning :: This API experimental and subject to change.

        Parameters
        ----------
        molecule : Molecule
            Molecule for which partial charges are to be computed


        Raises
        ------
        ValueError if the requested charge method could not be handled

        """

        import os
        from simtk import unit


        # Find the path to antechamber
        # TODO: How should we implement find_executable?
        ANTECHAMBER_PATH = find_executable("antechamber")
        if ANTECHAMBER_PATH is None:
            raise (IOError("Antechamber not found, cannot run "
                           "AmberToolsToolkitWrapper.compute_partial_charges_am1bcc()"))

        if len(molecule._conformers) == 0:
            raise ValueError(
                "No conformers present in molecule submitted for partial charge calculation. Consider "
                "loading the molecule from a file with geometry already present or running "
                "molecule.generate_conformers() before calling molecule.compute_partial_charges"
            )
        if len(molecule._conformers) > 1:
            logger.warning("In AmberToolsToolkitwrapper.computer_partial_charges_am1bcc: "
                           "Molecule '{}' has more than one conformer, but this function "
                           "will only generate charges for the first one.".format(molecule.name))


        # Compute charges
        from openforcefield.utils import temporary_directory, temporary_cd
        with temporary_directory() as tmpdir:
            with temporary_cd(tmpdir):
                net_charge = molecule.total_charge
                # Write out molecule in SDF format
                ## TODO: How should we handle multiple conformers?
                self._rdkit_toolkit_wrapper.to_file(
                    molecule, 'molecule.sdf', file_format='sdf')
                #os.system('ls')
                #os.system('cat molecule.sdf')
                # Compute desired charges
                # TODO: Add error handling if antechamber chokes
                # TODO: Add something cleaner than os.system
                os.system(
                    "antechamber -i molecule.sdf -fi sdf -o charged.mol2 -fo mol2 -pf "
                    "yes -c bcc -nc {}".format(net_charge))
                #os.system('cat charged.mol2')

                # Write out just charges
                os.system(
                    "antechamber -i charged.mol2 -fi mol2 -o charges2.mol2 -fo mol2 -c wc "
                    "-cf charges.txt -pf yes")
                #os.system('cat charges.txt')
                # Check to ensure charges were actually produced
                if not os.path.exists('charges.txt'):
                    # TODO: copy files into local directory to aid debugging?
                    raise Exception(
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
        return charges


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
            toolkit_wrapper = toolkit_wrapper()

        # Raise exception if not available.
        if not toolkit_wrapper.is_available():
            msg = "Unable to load toolkit {}.".format(toolkit_wrapper)
            if exception_if_unavailable:
                raise ToolkitUnavailableException(msg)
            else:
                logger.warning(msg)
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
        Resolve the requested method name by checking all registered toolkits in order of precedence for one that provides the requested method.

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

        value_errors = list()
        for toolkit in self._toolkits:
            if hasattr(toolkit, method_name):
                method = getattr(toolkit, method_name)
                try:
                    return method(*args, **kwargs)
                except NotImplementedError:
                    pass
                except ValueError as value_error:
                    value_errors.append((toolkit, value_error))

        # No toolkit was found to provide the requested capability
        # TODO: Can we help developers by providing a check for typos in expected method names?
        msg = 'No registered toolkits can provide the capability "{}".\n'.format(
            method_name)
        msg += 'Available toolkits are: {}\n'.format(self.registered_toolkits)
        # Append information about toolkits that implemented the method, but could not handle the provided parameters
        for toolkit, value_error in value_errors:
            msg += ' {} : {}\n'.format(toolkit, value_error)
        raise NotImplementedError(msg)


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
