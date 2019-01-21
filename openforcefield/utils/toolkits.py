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
#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import importlib
from functools import wraps
from openforcefield.utils.utils import inherit_docstrings
from openforcefield.utils import all_subclasses
from openforcefield.typing.chemistry.environment import SMIRKSParsingError
from distutils.spawn import find_executable
from simtk import unit
import numpy as np

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


class MissingPackageError(Exception):
    """This function requires a package that is not installed."""
    pass


class ToolkitUnavailableException(Exception):
    """The requested toolkit is unavailable."""
    # TODO: Allow toolkit to be specified and used in formatting/printing exception.
    pass


class InvalidToolkitError(Exception):
    """A non-toolkit object was received when a toolkit object was expected"""
    def __init__(self, msg):
        super().__init__(self, msg)
        self.msg = msg


class UndefinedStereochemistryError(Exception):
    """A molecule was attempted to be loaded with undefined stereochemistry"""
    def __init__(self, msg):
        super().__init__(self, msg)
        self.msg = msg

#=============================================================================================
# TOOLKIT UTILITY DECORATORS
#=============================================================================================

#=============================================================================================
# UTILITY FUNCTIONS
#=============================================================================================

#=============================================================================================
# CHEMINFORMATICS TOOLKIT WRAPPERS
#=============================================================================================


class ToolkitWrapper(object):
    """
    Toolkit wrapper base class.

    .. warning :: This API experimental and subject to change.
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
                if not cls.toolkit_is_available():
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
    def toolkit_is_available():
        """
        Check whether the corresponding toolkit can be imported

        .. note :: This method call may be expensive.

        Returns
        -------
        is_installed : bool
            True if corresponding toolkit is installed, False otherwise.

        """
        return NotImplementedError

    @classmethod
    def is_available(cls):
        """
        Check whether this toolkit wrapper is available for use because the underlying toolkit can be found.

        .. note :: This method caches the result of any costly checks for file paths or module imports.

        Parameters
        ----------
        is_available : bool
            True if toolkit is available for use, False otherwise

        """
        return NotImplementedError

    def from_file(self,
                  filename,
                  file_format,
                  exception_if_undefined_stereo=True):
        """
        Return an openforcefield.topology.Molecule from a file using this toolkit.
        
        Parameters
        ----------
        filename : str
            The file to read the molecule from
        file_format : str
            Format specifier, usually file suffix (eg. 'MOL2', 'SMI')
            Note that not all toolkits support all formats. Check ToolkitWrapper.toolkit_file_read_formats for details.
        exception_if_undefined_stereo : bool, default=True
            If true, raises an exception if any molecules contain undefined stereochemistry. If false, the function
            skips loading the molecule.
        Returns
        -------
        molecules : Molecule or list of Molecules
            a list of Molecule objects is returned.

        """
        return NotImplementedError

    def from_file_obj(self,
                      file_obj,
                      file_format,
                      exception_if_undefined_stereo=True):
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
        exception_if_undefined_stereo : bool, default=True
            If true, raises an exception if any molecules contain undefined stereochemistry. If false, the function
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
    def toolkit_is_available(
            oetools=('oechem', 'oequacpac', 'oeiupac', 'oeomega')):
        """
        Check if a given OpenEye toolkit component (or set of components) is installed and Licensed

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
            oetools = (oetools, )
        tool_set = set(oetools)
        valid_tool_set = set(tool_keys)
        if tool_set & valid_tool_set == set():
            # Check for empty set intersection
            raise ValueError(
                "Expected OpenEye tools to have at least of the following {}, "
                "but instead got {}".format(tool_keys, oetools))
        try:
            for tool in oetools:
                if tool in tool_keys:
                    # Try loading the module
                    try:
                        module = importlib.import_module('openeye', tool)
                    except SystemError:  # Python 3.4 relative import fix
                        module = importlib.import_module('openeye.' + tool)
                    # Check that we have the license
                    if not getattr(module, tools_license[tool])():
                        raise ImportError
        except ImportError:
            return False
        return True

    @classmethod
    def is_available(cls):
        """
        Check whether this toolkit wrapper is available for use because the underlying toolkit can be found.

        .. note :: This method caches the result of any costly checks for file paths or module imports.

        Parameters
        ----------
        is_available : bool
            True if toolkit wrapper is available for use, False otherwise

        """
        if cls._is_available is None:
            cls._is_available = cls.toolkit_is_available()
                #oetools=('oechem', 'oequacpac')
        return cls._is_available

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
        Molecule or False
            An openforcefield.topology.molecule Molecule, or False if loading was unsuccessful
        """
        # TODO: Add tests for the from_object functions
        from openeye import oechem
        if isinstance(object, oechem.OEMolBase):
            mol = self.from_openeye(object)
            return mol
        else:
            return False

    def from_file(self,
                  filename,
                  file_format,
                  exception_if_undefined_stereo=True):
        """
        Return an openforcefield.topology.Molecule from a file using this toolkit.

        .. warning :: This API is experimental and subject to change.
        
        Parameters
        ----------
        filename : str
            The file to read the molecule from
        file_format : str
            Format specifier, usually file suffix (eg. 'MOL2', 'SMI')
            Note that not all toolkits support all formats. Check ToolkitWrapper.toolkit_file_read_formats for details.
        exception_if_undefined_stereo : bool, default=True
            If true, raises an exception if oemol contains undefined stereochemistry. If false, the function skips
            loading the molecule.

        Returns
        -------
        molecules : list of Molecules
            a list of Molecule objects is returned.

        """
        from openforcefield.topology import Molecule
        from openeye import oechem
        mols = list()
        oemol = oechem.OEMol()
        ifs = oechem.oemolistream(filename)
        while oechem.OEReadMolecule(ifs, oemol):
            oechem.OEPerceiveChiral(oemol)
            oechem.OEAssignAromaticFlags(oemol, oechem.OEAroModel_MDL)
            oechem.OE3DToInternalStereo(oemol)
            mol = Molecule.from_openeye(
                oemol,
                exception_if_undefined_stereo=exception_if_undefined_stereo)
            mols.append(mol)
        return mols

    def from_file_obj(self,
                      file_obj,
                      file_format,
                      exception_if_undefined_stereo=True):
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
        exception_if_undefined_stereo : bool, default=True
            If true, raises an exception if oemol contains undefined stereochemistry. If false, the function skips
            loading the molecule.

        Returns
        -------
        molecules : Molecule or list of Molecules
            a list of Molecule objects is

        """
        from openforcefield.topology import Molecule
        from openeye import oechem
        mols = list()
        oemol = oechem.OEMol()
        file_data = file_obj.read()
        ifs = oechem.oemolistream()
        ifs.openstring(file_data)
        oeformat = getattr(oechem, 'OEFormat_' + file_format)
        ifs.SetFormat(oeformat)
        while oechem.OEReadMolecule(ifs, oemol):
            oechem.OEPerceiveChiral(oemol)
            oechem.OEAssignAromaticFlags(oemol, oechem.OEAroModel_MDL)
            oechem.OE3DToInternalStereo(oemol)
            mol = Molecule.from_openeye(
                oemol,
                exception_if_undefined_stereo=exception_if_undefined_stereo)
            mols.append(mol)
        return mols

    def to_file_obj(self, molecule, file_obj, outfile_format):
        """
        Writes an OpenFF Molecule to a file-like object

        Parameters
        ----------
        molecule : an OpenFF Molecule
            The molecule to write
        file_obj
            The file-like object to write to
        outfile_format
            The format for writing the molecule data

        """
        from openeye import oechem
        from openforcefield.utils import temporary_directory, temporary_cd

        oemol = self.to_openeye(molecule)

        # TODO: This is inefficiently implemented. Is there any way to attach a file-like object to an oemolstream?
        with temporary_directory() as tmpdir:
            with temporary_cd(tmpdir):
                outfile = 'temp_molecule.' + outfile_format
                ofs = oechem.oemolostream(outfile)
                openeye_format = getattr(oechem, 'OEFormat_' + outfile_format)
                ofs.SetFormat(outfile_format)
                oechem.OEWriteMolecule(ofs, oemol)
                ofs.close()
                file_data = open(outfile).read()
        file_obj.write(file_data)

    def to_file(self, molecule, outfile, outfile_format):
        """
        Writes an OpenFF Molecule to a file-like object

        Parameters
        ----------
        molecule : an OpenFF Molecule
            The molecule to write
        outfile
            The filename to write to
        outfile_format
            The format for writing the molecule data

        """
        from openeye import oechem
        oemol = self.to_openeye(molecule)
        ofs = oechem.oemolostream(outfile)
        openeye_format = getattr(oechem, 'OEFormat_' + outfile_format)
        ofs.SetFormat(openeye_format)
        oechem.OEWriteMolecule(ofs, oemol)
        ofs.close()

    @staticmethod
    def _openeye_cip_atom_stereochemistry(oemol, oeatom):
        """
        .. warning :: This API experimental and subject to change.

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
        .. warning :: This API experimental and subject to change.

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
    def from_openeye(oemol, exception_if_undefined_stereo=True):
        """
        Create a Molecule from an OpenEye molecule.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        oemol : openeye.oechem.OEMol
            An OpenEye molecule
        exception_if_undefined_stereo : bool, default=True
            If true, raises an exception if oemol contains undefined stereochemistry. If false, the function skips
            loading the molecule.

        Returns
        -------
        molecule : openforcefield.topology.Molecule
            An openforcefield molecule

        Examples
        --------

        Create a Molecule from an OpenEye OEMol

        >>> toolkit_wrapper = OpenEyeToolkitWrapper()
        >>> molecule = toolkit_wrapper.from_openeye(oemol)

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

            msg = "Unable to make OFFMol from OEMol: OEMol has unspecified stereochemistry. " \
                  "oemol.GetTitle(): {}\n".format(oemol.GetTitle())
            if len(problematic_atoms) != 0:
                msg += "Problematic atoms are:\n"
                for problematic_atom in problematic_atoms:
                    msg += describe_oeatom(problematic_atom) + '\n'
            if len(problematic_bonds) != 0:
                msg += "Problematic bonds are: {}\n".format(problematic_bonds)
            if exception_if_undefined_stereo:
                raise UndefinedStereochemistryError(msg)
            else:
                print(msg)
                return

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
                if (positions == 0 * unit.angstrom).all():
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

        molecule.set_partial_charges(partial_charges)

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
        # TODO: Add exception_if_undefined_stereo to this function, and pass to from_openeye?
        molecule = self.from_openeye(oemol)
        return molecule

    def generate_conformers(self, molecule, clear_existing=True):
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
            The molecule to generate conformers for
        clear_existing : bool, default=True
            Whether to overwrite existing conformers for the molecule
        
        """
        from openeye import oeomega
        oemol = self.to_openeye(molecule)
        omega = oeomega.OEOmega()
        omega.SetMaxConfs(800)
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

    def compute_partial_charges(self, molecule, charge_model="am1bcc"):
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
        from openeye import oequacpac
        import numpy as np

        if molecule.n_conformers == 0:
            raise Exception(
                "No conformers present in molecule submitted for partial charge calculation. Consider "
                "loading the molecule from a file with geometry already present or running "
                "molecule.generate_conformers() before calling molecule.compute_partial_charges"
            )
        oemol = molecule.to_openeye()

        ## This seems like a big decision. Implemented a simple solution here. Not to be considered final.
        ## Some discussion at https://github.com/openforcefield/openforcefield/pull/86#issuecomment-350111236

        if charge_model is None:
            charge_model = "am1bcc"

        if charge_model == "noop":
            result = oequacpac.OEAssignCharges(oemol,
                                               oequacpac.OEChargeEngineNoOp())
        elif charge_model == "mmff" or charge_model == "mmff94":
            result = oequacpac.OEAssignCharges(oemol,
                                               oequacpac.OEMMFF94Charges())
        elif charge_model == "am1bcc":
            result = oequacpac.OEAssignCharges(oemol,
                                               oequacpac.OEAM1BCCCharges())
        elif charge_model == "am1bccnosymspt":
            optimize = True
            symmetrize = True
            result = oequacpac.OEAssignCharges(
                oemol, oequacpac.OEAM1BCCCharges(not optimize, not symmetrize))
        elif charge_model == "amber" or charge_model == "amberff94":
            result = oequacpac.OEAssignCharges(oemol,
                                               oequacpac.OEAmberFF94Charges())
        elif charge_model == "am1bccelf10":
            result = oequacpac.OEAssignCharges(
                oemol, oequacpac.OEAM1BCCELF10Charges())
        else:
            raise ValueError('charge_model {} unknown'.format(charge_model))

        if result is False:
            raise Exception('Unable to assign charges')

        # Extract and return charges
        ## TODO: Behavior when given multiple conformations?
        ## TODO: Make sure atom mapping remains constant

        charges = unit.Quantity(
            np.zeros([oemol.NumAtoms()], np.float64), unit.elementary_charge)
        for index, atom in enumerate(oemol.GetAtoms()):
            charge = atom.GetPartialCharge()
            charge = charge * unit.elementary_charge
            charges[index] = charge

        if ((charges / unit.elementary_charge) == 0.
            ).all() and not (charge_model == 'noop'):
            # TODO: These will be 0 if the charging failed. What behavior do we want in that case?
            raise Exception(
                "Partial charge calculation failed. Charges from compute_partial_charges() are all 0."
            )
        molecule.set_partial_charges(charges)

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

        .. warning :: This API experimental and subject to change.

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
    """
    _toolkit_name = 'The RDKit'
    _toolkit_installation_instructions = 'A conda-installable version of the free and open source RDKit cheminformatics ' \
                                         'toolkit can be found at: https://anaconda.org/rdkit/rdkit'
    _toolkit_file_read_formats = ['SDF', 'MOL', 'SMI']  #TODO: Add TDT support
    _toolkit_file_write_formats = ['SDF', 'MOL', 'SMI', 'PDB']

    @staticmethod
    def toolkit_is_available():
        """
        Check whether the RDKit toolkit can be imported

        Returns
        -------
        is_installed : bool
            True if RDKit is installed, False otherwise.

        """
        try:
            module = importlib.import_module('rdkit', 'Chem')
            return True
        except ImportError:
            return False

    @classmethod
    def is_available(cls):
        """
        Check whether toolkit is available for use.

        Parameters
        ----------
        is_available : bool
            True if toolkit is available for use, False otherwise

        """
        if cls._is_available is None:
            cls._is_available = cls.toolkit_is_available()
        return cls._is_available

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
            An openforcefield.topology.molecule Molecule, or False if loading was unsuccessful
        """
        # TODO: Add tests for the from_object functions
        from rdkit import Chem
        if isinstance(object, Chem.rdchem.Mol):
            mol = self.from_rdkit(object)
            return mol
        else:
            return False

    def from_file(self,
                  filename,
                  file_format,
                  exception_if_undefined_stereo=True):
        """
        Create an openforcefield.topology.Molecule from a file using this toolkit.

        .. warning :: This API is experimental and subject to change.


        Parameters
        ----------
        filename : str
            The file to read the molecule from
        file_format : str
            Format specifier, usually file suffix (eg. 'MOL2', 'SMI')
            Note that not all toolkits support all formats. Check ToolkitWrapper.toolkit_file_read_formats for details.
        exception_if_undefined_stereo : bool, default=True
            If true, raises an exception if oemol contains undefined stereochemistry. If false, the function skips
            loading the molecule.

        Returns
        -------
        molecules : iterable of Molecules
            a list of Molecule objects is returned.

        """
        from openforcefield.topology import Molecule
        from rdkit import Chem
        mols = list()
        if (file_format == 'MOL') or (file_format == 'SDF'):
            for rdmol in Chem.SupplierFromFilename(filename, removeHs=False, sanitize=False, strictParsing=True):
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
                    exception_if_undefined_stereo=exception_if_undefined_stereo
                )

                mols.append(mol)
        elif (file_format == 'SMI'):
            # TODO: We have to do some special stuff when we import SMILES (currently
            # just adding H's, but could get fancier in the future). It might be
            # worthwhile to parse the SMILES file ourselves and pass each SMILES
            # through the from_smiles function instead
            for rdmol in Chem.SmilesMolSupplier(filename):
                rdmol = Chem.AddHs(rdmol)
                mol = Molecule.from_rdkit(rdmol)
                mols.append(mol)

        elif (file_format == 'PDB'):
            raise Exception(
                "RDKit can not safely read PDBs on their own. Information about bond order and aromaticity "
                "is likely to be lost.")
            # TODO: See if we can implement PDB+mol/smi combinations to get complete bond information.
            # https://github.com/openforcefield/openforcefield/issues/121
            rdmol = Chem.MolFromPDBFile(filename, removeHs=False)
            mol = Molecule.from_rdkit(rdmol)
            mols.append(mol)
            # TODO: Add SMI, TDT(?) support

        return mols

    def from_file_obj(self,
                      file_obj,
                      file_format,
                      exception_if_undefined_stereo=True):
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
        exception_if_undefined_stereo : bool, default=True
            If true, raises an exception if oemol contains undefined stereochemistry. If false, the function skips
            loading the molecule.

        Returns
        -------
        molecules : Molecule or list of Molecules
            a list of Molecule objects is returned.

        """
        from openforcefield.topology import Molecule
        from rdkit import Chem
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
            rdmol = Chem.MolFromPDBBlock(file_data)
            mol = Molecule.from_rdkit(rdmol)
            mols.append(mol)
        # TODO: TDT file support
        return mols

    def to_file_obj(self, molecule, file_obj, outfile_format):
        """
        Writes an OpenFF Molecule to a file-like object

        Parameters
        ----------
        molecule : an OpenFF Molecule
            The molecule to write
        file_obj
            The file-like object to write to
        outfile_format
            The format for writing the molecule data

        Returns
        -------

        """
        from rdkit import Chem
        outfile_format = outfile_format.upper()
        rdmol = self.to_rdkit(molecule)
        rdkit_writers = {
            'SDF': Chem.SDWriter,
            'PDB': Chem.PDBWriter,
            'SMI': Chem.SmilesWriter,
            'TDT': Chem.TDTWriter
        }
        writer = rdkit_writers[outfile_format](file_obj)
        writer.write(rdmol)
        writer.close()

    def to_file(self, molecule, outfile, outfile_format):
        """
        Writes an OpenFF Molecule to a file-like object

        Parameters
        ----------
        molecule : an OpenFF Molecule
            The molecule to write
        outfile
            The filename to write to
        outfile_format
            The format for writing the molecule data

        Returns
        ------

        """
        from rdkit import Chem
        outfile_format = outfile_format.upper()
        with open(outfile, 'w') as file_obj:
            rdmol = self.to_rdkit(molecule)
            rdkit_writers = {
                'SDF': Chem.SDWriter,
                'PDB': Chem.PDBWriter,
                'SMI': Chem.SmilesWriter,
                'TDT': Chem.TDTWriter
            }
            writer = rdkit_writers[outfile_format](file_obj)
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
        hydrogens_are_explicit : bool, default = False
            If False, RDKit will perform hydrogen addition using Chem.AddHs

        Returns
        -------
        molecule : openforcefield.topology.Molecule
            An openforcefield-style molecule.
        """
        from openforcefield.topology.molecule import Molecule
        # inherits base class docstring
        from rdkit import Chem
        from rdkit.Chem import EnumerateStereoisomers

        rdmol = Chem.MolFromSmiles(smiles, sanitize=False)
        Chem.SanitizeMol(rdmol, Chem.SANITIZE_ALL ^ Chem.SANITIZE_ADJUSTHS ^ Chem.SANITIZE_SETAROMATICITY)

        # Adding H's can hide undefined bond stereochemistry, so we have to test for undefined stereo here
        unspec_stereo = False
        rdmol_copy = Chem.Mol(rdmol)
        enumsi_opt = EnumerateStereoisomers.StereoEnumerationOptions(
            maxIsomers=2, onlyUnassigned=True)
        stereoisomers = [
            isomer
            for isomer in Chem.EnumerateStereoisomers.EnumerateStereoisomers(
                rdmol_copy, enumsi_opt)
        ]
        if len(stereoisomers) != 1:
            unspec_stereo = True

        if unspec_stereo:
            raise Exception(
                "Unable to make OFFMol from SMILES: SMILES has unspecified stereochemistry: {}"
                .format(smiles))

        # Add explicit hydrogens if they aren't there already
        if not (hydrogens_are_explicit):
            rdmol = Chem.AddHs(rdmol)


        # TODO: Add exception_if_undefined_stereo to this function, and pass to from_rdkit?
        molecule = Molecule.from_rdkit(rdmol)

        return molecule

    def generate_conformers(self, molecule, clear_existing=True):
        """
        Generate molecule conformers using RDKit. 

        .. warning :: This API is experimental and subject to change.

        .. todo ::
        
           * which parameters should we expose? (or can we implement a general system with **kwargs?)
           * will the coordinates be returned in the OpenFF Molecule's own indexing system? Or is there a chance that they'll get reindexed when we convert the input into an RDMol?
        
        Parameters
        ---------
        molecule : a :class:`Molecule` 
            The molecule to generate conformers for
        clear_existing : bool, default=True
            Whether to overwrite existing conformers for the molecule
        
        
        """
        from rdkit.Chem import AllChem
        rdmol = self.to_rdkit(molecule)
        # TODO: This generates way more conformations than omega, given the same nConfs and RMS threshold. Is there some way to set an energy cutoff as well?
        AllChem.EmbedMultipleConfs(
            rdmol,
            numConfs=800,
            pruneRmsThresh=1.0,
            randomSeed=1,
            #params=AllChem.ETKDG()
        )
        molecule2 = self.from_rdkit(rdmol)

        if clear_existing:
            molecule._conformers = list()

        for conformer in molecule2._conformers:
            molecule._add_conformer(conformer)

    def from_rdkit(self, rdmol, exception_if_undefined_stereo=True):
        """
        Create a Molecule from an RDKit molecule.

        Requires the RDKit to be installed.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        rdmol : rkit.RDMol
            An RDKit molecule
        exception_if_undefined_stereo : bool, default=True
            If true, raises an exception if rdmol contains undefined stereochemistry. If false, the function skips
            loading the molecule.

        Returns
        -------
        molecule : openforcefield.Molecule
            An openforcefield molecule

        Examples
        --------

        Create a molecule from an RDKit molecule

        >>> molecule = Molecule.from_rdkit(rdmol)

        """
        from rdkit import Chem
        from openforcefield.topology.molecule import Molecule

        # Check for undefined stereochemistry
        from rdkit.Chem import EnumerateStereoisomers
        # TODO: Does this work for molecules with 3D geometry?
        unspec_stereo = False
        # Use a copy of the input, in case EnumerateStereochemstry changes anything in-place
        rdmol_copy = Chem.Mol(rdmol)
        enumsi_opt = EnumerateStereoisomers.StereoEnumerationOptions(
            maxIsomers=2, onlyUnassigned=True)
        try:
            stereoisomers = [
                isomer
                for isomer in Chem.EnumerateStereoisomers.EnumerateStereoisomers(
                    rdmol_copy, enumsi_opt)
            ]
        except RuntimeError as e:
            msg = "Unable to check stereochemistry for {}. Original error:\n".format(rdmol.GetProp('_Name'))
            msg += str(e)
            if exception_if_undefined_stereo:
                raise UndefinedStereochemistryError(msg)
            else:
                stereoisomers = []
                print(msg)

        # TODO: This will catch undefined tetrahedral centers, but not bond stereochemistry. How can we check for that?
        if len(stereoisomers) != 1:
            unspec_stereo = True

        if unspec_stereo:
            msg = "RDMol has unspecified stereochemistry\n"
            msg += "RDMol name: " + rdmol.GetProp("_Name")
            if exception_if_undefined_stereo:
                raise UndefinedStereochemistryError(
                    "Unable to make OFFMol from RDMol: " + msg
                )
            else:
                print(
                    "WARNING: " + msg
                )
                # TODO: Can we find a way to print more about the error here?
        # Create a new openforcefield Molecule
        mol = Molecule()

        # These checks cause rdkit to choke on one member of our test set: ZINC16448882
        # http://zinc.docking.org/substance/16448882
        # This has a pentavalent nitrogen, which I think is really resonance-stabilized.
        # I think we should allow this as input, since a fractional bond order calculation will probably sort it out.

        #Chem.SanitizeMol(rdmol, Chem.SANITIZE_ALL ^ Chem.SANITIZE_SETAROMATICITY ^ Chem.SANITIZE_ADJUSTHS)
        #Chem.SetAromaticity(rdmol, Chem.AromaticityModel.AROMATICITY_MDL)


        # If RDMol has a title save it
        if rdmol.HasProp("_Name"):
            #raise Exception('{}'.format(rdmol.GetProp('name')))
            mol.name = rdmol.GetProp("_Name")
        else:
            mol.name = ""

        # Store all properties
        # TODO: Should there be an API point for storing properties?
        properties = rdmol.GetPropsAsDict()
        mol._properties = properties

        # We store bond orders as integers regardless of aromaticity.
        # In order to properly extract these, we need to have the "Kekulized" version of the rdkit mol
        kekul_mol = Chem.Mol(rdmol)
        Chem.Kekulize(kekul_mol, clearAromaticFlags=False)#True)

        # setting chirality in openeye requires using neighbor atoms
        # therefore we can't do it until after the atoms and bonds are all added
        chiral_atoms = dict()  # {rd_idx: openeye chirality}
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
                if stereo_code == 'S':
                    stereochemistry = 'S'

            atom_index = mol.add_atom(
                atomic_number,
                formal_charge,
                is_aromatic,
                name=name,
                stereochemistry=stereochemistry)
            map_atoms[rd_idx] = atom_index

        # Similar to chirality, stereochemistry of bonds in OE is set relative to their neighbors
        stereo_bonds = list()
        # stereo_bonds stores tuples in the form (oe_bond, rd_idx1, rd_idx2, OE stereo specification)
        # where rd_idx1 and 2 are the atoms on the outside of the bond
        # i.e. Cl and F in the example above
        aro_bond = 0
        for rdb in rdmol.GetBonds():
            rdb_idx = rdb.GetIdx()
            a1 = rdb.GetBeginAtomIdx()
            a2 = rdb.GetEndAtomIdx()

            # Determine bond aromaticity and Kekulized bond order
            is_aromatic = False
            order = rdb.GetBondTypeAsDouble()
            if order == 1.5:
                # get the bond order for this bond in the kekulized molecule
                order = kekul_mol.GetBondWithIdx(
                    rdb.GetIdx()).GetBondTypeAsDouble()
                is_aromatic = True
            # Convert floating-point bond order to integral bond order
            order = int(order)

            # create a new bond
            bond_index = mol.add_bond(map_atoms[a1], map_atoms[a2], order,
                                      is_aromatic)
            map_bonds[rdb_idx] = bond_index

        # Now fill in the cached (structure-dependent) properties. We have to have the 2D structure of the molecule
        # in place first, because each call to add_atom and add_bond invalidates all cached properties
        for rdb in rdmol.GetBonds():
            rdb_idx = rdb.GetIdx()
            offb_idx = map_bonds[rdb_idx]
            offb = mol.bonds[offb_idx]
            # determine if stereochemistry is needed
            stereochemistry = None
            tag = rdb.GetStereo()
            if tag == Chem.BondStereo.STEREOZ:
                stereochemistry = 'Z'
            elif tag == Chem.BondStereo.STEREOE:
                stereochemistry = 'E'
            elif tag == Chem.BondStereo.STEREOTRANS or tag == Chem.BondStereo.STEREOCIS:
                raise Exception(
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
                n_atoms = mol.n_atoms
                # TODO: Will this always be angstrom when loading from RDKit?
                positions = unit.Quantity(
                    np.zeros((n_atoms, 3)), unit.angstrom)
                for rd_idx, off_idx in map_atoms.items():
                    atom_coords = conf.GetPositions()[rd_idx, :] * unit.angstrom
                    positions[off_idx, :] = atom_coords
                mol.add_conformer(positions)

        partial_charges = unit.Quantity(
            np.zeros(mol.n_atoms, dtype=np.float), unit=unit.elementary_charge)

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

            mol.set_partial_charges(partial_charges)
        return mol

    @staticmethod
    def to_rdkit(molecule, aromaticity_model=DEFAULT_AROMATICITY_MODEL):
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

        >>> rdmol = molecule.to_rdkit()

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

        # atom map lets you find atoms again
        map_atoms = dict()  # { molecule index : rdkit index }
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

            map_atoms[index] = rd_index

        for bond in molecule.bonds:
            rdatom1 = map_atoms[bond.atom1.molecule_atom_index]
            rdatom2 = map_atoms[bond.atom2.molecule_atom_index]
            rdmol.AddBond(rdatom1, rdatom2)
            rdbond = rdmol.GetBondBetweenAtoms(rdatom1, rdatom2)
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


        #Assign atom stereochemsitry
        for index, atom in enumerate(molecule.atoms):
            rdatom = rdmol.GetAtomWithIdx(map_atoms[index])
            if atom.stereochemistry:
                if atom.stereochemistry == "R":
                    desired_rdk_stereo_code = "R" # Yes, it's just a string
                if atom.stereochemistry == "S":
                    desired_rdk_stereo_code = "S"
                # Let's randomly assign this atom's stereo to CW
                rdatom.SetChiralTag(Chem.CHI_TETRAHEDRAL_CW)
                # We need to do force and cleanIt to recalculate CIP stereo
                Chem.AssignStereochemistry(rdmol, force=True, cleanIt=True)
                # If our random initial assignment worked, then we're set
                if rdatom.GetProp("_CIPCode") == desired_rdk_stereo_code:
                    continue
                # Otherwise, set it to CCW
                rdatom.SetChiralTag(Chem.CHI_TETRAHEDRAL_CCW)
                # We need to do force and cleanIt to recalculate CIP stereo
                Chem.AssignStereochemistry(rdmol, force=True, cleanIt=True)
                # Hopefully this worked, otherwise something's wrong
                if rdatom.GetProp("_CIPCode") == desired_rdk_stereo_code:
                    continue
                else:
                    raise Exception("Unknown atom stereochemistry encountered in "
                                    "to_rdkit. Desired stereochemistry: {}. Set stereochemistry {}".format(atom.stereochemistry,
                                                                                                           rdatom.GetProp("_CIPCode")))


        # Assign bond stereochemistry
        for bond in molecule.bonds:
            if bond.stereochemistry:
                # Determine neighbors
                # TODO: This API needs to be created
                n1 = [
                    n.molecule_atom_index for n in bond.atom1.bonded_atoms
                    if n != bond.atom2
                ][0]
                n2 = [
                    n.molecule_atom_index for n in bond.atom2.bonded_atoms
                    if n != bond.atom1
                ][0]
                # Get rdmol bonds
                bond_atom1_index = molecule.atoms.index(bond.atom1)
                bond_atom2_index = molecule.atoms.index(bond.atom2)
                bond1 = rdmol.GetBondBetweenAtoms(map_atoms[n1],
                                                  map_atoms[bond.atom1_index])
                bond2 = rdmol.GetBondBetweenAtoms(map_atoms[bond_atom1_index],
                                                  map_atoms[bond.atom2_index])
                bond3 = rdmol.GetBondBetweenAtoms(map_atoms[bond_atom2_index],
                                                  map_atoms[n2])
                # Set arbitrary stereochemistry
                # Since this is relative, the first bond always goes up
                # as explained above these names come from SMILES slashes so UP/UP is Trans and Up/Down is cis
                bond1.SetBondDir(Chem.BondDir.ENDUPRIGHT)
                bond3.SetBondDir(Chem.BondDir.ENDDOWNRIGHT)
                # Flip the stereochemistry if it is incorrect
                # TODO: Clean up _CIPCode atom and bond properties
                Chem.AssignStereochemistry(rdmol, cleanIt=True, force=True)
                if bond.stereochemistry == 'E':
                    desired_rdk_stereo_code = Chem.rdchem.BondStereo.STEREOE
                elif bond.stereochemistry == 'Z':
                    desired_rdk_stereo_code = Chem.rdchem.BondStereo.STEREOZ
                else:
                    raise Exception(
                        "Unknown bond stereochemistry encountered in "
                        "to_rdkit : {}".format(bond.stereochemistry))

                if bond2.GetStereo() != desired_rdk_stereo_code:
                    # Flip it
                    bond3.SetBondDir(Chem.BondDir.ENDUPRIGHT)
                    # Validate we have the right stereochemistry as a sanity check
                    Chem.AssignStereochemistry(rdmol, cleanIt=True, force=True)
                    #if rdmol.GetProp('_CIPCode') != bond.stereochemistry:
                    if bond2.GetStereo() != desired_rdk_stereo_code:
                        raise Exception(
                            'Programming error with assumptions about RDKit stereochemistry model'
                        )

        # Set coordinates if we have them
        if molecule._conformers:
            for conformer in molecule._conformers:
                rdmol_conformer = Chem.Conformer()
                for index, rd_idx in map_atoms.items():
                    (x, y, z) = conformer[index, :] / unit.angstrom
                    rdmol_conformer.SetAtomPosition(rd_idx,
                                                    Geometry.Point3D(x, y, z))
                rdmol.AddConformer(rdmol_conformer)

        # Retain charges, if present
        if not (molecule._partial_charges is None):

            rdk_indexed_charges = np.zeros((molecule.n_atoms), dtype=np.float)
            for off_idx, charge in enumerate(molecule._partial_charges):
                rdk_idx = map_atoms[off_idx]
                charge_unitless = charge / unit.elementary_charge
                rdk_indexed_charges[rdk_idx] = charge_unitless
            for rdk_idx, rdk_atom in enumerate(rdmol.GetAtoms()):
                rdk_atom.SetDoubleProp('partial_charge',
                                       rdk_indexed_charges[rdk_idx])

        # Cleanup the rdmol
        rdmol.UpdatePropertyCache(strict=False)
        Chem.GetSSSR(rdmol)
        # I added AssignStereochemistry which takes the directions of the bond set
        # and assigns the stereochemistry tags on the double bonds
        Chem.AssignStereochemistry(rdmol, force=False)

        # Return non-editable version
        return Chem.Mol(rdmol)

    @staticmethod
    def _find_smarts_matches(rdmol, smirks,
                             aromaticity_model='OEAroModel_MDL'):
        """Find all sets of atoms in the provided RDKit molecule that match the provided SMARTS string.

        .. warning :: This API experimental and subject to change.

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
            raise SMIRKSParsingError(
                'RDKit could not parse the SMIRKS string "{}"'.format(smirks))

        # Create atom mapping for query molecule
        idx_map = dict()
        for atom in qmol.GetAtoms():
            smirks_index = atom.GetAtomMapNum()
            if smirks_index != 0:
                idx_map[smirks_index - 1] = atom.GetIdx()
        map_list = [idx_map[x] for x in sorted(idx_map)]

        # Perform matching
        # TODO: The MoleculeImage mapping should preserve ordering of template molecule for equivalent atoms
        #       and speed matching for larger molecules.
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


class AmberToolsToolkitWrapper(ToolkitWrapper):
    """
    AmberTools toolkit wrapper

    """
    _toolkit_name = 'AmberTools'
    _toolkit_installation_instructions = 'The AmberTools toolkit (free and open source) can be found at ' \
                                         'https://anaconda.org/omnia/ambertools'
    _toolkit_file_read_formats = []
    _toolkit_file_write_formats = []

    @staticmethod
    def toolkit_is_available():
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

    @classmethod
    def is_available(cls):
        """
        Check whether this toolkit wrapper is available for use because the underlying toolkit can be found.

        .. note :: This method caches the result of any costly checks for file paths or module imports.

        Parameters
        ----------
        is_available : bool
            True if toolkit wrapper is available for use, False otherwise

        """
        if cls._is_available is None:
            cls._is_available = cls.toolkit_is_available()
        return cls._is_available

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
        import os
        from simtk import unit

        if charge_model is None:
            charge_model = 'bcc'

        # Check that the requested charge method is supported
        # Needs to be fixed: 'cm1', 'cm2',
        SUPPORTED_ANTECHAMBER_CHARGE_MODELS = ['gas', 'mul', 'bcc']
        if charge_model not in SUPPORTED_ANTECHAMBER_CHARGE_MODELS:
            raise ValueError(
                'Requested charge method {} not among supported charge '
                'methods {}'.format(charge_model,
                                    SUPPORTED_ANTECHAMBER_CHARGE_MODELS))

        # Find the path to antechamber
        # TODO: How should we implement find_executable?
        ANTECHAMBER_PATH = find_executable("antechamber")
        if ANTECHAMBER_PATH is None:
            raise (IOError("Antechamber not found, cannot run charge_mol()"))

        if len(molecule._conformers) == 0:
            raise Exception(
                "No conformers present in molecule submitted for partial charge calculation. Consider "
                "loading the molecule from a file with geometry already present or running "
                "molecule.generate_conformers() before calling molecule.compute_partial_charges"
            )


        # Compute charges
        from openforcefield.utils import temporary_directory, temporary_cd
        with temporary_directory() as tmpdir:
            with temporary_cd(tmpdir):
                net_charge = molecule.total_charge
                # Write out molecule in SDF format
                ## TODO: How should we handle multiple conformers?
                self._rdkit_toolkit_wrapper.to_file(
                    molecule, 'molecule.sdf', outfile_format='sdf')
                #os.system('ls')
                #os.system('cat molecule.sdf')
                # Compute desired charges
                # TODO: Add error handling if antechamber chokes
                # TODO: Add something cleaner than os.system
                os.system(
                    "antechamber -i molecule.sdf -fi sdf -o charged.mol2 -fo mol2 -pf "
                    "yes -c {} -nc {}".format(charge_model, net_charge))
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

        molecule.set_partial_charges(charges)


#=============================================================================================
# Toolkit registry
#=============================================================================================


class ToolkitRegistry(object):
    """
    Registry for ToolkitWrapper objects

    Examples
    --------

    Register toolkits in a specified order, skipping if unavailable

    >>> from openforcefield.utils.toolkits import ToolkitRegistry
    >>> toolkit_registry = ToolkitRegistry()
    >>> toolkit_precedence = [OpenEyeToolkitWrapper, RDKitToolkitWrapper, AmberToolsToolkitWrapper]
    >>> [ toolkit_registry.register(toolkit) for toolkit in toolkit_precedence if toolkit.is_available() ]

    Register specified toolkits, raising an exception if one is unavailable

    >>> toolkit_registry = ToolkitRegistry()
    >>> toolkits = [OpenEyeToolkitWrapper, AmberToolsToolkitWrapper]
    >>> for toolkit in toolkits: toolkit_registry.register(toolkit)

    Register all available toolkits in arbitrary order

    >>> from openforcefield.utils import all_subclasses
    >>> toolkits = all_subclasses(ToolkitWrapper)
    >>> [ toolkit_registry.register(toolkit) for toolkit in toolkits if toolkit.is_available() ]

    Retrieve the global singleton toolkit registry, which is created when this module is imported from all available
    toolkits:

    >>> from openforcefield.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY as toolkit_registry
    >>> print(toolkit_registry.registered_toolkits())

    """

    def __init__(self,
                 register_imported_toolkit_wrappers=False,
                 toolkit_precedence=None,
                 exception_if_unavailable=True):
        """
        Create an empty toolkit registry.

        .. warning :: This API is experimental and subject to change.

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

        if toolkit_precedence == None:
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
                         toolkit_wrapper_class,
                         exception_if_unavailable=True):
        """
        Register the provided toolkit wrapper class, instantiating an object of it.

        .. warning :: This API is experimental and subject to change.

        .. todo ::

           This method should raise an exception if the toolkit is unavailable, unless an optional argument
           is specified that silently avoids registration of toolkits that are unavailable.

        Parameters
        ----------
        toolkit_wrapper_class : subclass of ToolkitWrapper
            The class of the toolkit wrapper to register.
        exception_if_unavailable : bool, optional, default=True
            If True, an exception will be raised if the toolkit is unavailable

        """
        # TODO: Instantiate class if class, or just add if already instantiated.
        try:
            toolkit_wrapper = toolkit_wrapper_class()
            if not(toolkit_wrapper.is_available()):
                raise ToolkitUnavailableException()
            self._toolkits.append(toolkit_wrapper)
        except ToolkitUnavailableException as e:
            if exception_if_unavailable:
                raise e
            print("Unable to load toolkit {}.".format(toolkit_wrapper))

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

        .. warning :: This API is experimental and subject to change.

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

        >>> molecule = Molecule.from_smiles('Cc1ccccc1')
        >>> toolkit_registry = ToolkitRegistry(register_imported_toolkit_wrappers=True)
        >>> method = toolkit_registry.resolve('to_smiles')
        >>> smiles = method(molecule)

        .. todo :: Is there a better way to figure out which toolkits implement given methods by introspection?

        """
        for toolkit in self._toolkits:
            if hasattr(toolkit, method_name):
                method = getattr(toolkit, method_name)
                try:
                    return method
                except NotImplementedError as e:
                    pass

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

        .. warning :: This API is experimental and subject to change.

        ``*args`` and ``**kwargs`` are passed to the desired method, and return values of the method are returned

        This is a convenient shorthand for

        >>> toolkit_registry.resolve_method(method_name)(*args, **kwargs)

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

        >>> molecule = Molecule.from_smiles('Cc1ccccc1')
        >>> toolkit_registry = ToolkitRegistry(register_imported_toolkit_wrappers=True)
        >>> smiles = toolkit_registry.call('to_smiles', molecule)

        """
        # TODO: catch ValueError and compile list of methods that exist but rejected the specific parameters because they did not implement the requested methods

        value_errors = list()
        for toolkit in self._toolkits:
            if hasattr(toolkit, method_name):
                method = getattr(toolkit, method_name)
                #return method(*args, **kwargs)
                try:
                    return method(*args, **kwargs)
                except NotImplementedError as e:
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
        tk.toolkit_is_available()
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
