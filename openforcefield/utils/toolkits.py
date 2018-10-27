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
from distutils.spawn import find_executable
from simtk import unit
import numpy as np

#=============================================================================================
# SUPPORTED MODELS
#
# TODO: We may no longer need these since we now require SMIRNOFF to specify these models explicitly.
#=============================================================================================

DEFAULT_AROMATICITY_MODEL = 'OEAroModel_MDL' # TODO: Is there a more specific name and reference for the aromaticity model?
ALLOWED_AROMATICITY_MODELS = ['OEAroModel_MDL']

DEFAULT_FRACTIONAL_BONDORDER_MODEL = 'Wiberg' # TODO: Is there a more specific name and reference for the fractional bond order models?
ALLOWED_FRACTIONAL_BONDORDER_MODELS = ['Wiberg']

DEFAULT_CHARGE_MODEL = 'AM1-BCC' # TODO: Should this be `AM1-BCC`, or should we encode BCCs explicitly via AM1-CM2 preprocessing?
ALLOWED_CHARGE_MODELS = ['AM1-BCC'] # TODO: Which models do we want to support?



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



#=============================================================================================
# TOOLKIT UTILITY DECORATORS
#=============================================================================================



# TODO : Wrap toolkits in a much more modular way to make it easier to query their capabilities
## From Jeff: Maybe we just put these in the toolkit definitions themselves
#SUPPORTED_FILE_FORMATS = dict()
#SUPPORTED_FILE_FORMATS['OpenEye Toolkit'] = ['CAN', 'CDX', 'CSV', 'FASTA', 'INCHI', 'INCHIKEY', 'ISM', 'MDL', 'MF', 'MMOD', 'MOL2', 'MOL2H', 'MOPAC',
#                                     'OEB', 'PDB', 'RDF', 'SDF', 'SKC', 'SLN', 'SMI', 'USM', 'XYC']
#SUPPORTED_FILE_FORMATS['The RDKit'] = ['SDF', 'PDB', 'SMI', 'TDT'] # Don't put MOL2 in here -- RDKit can only handle corina format and fails on SYBYL
#SUPPORTED_FILE_FORMATS['AmberTools'] = ['MOL2']

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
    _is_available = None # True if toolkit is available
    _toolkit_name = None # Name of the toolkit
    _toolkit_installation_instructions = None # Installation instructions for the toolkit
    _toolkit_file_read_formats = None # The file types that this toolkit can read
    _toolkit_file_write_formats = None # The file types that this toolkit can write
    
    #@staticmethod
    ## From Jeff: This is confusing, but I changed things to make it run.
    ## Did I actually break it?
    # TODO: Right now, to access the class definition, I have to make this a classmethod
    # and thereby call it with () on the outermost decorator. Is this wasting time? Are we caching
    # the is_available results?
    @classmethod
    def requires_toolkit(cls): #remember cls is a ToolkitWrapper subclass here
        def decorator(func):
            @wraps(func)
            def wrapped_function(*args, **kwargs):
                if not cls.toolkit_is_available():
                    msg = 'This function requires the {} toolkit'.format(cls._toolkit_name)
                    raise LicenseError(msg)
                value = func(*args, **kwargs)
                return value
            return wrapped_function
        return decorator

    
    #def requires_toolkit(func):
    #    @wraps(func)
    #    def wrapper_decorator(*args, **kwargs):
    #        if not toolkit_is_available():
    #            msg = 'This function requires the {} toolkit'.format(_toolkit_name)
    #            raise LicenseError(msg)
    #        value = func(*args, **kwargs)
    #        return value
    #    return wrapped_function
    

    @property
    def toolkit_name(self):
        """
        The name of the toolkit wrapped by this class.
        """
        return self._toolkit_name

    @classmethod
    @property
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


    def from_file(self, filename, file_format):
        """
        Return an openforcefield.topology.Molecule from a file using this toolkit.
        
        Parameters
        ----------
        filename : str
            The file to read the molecule from
        file_format : str
            Format specifier, usually file suffix (eg. 'MOL2', 'SMI')
            Note that not all toolkits support all formats. Check ToolkitWrapper.toolkit_file_read_formats for details.
        
        Returns
        -------
        molecules : Molecule or list of Molecules
            a list of Molecule objects is returned.

        """
        return NotImplementedError
     
    def from_file_obj(self, file_obj, file_format):
        """
        Return an openforcefield.topology.Molecule from a file-like object (an object with a ".read()" method using this toolkit.
        
        Parameters
        ----------
        file_obj : file-like object
            The file-like object to read the molecule from
        file_format : str
            Format specifier, usually file suffix (eg. 'MOL2', 'SMI')
            Note that not all toolkits support all formats. Check ToolkitWrapper.toolkit_file_read_formats for details.
        
        Returns
        -------
        molecules : Molecule or list of Molecules
            a list of Molecule objects is returned.
        """
        return NotImplementedError
     
     
    #@staticmethod
    def to_smiles(self, molecule):
        """
        Return a canonical isomeric SMILES representation of the current molecule

        .. warning :: This API experimental and subject to change.

        .. todo :: Is this needed at the base class level?

        Parameters
        ----------
        molecule : Molecule
            The molecule for which canonical isomeric SMILES is to be computed

        Returns
        -------
        smiles : str
            Canonical isomeric explicit-hydrogen SMILES

        Examples
        --------

        >>> molecule = Molecule.from_smiles('Cc1ccccc1')
        >>> smiles = toolkit_wrapper.to_smiles(molecule)

        """
        raise NotImplementedError

    
    
    def from_smiles(self, smiles):
        """
        Create a Molecule object from SMILES

        .. warning :: This API experimental and subject to change.

        .. todo :: Is this needed at the base class level?

        Parameters
        ----------
        smiles : str
            SMILES string specifying the molecule

        Returns
        -------
        molecule : Molecule
            The molecule for which canonical isomeric SMILES is to be computed

        Examples
        --------

        >>> molecule = toolkit_wrapper.from_smiles('Cc1ccccc1')

        .. todo :: How is ambiguous stereochemistry and protonation states to be handled?

        """
        raise NotImplementedError

    def compute_partial_charges(self, molecule, charge_model="bcc"):
        """
        Compute partial charges

        .. warning :: This API experimental and subject to change.

        .. todo ::

           * Do we want to also allow ESP/RESP charges?

        Parameters
        ----------
        molecule : Molecule
            Molecule for which partial charges are to be computed
        charge_model : str, optional, default='bcc'
            The charge model to use. One of ['gas', 'mul', 'cm1', 'cm2', 'bcc']

        Returns
        -------
        charges : numpy.array of shape (natoms) of type float
            The partial charges

        Raises
        ------
        ValueError if the requested charge method could not be handled

        Notes
        -----
        Currently only sdf file supported as input and mol2 as output
        https://github.com/choderalab/openmoltools/blob/master/openmoltools/packmol.py

        """
        raise NotImplementedError

@inherit_docstrings
class OpenEyeToolkitWrapper(ToolkitWrapper):
    """
    OpenEye toolkit wrapper
    """
    _toolkit_name = 'OpenEye Toolkit'
    _toolkit_installation_instructions = 'The OpenEye toolkit requires a (free for academics) license, and can be found at: https://docs.eyesopen.com/toolkits/python/quickstart-python/install.html'
    _toolkit_file_read_formats =  ['CAN', 'CDX', 'CSV', 'FASTA', 'INCHI', 'INCHIKEY', 'ISM',
                                   'MDL', 'MF', 'MMOD', 'MOL2', 'MOL2H', 'MOPAC',
                                   'OEB', 'PDB', 'RDF', 'SDF', 'SKC', 'SLN', 'SMI', 'USM', 'XYC']
    _toolkit_file_write_formats = ['CAN', 'CDX', 'CSV', 'FASTA', 'INCHI', 'INCHIKEY', 'ISM',
                                   'MDL', 'MF', 'MMOD', 'MOL2', 'MOL2H', 'MOPAC',
                                   'OEB', 'PDB', 'RDF', 'SDF', 'SKC', 'SLN', 'SMI', 'USM', 'XYC']


    @staticmethod
    def toolkit_is_available(oetools=('oechem', 'oequacpac', 'oeiupac', 'oeomega')):
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
            cls._is_available = cls.toolkit_is_available(oetools=('oechem', 'oequacpac'))
        return cls._is_available

    def from_file(self, filename, file_format):
        """
        Return an openforcefield.topology.Molecule from a file using this toolkit.
        
        Parameters
        ----------
        filename : str
            The file to read the molecule from
        file_format : str
            Format specifier, usually file suffix (eg. 'MOL2', 'SMI')
            Note that not all toolkits support all formats. Check ToolkitWrapper.toolkit_file_read_formats for details.
        
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
            mol = Molecule.from_openeye(oemol)
            mols.append(mol)
        return mols
     
    def from_file_obj(self, file_obj, file_format):
        """
        Return an openforcefield.topology.Molecule from a file-like object (an object with a ".read()" method using this toolkit.
        
        Parameters
        ----------
        file_obj : file-like object
            The file-like object to read the molecule from
        file_format : str
            Format specifier, usually file suffix (eg. 'MOL2', 'SMI')
            Note that not all toolkits support all formats. Check ToolkitWrapper.toolkit_file_read_formats for details.
        
        Returns
        -------
        molecules : Molecule or list of Molecules
            a list of Molecule objects is returned.

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
            mol = Molecule.from_openeye(oemol)
            mols.append(mol)
        return mols
    
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

        cip = oechem.OEPerceiveCIPStereo(mol, bond)

        if cip == oechem.OECIPBondStereo_E:
            return 'E'
        elif cip == oechem.OECIPBondStereo_Z:
            return 'Z'
        elif cip == oechem.OECIPBondStereo_NotStereo:
            return None

    @staticmethod
    def from_openeye(oemol):
        """
        Create a Molecule from an OpenEye molecule.

        Parameters
        ----------
        oemol : openeye.oechem.OEMol
            An OpenEye molecule

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
        #from openforcefield.utils.toolkits.OpenEyeToolkitWrapper import _openeye_cip_atom_stereochemistry, openeye_cip_bond_stereochemistry

        # TODO: Decide if this is where we want to add explicit hydrogens
        result = oechem.OEAddExplicitHydrogens(oemol)
        if result == False:
            raise Exception("Addition of explicit hydrogens failed in from_openeye")
        
        # TODO: What other information should we preserve besides name?
        # TODO: How should we preserve the name?
        
        molecule = Molecule()
        molecule._name = oemol.GetTitle()
        


        # Copy any attached SD tag information
        # TODO: Should we use an API for this?
        molecule._properties = dict()
        for dp in oechem.OEGetSDDataPairs(oemol):
            molecule._properties[dp.GetTag()] = dp.GetValue()

        map_atoms = dict() # {oemol_idx: molecule_idx}
        for oeatom in oemol.GetAtoms():
            oe_idx = oeatom.GetIdx()
            atomic_number = oeatom.GetAtomicNum()
            formal_charge = oeatom.GetFormalCharge()
            is_aromatic = oeatom.IsAromatic()
            stereochemistry = OpenEyeToolkitWrapper._openeye_cip_atom_stereochemistry(oemol, oeatom)
            #stereochemistry = self._openeye_cip_atom_stereochemistry(oemol, oeatom)
            atom_index = molecule.add_atom(atomic_number, formal_charge, is_aromatic, stereochemistry=stereochemistry)
            map_atoms[oe_idx] = atom_index # store for mapping oeatom to molecule atom indices below

        for oebond in oemol.GetBonds():
            atom1_index = map_atoms[oebond.GetBgnIdx()]
            atom2_index = map_atoms[oebond.GetEndIdx()]
            bond_order = oebond.GetOrder()
            is_aromatic = oebond.IsAromatic()
            stereochemistry = OpenEyeToolkitWrapper._openeye_cip_bond_stereochemistry(oemol, oebond)
            #stereochemistry = self._openeye_cip_bond_stereochemistry(oemol, oebond)
            molecule.add_bond(atom1_index, atom2_index, bond_order, is_aromatic=is_aromatic, stereochemistry=stereochemistry)

        # TODO: Copy conformations, if present
        # TODO: Come up with some scheme to know when to import coordinates
        # From SMILES: no
        # From MOL2: maybe
        # From other: maybe
        if hasattr(oemol,'GetConfs'):
            for conf in oemol.GetConfs():
                n_atoms = molecule.n_atoms
                positions = unit.Quantity(np.zeros([n_atoms, 3], np.float), unit.angstrom)
                for oe_id in conf.GetCoords().keys():
                    off_atom_coords = unit.Quantity(conf.GetCoords()[oe_id], unit.angstrom)
                    off_atom_index = map_atoms[oe_id]
                    positions[off_atom_index,:] = off_atom_coords
                if (positions == 0*unit.angstrom).all():
                    continue
                molecule.add_conformer(positions)
                
        ## TODO: Partial charges
                
        return molecule

    # TODO: We could make this a staticmethod. It seems to have formerly belonged to
    # the Molecule class
    @staticmethod
    def to_openeye(molecule, aromaticity_model=DEFAULT_AROMATICITY_MODEL):
        """
        Create an OpenEye molecule using the specified aromaticity model

        .. todo ::

           * Use stored conformer positions instead of an argument.
           * Should the aromaticity model be specified in some other way?

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
        #from openforcefield.utils.toolkits import openeye_cip_atom_stereochemistry, openeye_cip_bond_stereochemistry

        oemol = oechem.OEMol()
        map_atoms = {} # {molecule_index : rdkit_index}
        # Add atoms
        oemol_atoms = list() # list of corresponding oemol atoms
        for atom in molecule.atoms:
            oeatom = oemol.NewAtom(atom.atomic_number)
            oeatom.SetFormalCharge(atom.formal_charge)
            oeatom.SetAromatic(atom.is_aromatic)
            oemol_atoms.append(oeatom)
            map_atoms[atom.molecule_atom_index] = oeatom.GetIdx()

        # Add bonds
        oemol_bonds = list() # list of corresponding oemol bonds
        for bond in molecule.bonds:
            #atom1_index = molecule.atoms.index(bond.atom1)
            #atom2_index = molecule.atoms.index(bond.atom2)
            atom1_index = bond.atom1_index
            atom2_index = bond.atom2_index
            oebond = oemol.NewBond(oemol_atoms[atom1_index], oemol_atoms[atom2_index])
            oebond.SetOrder(bond.bond_order)
            oebond.SetAromatic(bond.is_aromatic)
            oemol_bonds.append(oebond)

        # Set atom stereochemistry now that all connectivity is in place
        for atom, oeatom in zip(molecule.atoms, oemol_atoms):
            if not atom.stereochemistry:
                continue

            # Set arbitrary initial stereochemistry
            neighs = [n for n in oeatom.GetAtoms()]
            oeatom.SetStereo(neighs, oechem.OEAtomStereo_Tetra, oechem.OEAtomStereo_Right)

            # Flip chirality if stereochemistry is incorrect
            oeatom_stereochemistry = OpenEyeToolkitWrapper._openeye_cip_atom_stereochemistry(oemol, oeatom)
            if oeatom_stereochemistry != atom.stereochemistry:
                # Flip the stereochemistry
                oea.SetStereo(neighs, oechem.OEAtomStereo_Tetra, oechem.OEAtomStereo_Left)
                # Verify it matches now as a sanity check
                oeatom_stereochemistry = OpenEyeToolkitWrapper._openeye_cip_atom_stereochemistry(oemol, oeatom)
                if oeatom_stereochemistry != atom.stereochemistry:
                    raise Exception('Programming error: OpenEye atom stereochemistry assumptions failed.')

        # Set bond stereochemistry
        for bond, oebond in zip(molecule.bonds, oemol_bonds):
            if not bond.stereochemistry:
                continue

            atom1_index = bond.molecule.atoms.index(bond.atom1)
            atom2_index = bond.molecule.atoms.index(bond.atom2)
            # Set arbitrary initial stereochemistry
            oeatom1, oeatom2 = oemol_atoms[atom1_index], oemol_atoms[atom2_index]
            oeatom1_neighbor = [n for n in oeatom1.GetAtoms()][0]
            oeatom2_neighbor = [n for n in oeatom2.GetAtoms()][0]
            oebond.SetStereo([oeatom1, oeatom2], oechem.OEBondStereo_CisTrans, oechem.OEBondStereo_Cis)

            # Flip stereochemistry if incorrect
            oebond_stereochemistry = OpenEyeToolkitWrapper._openeye_cip_bond_stereochemistry(oemol, oebond)
            if oebond_stereochemistry != bond.sterechemistry:
                # Flip the stereochemistry
                oebond.SetStereo([oeatom1, oeatom2], oechem.OEBondStereo_CisTrans, oechem.OEBondStereo_Trans)
                # Verify it matches now as a sanity check
                oebond_stereochemistry = OpenEyeToolkitWrapper._openeye_cip_bond_stereochemistry(oemol, oebond)
                if oebond_stereochemistry != bond.stereochemistry:
                    raise Exception('Programming error: OpenEye bond stereochemistry assumptions failed.')

        # TODO: Retain conformations, if present
        # TODO: Are atom indexing schemes preserved between OFF molecules and OE molecules?
        if molecule._conformers != None:
            oemol.DeleteConfs()
            for conf in molecule._conformers:
                # OE needs a 1 x (3*n_Atoms) double array as input
                flat_coords = np.zeros((oemol.NumAtoms()*3), dtype=np.float32)
                for index, oe_idx in map_atoms.items():
                    (x,y,z) = conf[index,:] / unit.angstrom
                    flat_coords[(3*oe_idx)] = x
                    flat_coords[(3*oe_idx)+1] = y
                    flat_coords[(3*oe_idx)+2] = z
                    
                # TODO: Do we need to do these internal unit checks?
                # TODO: Is there any risk that the atom indexing systems will change?
                #flat_coords = (conf.in_units_of(unit.angstrom) / unit.angstrom).flatten()
                oecoords = oechem.OEFloatArray(flat_coords)
                oemol.NewConf(oecoords)
        
        # TODO: Retain name and properties, if present
        
        # Clean Up phase
        # The only feature of a molecule that wasn't perceived above seemed to be ring connectivity, better to run it
        # here then for someone to inquire about ring sizes and get 0 when it shouldn't be
        oechem.OEFindRingAtomsAndBonds(oemol)

        return oemol

    @staticmethod
    def to_smiles(molecule):
        # inherits base class docstring
        from openeye import oechem
        oemol = OpenEyeToolkitWrapper.to_openeye(molecule)
        smiles = oechem.OECreateSmiString(oemol,
                                          oechem.OESMILESFlag_DEFAULT |
                                          oechem.OESMILESFlag_Hydrogens |
                                          oechem.OESMILESFlag_Isotopes)
        return smiles

    
    def from_smiles(self, smiles):
        from openeye import oechem
        oemol = oechem.OEGraphMol()
        oechem.OESmilesToMol(oemol, smiles)
        molecule = self.from_openeye(oemol)
        return molecule
    
    def generate_conformers(self, molecule, clear_existing=True):
        """
        Generate molecule conformers using OpenEye Omega. 

        .. todo ::
        
           * which parameters should we expose? (or can we implement a general system with **kwargs?)
           * will the coordinates be returned in the OpenFF Molecule's own indexing system? Or is there a chance that they'll get reindexed when we convert the input into an OEmol?
        
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
        omega.SetEnergyWindow(15.0) #unit?
        omega.SetRMSThreshold(1.0)
        #Don't generate random stereoisomer if not specified
        omega.SetStrictStereo(True) 
        status = omega(oemol)

        molecule2 = self.from_openeye(oemol)

        if clear_existing:
            molecule._conformers = list()
        
        for conformer in molecule2._conformers:
            molecule.add_conformer(conformer)
        
    def compute_partial_charges(self, molecule, charge_model="am1bcc"):
        """
        Compute partial charges with OpenEye quacpac

        .. warning :: This API experimental and subject to change.

        .. todo ::

           * Should the default be ELF?
           * Can we expose more charge models?

        Parameters
        ----------
        molecule : Molecule
            Molecule for which partial charges are to be computed
        charge_model : str, optional, default='bcc'
            The charge model to use. One of ['noop', 'mmff', 'mmff94', 'am1bcc', 'am1bccnosymspt', 'amber', 'amberff94', 'am1bccelf10']

        Returns
        -------
        charges : numpy.array of shape (natoms) of type float
            The partial charges

        """
        from openeye import oequacpac
        from openeye import oeomega
        from openeye import oechem
        import numpy as np
        oemol = molecule.to_openeye()
        
        ## This seems like a big decision. Implemented a simple solution here. Not to be considered final.
        ## Some discussion at https://github.com/openforcefield/openforcefield/pull/86#issuecomment-350111236
        ## The following code is taken from the just-openeye version of the openforcefield repo https://github.com/openforcefield/openforcefield/blob/65f6b45954bde02c6cec1059661635c53a7f4e35/openforcefield/typing/engines/smirnoff/forcefield.py#L850
        
        if charge_model == "noop":
            result = oequacpac.OEAssignCharges(oemol, oequacpac.OEChargeEngineNoOp())
        elif charge_model == "mmff" or charge_model == "mmff94":
            result = oequacpac.OEAssignCharges(oemol, oequacpac.OEMMFF94Charges())
        elif charge_model == "am1bcc":
            result = oequacpac.OEAssignCharges(oemol, oequacpac.OEAM1BCCCharges())
        elif charge_model == "am1bccnosymspt":
            optimize = True
            symmetrize = True
            result = oequacpac.OEAssignCharges(oemol, oequacpac.OEAM1BCCCharges(not optimize, not symmetrize))
        elif charge_model == "amber" or charge_model == "amberff94":
            result = oequacpac.OEAssignCharges(oemol, oequacpac.OEAmberFF94Charges())
        elif charge_model == "am1bccelf10":
            result = oequacpac.OEAssignCharges(oemol, oequacpac.OEAM1BCCELF10Charges())
        else:
            raise ValueError('charge_model {} unknown'.format(charge_model))

        if result is False:
            raise Exception('Unable to assign charges')

        # Extract and return charges
        ## TODO: Make sure this can handle multiple conformations
        ## TODO: Make sure atom mapping remains constant
        charges = unit.Quantity(np.zeros([oemol.NumAtoms()], np.float64),
                                unit.elementary_charge)
        for index, atom in enumerate(oemol.GetAtoms()):
            charge = atom.GetPartialCharge()
            charge = charge * unit.elementary_charge
            charges[index] = charge
        molecule.set_partial_charges(charges)
        

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
                    oearomodel = getattr(oechem, 'OEAroModel_' + aromaticity_model)
                else:
                    raise ValueError("Error: provided aromaticity model not recognized by oechem.")
            else:
                raise ValueError("Error: provided aromaticity model must be a string.")

            # If aromaticity model was provided, prepare molecule
            oechem.OEClearAromaticFlags(mol)
            oechem.OEAssignAromaticFlags(mol, oearomodel)
            # Avoid running OEPrepareSearch or we lose desired aromaticity, so instead:
            oechem.OEAssignHybridization(mol)

        # Build list of matches
        # TODO: The MoleculeImage mapping should preserve ordering of template molecule for equivalent atoms
        #       and speed matching for larger molecules.
        unique = False # We require all matches, not just one of each kind
        substructure_search = oechem.OESubSearch(qmol)
        matches = list()
        for match in substructure_search.Match(mol, unique):
            # Compile list of atom indices that match the pattern tags
            atom_indices = dict()
            for matched_atom in match.GetAtoms():
                if matched_atom.pattern.GetMapIdx() != 0:
                    atom_indices[matched_atom.pattern.GetMapIdx()-1] = matched_atom.target.GetIdx()
            # Compress into list
            atom_indices = [ atom_indices[index] for index in range(len(atom_indices)) ]
            # Convert to tuple
            matches.append( tuple(atom_indices) )

        return matches

    def find_smarts_matches(self, molecule, smarts, aromaticity_model='OEAroModel_MDL'):
        """
        Find all SMARTS matches for the specified molecule, using the specified aromaticity model.

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
        #oemol = self.to_oemol(molecule, aromaticity_model=aromaticity_model)
        oemol = molecule.to_openeye()
        return self._find_smarts_matches(oemol, smarts)

class RDKitToolkitWrapper(ToolkitWrapper):
    """
    RDKit toolkit wrapper
    """
    _toolkit_name = 'The RDKit'
    _toolkit_installation_instructions = 'A conda-installable version of the free and open source RDKit cheminformatics toolkit can be found at: https://anaconda.org/rdkit/rdkit'
    _toolkit_file_read_formats = ['SDF', 'MOL', 'SMI'] #TODO: Add TDT support
    _toolkit_file_write_formats =['SDF', 'MOL', 'SMI', 'PDB'] 

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


    def from_file(self, filename, file_format):
        """
        Return an openforcefield.topology.Molecule from a file using this toolkit.
        
        Parameters
        ----------
        filename : str
            The file to read the molecule from
        file_format : str
            Format specifier, usually file suffix (eg. 'MOL2', 'SMI')
            Note that not all toolkits support all formats. Check ToolkitWrapper.toolkit_file_read_formats for details.
        
        Returns
        -------
        molecules : list of Molecules
            a list of Molecule objects is returned.

        """
        from openforcefield.topology import Molecule
        from rdkit import Chem
        mols = list()
        if (file_format == 'MOL') or (file_format == 'SDF'):
            for rdmol in Chem.SupplierFromFilename(filename, removeHs=False):
                mol = Molecule.from_rdkit(rdmol)
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
            raise Exception("RDKit can not safely read PDBs on their own. Information about bond order and aromaticity is likely to be lost.")
            # TODO: See if we can implement PDB+mol/smi combinations to get complete bond information.
            # https://github.com/openforcefield/openforcefield/issues/121
            rdmol = Chem.MolFromPDBFile(filename, removeHs=False)
            mol = Molecule.from_rdkit(rdmol)
            mols.append(mol)
            # TODO: Add SMI, TDT(?) support
            
        return mols
    def from_file_obj(self, file_obj, file_format):
        """
        Return an openforcefield.topology.Molecule from a file-like object (an object with a ".read()" method using this toolkit.
        
        Parameters
        ----------
        file_obj : file-like object
            The file-like object to read the molecule from
        file_format : str
            Format specifier, usually file suffix (eg. 'MOL2', 'SMI')
            Note that not all toolkits support all formats. Check ToolkitWrapper.toolkit_file_read_formats for details.
        
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
            raise Exception("RDKit can not safely read PDBs on their own. Information about bond order and aromaticity is likely to be lost.")
            # TODO: See if we can implement PDB+mol/smi combinations to get complete bond information.
             # https://github.com/openforcefield/openforcefield/issues/121
            rdmol = Chem.MolFromPDBBlock(file_data)
            mol = Molecule.from_rdkit(rdmol)
            mols.append(mol)
        # TODO: TDT file support
        return mols




    
    @classmethod
    def to_smiles(cls, molecule):
        # inherits base class docstring
        from rdkit import Chem
        rdmol = cls.to_rdkit(molecule)
        #rdmol = Chem.RemoveHs(rdmol)
        return Chem.MolToSmiles(rdmol, isomericSmiles=True, allHsExplicit=True)

    def from_smiles(self, smiles):
        from openforcefield.topology.molecule import Molecule
        # inherits base class docstring
        from rdkit import Chem
        rdmol = Chem.MolFromSmiles(smiles)
        # Add explicit hydrogens if they aren't there already
        rdmol = Chem.AddHs(rdmol)

        molecule = Molecule.from_rdkit(rdmol)
        
        return molecule
    
    
    def generate_conformers(self, molecule, clear_existing=True):
        """
        Generate molecule conformers using RDKit. 

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
        AllChem.EmbedMultipleConfs(rdmol,
                                   numConfs=800,
                                   pruneRmsThresh=1.0,
                                   randomSeed=1,
                                   #params=AllChem.ETKDG()
                                  )
        molecule2 = self.from_rdkit(rdmol)
        
        if clear_existing:
            molecule._conformers = list()
        
        for conformer in molecule2._conformers:
            molecule.add_conformer(conformer)
        

        
    
    
    def from_rdkit(self, rdmol):
        """
        Create a Molecule from an RDKit molecule.

        Requires the RDKit to be installed.

        Parameters
        ----------
        rdmol : rkit.RDMol
            An RDKit molecule

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

        # Create a new openforcefield Molecule
        mol = Molecule()
        
        # These checks cause rdkit to choke on one member of our test set: ZINC16448882
        # http://zinc.docking.org/substance/16448882
        # This has a pentavalent nitrogen, which I think is really resonance-stabilized.
        # I think we should allow this as input, since a fractional bond order calculation will probably sort it out.
        #Chem.SanitizeMol(rdmol, Chem.SANITIZE_ALL^Chem.SANITIZE_SETAROMATICITY)
        #Chem.SetAromaticity(rdmol, Chem.AromaticityModel.AROMATICITY_MDL)
        # If RDMol has a title save it
        if rdmol.HasProp("_Name"):
            mol.name == rdmol.GetProp("_Name")
        else:
            mol.name = None

        # Store all properties
        # TODO: Should Title or _Name be a special property?
        # TODO: Should there be an API point for storing properties?
        properties = rdmol.GetPropsAsDict()
        mol._properties = properties

        # We store bond orders as integers regardless of aromaticity.
        # In order to properly extract these, we need to have the "Kekulized" version of the rdkit mol
        kekul_mol = Chem.Mol(rdmol)
        Chem.Kekulize(kekul_mol, True)

        # setting chirality in openeye requires using neighbor atoms
        # therefore we can't do it until after the atoms and bonds are all added
        chiral_atoms = dict() # {rd_idx: openeye chirality}
        map_atoms = {}
        for rda in rdmol.GetAtoms():
            rd_idx = rda.GetIdx()

            # create a new atom
            #atomic_number = oemol.NewAtom(rda.GetAtomicNum())
            atomic_number = rda.GetAtomicNum()
            formal_charge = rda.GetFormalCharge()
            is_aromatic = rda.GetIsAromatic()

            # If chiral, store the chirality to be set later
            stereochemistry = None
            tag = rda.GetChiralTag()
            if tag == Chem.CHI_TETRAHEDRAL_CCW:
                stereochemistry = 'R'
            if tag == Chem.CHI_TETRAHEDRAL_CW:
                stereochemistry = 'S'
            
            atom_index = mol.add_atom(atomic_number, formal_charge, is_aromatic, stereochemistry=stereochemistry)
            map_atoms[rd_idx] = atom_index

        # Similar to chirality, stereochemistry of bonds in OE is set relative to their neighbors
        stereo_bonds = list()
        # stereo_bonds stores tuples in the form (oe_bond, rd_idx1, rd_idx2, OE stereo specification)
        # where rd_idx1 and 2 are the atoms on the outside of the bond
        # i.e. Cl and F in the example above
        aro_bond = 0
        for rdb in rdmol.GetBonds():
            a1 = rdb.GetBeginAtomIdx()
            a2 = rdb.GetEndAtomIdx()

            # Determine bond aromaticity and Kekulized bond order
            is_aromatic = False
            order = rdb.GetBondTypeAsDouble()
            if order == 1.5:
                # get the bond order for this bond in the kekulized molecule
                order = kekul_mol.GetBondWithIdx(rdb.GetIdx()).GetBondTypeAsDouble()
                is_aromatic = True
            # Convert floating-point bond order to integral bond order
            order = int(order)

            # determine if stereochemistry is needed
            stereochemistry = None
            tag = rdb.GetStereo()
            if tag == Chem.BondStereo.STEREOCIS or tag == Chem.BondStereo.STEREOZ:
                stereochemistry = 'Z'
            if tag == Chem.BondStereo.STEREOTRANS or tag == Chem.BondStereo.STEREOE:
                stereochemistry = 'E'

            # create a new bond
            bond_index = mol.add_bond(map_atoms[a1], map_atoms[a2], order, is_aromatic, stereochemistry=stereochemistry)

        # TODO: Save conformer(s), if present
        # If the rdmol has a conformer, store its coordinates
        if len(rdmol.GetConformers()) != 0:
            for conf in rdmol.GetConformers():
                n_atoms = mol.n_atoms
                # TODO: Will this always be angstrom when loading from RDKit?
                positions = unit.Quantity(np.zeros((n_atoms,3)), unit.angstrom)
                for rd_idx, off_idx in map_atoms.items():
                    atom_coords = conf.GetPositions()[rd_idx,:] * unit.angstrom
                    positions[off_idx,:] = atom_coords
                mol.add_conformer(positions)
            

        return mol

    @staticmethod
    def to_rdkit(molecule, aromaticity_model=DEFAULT_AROMATICITY_MODEL):
        """
        Create an RDKit molecule

        Requires the RDKit to be installed.

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
        if not(molecule.name == None):
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

        _bondtypes = {1: Chem.BondType.SINGLE,
                      1.5: Chem.BondType.AROMATIC,
                      2: Chem.BondType.DOUBLE,
                      3: Chem.BondType.TRIPLE,
                      4: Chem.BondType.QUADRUPLE,
                      5: Chem.BondType.QUINTUPLE,
                      6: Chem.BondType.HEXTUPLE,
                      7: Chem.BondType.ONEANDAHALF,}

        # atom map lets you find atoms again
        map_atoms = dict() # { molecule index : rdkit index }
        for index, atom in enumerate(molecule.atoms):
            rdatom = Chem.Atom(atom.atomic_number)
            rdatom.SetFormalCharge(atom.formal_charge)
            rdatom.SetIsAromatic(atom.is_aromatic)

            if atom.stereochemistry == 'S':
                rdatom.SetChiralTag(Chem.CHI_TETRAHEDRAL_CW)
            elif atom.stereochemistry == 'R':
                rdatom.SetChiralTag(Chem.CHI_TETRAHEDRAL_CCW)
                
            rd_index = rdmol.AddAtom(rdatom)
                
            map_atoms[index] = rd_index

        for bond in molecule.bonds:
            ## TODO: Horribly inefficient lookup. Should get a better atom index method
            rdatom1 = map_atoms[molecule.atoms.index(bond.atom1)]
            rdatom2 = map_atoms[molecule.atoms.index(bond.atom2)]
            rdmol.AddBond(rdatom1, rdatom2)
            rdbond = rdmol.GetBondBetweenAtoms(rdatom1, rdatom2)

            # Assign bond type, which is based on order unless it is aromatic
            if bond.is_aromatic:
                rdbond.SetBondType(_bondtypes[1.5])
                rdbond.SetIsAromatic(True)
            else:
                rdbond.SetBondType(_bondtypes[bond.bond_order])
                rdbond.SetIsAromatic(False)

        # Assign bond stereochemistry
        for bond in molecule.bonds:
            if bond.stereochemistry:
                # Determine neighbors
                # TODO: This API needs to be created
                n1 = [n.index for n in bond.atom1.bonded_atoms if n != bond.atom2][0]
                n2 = [n.index for n in bond.atom2.bonded_atoms if n != bond.atom1][0]
                # Get rdmol bonds
                bond_atom1_index = molecule.atoms.index(bond.atom1)
                bond_atom2_index = molecule.atoms.index(bond.atom2)
                bond1 = rdmol.GetBondBetweenAtoms(map_atoms[n1], map_atoms[bond.atom1_index])
                bond2 = rdmol.GetBondBetweenAtoms(map_atoms[bond_atom1_index], map_atoms[bon_.atom2_index])
                bond3 = rdmol.GetBondBetweenAtoms(map_atoms[bond_atom2_index], map_atoms[n2])
                # Set arbitrary stereochemistry
                # Since this is relative, the first bond always goes up
                # as explained above these names come from SMILES slashes so UP/UP is Trans and Up/Down is cis
                bond1.SetBondDir(Chem.BondDir.ENDUPRIGHT)
                bond3.SetBondDir(Chem.BondDir.ENDDOWNRIGHT)
                # Flip the stereochemistry if it is incorrect
                # TODO: Clean up _CIPCode atom and bond properties
                Chem.AssignStereochemistry(rdmol, cleanIt=True, force=True)
                if rdmol.GetProp('_CIPCode') != bond.stereochemistry:
                    # Flip it
                    bond3.SetBondDir(Chem.BondDir.ENDUPRIGHT)
                    # Validate we have the right stereochemistry as a sanity check
                    Chem.AssignStereochemistry(rdmol, cleanIt=True, force=True)
                    if rdmol.GetProp('_CIPCode') != bond.stereochemistry:
                        raise Exception('Programming error with assumptions about RDKit stereochemistry model')

        # Set coordinates if we have them
        if molecule._conformers:
            for conformer in molecule._conformers:
                rdmol_conformer = Chem.Conformer()
                for index, rd_idx in map_atoms.items():
                    (x,y,z) = conformer[index,:] / unit.angstrom
                    rdmol_conformer.SetAtomPosition(rd_idx, Geometry.Point3D(x,y,z))
                rdmol.AddConformer(rdmol_conformer)

        # Cleanup the rdmol
        # Note I copied UpdatePropertyCache and GetSSSR from Shuzhe's code to convert oemol to rdmol here:
        rdmol.UpdatePropertyCache(strict=False)
        Chem.GetSSSR(rdmol)
        # I added AssignStereochemistry which takes the directions of the bond set
        # and assigns the stereochemistry tags on the double bonds
        Chem.AssignStereochemistry(rdmol, force=False)

        # Return non-editable version
        return Chem.Mol(rdmol)

    @staticmethod
    def _find_smarts_matches(rdmol, smirks, aromaticity_model='OEAroModel_MDL'):
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
            Chem.SanitizeMol(mol, Chem.SANITIZE_ALL^Chem.SANITIZE_SETAROMATICITY)
            Chem.SetAromaticity(mol, Chem.AromaticityModel.AROMATICITY_MDL)
        else:
            # Only the OEAroModel_MDL is supported for now
            raise ValueError('Unknown aromaticity model: {}'.aromaticity_models)

        # Set up query.
        qmol = Chem.MolFromSmarts(smirks)   #cannot catch the error
        if qmol is None:
            raise SMIRKSParsingError('RDKit could not parse the SMIRKS string "{}"'.format(smirks))

        # Create atom mapping for query molecule
        index_map = dict()
        for atom in qmol.GetAtoms():
             smirks_index = atom.GetAtomMapNum()
             if smirks_index != 0:
                ind_map[smirks_index - 1] = atom.GetIdx()
        map_list = [ index_map[x] for x in sorted(index_map) ]

        # Perform matching
        # TODO: The MoleculeImage mapping should preserve ordering of template molecule for equivalent atoms
        #       and speed matching for larger molecules.
        matches = list()
        for match in rdmol.GetSubstructMatches(qmol, uniquify=False):
            mas = [ match[x] for x in map_list ]
            matches.append(tuple(mas))

        return matches

    def find_smarts_matches(self, molecule, smarts, aromaticity_model='OEAroModel_MDL'):
        """
        Find all SMARTS matches for the specified molecule, using the specified aromaticity model.

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
        rdmol = self.to_rdmol(molecule, aromaticity_model=aromaticity_model)
        return _find_smarts_matches(rdmol, smarts, aromaticity_model='OEAroModel_MDL')

class AmberToolsToolkitWrapper(ToolkitWrapper):
    """
    AmberTools toolkit wrapper

    """
    _toolkit_name = 'AmberTools'
    _toolkit_installation_instructions = 'The AmberTools toolkit (free and open source) can be found at https://anaconda.org/omnia/ambertools'
    _toolkit_file_read_formats = [] 
    _toolkit_file_write_formats = []

    @staticmethod
    def toolkit_is_available():
        """
        Check whether the AmberTools toolkit is installed

        Returns
        -------
        is_installed : bool
            True if RDKit is installed, False otherwise.

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
        pass

    def compute_partial_charges(self, molecule, charge_model="bcc"):
        """
        Compute partial charges with AmberTools using antechamber/sqm

        .. warning :: This API experimental and subject to change.

        .. todo ::

           * Do we want to also allow ESP/RESP charges?

        Parameters
        ----------
        molecule : Molecule
            Molecule for which partial charges are to be computed
        charge_model : str, optional, default='bcc'
            The charge model to use. One of ['gas', 'mul', 'cm1', 'cm2', 'bcc']

        Returns
        -------
        charges : numpy.array of shape (natoms) of type float
            The partial charges

        Raises
        ------
        ValueError if the requested charge method could not be handled

        Notes
        -----
        Currently only sdf file supported as input and mol2 as output
        https://github.com/choderalab/openmoltools/blob/master/openmoltools/packmol.py

        """
        import os
        # Check that the requested charge method is supported
        SUPPORTED_ANTECHAMBER_CHARGE_MODELS = ['gas', 'mul', 'cm1', 'cm2', 'bcc']
        if charge_model not in SUPPORTED_ANTECHAMBER_CHARGE_MODELS:
            raise ValueError('Requested charge method {} not among supported charge methods {}'.format(charge_model, SUPPORTED_ANTECHAMBER_CHARGE_MODELS))

        # Find the path to antechamber
        # TODO: How should we implement find_executable?
        ANTECHAMBER_PATH = find_executable("antechamber")
        if ANTECHAMBER_PATH is None:
            raise(IOError("Antechamber not found, cannot run charge_mol()"))

        # Compute charges
        from openforcefield.utils import temporary_directory, temporary_cd
        with temporary_directory() as tmpdir:
            with temporary_cd(tmpdir):
                net_charge = molecule.total_charge
                # Write out molecule in SDF format
                ## TODO: Where will this get coordinates from?
                molecule.to_file('molecule.sdf', outfile_format='sdf')
                os.system('ls')
                os.system('cat molecule.sdf')
                # Compute desired charges
                # TODO: Add error handling if antechamber chokes
                # TODO: Add something cleaner than os.system
                os.system("antechamber -i molecule.sdf -fi sdf -o charged.mol2 -fo mol2 -pf yes -c {} -nc {}".format(charge_model, net_charge))
                # Write out just charges
                os.system("antechamber -i charges.mol2 -fi mol2 -o charges2.mol2 -fo mol2 -c wc -cf charges.txt -pf yes")
                # Read the charges
                with open('charges.txt', 'r') as infile:
                    contents = infile.read()
                text_charges = contents.split()
                charges = np.zeros([self.n_atoms], np.float64)
                for index, token in enumerate(text_charges):
                    charges[index] = float(token)

        # TODO: Read the charges back into the molecule?
        return charges

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

    Retrieve the global singleton toolkit registry, which is created when this module is imported from all available toolkits:

    >>> from openforcefield.utils.toolkits import DEFAULT_TOOLKIT_REGISTRY as toolkit_registry
    >>> print(toolkit_registry.registered_toolkits())

    """
    def __init__(self, register_imported_toolkit_wrappers=False, toolkit_precedence=None):
        """
        Create an empty toolkit registry.

        .. warning :: This API experimental and subject to change.

        Parameters
        ----------
        register_imported_toolkit_wrappers : bool, optional, default=False
            If True, will attempt to register all imported ToolkitWrapper subclasses that can be found, in no particular order.
        toolkit_precedence : list, optional, defailt=None
            List of toolkit wrapper classes, in order of desired precedence when performing molecule operations. If None, defaults to [OpenEyeToolkitWrapper, RDKitToolkitWrapper, AmberToolsToolkitWrapper]. 
        """

        self._toolkits = list()

        if toolkit_precedence == None:
            toolkit_precedence = [OpenEyeToolkitWrapper, RDKitToolkitWrapper, AmberToolsToolkitWrapper]

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
            self.register_toolkit(toolkit)

    @property
    def registered_toolkits(self):
        """
        List registered toolkits.

        .. warning :: This API experimental and subject to change.

        .. todo :: Should this return a generator? Deep copies? Classes? Toolkit names?

        """
        return list(self._toolkits)

    def register_toolkit(self, toolkit_wrapper_class, exception_if_unavailable=True):
        """
        Register the provided toolkit wrapper.

        .. warning :: This API experimental and subject to change.

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
            self._toolkits.append(toolkit_wrapper)
        except ToolkitUnavailableException as e:
            if exception_if_unavailable:
                raise e

    # TODO: Can we automatically resolve calls to methods that are not explicitly defined using some Python magic?

    def resolve(self, method_name):
        """
        Resolve the requested method name by checking all registered toolkits in order of precedence for one that provides the requested method.

        .. warning :: This API experimental and subject to change.

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
        msg = 'No registered toolkits can provide the capability "{}".\n'.format(method_name)
        msg += 'Available toolkits are: {}\n'.format(self.registered_toolkits)
        raise NotImplementedError(msg)

    # TODO: Can we instead register available methods directly with `ToolkitRegistry`, so we can just use `ToolkitRegistry.method()`?
    def call(self, method_name, *args, **kwargs):
        """
        Execute the requested method by attempting to use all registered toolkits in order of precedence.

        .. warning :: This API experimental and subject to change.

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
                    value_errors.append( (toolkit, value_error) )

        # No toolkit was found to provide the requested capability
        # TODO: Can we help developers by providing a check for typos in expected method names?
        msg = 'No registered toolkits can provide the capability "{}".\n'.format(method_name)
        msg += 'Available toolkits are: {}\n'.format(self.registered_toolkits)
        # Append information about toolkits that implemented the method, but could not handle the provided parameters
        for toolkit, value_error in value_errors:
            msg += ' {} : {}\n'.format(toolkit, value_error)
        raise NotImplementedError(msg)

#=============================================================================================
# GLOBAL TOOLKIT REGISTRY
#=============================================================================================

# Create global toolkit registry, where all available toolkits are registered

## From jeff: The commented-out functionality has been moved into the ToolkitRegistry constructor

GLOBAL_TOOLKIT_REGISTRY = ToolkitRegistry(register_imported_toolkit_wrappers=True)
#for toolkit in all_subclasses(ToolkitWrapper):
#    GLOBAL_TOOLKIT_REGISTRY.register_toolkit(toolkit, exception_if_unavailable=False)

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

if len(GLOBAL_TOOLKIT_REGISTRY.registered_toolkits) == 0:
    msg = 'WARNING: No cheminfomatics toolkits are available.\n'
    msg += 'Please install at least one of the following toolkits:\n'
    for wrapper in all_subclasses(ToolkitWrapper):
        if wrapper.toolkit_name is not None:
            msg += '{} : {}\n'.format(wrapper.toolkit_name, wrapper.installation_instructions)
    print(msg)
