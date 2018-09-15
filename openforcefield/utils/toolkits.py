#!/usr/bin/env python

"""
Wrapper classes for providing a minimal consistent interface to cheminformatics toolkits

Currently supported toolkits:

* The `OpenEye Toolkit <https://docs.eyesopen.com/toolkits/python/quickstart-python/index.html>`_
* The `RDKit <http://www.rdkit.org/>`_
* `AmberTools <http://ambermd.org/AmberTools.php>`_

.. todo::

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

#=============================================================================================
# CHEMINFORMATICS TOOLKITS
#=============================================================================================

DEFAULT_AROMATICITY_MODEL = 'OEAroModel_MDL' # TODO: Is there a more specific name and reference for the aromaticity model?
ALLOWED_AROMATICITY_MODELS = ['OEAroModel_MDL']

DEFAULT_FRACTIONAL_BONDORDER_MODEL = 'Wiberg' # TODO: Is there a more specific name and reference for the fractional bond order models?
ALLOWED_FRACTIONAL_BONDORDER_MODELS = ['Wiberg']

DEFAULT_CHARGE_MODEL = 'AM1-BCC' # TODO: Should this be `AM1-BCC`, or should we encode BCCs explicitly via AM1-CM2 preprocessing?
ALLOWED_CHARGE_MODELS = ['AM1-BCC'] # TODO: Which models do we want to support?

# Control the precedence order in which cheminformatics toolkits are used
TOOLKIT_PRECEDENCE = ['openeye', 'rdkit']

# List of supported toolkits and messages indicating how they can be installed
SUPPORTED_TOOLKITS = {
    'rdkit' : 'A conda-installable version of the free and open source RDKit cheminformatics toolkit can be found at: https://anaconda.org/rdkit/rdkit',
    'openeye' : 'The OpenEye toolkit requires a (free for academics) license, and can be found at: https://docs.eyesopen.com/toolkits/python/quickstart-python/install.html',
    'ambertools' : 'The AmberTools toolkit (free and open source) can be found at https://anaconda.org/omnia/ambertools'
}

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

# TODO : Wrap toolkits in a much more modular way to make it easier to query their capabilities
SUPPORTED_FILE_FORMATS = dict()
SUPPORTED_FILE_FORMATS['openeye'] = ['CAN', 'CDX', 'CSV', 'FASTA', 'INCHI', 'INCHIKEY', 'ISM', 'MDL', 'MF', 'MMOD', 'MOL2', 'MOL2H', 'MOPAC',
                                     'OEB', 'PDB', 'RDF', 'SDF', 'SKC', 'SLN', 'SMI', 'USM', 'XYC']
SUPPORTED_FILE_FORMATS['rdkit'] = ['SDF', 'PDB', 'SMI', 'TDT']

#=============================================================================================
# UTILITY FUNCTIONS
#=============================================================================================

from inspect import getmembers, isfunction

def inherit_docstrings(cls):
    for name, func in getmembers(cls, isfunction):
        if func.__doc__: continue
        for parent in cls.__mro__[1:]:
            if hasattr(parent, name):
                func.__doc__ = getattr(parent, name).__doc__
    return cls

#=============================================================================================
# CHEMINFORMATICS TOOLKIT WRAPPERS
#=============================================================================================

class ToolkitWrapper(object):
    """
    Toolkit wrapper base class.

    .. warning :: This API experimental and subject to change.
    """

    def compute_partial_charges(self, molecule, charge_model="bcc"):
        """
        Compute partial charges

        .. warning :: This API experimental and subject to change.

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

    @staticmethod
    def to_smiles(molecule):
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

    @staticmethod
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

        """
        raise NotImplementedError

@inherit_docstrings
class OpenEyeToolkitWrapper(object):
    """
    OpenEye toolkit wrapper
    """
    @staticmethod
    @requires_openeye('oechem')
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
    @requires_openeye('oechem')
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
    @requires_openeye('oechem') # TODO: Is this still needed?
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

        molecule = Molecule()

        # TODO: What other information should we preserve besides name?
        # TODO: How should we preserve the name?
        molecule.name = oemol.GetTitle()

        # Copy any attached SD tag information
        # TODO: Should we use an API for this?
        molecule.properties = dict()
        for dp in oechem.OEGetSDDataPairs(oemol):
            molecule[dp.GetTag()] = dp.GetValue()

        map_atoms = dict() # {oemol_idx: molecule_idx}
        for oeatom in oemol.GetAtoms():
            oe_idx = oeatom.GetIdx()
            atomic_number = oeatom.GetAtomicNum()
            formal_charge = oeatom.GetFormalCharge()
            is_aromatic = oeatom.IsAromatic()
            stereochemistry = openeye_cip_atom_stereochemistry(oemol, oeatom)
            atom_index = molecule.add_atom(atomic_number=atomic_number, formal_charge=formal_charge, is_aromatic=is_aromatic, stereochemistry=stereochemistry)
            map_atoms[oe_idx] = atom_index # store for mapping oeatom to molecule atom indices below

        for oebond in oemol.GetBonds():
            atom1_index = map_atoms[oebond.GetBgnIdx()]
            atom2_index = map_atoms[oebond.GetEndIdx()]
            order = oeb.GetOrder()
            is_aromatic = oeb.IsAromatic()
            stereochemistry = openeye_cip_bond_stereochemistry(oemol, oebond)
            molecule.add_bond(atom1_index, atom2_index, order=order, is_aromatic=is_aromatic, stereochemistry=stereochemistry)

        # TODO: Copy conformations, if present

        return molecule

    @staticmethod
    @requires_openeye('oechem') # TODO: Is this still needed?
    def to_openeye(molecule, aromaticity_model=DEFAULT_AROMATICITY_MODEL):
        """
        Create an OpenEye molecule using the specified aromaticity model

        .. todo ::

           * Use stored conformer positions instead of an argument.
           * Should the aromaticity model be specified in some other way?

        Parameters
        ----------
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
        from openforcefield.utils.toolkits import openeye_cip_atom_stereochemistry, openeye_cip_bond_stereochemistry

        oemol = oechem.OEMol()

        # Add atoms
        oemol_atoms = list() # list of corresponding oemol atoms
        for atom in self.atoms:
            oeatom = oemol.NewAtom(atom.atomic_number)
            oeatom.SetFormalCharge(atom.formal_charge)
            oeatom.SetAromatic(atom.is_aromatic)
            oemol_atoms.append(oeatom)

        # Add bonds
        oemol_bonds = list() # list of corresponding oemol bonds
        for bond in self.bonds:
            oebond = oemol.NewBond(oemol_atoms[bond.atom1_index], oemol_atoms[bond.atom2_index])
            newbond.SetOrder(bond.order)
            newbond.SetAromatic(bond.is_aromatic)
            oemol_bonds.append(oebond)

        # Set atom stereochemistry now that all connectivity is in place
        for atom, oeatom in zip(self.atoms, oemol_atoms):
            if not atom.stereochemistry:
                continue

            # Set arbitrary initial stereochemistry
            neighs = [n for n in oeatom.GetAtoms()]
            oeatom.SetStereo(neighs, oechem.OEAtomStereo_Tetra, oechem.OEAtomStereo_Right)

            # Flip chirality if stereochemistry is incorrect
            oeatom_stereochemistry = _openeye_cip_atom_stereochemistry(oemol, oeatom)
            if oeatom_stereochemistry != atom.sterechemistry:
                # Flip the stereochemistry
                oea.SetStereo(neighs, oechem.OEAtomStereo_Tetra, oechem.OEAtomStereo_Left)
                # Verify it matches now as a sanity check
                oeatom_stereochemistry = _openeye_cip_atom_stereochemistry(oemol, oeatom)
                if oeatom_stereochemistry != atom.stereochemistry:
                    raise Exception('Programming error: OpenEye atom stereochemistry assumptions failed.')

        # Set bond stereochemistry
        for bond, oebond in zip(self.atoms, oemol_bonds):
            if not bond.stereochemistry:
                continue

            # Set arbitrary initial stereochemistry
            oeatom1, oeatom2 = oemol_atoms[bond.atom1_index], oemol_atoms[bond.atom2_index]
            oeatom1_neighbor = [n for n in oeatom1.GetAtoms()][0]
            oeatom2_neighbor = [n for n in oeatom2.GetAtoms()][0]
            oebond.SetStereo([oeatom1, oeatom2], oechem.OEBondStereo_CisTrans, oechem.OEBondStereo_Cis)

            # Flip stereochemistry if incorrect
            oebond_stereochemistry = _openeye_cip_bond_stereochemistry(oemol, oebond)
            if oebond_stereochemistry != bond.sterechemistry:
                # Flip the stereochemistry
                oebond.SetStereo([oeatom1, oeatom2], oechem.OEBondStereo_CisTrans, oechem.OEBondStereo_Trans)
                # Verify it matches now as a sanity check
                oebond_stereochemistry = _openeye_cip_bond_stereochemistry(oemol, oebond)
                if oebond_stereochemistry != bond.stereochemistry:
                    raise Exception('Programming error: OpenEye bond stereochemistry assumptions failed.')

        # TODO: Retain conformations, if present

        # TODO: Retain name and properties, if present

        # Clean Up phase
        # The only feature of a molecule that wasn't perceived above seemed to be ring connectivity, better to run it
        # here then for someone to inquire about ring sizes and get 0 when it shouldn't be
        oechem.OEFindRingAtomsAndBonds(oemol)

        return oemol

    def to_smiles(self, molecule):
        # inherits base class docstring
        from openeye import oechem
        oemol = self.to_openeye(molecule)
        return oechem.OEMolToSmiles(oemol)

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
        oemol = molecule.to_openeye()

        result = False
        if name == "noop":
            result = oequacpac.OEAssignCharges(oemol, oequacpac.OEChargeEngineNoOp())
        elif name == "mmff" or name == "mmff94":
            result = oequacpac.OEAssignCharges(oemol, oequacpac.OEMMFF94Charges())
        elif name == "am1bcc":
            result = oequacpac.OEAssignCharges(oemol, oequacpac.OEAM1BCCCharges())
        elif name == "am1bccnosymspt":
            optimize = True
            symmetrize = True
            result = oequacpac.OEAssignCharges(mol, oequacpac.OEAM1BCCCharges(not optimize, not symmetrize))
        elif name == "amber" or name == "amberff94":
            result = oequacpac.OEAssignCharges(mol, oequacpac.OEAmberFF94Charges())
        elif name == "am1bccelf10":
            result = oequacpac.OEAssignCharges(mol, oequacpac.OEAM1BCCELF10Charges())
        else:
            raise ValueError('charge_model {} unknown'.format(charge_model))

        if result is False:
            raise Exception('Unable to assign charges')

        # Extract and return charges
        charges = np.zeros([oemol.NumAtoms()], np.float64)
        for index, atom in enumerate(oemol.GetAtoms()):
            charges[index] = atom.GetPartialCharge()
        return charges

class RDKitToolkitWrapper(ToolkitWrapper):
    """
    RDKit toolkit wrapper
    """
    def to_smiles(self, molecule):
        # inherits base class docstring
        from rdkit import Chem
        rdmol = self.to_rdkit(molecule)
        return Chem.MolToSmiles(rdmol, isomericSmiles=True)

    @staticmethod
    @requires_rdkit()
    def from_rdkit(rdmol):
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

        # Create a new openforcefield Molecule
        mol = Molecule()

        # If RDMol has a title save it
        if rdmol.HasProp("_Name"):
            self.name == rdmol.GetProp("_Name")

        # Store all properties
        # TODO: Should Title or _Name be a special property?
        # TODO: Should there be an API point for storing properties?
        properties = rdmol.GetPropsAsDict()
        mol.properties = properties

        # We store bond orders as integers regardless of aromaticity.
        # In order to properly extract these, we need to have the "Kekulized" version of the rdkit mol
        kekul_mol = Chem.Mol(rdmol)
        Chem.Kekulize(kekul_mol, True)

        # setting chirality in openeye requires using neighbor atoms
        # therefore we can't do it until after the atoms and bonds are all added
        chiral_atoms = dict() # {rd_idx: openeye chirality}
        for rda in rdmol.GetAtoms():
            rd_idx = rda.GetIdx()

            # create a new atom
            atomic_number = oemol.NewAtom(rda.GetAtomicNum())
            formal_charge = rda.GetFormalCharge()
            is_aromatic = rda.GetIsAromatic()

            # If chiral, store the chirality to be set later
            stereochemistry = None
            tag = rda.GetChiralTag()
            if tag == Chem.CHI_TETRAHEDRAL_CCW:
                stereochemistry = 'R'
            if tag == Chem.CHI_TETRAHEDRAL_CW:
                stereochemistry = 'S'

            atom_index = mol.add_atom(atomic_number=atomic_number, formal_charge=formal_charge, is_aromatic=is_aromatic, stereochemistry=stereochemistry)
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
            bond_index = mol.add_bond(map_atoms[a1], map_atoms[a2], is_aromatic=is_aromatic, order=order, stereochemistry=stereochemistry)

        # TODO: Save conformer(s), if present
        # If the rdmol has a conformer, store its coordinates
        # TODO: Note, this currently only adds the first conformer, it will need to be adjusted if the
        # you wanted to convert multiple sets of coordinates
        if rdmol.GetConformers():
            conf = rdmol.GetConformer()
            # TODO: Store conformers
            #for rd_idx, oeatom in map_atoms.items():
            #    coords = conf.GetAtomPosition(rd_idx)
            #    oemol.SetCoords(oeatom, oechem.OEFloatArray(coords))

        return mol

    @staticmethod
    @requires_rdkit()
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
        rdmol.SetProp('_Name', self.name)

        # TODO: Set other properties
        for name, value in self.properties.items():
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
        for index, atom in enumerate(self.atoms):
            rdatom = Chem.Atom(atom.atomic_number)
            rdatom.SetFormalCharge(atom.formal_charge)
            rdatom.SetIsAromatic(atom.is_aromatic)

            if atom.stereochemistry == 'S':
                rdatom.SetChiralTag(Chem.CHI_TETRAHEDRAL_CW)
            elif atom.stereochemistry == 'R':
                rdatom.SetChiralTag(Chem.CHI_TETRAHEDRAL_CCW)

            map_atoms[oe_idx] = rdmol.AddAtom(rdatom)

        for bond in self.bonds:
            rdatom1 = map_atoms[bond.atom1_index]
            rdatom2 = map_atoms[bond.atom2_index]
            rdmol.AddBond(rdatom1, rdatom2)
            rdbond = rdmol.GetBondBetweenAtoms(rdatom1, rdatom2)

            # Assign bond type, which is based on order unless it is aromatic
            if bond.is_aromatic:
                rdbond.SetBondType(_bondtypes[1.5])
                rdbond.SetIsAromatic(True)
            else:
                rdbond.SetBondType(_bondtypes[bond.order])
                rdbond.SetIsAromatic(False)

        # Assign bond stereochemistry
        for bond in self.bonds:
            if bond.stereochemistry:
                # Determine neighbors
                # TODO: This API needs to be created
                n1 = [n.index for n in bond.atom1.bonded_atoms if n != bond.atom2][0]
                n2 = [n.index for n in bond.atom2.bonded_atoms if n != bond.atom1][0]
                # Get rdmol bonds
                bond1 = rdmol.GetBondBetweenAtoms(map_atoms[n1], map_atoms[bond.atom1_index])
                bond2 = rdmol.GetBondBetweenAtoms(map_atoms[bond.atom1_index], map_atoms[bond.atom2_index])
                bond3 = rdmol.GetBondBetweenAtoms(map_atoms[bond.atom2_index], map_atoms[n2])
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
        # TODO: Fix this once conformer API is defined
        if self._conformers:
            for conformer in self._conformers:
                rdmol_conformer = Chem.Conformer()
                for index, rd_idx in map_atoms.items():
                    (x,y,z) = conformer[index,:]
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
        return rdkit.Mol(rdmol)

class AmberToolsWrapper(ToolkitWrapper):
    """
    AmberTools toolkit wrapper

    """
    def __init__(self):
        # TODO: Find AMBERHOME or executable home
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
        # Check that the requested charge method is supported
        SUPPORTED_ANTECHAMBER_CHARGE_METHODS = ['gas', 'mul', 'cm1', 'cm2', 'bcc']
        if charge_method not in SUPPORTED_ANTECHAMBER_CHARGE_METHODS:
            raise ValueError('Requested charge method {} not among supported charge methods {}'.format(charge_method, SUPPORTED_ANTECHAMBER_CHARGE_METHODS))

        # Find the path to antechamber
        ANTECHAMBER_PATH = find_executable("antechamber")
        if ANTECHAMBER_PATH is None:
            raise(IOError("Antechamber not found, cannot run charge_mol()"))

        # Compute charges
        from openmmtools.utils import temporary_directory, temporary_cd
        with temporary_directory() as tmpdir:
            with temporary_cd(tmpdir):
                net_charge = molecule.net_charge()
                # Write out molecule in SDF format
                molecule.to_file('molecule.sdf', format='SDF')
                # Compute desired charges
                # TODO: Add error handling if antechamber chokes
                os.system("antechamber -i molecule.sdf -fi sdf -o charged.mol2 -fo mol2 -pf yes -c {} -nc {}".format(charge_method, net_charge))
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

class ToolkitRegistry(object):
    """
    Registry for ToolkitWrapper objects

    Examples
    --------

    Register toolkits in a specified order, skipping if unavailable

    >>> toolkit_registry = ToolkitRegistry()
    >>> toolkit_precedence = [OpenEyeToolkitWrapper, RDKitToolkitWrapper, AmberToolsWrapper]
    >>> for toolkit in toolkit_precedence: toolkit_registry.register(toolkit, exception_if_unavailable=False)

    Register specified toolkits, raising an exception if one is unavailable

    >>> toolkit_registry = ToolkitRegistry()
    >>> toolkits = [OpenEyeToolkitWrapper, AmberToolsWrapper]
    >>> for toolkit in toolkits: toolkit_registry.register(toolkit)

    Register all available toolkits in arbitrary order

    >>> from openforcefield.utils import all_subclasses
    >>> toolkits = all_subclasses(ToolkitWrapper)
    >>> for toolkit in toolkits: toolkit_registry.register_toolkit(toolkit, exception_if_unavailable)

    """
    def __init__(self, register_imported_toolkit_wrappers=False):
        """
        Create an empty toolkit registry.

        .. warning :: This API experimental and subject to change.

        Parameters
        ----------
        register_imported_toolkit_wrappers : bool, optional, default=False
            If True, will attempt to register all imported ToolkitWrapper subclasses that can be found.

        """
        self._toolkits = list()

        if register_imported_toolkit_wrappers:
            toolkits = all_subclasses(ToolkitWrapper)
            for toolkit in toolkit_precedence:
                toolkit_registry.register_toolkit(toolkit)

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
                return method(*args, **kwargs)
                try:
                    return method
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
