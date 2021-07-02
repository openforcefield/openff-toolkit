"""
Wrapper class providing a minimal consistent interface to
the `OpenEye Toolkit <https://docs.eyesopen.com/toolkits/python/quickstart-python/index.html>`_
"""

__all__ = ("OpenEyeToolkitWrapper",)


# =============================================================================================
# IMPORTS
# =============================================================================================

import importlib
import logging
import re
import pathlib
import tempfile
from collections import defaultdict
from functools import wraps
from typing import TYPE_CHECKING, List, Optional, Tuple

import numpy as np
from simtk import unit

if TYPE_CHECKING:
    from openforcefield.topology.molecule import Molecule

from . import base_wrapper
from .constants import DEFAULT_AROMATICITY_MODEL
from .exceptions import (
    ChargeCalculationError,
    ChargeMethodUnavailableError,
    GAFFAtomTypeWarning,
    InvalidIUPACNameError,
    LicenseError,
    ToolkitUnavailableException,
    UndefinedStereochemistryError,
    ParseError,
)
from .utils import inherit_docstrings

# =============================================================================================
# CONFIGURE LOGGER
# =============================================================================================

logger = logging.getLogger(__name__)


# =============================================================================================
# IMPLEMENTATION
# =============================================================================================

def get_oeformat(file_format):
    from openeye import oechem
    
    file_format = file_format.upper()
    # XXX This is what RDKit does. Should be supported here too?
    if file_format == "MOL":
        file_format = "SDF"
        
    oeformat = getattr(oechem, "OEFormat_" + file_format, None)
    if oeformat is None:
        raise ValueError(f"Unsupported file format: {file_format}")
    return oeformat

@inherit_docstrings
class OpenEyeToolkitWrapper(base_wrapper.ToolkitWrapper):
    """
    OpenEye toolkit wrapper

    .. warning :: This API is experimental and subject to change.
    """

    _toolkit_name = "OpenEye Toolkit"
    _toolkit_installation_instructions = (
        "The OpenEye toolkit requires a (free for academics) license, and can be "
        "found at: "
        "https://docs.eyesopen.com/toolkits/python/quickstart-python/install.html"
    )
    # This could belong to ToolkitWrapper, although it seems strange
    # to carry that data for open-source toolkits
    _is_licensed = None
    # Only for OpenEye is there potentially a difference between
    # being available and installed
    _is_installed = None
    _license_functions = {
        "oechem": "OEChemIsLicensed",
        "oequacpac": "OEQuacPacIsLicensed",
        "oeiupac": "OEIUPACIsLicensed",
        "oeomega": "OEOmegaIsLicensed",
    }

    def __init__(self):

        self._toolkit_file_read_formats = [
            "CAN",
            "CDX",
            "CSV",
            "FASTA",
            "INCHI",
            "INCHIKEY",
            "ISM",
            "MDL",
            "MF",
            "MMOD",
            "MOL2",
            "MOL2H",
            "MOPAC",
            "OEB",
            "PDB",
            "RDF",
            "SDF",
            "SKC",
            "SLN",
            "SMI",
            "USM",
            "XYC",
        ]
        self._toolkit_file_write_formats = [
            "CAN",
            "CDX",
            "CSV",
            "FASTA",
            "INCHI",
            "INCHIKEY",
            "ISM",
            "MDL",
            "MF",
            "MMOD",
            "MOL2",
            "MOL2H",
            "MOPAC",
            "OEB",
            "PDB",
            "RDF",
            "SDF",
            "SKC",
            "SLN",
            "SMI",
            "USM",
            "XYC",
        ]

        # check if the toolkit can be loaded
        if not self.is_available():
            msg = (
                f"The required toolkit {self._toolkit_name} is not "
                f"available. {self._toolkit_installation_instructions}"
            )
            if self._is_installed is False:
                raise ToolkitUnavailableException(msg)
            if self._is_licensed is False:
                raise LicenseError(msg)

        from openeye import __version__ as openeye_version

        self._toolkit_version = openeye_version

    @classmethod
    def _check_licenses(cls):
        """Check license of all known OpenEye tools. Returns True if any are found
        to be licensed, False if any are not."""
        for (tool, license_func) in cls._license_functions.items():
            try:
                module = importlib.import_module("openeye." + tool)
            except (ImportError, ModuleNotFoundError):
                continue
            else:
                if getattr(module, license_func)():
                    return True
        return False

    @classmethod
    def is_available(cls):
        """
        Check if the given OpenEye toolkit components are available.

        If the OpenEye toolkit is not installed or no license is found
        for at least one the required toolkits , ``False`` is returned.

        Returns
        -------
        all_installed : bool
            ``True`` if all required OpenEye tools are installed and licensed,
            ``False`` otherwise

        """
        if cls._is_available is None:
            if cls._is_licensed is None:
                cls._is_licensed = cls._check_licenses()
            if cls._is_installed is None:
                for tool in cls._license_functions.keys():
                    cls._is_installed = True
                    try:
                        importlib.import_module("openeye." + tool)
                    except (ImportError, ModuleNotFoundError):
                        cls._is_installed = False
            cls._is_available = cls._is_installed and cls._is_licensed
        return cls._is_available

    def from_object(self, obj, allow_undefined_stereo=False, _cls=None):
        """
        Convert an OEMol (or OEMol-derived object) into an openff.toolkit.topology.molecule

        Parameters
        ----------
        obj : A molecule-like object
            An object to by type-checked.
        allow_undefined_stereo : bool, default=False
            Whether to accept molecules with undefined stereocenters. If False,
            an exception will be raised if a molecule with undefined stereochemistry
            is passed into this function.
        _cls : class
            Molecule constructor
        Returns
        -------
        Molecule
            An openff.toolkit.topology.molecule Molecule.

        Raises
        ------
        NotImplementedError
            If the object could not be converted into a Molecule.
        """
        # TODO: Add tests for the from_object functions
        from openeye import oechem

        if _cls is None:
            from openff.toolkit.topology.molecule import Molecule

            _cls = Molecule

        if isinstance(obj, oechem.OEMolBase):
            return self.from_openeye(
                oemol=obj, allow_undefined_stereo=allow_undefined_stereo, _cls=_cls
            )
        raise NotImplementedError(
            "Cannot create Molecule from {} object".format(type(obj))
        )

    def from_file(
        self, file_path, file_format, allow_undefined_stereo=False, _cls=None
    ):
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
            If false, raises an exception if oemol contains undefined stereochemistry.
        _cls : class
            Molecule constructor

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

        >>> from openff.toolkit.utils import get_data_file_path
        >>> mol2_file_path = get_data_file_path('molecules/cyclohexane.mol2')
        >>> toolkit = OpenEyeToolkitWrapper()
        >>> molecule = toolkit.from_file(mol2_file_path, file_format='mol2')

        """
        from openeye import oechem

        oeformat = get_oeformat(file_format)
        ifs = oechem.oemolistream(file_path)
        if not ifs.IsValid():
            # Get Python to report ane error message, if possible.
            # This can distinguish between FileNotFound, IsADirectoryError, etc.
            open(file_path).close()
            # If that worked, then who knows. Fail anyway.
            raise OSError("Unable to open file")
        
        ifs.SetFormat(oeformat)
        
        return self._read_oemolistream_molecules(
            ifs, allow_undefined_stereo, file_path=file_path, _cls=_cls
        )

    def from_file_obj(
        self, file_obj, file_format, allow_undefined_stereo=False, _cls=None
    ):
        """
        Return an openff.toolkit.topology.Molecule from a file-like object (an object with a ".read()" method using
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
        _cls : class
            Molecule constructor

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
        oeformat = get_oeformat(file_format)
        ifs.SetFormat(oeformat)

        return self._read_oemolistream_molecules(ifs, allow_undefined_stereo, _cls=_cls)

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
        # This function requires a text-mode file_obj.
        try:
            file_obj.write("")
        except TypeError:
            # Switch to a ValueError and use a more informative exception
            # message to match RDKit.
            raise ValueError("Need a text mode file object like StringIO or a file opened with mode 't'") from None

        with tempfile.TemporaryDirectory() as tmpdir:
            path = pathlib.Path(tmpdir, f"input.{file_format.lower()}")
            self.to_file(molecule, str(path), file_format)
            file_data = path.read_text()
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
        openeye_format = get_oeformat(file_format)
        ofs.SetFormat(openeye_format)
        
        if openeye_format == oechem.OEFormat_SMI:
            ofs.SetFlavor(openeye_format, self._get_smiles_flavor(isomeric=True, explicit_hydrogens=True))

        # OFFTK strictly treats SDF as a single-conformer format.
        # We need to override OETK's behavior here if the user is saving a multiconformer molecule.

        # Remove all but the first conformer when writing to SDF as we only support single conformer format
        if (file_format.lower() == "sdf") and oemol.NumConfs() > 1:
            conf1 = [conf for conf in oemol.GetConfs()][0]
            flat_coords = list()
            for idx, coord in conf1.GetCoords().items():
                flat_coords.extend(coord)
            oemol.DeleteConfs()
            oecoords = oechem.OEFloatArray(flat_coords)
            oemol.NewConf(oecoords)
        # We're standardizing on putting partial charges into SDFs under the `atom.dprop.PartialCharge` property
        if (file_format.lower() == "sdf") and (molecule.partial_charges is not None):
            partial_charges_list = [
                oeatom.GetPartialCharge() for oeatom in oemol.GetAtoms()
            ]
            partial_charges_str = " ".join([f"{val:f}" for val in partial_charges_list])
            # TODO: "dprop" means "double precision" -- Is there any way to make Python more accurately
            #  describe/infer the proper data type?
            oechem.OESetSDData(oemol, "atom.dprop.PartialCharge", partial_charges_str)

        # If the file format is "pdb" using OEWriteMolecule() rearranges the atoms (hydrogens are pushed to the bottom)
        # Issue #475 (https://github.com/openforcefield/openff-toolkit/issues/475)
        # dfhahn's workaround: Using OEWritePDBFile does not alter the atom arrangement
        if file_format.lower() == "pdb":
            if oemol.NumConfs() > 1:
                for conf in oemol.GetConfs():
                    oechem.OEWritePDBFile(ofs, conf, oechem.OEOFlavor_PDB_BONDS)
            else:
                oechem.OEWritePDBFile(ofs, oemol, oechem.OEOFlavor_PDB_BONDS)
        else:
            oechem.OEWriteMolecule(ofs, oemol)
        ofs.close()

    @staticmethod
    def _turn_oemolbase_sd_charges_into_partial_charges(oemol):
        """
        Process an OEMolBase object and check to see whether it has an SD data pair
        where the tag is "atom.dprop.PartialCharge", indicating that it has a list of
        atomic partial charges. If so, apply those charges to the OEAtoms in the OEMolBase,
        and delete the SD data pair.

        Parameters
        ----------
        oemol : openeye.oechem.OEMolBase
            The molecule to process

        Returns
        -------
        charges_are_present : bool
            Whether charges are present in the SD file. This is necessary because OEAtoms
            have a default partial charge of 0.0, which makes truly zero-charge molecules
            (eg "N2", "Ar"...) indistinguishable from molecules for which partial charges
            have not been assigned. The OFF Toolkit allows this distinction with
            mol.partial_charges=None. In order to complete roundtrips within the OFFMol
            spec, we must interpret the presence or absence of this tag as a proxy for
            mol.partial_charges=None.
        """
        from openeye import oechem

        for dp in oechem.OEGetSDDataPairs(oemol):
            if dp.GetTag() == "atom.dprop.PartialCharge":
                charges_str = oechem.OEGetSDData(oemol, "atom.dprop.PartialCharge")
                charges_unitless = [float(i) for i in charges_str.split()]
                assert len(charges_unitless) == oemol.NumAtoms()
                for charge, oeatom in zip(charges_unitless, oemol.GetAtoms()):
                    oeatom.SetPartialCharge(charge)
                oechem.OEDeleteSDData(oemol, "atom.dprop.PartialCharge")
                return True
        return False

    def _read_oemolistream_molecules(
        self, oemolistream, allow_undefined_stereo, file_path=None, _cls=None
    ):
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
        _cls : class
            Molecule constructor

        Returns
        -------
        molecules : List[Molecule]
            The list of Molecule objects in the stream.

        """
        from openeye import oechem

        mols = list()
        oemol = oechem.OEMol()
        while oechem.OEReadMolecule(oemolistream, oemol):
            oechem.OEPerceiveChiral(oemol)
            oechem.OEAssignAromaticFlags(oemol, oechem.OEAroModel_MDL)
            oechem.OE3DToInternalStereo(oemol)

            # If this is either a multi-conformer or multi-molecule SD file, check to see if there are partial charges
            if (oemolistream.GetFormat() == oechem.OEFormat_SDF) and hasattr(
                oemol, "GetConfs"
            ):
                # The openFF toolkit treats each conformer in a "multiconformer" SDF as
                # a separate molecule.
                # https://github.com/openforcefield/openff-toolkit/issues/202
                # Note that there is ambiguity about how SD data and "multiconformer" SD files should be stored.
                # As a result, we have to do some weird stuff below, as discussed in
                # https://docs.eyesopen.com/toolkits/python/oechemtk/oemol.html#dude-where-s-my-sd-data

                # Jeff: I was unable to find a way to distinguish whether a SDF was multiconformer or not.
                # The logic below should handle either single- or multi-conformer SDFs.
                for conf in oemol.GetConfIter():
                    # First, we turn "conf" into an OEMCMol (OE multiconformer mol), since OTHER file formats
                    # really are multiconformer, and we will eventually feed this into the `from_openeye` function,
                    # which is made to ingest multiconformer mols.
                    this_conf_oemcmol = conf.GetMCMol()

                    # Then, we take any SD data pairs that were on the oemol, and copy them on to "this_conf_oemcmol".
                    # These SD pairs will be populated if we're dealing with a single-conformer SDF.
                    for dp in oechem.OEGetSDDataPairs(oemol):
                        oechem.OESetSDData(
                            this_conf_oemcmol, dp.GetTag(), dp.GetValue()
                        )
                    # On the other hand, these SD pairs will be populated if we're dealing with a MULTI-conformer SDF.
                    for dp in oechem.OEGetSDDataPairs(conf):
                        oechem.OESetSDData(
                            this_conf_oemcmol, dp.GetTag(), dp.GetValue()
                        )
                    # This function fishes out the special SD data tag we use for partial charge
                    # ("atom.dprop.PartialCharge"), and applies those as OETK-supported partial charges on the OEAtoms
                    has_charges = self._turn_oemolbase_sd_charges_into_partial_charges(
                        this_conf_oemcmol
                    )

                    # Finally, we feed the molecule into `from_openeye`, where it converted into an OFFMol
                    mol = self.from_openeye(
                        this_conf_oemcmol,
                        allow_undefined_stereo=allow_undefined_stereo,
                        _cls=_cls,
                    )

                    # If the molecule didn't even have the `PartialCharges` tag, we set it from zeroes to None here.
                    if not (has_charges):
                        mol.partial_charges = None
                    mols.append(mol)

            else:
                # In case this is being read from a SINGLE-molecule SD file, convert the SD field where we
                # stash partial charges into actual per-atom partial charges
                self._turn_oemolbase_sd_charges_into_partial_charges(oemol)
                mol = self.from_openeye(
                    oemol, allow_undefined_stereo=allow_undefined_stereo, _cls=_cls
                )
                mols.append(mol)

            # Check if this is an AMBER-produced mol2 file, which we can not load because they use GAFF atom types.
            if oemolistream.GetFormat() == oechem.OEFormat_MOL2:
                self._check_mol2_gaff_atom_type(mol, file_path)

        return mols

    def enumerate_protomers(self, molecule, max_states=10):
        """
        Enumerate the formal charges of a molecule to generate different protomoers.

        Parameters
        ----------
        molecule: openff.toolkit.topology.Molecule
            The molecule whose state we should enumerate

        max_states: int optional, default=10,
            The maximum number of protomer states to be returned.

        Returns
        -------
        molecules: List[openff.toolkit.topology.Molecule],
            A list of the protomers of the input molecules not including the input.
        """

        from openeye import oequacpac

        options = oequacpac.OEFormalChargeOptions()
        # add one as the input is included
        options.SetMaxCount(max_states + 1)

        molecules = []

        oemol = self.to_openeye(molecule=molecule)
        for protomer in oequacpac.OEEnumerateFormalCharges(oemol, options):

            mol = self.from_openeye(
                protomer, allow_undefined_stereo=True, _cls=molecule.__class__
            )

            if mol != molecule:
                molecules.append(mol)

        return molecules

    def enumerate_stereoisomers(
        self, molecule, undefined_only=False, max_isomers=20, rationalise=True
    ):
        """
        Enumerate the stereocenters and bonds of the current molecule.

        Parameters
        ----------
        molecule: openff.toolkit.topology.Molecule
            The molecule whose state we should enumerate

        undefined_only: bool optional, default=False
            If we should enumerate all stereocenters and bonds or only those with undefined stereochemistry

        max_isomers: int optional, default=20
            The maximum amount of molecules that should be returned

        rationalise: bool optional, default=True
            If we should try to build and rationalise the molecule to ensure it can exist


        Returns
        --------
        molecules: List[openff.toolkit.topology.Molecule]
            A list of openff.toolkit.topology.Molecule instances

        """
        from openeye import oechem, oeomega

        oemol = self.to_openeye(molecule=molecule)

        # arguments for this function can be found here
        # <https://docs.eyesopen.com/toolkits/python/omegatk/OEConfGenFunctions/OEFlipper.html?highlight=stereoisomers>

        molecules = []
        for isomer in oeomega.OEFlipper(oemol, 200, not undefined_only, True, False):

            if rationalise:
                # try and determine if the molecule is reasonable by generating a conformer with
                # strict stereo, like embedding in rdkit
                omega = oeomega.OEOmega()
                omega.SetMaxConfs(1)
                omega.SetCanonOrder(False)
                # Don't generate random stereoisomer if not specified
                omega.SetStrictStereo(True)
                mol = oechem.OEMol(isomer)
                status = omega(mol)
                if status:
                    isomol = self.from_openeye(mol, _cls=molecule.__class__)
                    if isomol != molecule:
                        molecules.append(isomol)

            else:
                isomol = self.from_openeye(isomer, _cls=molecule.__class__)
                if isomol != molecule:
                    molecules.append(isomol)

        return molecules[:max_isomers]

    def enumerate_tautomers(self, molecule, max_states=20):
        """
        Enumerate the possible tautomers of the current molecule

        Parameters
        ----------
        molecule: openff.toolkit.topology.Molecule
            The molecule whose state we should enumerate

        max_states: int optional, default=20
            The maximum amount of molecules that should be returned

        Returns
        -------
        molecules: List[openff.toolkit.topology.Molecule]
            A list of openff.toolkit.topology.Molecule instances excluding the input molecule.
        """
        from openeye import oequacpac

        oemol = self.to_openeye(molecule=molecule)

        tautomers = []

        # set the options
        tautomer_options = oequacpac.OETautomerOptions()
        tautomer_options.SetApplyWarts(False)
        tautomer_options.SetMaxTautomersGenerated(max_states + 1)
        tautomer_options.SetSaveStereo(True)
        # this aligns the outputs of rdkit and openeye for the example cases
        tautomer_options.SetCarbonHybridization(False)

        for tautomer in oequacpac.OEEnumerateTautomers(oemol, tautomer_options):
            # remove the input tautomer from the output
            taut = self.from_openeye(
                tautomer, allow_undefined_stereo=True, _cls=molecule.__class__
            )
            if taut != molecule:
                tautomers.append(
                    self.from_openeye(
                        tautomer, allow_undefined_stereo=True, _cls=molecule.__class__
                    )
                )

        return tautomers

    @staticmethod
    def _check_mol2_gaff_atom_type(molecule, file_path=None):
        """Attempts to detect the presence of GAFF atom types in a molecule loaded from a mol2 file.

        For now, this raises a ``GAFFAtomTypeWarning`` if the molecule
        include Osmium and Holmium atoms, which have GAFF types OS and
        HO respectively.

        Parameters
        ----------
        molecule : openff.toolkit.topology.molecule.Molecule
            The loaded molecule.
        file_path : str, optional
            The path to the mol2 file. This is used exclusively to make
            the error message more meaningful.

        """
        # Handle default.
        if file_path is None:
            file_path = ""
        else:
            # Append a ':' character that will separate the file
            # path from the molecule string representation.
            file_path = file_path + ":"
        # atomic_number: (GAFF_type, element_name)
        warning_atomic_numbers = {76: ("OS", "Osmium"), 67: ("HO", "Holmium")}

        for atom in molecule.atoms:
            try:
                atom_type, element_name = warning_atomic_numbers[atom.atomic_number]
            except KeyError:
                pass
            else:
                import warnings

                warn_msg = (
                    f'OpenEye interpreted the type "{atom_type}" in {file_path}{molecule.name}'
                    f" as {element_name}. Does your mol2 file uses Tripos SYBYL atom types?"
                    " Other atom types such as GAFF are not supported."
                )
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
            return "S"
        elif cip == oechem.OECIPAtomStereo_R:
            return "R"
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
            return "E"
        elif cip == oechem.OECIPBondStereo_Z:
            return "Z"
        elif cip == oechem.OECIPBondStereo_NotStereo:
            return None

    @staticmethod
    def from_openeye(oemol, allow_undefined_stereo=False, _cls=None):
        """
        Create a Molecule from an OpenEye molecule. If the OpenEye molecule has
        implicit hydrogens, this function will make them explicit.

        ``OEAtom`` s have a different set of allowed value for partial charges than
        ``openff.toolkit.topology.Molecule`` s. In the OpenEye toolkits, partial charges
        are stored on individual ``OEAtom`` s, and their values are initialized to ``0.0``.
        In the Open Force Field Toolkit, an ``openff.toolkit.topology.Molecule``'s
        ``partial_charges`` attribute is initialized to ``None`` and can be set to a
        ``simtk.unit.Quantity``-wrapped numpy array with units of
        elementary charge. The Open Force
        Field Toolkit considers an ``OEMol`` where every ``OEAtom`` has a partial
        charge of ``float('nan')`` to be equivalent to an Open Force Field Toolkit `Molecule`'s
        ``partial_charges = None``.
        This assumption is made in both ``to_openeye`` and ``from_openeye``.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        oemol : openeye.oechem.OEMol
            An OpenEye molecule
        allow_undefined_stereo : bool, default=False
            If false, raises an exception if oemol contains undefined stereochemistry.
        _cls : class
            Molecule constructor

        Returns
        -------
        molecule : openff.toolkit.topology.Molecule
            An OpenFF molecule

        Examples
        --------

        Create a Molecule from an OpenEye OEMol

        >>> from openeye import oechem
        >>> from openff.toolkit.tests.utils import get_data_file_path
        >>> ifs = oechem.oemolistream(get_data_file_path('systems/monomers/ethanol.mol2'))
        >>> oemols = list(ifs.GetOEGraphMols())

        >>> toolkit_wrapper = OpenEyeToolkitWrapper()
        >>> molecule = toolkit_wrapper.from_openeye(oemols[0])

        """
        import math

        from openeye import oechem

        oemol = oechem.OEMol(oemol)

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
        if unspec_chiral or unspec_db:

            def oeatom_to_str(oeatom):
                return "atomic num: {}, name: {}, idx: {}, aromatic: {}, chiral: {}".format(
                    oeatom.GetAtomicNum(),
                    oeatom.GetName(),
                    oeatom.GetIdx(),
                    oeatom.IsAromatic(),
                    oeatom.IsChiral(),
                )

            def oebond_to_str(oebond):
                return "order: {}, chiral: {}".format(
                    oebond.GetOrder(), oebond.IsChiral()
                )

            def describe_oeatom(oeatom):
                description = "Atom {} with bonds:".format(oeatom_to_str(oeatom))
                for oebond in oeatom.GetBonds():
                    description += "\nbond {} to atom {}".format(
                        oebond_to_str(oebond), oeatom_to_str(oebond.GetNbr(oeatom))
                    )
                return description

            msg = (
                "OEMol has unspecified stereochemistry. "
                "oemol.GetTitle(): {}\n".format(oemol.GetTitle())
            )
            if len(problematic_atoms) != 0:
                msg += "Problematic atoms are:\n"
                for problematic_atom in problematic_atoms:
                    msg += describe_oeatom(problematic_atom) + "\n"
            if len(problematic_bonds) != 0:
                msg += "Problematic bonds are: {}\n".format(problematic_bonds)
            if allow_undefined_stereo:
                msg = "Warning (not error because allow_undefined_stereo=True): " + msg
                logger.warning(msg)
            else:
                msg = "Unable to make OFFMol from OEMol: " + msg
                raise UndefinedStereochemistryError(msg)

        if _cls is None:
            from openff.toolkit.topology.molecule import Molecule

            _cls = Molecule

        molecule = _cls()
        molecule.name = oemol.GetTitle()

        # Copy any attached SD tag information
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
                oemol, oeatom
            )
            # stereochemistry = self._openeye_cip_atom_stereochemistry(oemol, oeatom)
            name = ""
            if oeatom.HasData("name"):
                name = oeatom.GetData("name")
            atom_index = molecule._add_atom(
                atomic_number,
                formal_charge,
                is_aromatic,
                stereochemistry=stereochemistry,
                name=name,
            )
            map_atoms[
                oe_idx
            ] = atom_index  # store for mapping oeatom to molecule atom indices below
            atom_mapping[atom_index] = map_id

        # If we have a full / partial atom map add it to the molecule. Zeroes 0
        # indicates no mapping
        if {*atom_mapping.values()} != {0}:

            molecule._properties["atom_map"] = {
                idx: map_idx for idx, map_idx in atom_mapping.items() if map_idx != 0
            }

        for oebond in oemol.GetBonds():
            atom1_index = map_atoms[oebond.GetBgnIdx()]
            atom2_index = map_atoms[oebond.GetEndIdx()]
            bond_order = oebond.GetOrder()
            is_aromatic = oebond.IsAromatic()
            stereochemistry = OpenEyeToolkitWrapper._openeye_cip_bond_stereochemistry(
                oemol, oebond
            )
            if oebond.HasData("fractional_bond_order"):
                fractional_bond_order = oebond.GetData("fractional_bond_order")
            else:
                fractional_bond_order = None

            molecule._add_bond(
                atom1_index,
                atom2_index,
                bond_order,
                is_aromatic=is_aromatic,
                stereochemistry=stereochemistry,
                fractional_bond_order=fractional_bond_order,
            )

        # TODO: Copy conformations, if present
        # TODO: Come up with some scheme to know when to import coordinates
        # From SMILES: no
        # From MOL2: maybe
        # From other: maybe
        if hasattr(oemol, "GetConfs"):
            for conf in oemol.GetConfs():
                n_atoms = molecule.n_atoms
                positions = unit.Quantity(
                    np.zeros(shape=[n_atoms, 3], dtype=np.float64), unit.angstrom
                )
                for oe_id in conf.GetCoords().keys():
                    off_atom_coords = unit.Quantity(
                        conf.GetCoords()[oe_id], unit.angstrom
                    )
                    off_atom_index = map_atoms[oe_id]
                    positions[off_atom_index, :] = off_atom_coords
                if (positions == 0 * unit.angstrom).all() and n_atoms > 1:
                    continue
                molecule._add_conformer(positions)

        # Copy partial charges, if present
        partial_charges = unit.Quantity(
            np.zeros(shape=molecule.n_atoms, dtype=np.float64),
            unit=unit.elementary_charge,
        )

        # If all OEAtoms have a partial charge of NaN, then the OFFMol should
        # have its partial_charges attribute set to None
        any_partial_charge_is_not_nan = False
        for oe_atom in oemol.GetAtoms():
            oe_idx = oe_atom.GetIdx()
            off_idx = map_atoms[oe_idx]
            unitless_charge = oe_atom.GetPartialCharge()
            if not math.isnan(unitless_charge):
                any_partial_charge_is_not_nan = True
                # break
            charge = unitless_charge * unit.elementary_charge
            partial_charges[off_idx] = charge

        if any_partial_charge_is_not_nan:
            molecule.partial_charges = partial_charges
        else:
            molecule.partial_charges = None

        return molecule

    @staticmethod
    def to_openeye(molecule, aromaticity_model=DEFAULT_AROMATICITY_MODEL):
        r"""
        Create an OpenEye molecule using the specified aromaticity model

        ``OEAtom`` s have a different set of allowed value for partial
        charges than ``openff.toolkit.topology.Molecule``\ s. In the
        OpenEye toolkits, partial charges are stored on individual
        ``OEAtom``\ s, and their values are initialized to ``0.0``. In
        the Open Force Field Toolkit, an``openff.toolkit.topology.Molecule``'s
        ``partial_charges`` attribute is initialized to ``None`` and can
        be set to a ``simtk.unit.Quantity``-wrapped numpy array with
        units of elementary charge. The Open Force Field Toolkit
        considers an ``OEMol`` where every ``OEAtom`` has a partial
        charge of ``float('nan')`` to be equivalent to an Open Force
        Field Toolkit ``Molecule``'s ``partial_charges = None``. This
        assumption is made in both ``to_openeye`` and ``from_openeye``.

        .. todo ::

           * Should the aromaticity model be specified in some other way?

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        molecule : openff.toolkit.topology.molecule.Molecule object
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

        >>> from openff.toolkit.topology import Molecule
        >>> toolkit_wrapper = OpenEyeToolkitWrapper()
        >>> molecule = Molecule.from_smiles('CC')
        >>> oemol = toolkit_wrapper.to_openeye(molecule)

        """
        from openeye import oechem

        if hasattr(oechem, aromaticity_model):
            oe_aro_model = getattr(oechem, aromaticity_model)
        else:
            raise ValueError(
                "Error: provided aromaticity model not recognized by oechem."
            )

        oemol = oechem.OEMol()
        # if not(molecule.name is None):
        oemol.SetTitle(molecule.name)
        map_atoms = {}  # {off_idx : oe_idx}
        # Add atoms
        oemol_atoms = list()  # list of corresponding oemol atoms
        for atom in molecule.atoms:
            oeatom = oemol.NewAtom(atom.atomic_number)
            oeatom.SetFormalCharge(
                atom.formal_charge.value_in_unit(unit.elementary_charge)
            )  # simtk.unit.Quantity(1, unit.elementary_charge)
            # TODO: Do we want to provide _any_ pathway for Atom.is_aromatic to influence the OEMol?
            # oeatom.SetAromatic(atom.is_aromatic)
            oeatom.SetData("name", atom.name)
            oeatom.SetPartialCharge(float("nan"))
            oemol_atoms.append(oeatom)
            map_atoms[atom.molecule_atom_index] = oeatom.GetIdx()

        # Add bonds
        oemol_bonds = list()  # list of corresponding oemol bonds
        for bond in molecule.bonds:
            # atom1_index = molecule.atoms.index(bond.atom1)
            # atom2_index = molecule.atoms.index(bond.atom2)
            atom1_index = bond.atom1_index
            atom2_index = bond.atom2_index
            oebond = oemol.NewBond(oemol_atoms[atom1_index], oemol_atoms[atom2_index])
            oebond.SetOrder(bond.bond_order)
            # TODO: Do we want to provide _any_ pathway for Bond.is_aromatic to influence the OEMol?
            # oebond.SetAromatic(bond.is_aromatic)
            if not (bond.fractional_bond_order is None):
                oebond.SetData("fractional_bond_order", bond.fractional_bond_order)
            oemol_bonds.append(oebond)

        oechem.OEAssignAromaticFlags(oemol, oe_aro_model)

        # Set atom stereochemistry now that all connectivity is in place
        for atom, oeatom in zip(molecule.atoms, oemol_atoms):
            if not atom.stereochemistry:
                continue

            # Set arbitrary initial stereochemistry
            neighs = [n for n in oeatom.GetAtoms()]
            oeatom.SetStereo(
                neighs, oechem.OEAtomStereo_Tetra, oechem.OEAtomStereo_Right
            )

            # Flip chirality if stereochemistry isincorrect
            oeatom_stereochemistry = (
                OpenEyeToolkitWrapper._openeye_cip_atom_stereochemistry(oemol, oeatom)
            )
            if oeatom_stereochemistry != atom.stereochemistry:
                # Flip the stereochemistry
                oeatom.SetStereo(
                    neighs, oechem.OEAtomStereo_Tetra, oechem.OEAtomStereo_Left
                )
                # Verify it matches now as a sanity check
                oeatom_stereochemistry = (
                    OpenEyeToolkitWrapper._openeye_cip_atom_stereochemistry(
                        oemol, oeatom
                    )
                )
                if oeatom_stereochemistry != atom.stereochemistry:
                    raise Exception(
                        "Programming error: OpenEye atom stereochemistry assumptions failed."
                    )

        # Set bond stereochemistry
        for bond, oebond in zip(molecule.bonds, oemol_bonds):
            if not bond.stereochemistry:
                continue

            atom1_index = bond.molecule.atoms.index(bond.atom1)
            atom2_index = bond.molecule.atoms.index(bond.atom2)
            # Set arbitrary initial stereochemistry
            oeatom1, oeatom2 = oemol_atoms[atom1_index], oemol_atoms[atom2_index]
            oeatom1_neighbor = [n for n in oeatom1.GetAtoms() if not n == oeatom2][0]
            oeatom2_neighbor = [n for n in oeatom2.GetAtoms() if not n == oeatom1][0]
            # oebond.SetStereo([oeatom1, oeatom2], oechem.OEBondStereo_CisTrans, oechem.OEBondStereo_Cis)
            oebond.SetStereo(
                [oeatom1_neighbor, oeatom2_neighbor],
                oechem.OEBondStereo_CisTrans,
                oechem.OEBondStereo_Cis,
            )

            # Flip stereochemistry if incorrect
            oebond_stereochemistry = (
                OpenEyeToolkitWrapper._openeye_cip_bond_stereochemistry(oemol, oebond)
            )
            if oebond_stereochemistry != bond.stereochemistry:
                # Flip the stereochemistry
                oebond.SetStereo(
                    [oeatom1_neighbor, oeatom2_neighbor],
                    oechem.OEBondStereo_CisTrans,
                    oechem.OEBondStereo_Trans,
                )
                # Verify it matches now as a sanity check
                oebond_stereochemistry = (
                    OpenEyeToolkitWrapper._openeye_cip_bond_stereochemistry(
                        oemol, oebond
                    )
                )
                if oebond_stereochemistry != bond.stereochemistry:
                    raise Exception(
                        "Programming error: OpenEye bond stereochemistry assumptions failed."
                    )

        # Retain conformations, if present
        if molecule.n_conformers != 0:
            oemol.DeleteConfs()
            for conf in molecule._conformers:
                # OE needs a 1 x (3*n_Atoms) double array as input
                flat_coords = np.zeros(shape=oemol.NumAtoms() * 3, dtype=np.float64)
                for index, oe_idx in map_atoms.items():
                    (x, y, z) = conf[index, :] / unit.angstrom
                    flat_coords[(3 * oe_idx)] = x
                    flat_coords[(3 * oe_idx) + 1] = y
                    flat_coords[(3 * oe_idx) + 2] = z

                oecoords = oechem.OEFloatArray(flat_coords)
                oemol.NewConf(oecoords)

        # Retain charges, if present. All atoms are initialized above with a partial charge of NaN.
        if molecule._partial_charges is not None:
            oe_indexed_charges = np.zeros(shape=molecule.n_atoms, dtype=np.float64)
            for off_idx, charge in enumerate(molecule._partial_charges):
                oe_idx = map_atoms[off_idx]
                charge_unitless = charge / unit.elementary_charge
                oe_indexed_charges[oe_idx] = charge_unitless
            # TODO: This loop below fails if we try to use an "enumerate"-style loop.
            #  It's worth investigating whether we make this assumption elsewhere in the codebase, since
            #  the OE docs may indicate that this sort of usage is a very bad thing to do.
            #  https://docs.eyesopen.com/toolkits/python/oechemtk/atombondindices.html#indices-for-molecule-lookup-considered-harmful
            # for oe_idx, oe_atom in enumerate(oemol.GetAtoms()):
            for oe_atom in oemol.GetAtoms():
                oe_idx = oe_atom.GetIdx()
                oe_atom.SetPartialCharge(oe_indexed_charges[oe_idx])

        # Retain properties, if present
        for key, value in molecule.properties.items():
            oechem.OESetSDData(oemol, str(key), str(value))

        # Clean Up phase
        # The only feature of a molecule that wasn't perceived above seemed to be ring connectivity, better to run it
        # here then for someone to inquire about ring sizes and get 0 when it shouldn't be
        oechem.OEFindRingAtomsAndBonds(oemol)

        return oemol

    def _get_smiles_flavor(self, isomeric, explicit_hydrogens):
        from openeye import oechem
        
        # this sets up the default settings following the old DEFAULT flag
        # more information on flags can be found here
        # <https://docs.eyesopen.com/toolkits/python/oechemtk/OEChemConstants/OESMILESFlag.html#OEChem::OESMILESFlag>
        smiles_options = (
            oechem.OESMILESFlag_Canonical
            | oechem.OESMILESFlag_Isotopes
            | oechem.OESMILESFlag_RGroups
        )

        # check if we want an isomeric smiles
        if isomeric:
            # add the atom and bond stereo flags
            smiles_options |= (
                oechem.OESMILESFlag_AtomStereo | oechem.OESMILESFlag_BondStereo
            )

        if explicit_hydrogens:
            # add the hydrogen flag
            smiles_options |= oechem.OESMILESFlag_Hydrogens

        return smiles_options
        
    def to_smiles(self, molecule, isomeric=True, explicit_hydrogens=True, mapped=False):
        """
        Uses the OpenEye toolkit to convert a Molecule into a SMILES string.
        A partially mapped smiles can also be generated for atoms of interest by supplying an `atom_map` to the
        properties dictionary.

        Parameters
        ----------
        molecule : An openff.toolkit.topology.Molecule
            The molecule to convert into a SMILES.
        isomeric: bool optional, default= True
            return an isomeric smiles
        explicit_hydrogens: bool optional, default=True
            return a smiles string containing all hydrogens explicitly
        mapped: bool optional, default=False
            return a explicit hydrogen mapped smiles, the atoms to be mapped can be controlled by supplying an
            atom map into the properties dictionary. If no mapping is passed all atoms will be mapped in order, else
            an atom map dictionary from the current atom index to the map id should be supplied with no duplicates.
            The map ids (values) should start from 0 or 1.

        Returns
        -------
        smiles : str
            The SMILES of the input molecule.
        """
        from openeye import oechem

        oemol = self.to_openeye(molecule)

        smiles_options = self._get_smiles_flavor(isomeric, explicit_hydrogens)
            
        if mapped:
            assert explicit_hydrogens is True, (
                "Mapped smiles require all hydrogens and "
                "stereochemsitry to be defined to retain order"
            )

            # if we only want to map specific atoms check for an atom map
            atom_map = molecule._properties.get("atom_map", None)
            if atom_map is not None:
                # make sure there are no repeated indices
                map_ids = set(atom_map.values())
                if len(map_ids) < len(atom_map):
                    atom_map = None
                elif 0 in atom_map.values():
                    # we need to increment the map index
                    for atom, map in atom_map.items():
                        atom_map[atom] = map + 1

            if atom_map is None:
                # now we need to add the atom map to the atoms
                for oeatom in oemol.GetAtoms():
                    oeatom.SetMapIdx(oeatom.GetIdx() + 1)
            else:
                for atom in oemol.GetAtoms():
                    try:
                        # try to set the atom map
                        map_idx = atom_map[atom.GetIdx()]
                        atom.SetMapIdx(map_idx)
                    except KeyError:
                        continue

            smiles_options |= oechem.OESMILESFlag_AtomMaps

        smiles = oechem.OECreateSmiString(oemol, smiles_options)
        return smiles

    def to_inchi(self, molecule, fixed_hydrogens=False):
        """
        Create an InChI string for the molecule using the RDKit Toolkit.
        InChI is a standardised representation that does not capture tautomers
        unless specified using the fixed hydrogen layer.

        For information on InChi see here https://iupac.org/who-we-are/divisions/division-details/inchi/

        Parameters
        ----------
        molecule : An openff.toolkit.topology.Molecule
            The molecule to convert into a SMILES.

        fixed_hydrogens: bool, default=False
            If a fixed hydrogen layer should be added to the InChI, if `True` this
            will produce a non standard specific InChI string of the molecule.

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
        InChIKey is a standardised representation that does not capture tautomers
        unless specified using the fixed hydrogen layer.

        For information on InChi see here https://iupac.org/who-we-are/divisions/division-details/inchi/

        Parameters
        ----------
        molecule : An openff.toolkit.topology.Molecule
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

    def to_iupac(self, molecule):
        """Generate IUPAC name from Molecule

        Parameters
        ----------
        molecule : An openff.toolkit.topology.Molecule
            The molecule to convert into a SMILES.

        Returns
        -------
        iupac_name : str
            IUPAC name of the molecule

        Examples
        --------

        >>> from openff.toolkit.topology import Molecule
        >>> from openff.toolkit.utils import get_data_file_path
        >>> sdf_filepath = get_data_file_path('molecules/ethanol.sdf')
        >>> molecule = Molecule(sdf_filepath)
        >>> toolkit = OpenEyeToolkitWrapper()
        >>> iupac_name = toolkit.to_iupac(molecule)

        """
        from openeye import oeiupac

        oemol = self.to_openeye(molecule)

        return oeiupac.OECreateIUPACName(oemol)

    def canonical_order_atoms(self, molecule):
        """
        Canonical order the atoms in the molecule using the OpenEye toolkit.

        Parameters
        ----------
        molecule: openff.toolkit.topology.Molecule
            The input molecule

         Returns
        -------
        molecule : openff.toolkit.topology.Molecule
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
            if (
                bond.GetBgn().GetAtomicNum() != oechem.OEElemNo_H
                and bond.GetEnd().GetAtomicNum() != oechem.OEElemNo_H
            ):
                vbnd.append(bond)
        oemol.OrderBonds(vbnd)

        oemol.Sweep()

        for bond in oemol.GetBonds():
            if bond.GetBgnIdx() > bond.GetEndIdx():
                bond.SwapEnds()

        return self.from_openeye(
            oemol, allow_undefined_stereo=True, _cls=molecule.__class__
        )

    def from_smiles(
        self,
        smiles,
        hydrogens_are_explicit=False,
        allow_undefined_stereo=False,
        _cls=None,
    ):
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
        _cls : class
            Molecule constructor
        Returns
        -------
        molecule : openff.toolkit.topology.Molecule
            An OpenFF style molecule.
        """
        from openeye import oechem

        oemol = oechem.OEGraphMol()
        if not oechem.OESmilesToMol(oemol, smiles):
            raise ParseError("Unable to parse the SMILES string")
        if not (hydrogens_are_explicit):
            result = oechem.OEAddExplicitHydrogens(oemol)
            if not result:
                raise ValueError(
                    "Addition of explicit hydrogens failed in from_openeye"
                )
        elif hydrogens_are_explicit and oechem.OEHasImplicitHydrogens(oemol):
            raise ValueError(
                f"'hydrogens_are_explicit' was specified as True, but OpenEye Toolkit interpreted "
                f"SMILES '{smiles}' as having implicit hydrogen. If this SMILES is intended to "
                f"express all explicit hydrogens in the molecule, then you should construct the "
                f"desired molecule as an OEMol (where oechem.OEHasImplicitHydrogens(oemol) returns "
                f"False), and then use Molecule.from_openeye() to create the desired OFFMol."
            )

        # Set partial charges to None, since they couldn't have been stored in a SMILES
        for atom in oemol.GetAtoms():
            atom.SetPartialCharge(float("nan"))

        molecule = self.from_openeye(
            oemol, _cls=_cls, allow_undefined_stereo=allow_undefined_stereo
        )
        return molecule

    def from_inchi(self, inchi, allow_undefined_stereo=False, _cls=None):
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

        _cls : class
            Molecule constructor

        Returns
        -------
        molecule : openff.toolkit.topology.Molecule
        """

        from openeye import oechem

        # This calls the same functions as OESmilesToMol
        oemol = oechem.OEGraphMol()
        oechem.OEInChIToMol(oemol, inchi)

        # try and catch InChI parsing fails
        # if there are no atoms don't build the molecule
        if oemol.NumAtoms() == 0:
            raise RuntimeError(
                "There was an issue parsing the InChI string, please check and try again."
            )

        molecule = self.from_openeye(
            oemol, allow_undefined_stereo=allow_undefined_stereo, _cls=_cls
        )

        return molecule

    def from_iupac(self, iupac_name, allow_undefined_stereo=False, _cls=None, **kwargs):
        """
        Construct a Molecule from an IUPAC name

        Parameters
        ----------
        iupac_name : str
            The IUPAC or common name of the molecule.
        allow_undefined_stereo : bool, default=False
            Whether to accept a molecule name with undefined stereochemistry. If False,
            an exception will be raised if a molecule name with undefined stereochemistry
            is passed into this function.
        _cls : class
            Molecule constructor

        Returns
        -------
        molecule : openff.toolkit.topology.Molecule

        """
        from openeye import oechem, oeiupac

        oemol = oechem.OEMol()
        parsing_result = oeiupac.OEParseIUPACName(oemol, iupac_name)
        if not parsing_result:
            raise InvalidIUPACNameError(
                f"OpenEye failed to parse {iupac_name} as a IUPAC name"
            )
        oechem.OETriposAtomNames(oemol)
        result = oechem.OEAddExplicitHydrogens(oemol)
        if not result:
            raise Exception("Addition of explicit hydrogens failed in from_iupac")

        molecule = self.from_openeye(
            oemol, allow_undefined_stereo=allow_undefined_stereo, _cls=_cls, **kwargs
        )

        return molecule

    def generate_conformers(
        self, molecule, n_conformers=1, rms_cutoff=None, clear_existing=True
    ):
        r"""
        Generate molecule conformers using OpenEye Omega.

        .. warning :: This API is experimental and subject to change.

        .. todo ::

            * which parameters should we expose? (or can we implement a general system with \*\*kwargs?)
            * will the coordinates be returned in the OpenFF Molecule's own indexing system? Or is there a chance that
              they'll get reindexed when we convert the input into an OEmol?

        Parameters
        ----------
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
        omega.SetEnergyWindow(15.0)  # unit?
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

        molecule2 = self.from_openeye(
            oemol, allow_undefined_stereo=True, _cls=molecule.__class__
        )

        if clear_existing:
            molecule._conformers = list()

        for conformer in molecule2._conformers:
            molecule._add_conformer(conformer)

    def apply_elf_conformer_selection(
        self,
        molecule: "Molecule",
        percentage: float = 2.0,
        limit: int = 10,
    ):
        """Applies the `ELF method
        <https://docs.eyesopen.com/toolkits/python/quacpactk/molchargetheory.html#elf-conformer-selection>`_
        to select a set of diverse
        conformers which have minimal electrostatically strongly interacting functional
        groups from a molecules conformers.

        Notes
        -----
        * The input molecule should have a large set of conformers already
          generated to select the ELF conformers from.
        * The selected conformers will be retained in the `molecule.conformers` list
          while unselected conformers will be discarded.

        See Also
        --------
        RDKitToolkitWrapper.apply_elf_conformer_selection

        Parameters
        ----------
        molecule
            The molecule which contains the set of conformers to select from.
        percentage
            The percentage of conformers with the lowest electrostatic interaction
            energies to greedily select from.
        limit
            The maximum number of conformers to select.
        """

        from openeye import oechem, oequacpac

        if molecule.n_conformers == 0:
            return

        oe_molecule = molecule.to_openeye()

        # Select a subset of the OMEGA generated conformers using the ELF10 method.
        oe_elf_options = oequacpac.OEELFOptions()
        oe_elf_options.SetElfLimit(limit)
        oe_elf_options.SetPercent(percentage)

        oe_elf = oequacpac.OEELF(oe_elf_options)

        output_stream = oechem.oeosstream()

        oechem.OEThrow.SetOutputStream(output_stream)
        oechem.OEThrow.Clear()

        status = oe_elf.Select(oe_molecule)

        oechem.OEThrow.SetOutputStream(oechem.oeerr)

        output_string = output_stream.str().decode("UTF-8")
        output_string = output_string.replace("Warning: ", "")
        output_string = re.sub("^: +", "", output_string, flags=re.MULTILINE)
        output_string = re.sub("\n$", "", output_string)

        # Check to make sure the call to OE was succesful, and re-route any
        # non-fatal warnings to the correct logger.
        if not status:
            raise RuntimeError("\n" + output_string)
        elif len(output_string) > 0:
            logger.warning(output_string)

        # Extract and store the ELF conformers on the input molecule.
        conformers = []

        for oe_conformer in oe_molecule.GetConfs():

            conformer = np.zeros((oe_molecule.NumAtoms(), 3))

            for atom_index, coordinates in oe_conformer.GetCoords().items():
                conformer[atom_index, :] = coordinates

            conformers.append(conformer * unit.angstrom)

        molecule._conformers = conformers

    def assign_partial_charges(
        self,
        molecule,
        partial_charge_method=None,
        use_conformers=None,
        strict_n_conformers=False,
        _cls=None,
    ):
        """
        Compute partial charges with OpenEye quacpac, and assign
        the new values to the partial_charges attribute.

        .. warning :: This API is experimental and subject to change.

        .. todo ::

           * Should the default be ELF?
           * Can we expose more charge models?


        Parameters
        ----------
        molecule : openff.toolkit.topology.Molecule
            Molecule for which partial charges are to be computed
        partial_charge_method : str, optional, default=None
            The charge model to use. One of ['amberff94', 'mmff', 'mmff94', `am1-mulliken`, 'am1bcc',
            'am1bccnosymspt', 'am1bccelf10']
            If None, 'am1-mulliken' will be used.
        use_conformers : iterable of simtk.unit.Quantity-wrapped numpy arrays, each with
            shape (n_atoms, 3) and dimension of distance. Optional, default = None
            Coordinates to use for partial charge calculation. If None, an appropriate number
            of conformers will be generated.
        strict_n_conformers : bool, default=False
            Whether to raise an exception if an invalid number of conformers is provided for the
            given charge method.
            If this is False and an invalid number of conformers is found, a warning will be raised.
        _cls : class
            Molecule constructor

        Raises
        ------
        ChargeMethodUnavailableError if the requested charge method can not be handled by this toolkit

        ChargeCalculationError if the charge method is supported by this toolkit, but fails
        """

        import numpy as np
        from openeye import oechem, oequacpac

        from openff.toolkit.topology import Molecule

        SUPPORTED_CHARGE_METHODS = {
            "am1bcc": {
                "oe_charge_method": oequacpac.OEAM1BCCCharges,
                "min_confs": 1,
                "max_confs": 1,
                "rec_confs": 1,
            },
            "am1-mulliken": {
                "oe_charge_method": oequacpac.OEAM1Charges,
                "min_confs": 1,
                "max_confs": 1,
                "rec_confs": 1,
            },
            "gasteiger": {
                "oe_charge_method": oequacpac.OEGasteigerCharges,
                "min_confs": 0,
                "max_confs": 0,
                "rec_confs": 0,
            },
            "mmff94": {
                "oe_charge_method": oequacpac.OEMMFF94Charges,
                "min_confs": 0,
                "max_confs": 0,
                "rec_confs": 0,
            },
            "am1bccnosymspt": {
                "oe_charge_method": oequacpac.OEAM1BCCCharges,
                "min_confs": 1,
                "max_confs": 1,
                "rec_confs": 1,
            },
            "am1elf10": {
                "oe_charge_method": oequacpac.OEELFCharges(
                    oequacpac.OEAM1Charges(optimize=True, symmetrize=True), 10
                ),
                "min_confs": 1,
                "max_confs": None,
                "rec_confs": 500,
            },
            "am1bccelf10": {
                "oe_charge_method": oequacpac.OEAM1BCCELF10Charges,
                "min_confs": 1,
                "max_confs": None,
                "rec_confs": 500,
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

        if _cls is None:
            _cls = Molecule

        # Make a temporary copy of the molecule, since we'll be messing with its conformers
        mol_copy = _cls(molecule)

        if use_conformers is None:
            if charge_method["rec_confs"] == 0:
                mol_copy._conformers = None
            else:
                self.generate_conformers(
                    mol_copy,
                    n_conformers=charge_method["rec_confs"],
                    rms_cutoff=0.25 * unit.angstrom,
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

        oemol = mol_copy.to_openeye()

        errfs = oechem.oeosstream()
        oechem.OEThrow.SetOutputStream(errfs)
        oechem.OEThrow.Clear()

        # The OpenFF toolkit has always supported a version of AM1BCC with no geometry optimization
        # or symmetry correction. So we include this keyword to provide a special configuration of quacpac
        # if requested.
        if partial_charge_method == "am1bccnosymspt":
            optimize = False
            symmetrize = False
            quacpac_status = oequacpac.OEAssignCharges(
                oemol, charge_method["oe_charge_method"](optimize, symmetrize)
            )
        else:
            oe_charge_method = charge_method["oe_charge_method"]

            if callable(oe_charge_method):
                oe_charge_method = oe_charge_method()

            quacpac_status = oequacpac.OEAssignCharges(oemol, oe_charge_method)

        oechem.OEThrow.SetOutputStream(oechem.oeerr)  # restoring to original state
        # This logic handles errors encountered in #34, which can occur when using ELF10 conformer selection
        if not quacpac_status:

            oe_charge_engine = (
                oequacpac.OEAM1Charges
                if partial_charge_method == "am1elf10"
                else oequacpac.OEAM1BCCCharges
            )

            if "SelectElfPop: issue with removing trans COOH conformers" in (
                errfs.str().decode("UTF-8")
            ):
                logger.warning(
                    f"Warning: charge assignment involving ELF10 conformer selection failed due to a "
                    f"known bug (toolkit issue #346). Downgrading to {oe_charge_engine.__name__} "
                    f"charge assignment for this molecule. More information "
                    f"is available at https://github.com/openforcefield/openff-toolkit/issues/346"
                )
                quacpac_status = oequacpac.OEAssignCharges(oemol, oe_charge_engine())

        if quacpac_status is False:
            raise ChargeCalculationError(
                f'Unable to assign charges: {errfs.str().decode("UTF-8")}'
            )

        # Extract and return charges
        # TODO: Make sure atom mapping remains constant

        charges = unit.Quantity(
            np.zeros(shape=oemol.NumAtoms(), dtype=np.float64), unit.elementary_charge
        )
        for oeatom in oemol.GetAtoms():
            index = oeatom.GetIdx()
            charge = oeatom.GetPartialCharge()
            charge = charge * unit.elementary_charge
            charges[index] = charge

        molecule.partial_charges = charges

    def compute_partial_charges_am1bcc(
        self, molecule, use_conformers=None, strict_n_conformers=False
    ):
        """
        Compute AM1BCC partial charges with OpenEye quacpac. This function will attempt to use
        the OEAM1BCCELF10 charge generation method, but may print a warning and fall back to
        normal OEAM1BCC if an error is encountered. This error is known to occur with some
        carboxylic acids, and is under investigation by OpenEye.


        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        molecule : Molecule
            Molecule for which partial charges are to be computed
        use_conformers : iterable of simtk.unit.Quantity-wrapped numpy arrays, each with
            shape (n_atoms, 3) and dimension of distance. Optional, default = None
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
            "Use assign_partial_charges(partial_charge_method='am1bccelf10') instead.",
            DeprecationWarning,
        )
        self.assign_partial_charges(
            molecule,
            partial_charge_method="am1bccelf10",
            use_conformers=use_conformers,
            strict_n_conformers=strict_n_conformers,
        )
        return molecule.partial_charges

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
            The charge model to use. One of ['am1-wiberg', 'am1-wiberg-elf10',
            'pm3-wiberg', 'pm3-wiberg-elf10']. If None, 'am1-wiberg' will be used.
        use_conformers : iterable of simtk.unit.Quantity(np.array) with shape (n_atoms, 3) and
            dimension of distance, optional, default=None
            The conformers to use for fractional bond order calculation. If None, an
            appropriate number of conformers will be generated by an available
            ToolkitWrapper. If the chosen ``bond_order_model`` is an ELF variant, the ELF
            conformer selection method will be applied to the provided conformers.
        _cls : class
            Molecule constructor
        """
        from openeye import oechem, oequacpac

        if _cls is None:
            from openff.toolkit.topology.molecule import Molecule

            _cls = Molecule

        # Make a copy since we'll be messing with this molecule's conformers
        temp_mol = _cls(molecule)

        if bond_order_model is None:
            bond_order_model = "am1-wiberg"

        is_elf_method = bond_order_model in ["am1-wiberg-elf10", "pm3-wiberg-elf10"]

        if use_conformers is None:
            temp_mol.generate_conformers(
                n_conformers=1 if not is_elf_method else 500,
                # 0.05 is the recommended RMS when generating a 'Dense' amount of
                # conformers using Omega: https://docs.eyesopen.com/toolkits/python/
                # omegatk/OEConfGenConstants/OEFragBuilderMode.html.
                rms_cutoff=None if not is_elf_method else 0.05 * unit.angstrom,
            )
        else:
            temp_mol._conformers = None
            for conformer in use_conformers:
                temp_mol._add_conformer(conformer)
        if temp_mol.n_conformers == 0:
            raise Exception(
                "No conformers present in molecule submitted for fractional bond order calculation. Consider "
                "loading the molecule from a file with geometry already present or running "
                "molecule.generate_conformers() before calling molecule.compute_wiberg_bond_orders()"
            )

        if is_elf_method:
            # Apply the ELF10 conformer selection method.
            temp_mol.apply_elf_conformer_selection()

        # Set the options to use when computing the WBOs. This is based on example at
        # https://docs.eyesopen.com/toolkits/python/quacpactk/examples_summary_wibergbondorders.html
        am1 = oequacpac.OEAM1()

        am1results = oequacpac.OEAM1Results()
        am1options = am1.GetOptions()

        if bond_order_model.startswith("am1-wiberg"):
            am1options.SetSemiMethod(oequacpac.OEMethodType_AM1)
        elif bond_order_model.startswith("pm3-wiberg"):
            # TODO: Make sure that modifying am1options actually works
            am1options.SetSemiMethod(oequacpac.OEMethodType_PM3)
        else:
            raise ValueError(
                f"Bond order model '{bond_order_model}' is not supported by "
                f"OpenEyeToolkitWrapper. Supported models are ['am1-wiberg', "
                f"'am1-wiberg-elf10', 'pm3-wiberg', 'pm3-wiberg-elf10']."
            )

        # Convert the conformers into OE friendly objects to make setting them one
        # at a time easier.
        oe_conformers = [
            oechem.OEFloatArray(conformer.value_in_unit(unit.angstrom).flatten())
            for conformer in temp_mol.conformers
        ]

        oemol = self.to_openeye(temp_mol)
        bond_orders = defaultdict(list)

        for oe_conformer in oe_conformers:

            oemol.DeleteConfs()
            oemol.NewConf(oe_conformer)

            status = am1.CalcAM1(am1results, oemol)

            if status is False:

                raise Exception(
                    "Unable to assign charges (in the process of calculating "
                    "fractional bond orders)"
                )

            for bond in oemol.GetBonds():

                bond_orders[bond.GetIdx()].append(
                    am1results.GetBondOrder(bond.GetBgnIdx(), bond.GetEndIdx())
                )

        # TODO: Will bonds always map back to the same index? Consider doing a
        #       topology mapping.
        for bond_idx, conformer_bond_orders in bond_orders.items():

            # Get bond order
            order = np.mean(conformer_bond_orders)

            mol_bond = molecule._bonds[bond_idx]
            mol_bond.fractional_bond_order = order

    def get_tagged_smarts_connectivity(self, smarts):
        """
        Returns a tuple of tuples indicating connectivity between tagged atoms in a SMARTS string. Does not
        return bond order.

        Parameters
        ----------
        smarts : str
            The tagged SMARTS to analyze

        Returns
        -------
        unique_tags : tuple of int
            A sorted tuple of all unique tagged atom map indices.
        tagged_atom_connectivity : tuple of tuples of int, shape n_tagged_bonds x 2
            A tuple of tuples, where each inner tuple is a pair of tagged atoms
            (tag_idx_1, tag_idx_2) which are bonded. The inner tuples are ordered
            smallest-to-largest, and the tuple of tuples is ordered lexically. The
            return value for an improper torsion would be ((1, 2), (2, 3), (2, 4)).

        Raises
        ------
        SMIRKSParsingError
            If OpenEye toolkit was unable to parse the provided smirks/tagged smarts
        """
        from openeye import oechem

        from openff.toolkit.typing.chemistry import SMIRKSParsingError

        qmol = oechem.OEQMol()
        status = oechem.OEParseSmarts(qmol, smarts)
        if not status:
            raise SMIRKSParsingError(
                f"OpenEye Toolkit was unable to parse SMIRKS {smarts}"
            )

        unique_tags = set()
        connections = set()
        for at1 in qmol.GetAtoms():
            if at1.GetMapIdx() == 0:
                continue
            unique_tags.add(at1.GetMapIdx())
            for at2 in at1.GetAtoms():
                if at2.GetMapIdx() == 0:
                    continue
                cxn_to_add = sorted([at1.GetMapIdx(), at2.GetMapIdx()])
                connections.add(tuple(cxn_to_add))
        connections = tuple(sorted(list(connections)))
        unique_tags = tuple(sorted(list(unique_tags)))
        return tuple(unique_tags), tuple(connections)

    @staticmethod
    def _find_smarts_matches(
        oemol, smarts, aromaticity_model=DEFAULT_AROMATICITY_MODEL
    ):
        """Find all sets of atoms in the provided OpenEye molecule that match the provided SMARTS string.

        Parameters
        ----------
        oemol : openeye.oechem.OEMol or similar
            oemol to process with the SMIRKS in order to find matches
        smarts : str
            SMARTS string with any number of sequentially tagged atoms.
            If there are N tagged atoms numbered 1..N, the resulting matches will be N-tuples of
            atoms that match the corresponding tagged atoms.
        aromaticity_model : str, optional, default=None
            OpenEye aromaticity model designation as a string, such as ``OEAroModel_MDL``.
            Molecule is prepared with this aromaticity model prior to querying.

        Returns
        -------
        matches : list of tuples of atoms indices within the ``oemol``
            matches[index] is an N-tuple of atom numbers from the ``oemol``
            Matches are returned in no guaranteed order.
            # TODO: What is returned if no matches are found? An empty list, or None?
            # TODO: Ensure that SMARTS numbers 1, 2, 3... are rendered into order of returned
            #    matches indexed by 0, 1, 2...

        .. notes ::

           * Raises ``LicenseError`` if valid OpenEye tools license is not found, rather than
               causing program to terminate
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

        # Apply aromaticity model
        if type(aromaticity_model) == str:
            # Check if the user has provided a manually-specified aromaticity_model
            if hasattr(oechem, aromaticity_model):
                oearomodel = getattr(oechem, aromaticity_model)
            else:
                raise ValueError(
                    "Error: provided aromaticity model not recognized by oechem."
                )
        else:
            raise ValueError("Error: provided aromaticity model must be a string.")

        # OEPrepareSearch will clobber our desired aromaticity model if we don't sync up mol
        # and qmol ahead of time.
        # Prepare molecule
        oechem.OEClearAromaticFlags(mol)
        oechem.OEAssignAromaticFlags(mol, oearomodel)

        # If aromaticity model was provided, prepare query molecule
        oechem.OEClearAromaticFlags(qmol)
        oechem.OEAssignAromaticFlags(qmol, oearomodel)
        oechem.OEAssignHybridization(mol)
        oechem.OEAssignHybridization(qmol)

        # Build list of matches
        # TODO: The MoleculeImage mapping should preserve ordering of template molecule for equivalent atoms
        #       and speed matching for larger molecules.
        unique = False  # We require all matches, not just one of each kind
        substructure_search = OESubSearch(qmol)
        substructure_search.SetMaxMatches(0)
        oechem.OEPrepareSearch(mol, substructure_search)
        matches = list()
        for match in substructure_search.Match(mol, unique):
            # Compile list of atom indices that match the pattern tags
            atom_indices = dict()
            for matched_atom in match.GetAtoms():
                if matched_atom.pattern.GetMapIdx() != 0:
                    atom_indices[
                        matched_atom.pattern.GetMapIdx() - 1
                    ] = matched_atom.target.GetIdx()
            # Compress into list
            atom_indices = [atom_indices[index] for index in range(len(atom_indices))]
            # Convert to tuple
            matches.append(tuple(atom_indices))
        return matches

    def find_smarts_matches(self, molecule, smarts, aromaticity_model="OEAroModel_MDL"):
        """
        Find all SMARTS matches for the specified molecule, using the specified aromaticity model.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        molecule : openff.toolkit.topology.Molecule
            The molecule for which all specified SMARTS matches are to be located
        smarts : str
            SMARTS string with optional SMIRKS-style atom tagging
        aromaticity_model : str, optional, default='OEAroModel_MDL'
            Molecule is prepared with this aromaticity model prior to querying.

        .. note :: Currently, the only supported ``aromaticity_model`` is ``OEAroModel_MDL``

        """
        oemol = self.to_openeye(molecule)
        return self._find_smarts_matches(
            oemol, smarts, aromaticity_model=aromaticity_model
        )


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
