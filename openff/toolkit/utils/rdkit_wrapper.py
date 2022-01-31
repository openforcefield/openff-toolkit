"""
Wrapper class providing a minimal consistent interface to the `RDKit <http://www.rdkit.org/>`.
"""

__all__ = ("RDKitToolkitWrapper",)

import copy
import importlib
import itertools
import logging
import tempfile
from typing import TYPE_CHECKING, List, Optional, Tuple

import numpy as np

try:
    from openmm import unit
except ImportError:
    from simtk import unit

if TYPE_CHECKING:
    from openff.toolkit.topology.molecule import Molecule

from openff.toolkit.utils import base_wrapper
from openff.toolkit.utils.constants import DEFAULT_AROMATICITY_MODEL
from openff.toolkit.utils.exceptions import (
    ChargeMethodUnavailableError,
    ConformerGenerationError,
    SMILESParseError,
    ToolkitUnavailableException,
    UndefinedStereochemistryError,
)

# =============================================================================================
# CONFIGURE LOGGER
# =============================================================================================

logger = logging.getLogger(__name__)

# =============================================================================================
# IMPLEMENTATION
# =============================================================================================


def normalize_file_format(file_format):
    return file_format.upper()


def _require_text_file_obj(file_obj):
    try:
        file_obj.write("")
    except TypeError:
        # Switch to a ValueError and use a more informative exception
        # message to match RDKit's toolkit writer.
        raise ValueError(
            "Need a text mode file object like StringIO or a file opened with mode 't'"
        ) from None


class RDKitToolkitWrapper(base_wrapper.ToolkitWrapper):
    """
    RDKit toolkit wrapper

    .. warning :: This API is experimental and subject to change.
    """

    _toolkit_name = "The RDKit"
    _toolkit_installation_instructions = (
        "A conda-installable version of the free and open source RDKit cheminformatics "
        "toolkit can be found at: https://anaconda.org/rdkit/rdkit"
    )

    def __init__(self):
        super().__init__()

        self._toolkit_file_read_formats = ["SDF", "MOL", "SMI"]  # TODO: Add TDT support

        if not self.is_available():
            raise ToolkitUnavailableException(
                f"The required toolkit {self._toolkit_name} is not "
                f"available. {self._toolkit_installation_instructions}"
            )
        else:
            from rdkit import __version__ as rdkit_version

            self._toolkit_version = rdkit_version

            from rdkit import Chem

            # we have to make sure the toolkit can be loaded before formatting this dict
            # Note any new file write formats should be added here only
            self._toolkit_file_write_formats = {
                "SDF": Chem.SDWriter,
                "MOL": Chem.SDWriter,
                "SMI": None,  # Special support to use to_smiles() instead of RDKit's SmilesWriter
                "PDB": Chem.PDBWriter,
                "TDT": Chem.TDTWriter,
            }

    @property
    def toolkit_file_write_formats(self):
        """
        List of file formats that this toolkit can write.
        """
        return list(self._toolkit_file_write_formats.keys())

    @classmethod
    def is_available(cls):
        """
        Check whether the RDKit toolkit can be imported

        Returns
        -------
        is_installed : bool
            True if RDKit is installed, False otherwise.

        """
        if cls._is_available is None:
            try:
                importlib.import_module("rdkit", "Chem")
            except ImportError:
                cls._is_available = False
            else:
                cls._is_available = True
        return cls._is_available

    def from_object(self, obj, allow_undefined_stereo=False, _cls=None):
        """
        If given an rdchem.Mol (or rdchem.Mol-derived object), this function will load it into an
        openff.toolkit.topology.molecule. Otherwise, it will return False.

        Parameters
        ----------
        obj : A rdchem.Mol-derived object
            An object to be type-checked and converted into a Molecule, if possible.
        allow_undefined_stereo : bool, default=False
            Whether to accept molecules with undefined stereocenters. If False,
            an exception will be raised if a molecule with undefined stereochemistry
            is passed into this function.
        _cls : class
            Molecule constructor

        Returns
        -------
        Molecule or False
            An openff.toolkit.topology.molecule Molecule.

        Raises
        ------
        NotImplementedError
            If the object could not be converted into a Molecule.
        """
        # TODO: Add tests for the from_object functions
        from rdkit import Chem

        if _cls is None:
            from openff.toolkit.topology.molecule import Molecule

            _cls = Molecule
        if isinstance(obj, Chem.rdchem.Mol):
            return _cls.from_rdkit(obj, allow_undefined_stereo=allow_undefined_stereo)
        raise NotImplementedError(
            "Cannot create Molecule from {} object".format(type(obj))
        )

    def from_pdb_and_smiles(
        self, file_path, smiles, allow_undefined_stereo=False, _cls=None
    ):
        """
        Create a Molecule from a pdb file and a SMILES string using RDKit.

        Requires RDKit to be installed.

        The molecule is created and sanitised based on the SMILES string, we then find a mapping
        between this molecule and one from the PDB based only on atomic number and connections.
        The SMILES molecule is then reindexed to match the PDB, the conformer is attached, and the
        molecule returned.

        Note that any stereochemistry in the molecule is set by the SMILES, and not the coordinates
        of the PDB.

        Parameters
        ----------
        file_path: str
            PDB file path
        smiles : str
            a valid smiles string for the pdb, used for stereochemistry, formal charges, and bond order
        allow_undefined_stereo : bool, default=False
            If false, raises an exception if SMILES contains undefined stereochemistry.
        _cls : class
            Molecule constructor

        Returns
        --------
        molecule : openff.toolkit.Molecule (or _cls() type)
            An OFFMol instance with ordering the same as used in the PDB file.

        Raises
        ------
        InvalidConformerError : if the SMILES and PDB molecules are not isomorphic.
        """

        from rdkit import Chem

        # Make the molecule from smiles
        offmol = self.from_smiles(
            smiles, allow_undefined_stereo=allow_undefined_stereo, _cls=_cls
        )

        # Make another molecule from the PDB. We squelch stereo errors here, since
        # RDKit's PDB loader doesn't attempt to perceive stereochemistry, bond order,
        # or formal charge (and we don't need those here).
        prev_log_level = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)
        pdbmol = self.from_rdkit(
            Chem.MolFromPDBFile(file_path, removeHs=False),
            allow_undefined_stereo=True,
            hydrogens_are_explicit=True,
            _cls=_cls,
        )
        logger.setLevel(prev_log_level)

        # check isomorphic and get the mapping if true the mapping will be
        # Dict[pdb_index: offmol_index] sorted by pdb_index
        isomorphic, mapping = _cls.are_isomorphic(
            pdbmol,
            offmol,
            return_atom_map=True,
            aromatic_matching=False,
            formal_charge_matching=False,
            bond_order_matching=False,
            atom_stereochemistry_matching=False,
            bond_stereochemistry_matching=False,
        )

        if mapping is not None:
            new_mol = offmol.remap(mapping)

            # the pdb conformer is in the correct order so just attach it here
            new_mol._add_conformer(pdbmol.conformers[0])

            return new_mol

        else:
            from openff.toolkit.topology.molecule import InvalidConformerError

            raise InvalidConformerError("The PDB and SMILES structures do not match.")

    def _process_sdf_supplier(self, sdf_supplier, allow_undefined_stereo, _cls):
        "Helper function to process RDKit molecules from an SDF input source"
        from rdkit import Chem

        for rdmol in sdf_supplier:
            if rdmol is None:
                continue

            # Sanitize the molecules (fails on nitro groups)
            try:
                Chem.SanitizeMol(
                    rdmol,
                    Chem.SANITIZE_ALL
                    ^ Chem.SANITIZE_SETAROMATICITY
                    ^ Chem.SANITIZE_ADJUSTHS,
                )
                Chem.AssignStereochemistryFrom3D(rdmol)
            except ValueError as e:
                logger.warning(rdmol.GetProp("_Name") + " " + str(e))
                continue
            Chem.SetAromaticity(rdmol, Chem.AromaticityModel.AROMATICITY_MDL)
            mol = self.from_rdkit(
                rdmol, allow_undefined_stereo=allow_undefined_stereo, _cls=_cls
            )
            yield mol

    def from_file(
        self, file_path, file_format, allow_undefined_stereo=False, _cls=None
    ):
        """
        Create an openff.toolkit.topology.Molecule from a file using this toolkit.



        Parameters
        ----------
        file_path : str
            The file to read the molecule from
        file_format : str
            Format specifier, usually file suffix (eg. 'MOL2', 'SMI')
            Note that not all toolkits support all formats. Check
            ToolkitWrapper.toolkit_file_read_formats for details.
        allow_undefined_stereo : bool, default=False
            If false, raises an exception if oemol contains undefined stereochemistry.
        _cls : class
            Molecule constructor
        Returns
        -------
        molecules : iterable of Molecules
            a list of Molecule objects is returned.

        """
        from rdkit import Chem

        file_format = normalize_file_format(file_format)

        mols = list()
        if (file_format == "MOL") or (file_format == "SDF"):
            sdf_supplier = Chem.ForwardSDMolSupplier(
                file_path, removeHs=False, sanitize=False, strictParsing=True
            )
            mols.extend(
                self._process_sdf_supplier(
                    sdf_supplier,
                    allow_undefined_stereo=allow_undefined_stereo,
                    _cls=_cls,
                )
            )

        elif file_format == "SMI":
            # TODO: We have to do some special stuff when we import SMILES (currently
            # just adding H's, but could get fancier in the future). It might be
            # worthwhile to parse the SMILES file ourselves and pass each SMILES
            # through the from_smiles function instead
            for rdmol in Chem.SmilesMolSupplier(file_path, titleLine=False):
                if rdmol is None:
                    # Skip any lines that could not be processed.
                    # This is consistent with the SDF reader and with
                    # the OpenEye toolkit wrapper.
                    continue
                rdmol = Chem.AddHs(rdmol)
                mol = self.from_rdkit(
                    rdmol, allow_undefined_stereo=allow_undefined_stereo, _cls=_cls
                )
                mols.append(mol)

        elif file_format == "PDB":
            raise Exception(
                "RDKit can not safely read PDBs on their own. Information about bond order "
                "and aromaticity is likely to be lost. To read a PDB using RDKit use "
                "Molecule.from_pdb_and_smiles()"
            )
            # TODO: See if we can implement PDB+mol/smi combinations to get complete bond information.
            #  testing to see if we can make a molecule from smiles and then use the PDB
            #  conformer as the geometry and just reorder the molecule
            # https://github.com/openforcefield/openff-toolkit/issues/121
            # rdmol = Chem.MolFromPDBFile(file_path, removeHs=False)
            # mol = Molecule.from_rdkit(rdmol, _cls=_cls)
            # mols.append(mol)
            # TODO: Add SMI, TDT(?) support

        else:
            raise ValueError(f"Unsupported file format: {file_format}")

        return mols

    def from_file_obj(
        self, file_obj, file_format, allow_undefined_stereo=False, _cls=None
    ):
        """
        Return an openff.toolkit.topology.Molecule from a file-like object using this toolkit.

        A file-like object is an object with a ".read()" method.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        file_obj : file-like object
            The file-like object to read the molecule from
        file_format : str
            Format specifier, usually file suffix (eg. 'MOL2', 'SMI')
            Note that not all toolkits support all formats. Check
            ToolkitWrapper.toolkit_file_read_formats for details.
        allow_undefined_stereo : bool, default=False
            If false, raises an exception if oemol contains undefined stereochemistry.
        _cls : class
            Molecule constructor
        Returns
        -------
        molecules : Molecule or list of Molecules
            a list of Molecule objects is returned.

        """
        from rdkit import Chem

        mols = []

        file_format = normalize_file_format(file_format)

        if (file_format == "MOL") or (file_format == "SDF"):
            # TODO: Iterate over all mols in file_data
            sdf_supplier = Chem.ForwardSDMolSupplier(
                file_obj, removeHs=False, sanitize=False, strictParsing=True
            )
            mols.extend(
                self._process_sdf_supplier(
                    sdf_supplier,
                    allow_undefined_stereo=allow_undefined_stereo,
                    _cls=_cls,
                )
            )

        elif file_format == "SMI":
            # There's no good way to create a SmilesMolSuppler from a string
            # other than to use a temporary file.
            with tempfile.NamedTemporaryFile(suffix=".smi") as tmpfile:
                content = file_obj.read()
                if isinstance(content, str):
                    # Backwards compatibility. Older versions of OpenFF supported
                    # file objects in "t"ext mode, but not file objects
                    # in "b"inary mode. Now we expect all input file objects
                    # to handle binary mode, but don't want to break older code.
                    tmpfile.write(content.encode("utf8"))
                else:
                    tmpfile.write(content)
                tmpfile.flush()
                return self.from_file(
                    tmpfile.name,
                    "SMI",
                    allow_undefined_stereo=allow_undefined_stereo,
                    _cls=_cls,
                )

        elif file_format == "PDB":
            raise Exception(
                "RDKit can not safely read PDBs on their own. Information about bond order and aromaticity "
                "is likely to be lost. To read a PDB using RDKit use Molecule.from_pdb_and_smiles()"
            )
            # TODO: See if we can implement PDB+mol/smi combinations to get complete bond information.
            # https://github.com/openforcefield/openff-toolkit/issues/121
            # file_data = file_obj.read()
            # rdmol = Chem.MolFromPDBBlock(file_data)
            # mol = Molecule.from_rdkit(rdmol, _cls=_cls)
            # mols.append(mol)

        else:
            raise ValueError(f"Unsupported file format: {file_format}")
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
        file_format = normalize_file_format(file_format)
        _require_text_file_obj(file_obj)

        if file_format == "SMI":
            # Special case for SMILES
            smiles = self.to_smiles(molecule)
            name = molecule.name
            if name is not None:
                output_line = f"{smiles} {name}\n"
            else:
                output_line = f"{smiles}\n"
            file_obj.write(output_line)
        else:
            try:
                writer_func = self._toolkit_file_write_formats[file_format]
            except KeyError:
                raise ValueError(f"Unsupported file format: {file_format})") from None
            rdmol = self.to_rdkit(molecule)
            writer = writer_func(file_obj)
            try:
                writer.write(rdmol)
            finally:
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

        # open a file object and pass to the object writer
        with open(file_path, "w") as file_obj:
            self.to_file_obj(
                molecule=molecule, file_obj=file_obj, file_format=file_format
            )

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
        from rdkit import Chem
        from rdkit.Chem.EnumerateStereoisomers import (  # type: ignore[import]
            EnumerateStereoisomers,
            StereoEnumerationOptions,
        )

        # create the molecule
        rdmol = self.to_rdkit(molecule=molecule)

        # in case any bonds/centers are missing stereo chem flag it here
        Chem.AssignStereochemistry(
            rdmol, cleanIt=True, force=True, flagPossibleStereoCenters=True
        )
        Chem.FindPotentialStereoBonds(rdmol)

        # set up the options
        stereo_opts = StereoEnumerationOptions(
            tryEmbedding=rationalise,
            onlyUnassigned=undefined_only,
            maxIsomers=max_isomers,
        )

        isomers = tuple(EnumerateStereoisomers(rdmol, options=stereo_opts))

        molecules = []
        for isomer in isomers:
            # isomer has CIS/TRANS tags so convert back to E/Z
            Chem.SetDoubleBondNeighborDirections(isomer)
            Chem.AssignStereochemistry(isomer, force=True, cleanIt=True)
            mol = self.from_rdkit(isomer, _cls=molecule.__class__)
            if mol != molecule:
                molecules.append(mol)

        return molecules

    def enumerate_tautomers(self, molecule, max_states=20):
        """
        Enumerate the possible tautomers of the current molecule.

        Parameters
        ----------
        molecule: openff.toolkit.topology.Molecule
            The molecule whose state we should enumerate

        max_states: int optional, default=20
            The maximum amount of molecules that should be returned

        Returns
        -------
        molecules: List[openff.toolkit.topology.Molecule]
            A list of openff.toolkit.topology.Molecule instances not including the input molecule.
        """

        from rdkit import Chem
        from rdkit.Chem.MolStandardize import rdMolStandardize  # type: ignore[import]

        enumerator = rdMolStandardize.TautomerEnumerator()
        enumerator.SetMaxTautomers(max_states)
        rdmol = Chem.RemoveHs(molecule.to_rdkit())

        tautomers = enumerator.Enumerate(rdmol)

        # make a list of OpenFF molecules excluding the input molecule
        molecules = []
        for taut in tautomers:
            taut_hs = Chem.AddHs(taut)
            mol = self.from_smiles(
                Chem.MolToSmiles(taut_hs), allow_undefined_stereo=True
            )
            if mol != molecule:
                molecules.append(mol)

        return molecules[:max_states]

    def canonical_order_atoms(self, molecule):
        """
        Canonical order the atoms in the molecule using the RDKit.

        Parameters
        ----------
        molecule: openff.toolkit.topology.Molecule
            The input molecule

         Returns
        -------
        molecule : openff.toolkit.topology.Molecule
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

    def to_smiles(self, molecule, isomeric=True, explicit_hydrogens=True, mapped=False):
        """
        Uses the RDKit toolkit to convert a Molecule into a SMILES string.
        A partially mapped smiles can also be generated for atoms of interest by supplying
        an `atom_map` to the properties dictionary.

        Parameters
        ----------
        molecule : An openff.toolkit.topology.Molecule
            The molecule to convert into a SMILES.
        isomeric: bool optional, default= True
            return an isomeric smiles
        explicit_hydrogens: bool optional, default=True
            return a smiles string containing all hydrogens explicitly
        mapped: bool optional, default=False
            return a explicit hydrogen mapped smiles, the atoms to be mapped can be controlled by
            supplying an atom map into the properties dictionary. If no mapping is passed all
            atoms will be mapped in order, else an atom map dictionary from the current atom
            index to the map id should be supplied with no duplicates. The map ids (values) should
            start from 0 or 1.

        Returns
        -------
        smiles : str
            The SMILES of the input molecule.
        """
        from rdkit import Chem

        rdmol = self.to_rdkit(molecule)

        if not explicit_hydrogens:
            # remove the hydrogens from the molecule
            rdmol = Chem.RemoveHs(rdmol)

        if mapped:
            assert explicit_hydrogens is True, (
                "Mapped smiles require all hydrogens and "
                "stereochemistry to be defined to retain order"
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
                # now we need to add the indexing to the rdmol to get it in the smiles
                for atom in rdmol.GetAtoms():
                    # the mapping must start from 1, as RDKit uses 0 to represent no mapping.
                    atom.SetAtomMapNum(atom.GetIdx() + 1)
            else:
                for atom in rdmol.GetAtoms():
                    try:
                        # try to set the atom map
                        map_idx = atom_map[atom.GetIdx()]
                        atom.SetAtomMapNum(map_idx)
                    except KeyError:
                        continue

        return Chem.MolToSmiles(
            rdmol, isomericSmiles=isomeric, allHsExplicit=explicit_hydrogens
        )

    def from_smiles(
        self,
        smiles,
        hydrogens_are_explicit=False,
        allow_undefined_stereo=False,
        _cls=None,
    ):
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
        _cls : class
            Molecule constructor

        Returns
        -------
        molecule : openff.toolkit.topology.Molecule
            An OpenFF style molecule.
        """
        from rdkit import Chem

        rdmol = Chem.MolFromSmiles(smiles, sanitize=False)
        if rdmol is None:
            raise SMILESParseError("Unable to parse the SMILES string")

        # strip the atom map from the molecule if it has one
        # so we don't affect the sterochemistry tags
        for atom in rdmol.GetAtoms():
            if atom.GetAtomMapNum() != 0:
                # set the map back to zero but hide the index in the atom prop data
                atom.SetProp("_map_idx", str(atom.GetAtomMapNum()))
                # set it back to zero
                atom.SetAtomMapNum(0)

        # Chem.SanitizeMol calls updatePropertyCache so we don't need to call it ourselves
        # https://www.rdkit.org/docs/cppapi/namespaceRDKit_1_1MolOps.html#a8d831787aaf2d65d9920c37b25b476f5
        Chem.SanitizeMol(
            rdmol,
            Chem.SANITIZE_ALL ^ Chem.SANITIZE_ADJUSTHS ^ Chem.SANITIZE_SETAROMATICITY,
        )
        Chem.SetAromaticity(rdmol, Chem.AromaticityModel.AROMATICITY_MDL)

        # Chem.MolFromSmiles adds bond directions (i.e. ENDDOWNRIGHT/ENDUPRIGHT), but
        # doesn't set bond.GetStereo(). We need to call AssignStereochemistry for that.
        Chem.AssignStereochemistry(rdmol)

        # Throw an exception/warning if there is unspecified stereochemistry.
        if not allow_undefined_stereo:
            self._detect_undefined_stereo(
                rdmol, err_msg_prefix="Unable to make OFFMol from SMILES: "
            )

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
                        f"Molecule.from_rdkit() to create the desired OFFMol."
                    )

        molecule = self.from_rdkit(
            rdmol,
            _cls=_cls,
            allow_undefined_stereo=allow_undefined_stereo,
            hydrogens_are_explicit=hydrogens_are_explicit,
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

        from rdkit import Chem

        # this seems to always remove the hydrogens
        rdmol = Chem.MolFromInchi(inchi, sanitize=False, removeHs=False)

        # try and catch an InChI parsing error
        if rdmol is None:
            raise RuntimeError(
                "There was an issue parsing the InChI string, please check and try again."
            )

        # process the molecule
        # TODO do we need this with inchi?
        rdmol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(
            rdmol,
            Chem.SANITIZE_ALL ^ Chem.SANITIZE_ADJUSTHS ^ Chem.SANITIZE_SETAROMATICITY,
        )
        Chem.SetAromaticity(rdmol, Chem.AromaticityModel.AROMATICITY_MDL)

        # add hydrogens back here
        rdmol = Chem.AddHs(rdmol)

        molecule = self.from_rdkit(
            rdmol, allow_undefined_stereo=allow_undefined_stereo, _cls=_cls
        )

        return molecule

    def generate_conformers(
        self, molecule, n_conformers=1, rms_cutoff=None, clear_existing=True, _cls=None
    ):
        r"""
        Generate molecule conformers using RDKit.

        .. warning :: This API is experimental and subject to change.

        .. todo ::

           * which parameters should we expose? (or can we implement a general system with \*\*kwargs?)
           * will the coordinates be returned in the OpenFF Molecule's own indexing system?
               Or is there a chance that they'll get reindexed when we convert the input into an RDMol?

        Parameters
        ----------
        molecule : a :class:`Molecule`
            The molecule to generate conformers for.
        n_conformers : int, default=1
            Maximum number of conformers to generate.
        rms_cutoff : openmm.Quantity-wrapped float, in units of distance, optional, default=None
            The minimum RMS value at which two conformers are considered redundant and one is deleted.
            If None, the cutoff is set to 1 Angstrom

        clear_existing : bool, default=True
            Whether to overwrite existing conformers for the molecule.
        _cls : class
            Molecule constructor

        """
        from rdkit.Chem import AllChem

        if rms_cutoff is None:
            rms_cutoff = 1.0 * unit.angstrom
        rdmol = self.to_rdkit(molecule)
        # TODO: This generates way more conformations than omega, given the same
        # nConfs and RMS threshold. Is there some way to set an energy cutoff as well?
        conformer_generation_status = AllChem.EmbedMultipleConfs(
            rdmol,
            numConfs=n_conformers,
            pruneRmsThresh=rms_cutoff.value_in_unit(unit.angstrom),
            randomSeed=1,
            # params=AllChem.ETKDG()
        )

        if not conformer_generation_status:
            raise ConformerGenerationError("RDKit conformer generation failed.")

        molecule2 = self.from_rdkit(
            rdmol, allow_undefined_stereo=True, _cls=molecule.__class__
        )

        if clear_existing:
            molecule._conformers = list()

        for conformer in molecule2._conformers:
            molecule._add_conformer(conformer)

    def assign_partial_charges(
        self,
        molecule,
        partial_charge_method=None,
        use_conformers=None,
        strict_n_conformers=False,
        normalize_partial_charges=True,
        _cls=None,
    ):
        """
        Compute partial charges with RDKit, and assign
        the new values to the partial_charges attribute.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        molecule : openff.toolkit.topology.Molecule
            Molecule for which partial charges are to be computed
        partial_charge_method : str, optional, default=None
            The charge model to use. One of ['mmff94']. If None, 'mmff94' will be used.

            * 'mmff94': Applies partial charges using the Merck Molecular Force Field
                        (MMFF). This method does not make use of conformers, and hence
                        ``use_conformers`` and ``strict_n_conformers`` will not impact
                        the partial charges produced.
        use_conformers : iterable of openmm.unit.Quantity-wrapped numpy arrays, each with
            shape (n_atoms, 3) and dimension of distance. Optional, default = None
            Coordinates to use for partial charge calculation. If None, an appropriate number of
            conformers will be generated.
        strict_n_conformers : bool, default=False
            Whether to raise an exception if an invalid number of conformers is provided for
            the given charge method.
            If this is False and an invalid number of conformers is found, a warning will be raised.
        normalize_partial_charges : bool, default=True
            Whether to offset partial charges so that they sum to the total formal charge of the molecule.
            This is used to prevent accumulation of rounding errors when the partial charge generation method has
            low precision.
        _cls : class
            Molecule constructor

        Raises
        ------
        ChargeMethodUnavailableError if the requested charge method can not be handled by this toolkit

        ChargeCalculationError if the charge method is supported by this toolkit, but fails
        """

        import numpy as np
        from rdkit.Chem import AllChem

        SUPPORTED_CHARGE_METHODS = {"mmff94"}

        if partial_charge_method is None:
            partial_charge_method = "mmff94"

        partial_charge_method = partial_charge_method.lower()

        if partial_charge_method not in SUPPORTED_CHARGE_METHODS:
            raise ChargeMethodUnavailableError(
                f"partial_charge_method '{partial_charge_method}' is not available from RDKitToolkitWrapper. "
                f"Available charge methods are {list(SUPPORTED_CHARGE_METHODS)} "
            )

        rdkit_molecule = molecule.to_rdkit()
        charges = None

        if partial_charge_method == "mmff94":

            mmff_properties = AllChem.MMFFGetMoleculeProperties(
                rdkit_molecule, "MMFF94"
            )
            charges = np.array(
                [
                    mmff_properties.GetMMFFPartialCharge(i)
                    for i in range(molecule.n_atoms)
                ]
            )

        molecule.partial_charges = charges * unit.elementary_charge

        if normalize_partial_charges:
            molecule._normalize_partial_charges()

    @classmethod
    def _elf_is_problematic_conformer(
        cls, molecule: "Molecule", conformer: unit.Quantity
    ) -> Tuple[bool, Optional[str]]:
        """A function which checks if a particular conformer is known to be problematic
        when computing ELF partial charges.

        Currently this includes conformers which:

        * contain a trans-COOH configuration. The trans conformer is discarded because
          it leads to strong electrostatic interactions when assigning charges, and these
          result in unreasonable charges. Downstream calculations have observed up to a
          4 log unit error in water-octanol logP calculations when using charges assigned
          from trans conformers.

        Returns
        -------
            A tuple of a bool stating whether the conformer is problematic and, if it
            is, a string message explaing why. If the conformer is not problematic, the
            second return value will be none.
        """
        from rdkit.Chem.rdMolTransforms import GetDihedralRad  # type: ignore[import]

        # Create a copy of the molecule which contains only this conformer.
        molecule_copy = copy.deepcopy(molecule)
        molecule_copy._conformers = [conformer]

        rdkit_molecule = molecule_copy.to_rdkit()

        # Check for trans-COOH configurations
        carboxylic_acid_matches = cls._find_smarts_matches(
            rdkit_molecule, "[#6X3:2](=[#8:1])(-[#8X2H1:3]-[#1:4])"
        )

        for match in carboxylic_acid_matches:

            dihedral_angle = GetDihedralRad(rdkit_molecule.GetConformer(0), *match)

            if dihedral_angle > np.pi / 2.0:
                # Discard the 'trans' conformer.
                return (
                    True,
                    "Molecules which contain COOH functional groups in a trans "
                    "configuration are discarded by the ELF method.",
                )

        return False, None

    @classmethod
    def _elf_prune_problematic_conformers(
        cls, molecule: "Molecule"
    ) -> List[unit.Quantity]:
        """A function which attempts to remove conformers which are known to be
        problematic when computing ELF partial charges.

        Currently this includes conformers which:

        * contain a trans-COOH configuration. These conformers ... TODO add reason.

        Notes
        -----
        * Problematic conformers are flagged by the
          ``RDKitToolkitWrapper._elf_is_problematic_conformer`` function.

        Returns
        -------
            The conformers to retain.
        """

        valid_conformers = []

        for i, conformer in enumerate(molecule.conformers):

            is_problematic, reason = cls._elf_is_problematic_conformer(
                molecule, conformer
            )

            if is_problematic:
                logger.warning(f"Discarding conformer {i}: {reason}")
            else:
                valid_conformers.append(conformer)

        return valid_conformers

    @classmethod
    def _elf_compute_electrostatic_energy(
        cls, molecule: "Molecule", conformer: unit.Quantity
    ) -> float:
        """Computes the 'electrostatic interaction energy' of a particular conformer
        of a molecule.

        The energy is computed as the sum of ``|q_i * q_j| * r_ij^-1`` over all pairs
        of atoms (i, j) excluding 1-2 and 1-3 terms, where q_i is the partial charge
        of atom i and r_ij the Euclidean distance between atoms i and j.

        Notes
        -----
        * The partial charges will be taken from the molecule directly.

        Parameters
        ----------
        molecule
            The molecule containing the partial charges.
        conformer
            The conformer to compute the energy of. This should be a unit wrapped
            numpy array with shape=(n_atoms, 3) with units compatible with angstroms.

        Returns
        -------
            The electrostatic interaction energy in units of [e^2 / Angstrom].
        """

        if molecule.partial_charges is None:
            raise ValueError("The molecule has no partial charges assigned.")

        partial_charges = np.abs(
            molecule.partial_charges.value_in_unit(unit.elementary_charge)
        ).reshape(-1, 1)

        # Build an exclusion list for 1-2 and 1-3 interactions.
        excluded_x, excluded_y = zip(
            *{
                *[(bond.atom1_index, bond.atom2_index) for bond in molecule.bonds],
                *[
                    (angle[0].molecule_atom_index, angle[-1].molecule_atom_index)
                    for angle in molecule.angles
                ],
            }
        )

        # Build the distance matrix between all pairs of atoms.
        coordinates = conformer.value_in_unit(unit.angstrom)

        # Equation 1: (a - b)^2 = a^2 - 2ab + b^2
        # distances_squared will eventually wind up as the squared distances
        # although it is first computed as the ab portion of Eq 1
        distances_squared = coordinates.dot(coordinates.T)
        # np.einsum is both faster than np.diag, and not read-only
        # we know that a^2 == b^2 == diag(ab)
        diag = np.einsum("ii->i", distances_squared)
        # we modify in-place so we can use the `diag` view
        # to make the diagonals 0
        distances_squared += distances_squared - diag - diag[..., np.newaxis]
        # Handle edge cases where the squared distance is slightly negative due to
        # precision issues
        diag[:] = -0.0  # this is somewhat faster than np.fill_diagonal
        distances = np.sqrt(-distances_squared)

        inverse_distances = np.reciprocal(
            distances, out=np.zeros_like(distances), where=~np.isclose(distances, 0.0)
        )

        # Multiply by the charge products.
        charge_products = partial_charges @ partial_charges.T

        charge_products[excluded_x, excluded_y] = 0.0
        charge_products[excluded_y, excluded_x] = 0.0

        interaction_energies = inverse_distances * charge_products

        return 0.5 * interaction_energies.sum()

    @classmethod
    def _elf_compute_rms_matrix(cls, molecule: "Molecule") -> np.ndarray:
        """Computes the symmetric RMS matrix of all conformers in a molecule taking
        only heavy atoms into account.

        Parameters
        ----------
        molecule
            The molecule containing the conformers.

        Returns
        -------
            The RMS matrix with shape=(n_conformers, n_conformers).
        """

        from rdkit import Chem
        from rdkit.Chem import AllChem

        rdkit_molecule: Chem.RWMol = Chem.RemoveHs(molecule.to_rdkit())

        n_conformers = len(molecule.conformers)

        # rdkit does not have conformer indices but conformer "ids"
        conformer_ids = [conf.GetId() for conf in rdkit_molecule.GetConformers()]

        # Compute the RMS matrix making sure to take into account any automorhism (e.g
        # a phenyl or nitro substituent flipped 180 degrees.
        rms_matrix = np.zeros((n_conformers, n_conformers))

        for i, j in itertools.combinations(np.arange(n_conformers), 2):

            rms_matrix[i, j] = AllChem.GetBestRMS(
                rdkit_molecule,
                rdkit_molecule,
                conformer_ids[i],
                conformer_ids[j],
            )

        rms_matrix += rms_matrix.T
        return rms_matrix

    @classmethod
    def _elf_select_diverse_conformers(
        cls,
        molecule: "Molecule",
        ranked_conformers: List[unit.Quantity],
        limit: int,
        rms_tolerance: unit.Quantity,
    ) -> List[unit.Quantity]:
        """Attempt to greedily select a specified number conformers which are maximally
        diverse.

        The conformer with the lowest electrostatic energy (the first conformer in the
        ``ranked_conformers`` list) is always chosen. After that selection proceeds by:

        a) selecting an un-selected conformer which is the most different from those
          already selected, and whose RMS compared to each selected conformer is
          greater than ``rms_tolerance``. Here most different means the conformer
          which has the largest sum of RMS with the selected conformers.

        b) repeating a) until either ``limit`` number of conformers have been selected,
           or there are no more distinct conformers to select from.

        Notes
        -----

        * As the selection is greedy there is no guarantee that the selected conformers
          will be the optimal distinct i.e. there may be other selections of conformers
          which are more distinct.

        Parameters
        ----------
        molecule
            The molecule object which matches the conformers to select from.
        ranked_conformers
            A list of conformers to select from, ranked by their electrostatic
            interaction energy (see ``_compute_electrostatic_energy``).
        limit
            The maximum number of conformers to select.
        rms_tolerance
            Conformers whose RMS is within this amount will be treated as identical and
            the duplicate discarded.

        Returns
        -------
            The select list of conformers.
        """

        # Compute the RMS between all pairs of conformers
        molecule = copy.deepcopy(molecule)
        molecule.conformers.clear()

        for conformer in ranked_conformers:
            molecule.add_conformer(conformer)

        rms_matrix = cls._elf_compute_rms_matrix(molecule)

        # Apply the greedy selection process.
        closed_list = np.zeros(limit).astype(int)
        closed_mask = np.zeros(rms_matrix.shape[0], dtype=bool)

        n_selected = 1

        for i in range(min(molecule.n_conformers, limit - 1)):

            distances = rms_matrix[closed_list[: i + 1], :].sum(axis=0)

            # Exclude already selected conformers or conformers which are too similar
            # to those already selected.
            closed_mask[
                np.any(
                    rms_matrix[closed_list[: i + 1], :]
                    < rms_tolerance.value_in_unit(unit.angstrom),
                    axis=0,
                )
            ] = True

            if np.all(closed_mask):
                # Stop of there are no more distinct conformers to select from.
                break

            distant_index = np.ma.array(distances, mask=closed_mask).argmax()
            closed_list[i + 1] = distant_index

            n_selected += 1

        return [ranked_conformers[i.item()] for i in closed_list[:n_selected]]

    def apply_elf_conformer_selection(
        self,
        molecule: "Molecule",
        percentage: float = 2.0,
        limit: int = 10,
        rms_tolerance: unit.Quantity = 0.05 * unit.angstrom,
    ):
        """Applies the `ELF method
        <https://docs.eyesopen.com/toolkits/python/quacpactk/molchargetheory.html#elf-conformer-selection>`_
        to select a set of diverse conformers which have minimal electrostatically
        strongly interacting functional groups from a molecules conformers.

        The diverse conformer selection is performed by the ``_elf_select_diverse_conformers``
        function, which attempts to greedily select conformers which are most distinct
        according to their RMS.

        Warnings
        --------
        * Although this function is inspired by the OpenEye ELF10 method, this
          implementation may yield slightly different conformers due to potential
          differences in this and the OE closed source implementation.

        Notes
        -----
        * The input molecule should have a large set of conformers already
          generated to select the ELF10 conformers from.
        * The selected conformers will be retained in the `molecule.conformers` list
          while unselected conformers will be discarded.
        * Only heavy atoms are included when using the RMS to select diverse conformers.

        See Also
        --------
        RDKitToolkitWrapper._elf_select_diverse_conformers

        Parameters
        ----------
        molecule
            The molecule which contains the set of conformers to select from.
        percentage
            The percentage of conformers with the lowest electrostatic interaction
            energies to greedily select from.
        limit
            The maximum number of conformers to select.
        rms_tolerance
            Conformers whose RMS is within this amount will be treated as identical and
            the duplicate discarded.
        """

        if molecule.n_conformers == 0:
            return

        # Copy the input molecule so we can directly perturb it within the method.
        molecule_copy = copy.deepcopy(molecule)

        # Prune any problematic conformers, such as trans-COOH configurations.
        conformers = self._elf_prune_problematic_conformers(molecule_copy)

        if len(conformers) == 0:

            raise ValueError(
                "There were no conformers to select from after discarding conformers "
                "which are known to be problematic when computing ELF partial charges. "
                "Make sure to generate a diverse array of conformers before calling the "
                "`RDKitToolkitWrapper.apply_elf_conformer_selection` method."
            )

        # Generate a set of absolute MMFF94 partial charges for the molecule and use
        # these to compute the electrostatic interaction energy of each conformer.
        self.assign_partial_charges(molecule_copy, "mmff94")

        conformer_energies = [
            (
                self._elf_compute_electrostatic_energy(molecule_copy, conformer),
                conformer,
            )
            for conformer in conformers
        ]

        # Rank the conformer energies and retain `percentage`% with the lowest energies.
        conformer_energies = sorted(conformer_energies, key=lambda x: x[0])
        cutoff_index = max(1, int(len(conformer_energies) * percentage / 100.0))

        low_energy_conformers = [
            conformer for _, conformer in conformer_energies[:cutoff_index]
        ]

        # Attempt to greedily select `limit` conformers which are maximally diverse.
        diverse_conformers = self._elf_select_diverse_conformers(
            molecule_copy, low_energy_conformers, limit, rms_tolerance
        )

        molecule._conformers = diverse_conformers

    def from_rdkit(
        self,
        rdmol,
        allow_undefined_stereo=False,
        hydrogens_are_explicit=False,
        _cls=None,
    ):
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
        hydrogens_are_explicit : bool, default=False
            If False, RDKit will perform hydrogen addition using Chem.AddHs
        _cls : class
            Molecule constructor

        Returns
        -------
        molecule : openff.toolkit.topology.Molecule
            An OpenFF molecule

        Examples
        --------

        Create a molecule from an RDKit molecule

        >>> from rdkit import Chem
        >>> from openff.toolkit.tests.utils import get_data_file_path
        >>> rdmol = Chem.MolFromMolFile(get_data_file_path('systems/monomers/ethanol.sdf'))

        >>> toolkit_wrapper = RDKitToolkitWrapper()
        >>> molecule = toolkit_wrapper.from_rdkit(rdmol)

        """
        from rdkit import Chem

        if _cls is None:
            from openff.toolkit.topology.molecule import Molecule

            _cls = Molecule

        # Make a copy of the RDKit Mol as we'll need to change it (e.g. assign stereo).
        rdmol = Chem.Mol(rdmol)

        if not hydrogens_are_explicit:
            rdmol = Chem.AddHs(rdmol, addCoords=True)

        # Sanitizing the molecule. We handle aromaticity and chirality manually.
        # This SanitizeMol(...) calls cleanUp, updatePropertyCache, symmetrizeSSSR,
        # assignRadicals, setConjugation, and setHybridization.
        Chem.SanitizeMol(
            rdmol,
            (
                Chem.SANITIZE_ALL
                ^ Chem.SANITIZE_SETAROMATICITY
                ^ Chem.SANITIZE_ADJUSTHS
                ^ Chem.SANITIZE_CLEANUPCHIRALITY
                ^ Chem.SANITIZE_KEKULIZE
            ),
        )
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
        self._detect_undefined_stereo(
            rdmol,
            raise_warning=allow_undefined_stereo,
            err_msg_prefix="Unable to make OFFMol from RDMol: ",
        )

        # Create a new OpenFF Molecule
        offmol = _cls()

        # If RDMol has a title save it
        if rdmol.HasProp("_Name"):
            # raise Exception('{}'.format(rdmol.GetProp('name')))
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
                map_id = int(rda.GetProp("_map_idx"))
            except KeyError:
                map_id = rda.GetAtomMapNum()

            # create a new atom
            # atomic_number = oemol.NewAtom(rda.GetAtomicNum())
            atomic_number = rda.GetAtomicNum()
            formal_charge = rda.GetFormalCharge() * unit.elementary_charge
            is_aromatic = rda.GetIsAromatic()
            if rda.HasProp("_Name"):
                name = rda.GetProp("_Name")
            else:
                # check for PDB names
                try:
                    name = rda.GetMonomerInfo().GetName().strip()
                except AttributeError:
                    name = ""

            # If chiral, store the chirality to be set later
            stereochemistry = None
            # tag = rda.GetChiralTag()
            if rda.HasProp("_CIPCode"):
                stereo_code = rda.GetProp("_CIPCode")
                # if tag == Chem.CHI_TETRAHEDRAL_CCW:
                if stereo_code == "R":
                    stereochemistry = "R"
                # if tag == Chem.CHI_TETRAHEDRAL_CW:
                elif stereo_code == "S":
                    stereochemistry = "S"
                else:
                    raise UndefinedStereochemistryError(
                        "In from_rdkit: Expected atom stereochemistry of R or S. "
                        "Got {} instead.".format(stereo_code)
                    )

            atom_index = offmol._add_atom(
                atomic_number,
                formal_charge,
                is_aromatic,
                name=name,
                stereochemistry=stereochemistry,
            )
            map_atoms[rd_idx] = atom_index
            atom_mapping[atom_index] = map_id

        # If we have a full / partial atom map add it to the molecule. Zeroes 0
        # indicates no mapping
        if {*atom_mapping.values()} != {0}:

            offmol._properties["atom_map"] = {
                idx: map_idx for idx, map_idx in atom_mapping.items() if map_idx != 0
            }

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
            bond_index = offmol._add_bond(
                map_atoms[a1], map_atoms[a2], order, is_aromatic
            )
            map_bonds[rdb_idx] = bond_index

        # Now fill in the cached (structure-dependent) properties. We have to have the 2D
        # structure of the molecule in place first, because each call to add_atom and
        # add_bond invalidates all cached properties
        for rdb in rdmol.GetBonds():
            rdb_idx = rdb.GetIdx()
            offb_idx = map_bonds[rdb_idx]
            offb = offmol.bonds[offb_idx]
            # determine if stereochemistry is needed
            # Note that RDKit has 6 possible values of bond stereo: CIS, TRANS, E, Z, ANY, or NONE
            # The logic below assumes that "ANY" and "NONE" mean the same thing.
            stereochemistry = None
            tag = rdb.GetStereo()
            if tag == Chem.BondStereo.STEREOZ:
                stereochemistry = "Z"
            elif tag == Chem.BondStereo.STEREOE:
                stereochemistry = "E"
            elif tag == Chem.BondStereo.STEREOTRANS or tag == Chem.BondStereo.STEREOCIS:
                raise ValueError(
                    "Expected RDKit bond stereochemistry of E or Z, got {} instead".format(
                        tag
                    )
                )
            offb._stereochemistry = stereochemistry
            fractional_bond_order = None
            if rdb.HasProp("fractional_bond_order"):
                fractional_bond_order = rdb.GetDoubleProp("fractional_bond_order")
            offb.fractional_bond_order = fractional_bond_order

        # TODO: Save conformer(s), if present
        # If the rdmol has a conformer, store its coordinates
        if len(rdmol.GetConformers()) != 0:
            for conf in rdmol.GetConformers():
                n_atoms = offmol.n_atoms
                # TODO: Will this always be angstrom when loading from RDKit?
                positions = unit.Quantity(np.zeros((n_atoms, 3)), unit.angstrom)
                for rd_idx, off_idx in map_atoms.items():
                    atom_coords = conf.GetPositions()[rd_idx, :] * unit.angstrom
                    positions[off_idx, :] = atom_coords
                offmol._add_conformer(positions)

        partial_charges = unit.Quantity(
            np.zeros(shape=offmol.n_atoms, dtype=np.float64),
            unit=unit.elementary_charge,
        )

        any_atom_has_partial_charge = False
        for rd_idx, rd_atom in enumerate(rdmol.GetAtoms()):
            off_idx = map_atoms[rd_idx]
            if rd_atom.HasProp("PartialCharge"):
                charge = rd_atom.GetDoubleProp("PartialCharge") * unit.elementary_charge
                partial_charges[off_idx] = charge
                any_atom_has_partial_charge = True
            else:
                # If some other atoms had partial charges but this one doesn't, raise an Exception
                if any_atom_has_partial_charge:
                    raise ValueError(
                        "Some atoms in rdmol have partial charges, but others do not."
                    )
        if any_atom_has_partial_charge:
            offmol.partial_charges = partial_charges
        else:
            offmol.partial_charges = None
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

        >>> from openff.toolkit.topology import Molecule
        >>> ethanol = Molecule.from_smiles('CCO')
        >>> rdmol = ethanol.to_rdkit()

        """
        from rdkit import Chem, Geometry

        # Create an editable RDKit molecule
        rdmol = Chem.RWMol()

        # Set name
        # TODO: What is the best practice for how this should be named?
        if not (molecule.name is None):
            rdmol.SetProp("_Name", molecule.name)

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
            rdatom.SetFormalCharge(
                atom.formal_charge.value_in_unit(unit.elementary_charge)
            )
            rdatom.SetIsAromatic(atom.is_aromatic)
            rdatom.SetProp("_Name", atom.name)

            # Stereo handling code moved to after bonds are added
            if atom.stereochemistry == "S":
                rdatom.SetChiralTag(Chem.CHI_TETRAHEDRAL_CW)
            elif atom.stereochemistry == "R":
                rdatom.SetChiralTag(Chem.CHI_TETRAHEDRAL_CCW)

            # Stop rdkit from adding implicit hydrogens
            rdatom.SetNoImplicit(True)

            rd_index = rdmol.AddAtom(rdatom)

            # Let's make sure al the atom indices in the two molecules
            # are the same, otherwise we need to create an atom map.
            assert index == atom.molecule_atom_index
            assert index == rd_index

        for bond in molecule.bonds:
            atom_indices = (
                bond.atom1.molecule_atom_index,
                bond.atom2.molecule_atom_index,
            )
            rdmol.AddBond(*atom_indices)
            rdbond = rdmol.GetBondBetweenAtoms(*atom_indices)
            if not (bond.fractional_bond_order is None):
                rdbond.SetDoubleProp(
                    "fractional_bond_order", bond.fractional_bond_order
                )
            # Assign bond type, which is based on order unless it is aromatic
            if bond.is_aromatic:
                rdbond.SetBondType(_bondtypes[1.5])
                rdbond.SetIsAromatic(True)
            else:
                rdbond.SetBondType(_bondtypes[bond.bond_order])
                rdbond.SetIsAromatic(False)

        Chem.SanitizeMol(
            rdmol,
            Chem.SANITIZE_ALL ^ Chem.SANITIZE_ADJUSTHS ^ Chem.SANITIZE_SETAROMATICITY,
        )

        # Fix for aromaticity being lost
        if aromaticity_model == "OEAroModel_MDL":
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
            if (
                rdatom.HasProp("_CIPCode")
                and rdatom.GetProp("_CIPCode") == atom.stereochemistry
            ):
                continue

            # Otherwise, set it to CCW.
            rdatom.SetChiralTag(Chem.CHI_TETRAHEDRAL_CCW)
            # We need to do force and cleanIt to recalculate CIP stereo.
            Chem.AssignStereochemistry(rdmol, force=True, cleanIt=True)
            # Hopefully this worked, otherwise something's wrong
            if (
                rdatom.HasProp("_CIPCode")
                and rdatom.GetProp("_CIPCode") == atom.stereochemistry
            ):
                continue

            # Keep track of undefined stereo atoms. We'll force stereochemistry
            # at the end to avoid the next AssignStereochemistry to overwrite.
            if not rdatom.HasProp("_CIPCode"):
                undefined_stereo_atoms[rdatom] = atom.stereochemistry
                continue

            # Something is wrong.
            err_msg = (
                "Unknown atom stereochemistry encountered in to_rdkit. "
                "Desired stereochemistry: {}. Set stereochemistry {}".format(
                    atom.stereochemistry, rdatom.GetProp("_CIPCode")
                )
            )
            raise RuntimeError(err_msg)

        # Copy bond stereo info from molecule to rdmol.
        cls._assign_rdmol_bonds_stereo(molecule, rdmol)

        # Set coordinates if we have them
        if molecule._conformers:
            for conformer in molecule._conformers:
                rdmol_conformer = Chem.Conformer()
                for atom_idx in range(molecule.n_atoms):
                    x, y, z = conformer[atom_idx, :].value_in_unit(unit.angstrom)
                    rdmol_conformer.SetAtomPosition(atom_idx, Geometry.Point3D(x, y, z))
                rdmol.AddConformer(rdmol_conformer, assignId=True)

        # Retain charges, if present
        if not (molecule._partial_charges is None):

            rdk_indexed_charges = np.zeros(shape=molecule.n_atoms, dtype=float)
            for atom_idx, charge in enumerate(molecule._partial_charges):
                charge_unitless = charge.value_in_unit(unit.elementary_charge)
                rdk_indexed_charges[atom_idx] = charge_unitless
            for atom_idx, rdk_atom in enumerate(rdmol.GetAtoms()):
                rdk_atom.SetDoubleProp("PartialCharge", rdk_indexed_charges[atom_idx])

            # Note: We could put this outside the "if" statement, which would result in all
            #     partial charges in the resulting file being set to "n/a" if they weren't
            #     set in the Open Force Field Toolkit ``Molecule``
            Chem.CreateAtomDoublePropertyList(rdmol, "PartialCharge")

        # Cleanup the rdmol
        rdmol.UpdatePropertyCache(strict=False)
        Chem.GetSSSR(rdmol)

        # Forcefully assign stereo information on the atoms that RDKit
        # can't figure out. This must be done last as calling AssignStereochemistry
        # again will delete these properties (see #196).
        for rdatom, stereochemistry in undefined_stereo_atoms.items():
            rdatom.SetProp("_CIPCode", stereochemistry)

        # Return non-editable version
        return Chem.Mol(rdmol)

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
            If a fixed hydrogen layer should be added to the InChI, if `True` this will produce a
            non standard specific InChI string of the molecule.

        Returns
        --------
        inchi: str
            The InChI string of the molecule.
        """

        from rdkit import Chem

        rdmol = self.to_rdkit(molecule)
        if fixed_hydrogens:
            inchi = Chem.MolToInchi(rdmol, options="-FixedH")
        else:
            inchi = Chem.MolToInchi(rdmol)
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
            If a fixed hydrogen layer should be added to the InChI, if `True` this will
            produce a non standard specific InChI string of the molecule.

        Returns
        --------
        inchi_key: str
            The InChIKey representation of the molecule.
        """

        from rdkit import Chem

        rdmol = self.to_rdkit(molecule)
        if fixed_hydrogens:
            inchi_key = Chem.MolToInchiKey(rdmol, options="-FixedH")
        else:
            inchi_key = Chem.MolToInchiKey(rdmol)
        return inchi_key

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
            A tuple of tuples, where each inner tuple is a pair of tagged atoms (tag_idx_1, tag_idx_2)
             which are bonded. The inner tuples are ordered smallest-to-largest, and the tuple of
             tuples is ordered lexically. So the return value for an improper torsion would be
             ((1, 2), (2, 3), (2, 4)).

        Raises
        ------
        SMIRKSParsingError
            If RDKit was unable to parse the provided smirks/tagged smarts
        """
        from rdkit import Chem

        from openff.toolkit.typing.chemistry import SMIRKSParsingError

        ss = Chem.MolFromSmarts(smarts)

        if ss is None:
            raise SMIRKSParsingError(f"RDKit was unable to parse SMIRKS {smarts}")

        unique_tags = set()
        connections = set()
        for at1 in ss.GetAtoms():
            if at1.GetAtomMapNum() == 0:
                continue
            unique_tags.add(at1.GetAtomMapNum())
            for at2 in at1.GetNeighbors():
                if at2.GetAtomMapNum() == 0:
                    continue
                cxn_to_add = sorted([at1.GetAtomMapNum(), at2.GetAtomMapNum()])
                connections.add(tuple(cxn_to_add))
        connections = tuple(sorted(list(connections)))
        unique_tags = tuple(sorted(list(unique_tags)))
        return unique_tags, connections

    @staticmethod
    def _find_smarts_matches(
        rdmol,
        smirks,
        aromaticity_model="OEAroModel_MDL",
        unique=False,
    ):
        """Find all sets of atoms in the provided RDKit molecule that match the provided SMARTS string.

        Parameters
        ----------
        rdmol : rdkit.Chem.Mol
            rdmol to process with the SMIRKS in order to find matches
        smarts : str
            SMARTS string with any number of sequentially tagged atoms.
            If there are N tagged atoms numbered 1..N, the resulting matches will be
            N-tuples of atoms that match the corresponding tagged atoms.
        aromaticity_model : str, optional, default='OEAroModel_MDL'
            OpenEye aromaticity model designation as a string, such as ``OEAroModel_MDL``.
            Molecule is prepared with this aromaticity model prior to querying.

        Returns
        -------
        matches : list of tuples of atoms indices within the ``rdmol``
            matches[index] is an N-tuple of atom numbers from the ``rdmol``
            Matches are returned in no guaranteed order.
            # TODO: What is returned if no matches are found? An empty list, or None?
            # TODO: Ensure that SMARTS numbers 1, 2, 3... are rendered into order of
            #    returned matches indexed by 0, 1, 2...

        .. notes ::

           * Raises ``ValueError`` if ``smarts`` query is malformed

        """
        from rdkit import Chem

        # This code is part of a possible performance optimization that hasn't been validated
        # for production use yet.
        def _match_smarts_with_heavy_atoms_first(rdmol, qmol, match_kwargs):
            for i, atom in enumerate(qmol.GetAtoms()):
                atom.SetIntProp("index", i)

            remove_params = Chem.rdmolops.RemoveHsParameters()
            remove_params.removeWithQuery = True
            heavy_query = Chem.RemoveHs(qmol, remove_params, sanitize=False)
            heavy_to_qmol = [
                atom.GetIntProp("index") for atom in heavy_query.GetAtoms()
            ]
            query_atoms = [Chem.Atom(i + 2) for i in range(len(heavy_to_qmol))]

            full_matches = set()

            for heavy_match in rdmol.GetSubstructMatches(heavy_query, **match_kwargs):
                rdmol_copy = Chem.RWMol(rdmol)
                qmol_copy = Chem.RWMol(qmol)
                # pin atoms by atom type
                for heavy_index, rdmol_index in enumerate(heavy_match):
                    qmol_index = heavy_to_qmol[heavy_index]
                    qmol_copy.ReplaceAtom(qmol_index, query_atoms[heavy_index])
                    rdmol_copy.ReplaceAtom(rdmol_index, query_atoms[heavy_index])

                rdmol_copy.UpdatePropertyCache(strict=False)
                qmol_copy.UpdatePropertyCache(strict=False)
                h_matches = rdmol_copy.GetSubstructMatches(qmol_copy, **match_kwargs)
                full_matches |= set(h_matches)

            return full_matches

        # Make a copy of the molecule
        rdmol = Chem.Mol(rdmol)
        # Use designated aromaticity model
        if aromaticity_model == "OEAroModel_MDL":
            Chem.SanitizeMol(rdmol, Chem.SANITIZE_ALL ^ Chem.SANITIZE_SETAROMATICITY)
            Chem.SetAromaticity(rdmol, Chem.AromaticityModel.AROMATICITY_MDL)
        else:
            # Only the OEAroModel_MDL is supported for now
            raise ValueError("Unknown aromaticity model: {}".aromaticity_models)

        # Set up query.
        qmol = Chem.MolFromSmarts(smirks)  # cannot catch the error
        if qmol is None:
            raise ValueError(
                'RDKit could not parse the SMIRKS string "{}"'.format(smirks)
            )

        # Create atom mapping for query molecule
        idx_map = dict()
        for atom in qmol.GetAtoms():
            smirks_index = atom.GetAtomMapNum()
            if smirks_index != 0:
                idx_map[smirks_index - 1] = atom.GetIdx()
        map_list = [idx_map[x] for x in sorted(idx_map)]

        # choose the largest unsigned int without overflow
        # since the C++ signature is a uint
        # TODO: max_matches = int(max_matches) if max_matches is not None else np.iinfo(np.uintc).max
        max_matches = np.iinfo(np.uintc).max
        match_kwargs = dict(uniquify=unique, maxMatches=max_matches, useChirality=True)
        # These variables are un-used, do they serve a purpose?
        # n_heavy = qmol.GetNumHeavyAtoms()
        # n_h = qmol.GetNumAtoms() - n_heavy
        # TODO: if match_heavy_first: full_matches = _match_smarts_with_heavy_atoms_first(...)
        full_matches = rdmol.GetSubstructMatches(qmol, **match_kwargs)

        matches = [tuple(match[x] for x in map_list) for match in full_matches]

        return matches

    def find_smarts_matches(
        self,
        molecule,
        smarts,
        aromaticity_model="OEAroModel_MDL",
        unique=False,
    ):
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
        rdmol = self.to_rdkit(molecule, aromaticity_model=aromaticity_model)
        return self._find_smarts_matches(
            rdmol,
            smarts,
            aromaticity_model="OEAroModel_MDL",
            unique=unique,
        )

    # --------------------------------
    # Stereochemistry RDKit utilities.
    # --------------------------------

    def find_rings(self, molecule) -> Tuple[Tuple[int]]:
        """Find the rings in a given molecule.

        .. note ::

            For systems containing some special cases of connected rings, this
            function may not be well-behaved and may report a different number
            rings than expected. Some problematic cases include networks of many
            (5+) rings or bicyclic moieties (i.e. norbornane).

        Parameters
        ----------
        molecule : openff.toolkit.topology.Molecule
            The molecule for which rings are to be found

        Returns
        -------
        rings : tuple of tuples of atom indices
            Nested tuples, each containing the indices of atoms in each ring

        """
        rdmol = molecule.to_rdkit()
        ring_info = rdmol.GetRingInfo()
        rings = ring_info.AtomRings()
        return rings

    def get_n_rings(self, molecule) -> int:
        """Return the number of rings in this molecule. See `find_rings` for more details."""
        return len(self.find_rings(molecule))

    def atom_is_in_ring(self, molecule: "Molecule", atom_index: int) -> bool:
        """Return whether or not an atom is in a ring.

        Parameters
        ----------
        molecule : openff.toolkit.topology.Molecule
            The molecule containing the atom of interest
        atom_index : int
            The index of the atom of interest

        Returns
        -------
        is_in_ring : bool
            Whether or not the atom of index `atom_index` is in a ring

        """
        rdmol = molecule.to_rdkit()
        is_in_ring = rdmol.GetAtomWithIdx(atom_index).IsInRing()

        return is_in_ring

    def bond_is_in_ring(self, molecule: "Molecule", bond_index: int) -> bool:
        """Return whether or not a bond is in a ring.

        Parameters
        ----------
        molecule : openff.toolkit.topology.Molecule
            The molecule containing the bond of interest
        bond_index : int
            The index of the bond of interest

        Returns
        -------
        is_in_ring : bool
            Whether or not the bond of index `bond_index` is in a ring

        """
        rdmol = molecule.to_rdkit()
        rdbond = rdmol.GetBondWithIdx(bond_index)
        is_in_ring = rdbond.IsInRing()

        return is_in_ring

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
            if atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED and atom.HasProp(
                "_ChiralityPossible"
            ):
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
        rdmol_copy = copy.deepcopy(rdmol)

        # Clear any previous assignments on the bonds, since FindPotentialStereo may not overwrite it
        for bond in rdmol_copy.GetBonds():
            bond.SetStereo(Chem.BondStereo.STEREONONE)

        # This function assigns Bond.GetStereo() == Bond.STEREOANY to bonds with
        # possible stereochemistry.
        Chem.FindPotentialStereoBonds(rdmol_copy, cleanIt=True)

        # Any TRULY stereogenic bonds in the molecule are now marked as STEREOANY in rdmol_copy.
        # Iterate through all the bonds, and for the ones where rdmol_copy is marked as STEREOANY,
        # ensure that they are cis/trans/E/Z (tested here be ensuring that they're NOT either
        # # of the other possible types (NONE or ANY))
        undefined_bond_indices = []
        for bond_idx, (orig_bond, repercieved_bond) in enumerate(
            zip(rdmol.GetBonds(), rdmol_copy.GetBonds())
        ):
            # print(repercieved_bond.GetStereo(), orig_bond.GetStereo())
            if (repercieved_bond.GetStereo() == Chem.BondStereo.STEREOANY) and (
                (orig_bond.GetStereo() == Chem.BondStereo.STEREOANY)
                or (orig_bond.GetStereo() == Chem.BondStereo.STEREONONE)
            ):
                undefined_bond_indices.append(bond_idx)
        return undefined_bond_indices

    @classmethod
    def _detect_undefined_stereo(cls, rdmol, err_msg_prefix="", raise_warning=False):
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
                msg += " - Atom {symbol} (index {index})\n".format(
                    symbol=rdmol.GetAtomWithIdx(undefined_atom_idx).GetSymbol(),
                    index=undefined_atom_idx,
                )

        # Details about undefined bond.
        if len(undefined_bond_indices) > 0:
            msg += "Bonds with undefined stereochemistry are:\n"
            for undefined_bond_idx in undefined_bond_indices:
                bond = rdmol.GetBondWithIdx(undefined_bond_idx)
                atom1, atom2 = bond.GetBeginAtom(), bond.GetEndAtom()
                msg += " - Bond {bindex} (atoms {aindex1}-{aindex2} of element ({symbol1}-{symbol2})\n".format(
                    bindex=undefined_bond_idx,
                    aindex1=atom1.GetIdx(),
                    aindex2=atom2.GetIdx(),
                    symbol1=atom1.GetSymbol(),
                    symbol2=atom2.GetSymbol(),
                )

        if msg is not None:
            if raise_warning:
                msg = "Warning (not error because allow_undefined_stereo=True): " + msg
                logger.warning(msg)
            else:
                msg = "Unable to make OFFMol from RDMol: " + msg
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
                raise RuntimeError("Cannot flip the bond direction consistently.")

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
                    if (
                        paired_rdbond.GetBeginAtomIdx(),
                        paired_rdbond.GetEndAtomIdx(),
                    ) != ignored:
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
            rdbond_atom_indices = (
                bond.atom1.molecule_atom_index,
                bond.atom2.molecule_atom_index,
            )
            stereo_rdbond = rdmol.GetBondBetweenAtoms(*rdbond_atom_indices)

            # Collect all neighboring rdbonds of atom1 and atom2.
            neighbor_rdbonds1 = [
                rdmol.GetBondBetweenAtoms(
                    n.molecule_atom_index, bond.atom1.molecule_atom_index
                )
                for n in bond.atom1.bonded_atoms
                if n != bond.atom2
            ]
            neighbor_rdbonds2 = [
                rdmol.GetBondBetweenAtoms(
                    bond.atom2.molecule_atom_index, n.molecule_atom_index
                )
                for n in bond.atom2.bonded_atoms
                if n != bond.atom1
            ]

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
            assert bond.stereochemistry in {"E", "Z"}
            if bond.stereochemistry == "E":
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
            neighbor_bond_indices = [
                (rdb.GetBeginAtomIdx(), rdb.GetEndAtomIdx()) for rdb in neighbor_rdbonds
            ]
            for i, bond_indices in enumerate(neighbor_bond_indices):
                try:
                    paired_bonds[bond_indices].append(neighbor_rdbonds[1 - i])
                except KeyError:
                    paired_bonds[bond_indices] = [neighbor_rdbonds[1 - i]]
