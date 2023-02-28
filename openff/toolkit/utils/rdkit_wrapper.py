"""
Wrapper class providing a minimal consistent interface to the `RDKit <http://www.rdkit.org/>`.
"""

__all__ = ("RDKitToolkitWrapper",)

import copy
import functools
import importlib
import itertools
import logging
import pathlib
import tempfile
from collections import defaultdict
from typing import TYPE_CHECKING, Dict, List, Optional, Tuple

import numpy as np
from cachetools import LRUCache, cached
from openff.units import Quantity, unit
from openff.units.elements import SYMBOLS

from openff.toolkit.utils import base_wrapper
from openff.toolkit.utils.constants import (
    ALLOWED_AROMATICITY_MODELS,
    DEFAULT_AROMATICITY_MODEL,
)
from openff.toolkit.utils.exceptions import (
    ChargeMethodUnavailableError,
    ConformerGenerationError,
    InvalidAromaticityModelError,
    NotAttachedToMoleculeError,
    RadicalsNotSupportedError,
    SMILESParseError,
    ToolkitUnavailableException,
    UnassignedChemistryInPDBError,
    UndefinedStereochemistryError,
)

if TYPE_CHECKING:
    from openff.toolkit.topology.molecule import Atom, Bond, Molecule

logger = logging.getLogger(__name__)


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
        "toolkit can be found at: https://anaconda.org/conda-forge/rdkit"
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
    def toolkit_file_write_formats(self) -> List[str]:
        """
        List of file formats that this toolkit can write.
        """
        return list(self._toolkit_file_write_formats.keys())

    @classmethod
    def is_available(cls) -> bool:
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

    def from_object(self, obj, allow_undefined_stereo: bool = False, _cls=None):
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
        self,
        file_path: str,
        smiles: str,
        allow_undefined_stereo: bool = False,
        _cls=None,
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
        # Dict[offmol_index, pdbmol_index] sorted by offmol index
        isomorphic, mapping = _cls.are_isomorphic(
            offmol,
            pdbmol,
            return_atom_map=True,
            aromatic_matching=False,
            formal_charge_matching=False,
            bond_order_matching=False,
            atom_stereochemistry_matching=False,
            bond_stereochemistry_matching=False,
        )

        if mapping is None:
            from openff.toolkit.topology.molecule import InvalidConformerError

            raise InvalidConformerError("The PDB and SMILES structures do not match.")

        new_mol = offmol.remap(mapping)

        # the pdb conformer is in the correct order so just attach it here
        new_mol._add_conformer(pdbmol.conformers[0])

        # Take residue info from PDB
        for pdbatom, newatom in zip(pdbmol.atoms, new_mol.atoms):
            newatom.metadata.update(pdbatom.metadata)
            newatom.name = pdbatom.name
        new_mol.add_default_hierarchy_schemes()

        return new_mol

    def _polymer_openmm_topology_to_offmol(self, omm_top, substructure_dictionary):
        rdkit_mol = self._polymer_openmm_topology_to_rdmol(
            omm_top, substructure_dictionary
        )
        offmol = self.from_rdkit(rdkit_mol, allow_undefined_stereo=True)
        return offmol

    def _polymer_openmm_topology_to_rdmol(
        self,
        omm_top,
        substructure_library,
    ):
        """
        Parameters
        ----------
        rdkit_mol : rdkit.Chem.Mol
            Currently invalid (bond orders and charge) Molecule
        substructure_library : dict{str:list[str, list[str]]}
            A dictionary of substructures. substructure_library[aa_name] = list[tagged SMARTS, list[atom_names]]
        toolkit_registry = ToolkitWrapper or ToolkitRegistry. Default = None
            Either a ToolkitRegistry, ToolkitWrapper

        Returns
        -------
        rdkit_mol : rdkit.Chem.Mol
            a copy of the original molecule with charges and bond order added

        Raises
        ------
        MissingChemistryFromPolymerError
            Raised when bonds or atoms in ``rdkit_mol`` are missing from the
            substructure library
        """
        from rdkit import Chem

        already_assigned_nodes = set()
        # TODO: We currently assume all single and modify a few
        # Therefore it's hard to know if we've missed any edges...
        # Notably assumes all bonds *between* fragments are single
        already_assigned_edges = set()

        # Keeping track of which atoms are matched where will help us with error
        # messages
        matches = defaultdict(list)

        rdkit_mol = self._get_connectivity_from_openmm_top(omm_top)
        mol = Chem.Mol(rdkit_mol)

        for res_name in substructure_library:
            # TODO: This is a hack for the moment since we don't have a more sophisticated way to resolve clashes
            # so it just does the biggest substructures first
            # NOTE: If this changes, MissingChemistryFromPolymerError needs to be updated too
            sorted_substructure_smarts = sorted(
                substructure_library[res_name], key=len, reverse=True
            )
            for substructure_smarts in sorted_substructure_smarts:
                # this is the molecule as defined in template
                ref = Chem.MolFromSmarts(substructure_smarts)
                # then create a looser definition for pattern matching...
                # be lax about double bonds and chirality
                fuzzy = self._fuzzy_query(ref)

                # It's important that we do the substructure search on `rdkit_mol`, but the chemical
                # info is added to `mol`. If we use the same rdkit molecule for search AND info addition,
                # then single bonds may no longer be present for subsequent overlapping matches.
                for match in rdkit_mol.GetSubstructMatches(fuzzy, maxMatches=0):
                    for i in match:
                        matches[i].append(res_name)

                    if any(m in already_assigned_nodes for m in match) and (
                        res_name not in ["PEPTIDE_BOND", "DISULFIDE"]
                    ):
                        continue
                    already_assigned_nodes.update(match)

                    for atom_i, j in zip(ref.GetAtoms(), match):
                        atom_j = mol.GetAtomWithIdx(j)
                        # copy over chirality
                        if atom_i.GetChiralTag():
                            mol.GetAtomWithIdx(j).SetChiralTag(atom_i.GetChiralTag())
                        atom_j.SetFormalCharge(atom_i.GetFormalCharge())

                    for b in ref.GetBonds():
                        x = match[b.GetBeginAtomIdx()]
                        y = match[b.GetEndAtomIdx()]
                        b2 = mol.GetBondBetweenAtoms(x, y)
                        b2.SetBondType(b.GetBondType())
                        already_assigned_edges.add(tuple(sorted([x, y])))

        unassigned_atoms = sorted(
            set(range(rdkit_mol.GetNumAtoms())) - already_assigned_nodes
        )
        all_bonds = set(
            [
                tuple(sorted([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()]))
                for bond in rdkit_mol.GetBonds()
            ]
        )
        unassigned_bonds = sorted(all_bonds - already_assigned_edges)

        if unassigned_atoms or unassigned_bonds:
            # Some advanced error reporting needs to interpret the substructure smarts to do things like
            # compare atom counts. Since OFFTK doesn't have a native class to hold fragments, we convert
            # the smarts into a sorted list of symbols to help with generating the error message.
            resname_to_symbols_and_atomnames = {}
            for resname, smarts_to_atom_names in substructure_library.items():
                resname_to_symbols_and_atomnames[resname] = list()
                for smarts, atom_names in smarts_to_atom_names.items():
                    ref = Chem.MolFromSmarts(smarts)
                    symbols = sorted(
                        [SYMBOLS[atom.GetAtomicNum()] for atom in ref.GetAtoms()]
                    )
                    resname_to_symbols_and_atomnames[resname].append(
                        (symbols, atom_names)
                    )

            raise UnassignedChemistryInPDBError(
                substructure_library=resname_to_symbols_and_atomnames,
                omm_top=omm_top,
                unassigned_atoms=unassigned_atoms,
                unassigned_bonds=unassigned_bonds,
                matches=matches,
            )

        return mol

    def _get_connectivity_from_openmm_top(self, omm_top):
        from rdkit import Chem

        # convert openmm topology to rdkit Molecule
        # all bonds initially SINGLE, all charge initially neutral
        rwmol = Chem.RWMol()
        for atom in omm_top.atoms():
            idx = rwmol.AddAtom(Chem.Atom(atom.element.atomic_number))
            res = Chem.AtomPDBResidueInfo()
            res.SetResidueName(atom.residue.name)
            res.SetResidueNumber(int(atom.residue.id))
            res.SetChainId(atom.residue.chain.id)
            rwatom = rwmol.GetAtomWithIdx(idx)
            rwatom.SetPDBResidueInfo(res)
        # we're fully explicit
        for atom in rwmol.GetAtoms():
            atom.SetNoImplicit(True)
        for bond in omm_top.bonds():
            rwmol.AddBond(bond[0].index, bond[1].index, Chem.BondType.SINGLE)

        # conf = Chem.Conformer()
        # for i, pos in enumerate(omm_top.getPositions()):
        #     conf.SetAtomPosition(i, list(pos.value_in_unit(pos.unit)))
        # rwmol.AddConformer(conf)

        return rwmol

    @staticmethod
    def _fuzzy_query(query):
        """return a copy of Query which is less specific:
        - ignore aromaticity and hybridization of atoms (i.e. [#6] not C)
        - ignore bond orders
        - ignore formal charges
        """
        from rdkit import Chem

        # it's tricky from the Python API to properly edit queries,
        # but you can do SetQuery on Atoms/Bonds to edit them quite powerfully
        generic = Chem.MolFromSmarts("**")
        generic_bond = generic.GetBondWithIdx(0)
        # N.B. This isn't likely to be an active
        generic_mol = (
            Chem.MolFromSmarts(  # TODO: optimisation, create this once somewhere
                "".join("[#{}]".format(i + 1) for i in range(112))
            )
        )

        fuzzy = Chem.Mol(query)
        for a in fuzzy.GetAtoms():
            a.SetFormalCharge(0)
            a.SetQuery(
                generic_mol.GetAtomWithIdx(a.GetAtomicNum() - 1)
            )  # i.e. H looks up atom 0 in our generic mol
            a.SetNoImplicit(True)
        for b in fuzzy.GetBonds():
            b.SetIsAromatic(False)
            b.SetBondType(Chem.rdchem.BondType.SINGLE)
            b.SetQuery(generic_bond)

        return fuzzy

    def _assign_aromaticity_and_stereo_from_3d(self, offmol):
        from rdkit import Chem

        rdmol = offmol.to_rdkit()
        Chem.SanitizeMol(
            rdmol,
            Chem.SANITIZE_ALL
            ^ Chem.SANITIZE_ADJUSTHS,  # ^ Chem.SANITIZE_SETAROMATICITY,
        )
        Chem.AssignStereochemistryFrom3D(rdmol)
        Chem.Kekulize(rdmol, clearAromaticFlags=True)
        Chem.SetAromaticity(rdmol, Chem.AromaticityModel.AROMATICITY_MDL)
        # To get HIS//TRP to get recognized as aromatic, we can use a different aromaticity model
        # Chem.SetAromaticity(rdmol, Chem.AromaticityModel.AROMATICITY_DEFAULT)

        offmol_w_stereo_and_aro = offmol.from_rdkit(
            rdmol, allow_undefined_stereo=True, hydrogens_are_explicit=True
        )
        return offmol_w_stereo_and_aro

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
        self,
        file_path: str,
        file_format: str,
        allow_undefined_stereo: bool = False,
        _cls=None,
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

        if isinstance(file_path, pathlib.Path):
            file_path: str = file_path.as_posix()

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
        self,
        file_obj,
        file_format: str,
        allow_undefined_stereo: bool = False,
        _cls=None,
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

    def to_file_obj(self, molecule: "Molecule", file_obj, file_format: str):
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

    def to_file(self, molecule: "Molecule", file_path: str, file_format: str):
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
        self,
        molecule: "Molecule",
        undefined_only: bool = False,
        max_isomers: int = 20,
        rationalise: bool = True,
    ) -> List["Molecule"]:
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

    def enumerate_tautomers(
        self, molecule: "Molecule", max_states: int = 20
    ) -> List["Molecule"]:
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

    def canonical_order_atoms(self, molecule: "Molecule") -> "Molecule":
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

    def to_smiles(
        self,
        molecule: "Molecule",
        isomeric: bool = True,
        explicit_hydrogens: bool = True,
        mapped: bool = False,
    ):
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
        smiles: str,
        hydrogens_are_explicit: bool = False,
        allow_undefined_stereo: bool = False,
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

        Raises
        ------
        RadicalsNotSupportedError : If any atoms in the RDKit molecule contain radical electrons.
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

    def from_inchi(self, inchi: str, allow_undefined_stereo: bool = False, _cls=None):
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
        self,
        molecule: "Molecule",
        n_conformers: int = 1,
        rms_cutoff: Optional[Quantity] = None,
        clear_existing: bool = True,
        _cls=None,
        make_carboxylic_acids_cis: bool = False,
    ):
        """
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
        rms_cutoff : unit-wrapped float, in units of distance, optional, default=None
            The minimum RMS value at which two conformers are considered redundant and one is deleted.
            If None, the cutoff is set to 1 Angstrom

        clear_existing : bool, default=True
            Whether to overwrite existing conformers for the molecule.
        _cls : class
            Molecule constructor
        make_carboxylic_acids_cis: bool, default=False
            Guarantee all conformers have exclusively cis carboxylic acid groups (COOH)
            by rotating the proton in any trans carboxylic acids 180 degrees around the C-O bond.

        """
        from rdkit.Chem import AllChem

        if rms_cutoff is None:
            rms_cutoff = unit.Quantity(1.0, unit.angstrom)
        rdmol = self.to_rdkit(molecule)
        # TODO: This generates way more conformations than omega, given the same
        # nConfs and RMS threshold. Is there some way to set an energy cutoff as well?
        first_conformer_generation_status = AllChem.EmbedMultipleConfs(
            rdmol,
            numConfs=n_conformers,
            pruneRmsThresh=rms_cutoff.m_as(unit.angstrom),
            randomSeed=1,
            # params=AllChem.ETKDG()
        )

        if not first_conformer_generation_status:
            # For some large molecules, conformer generation fails without `useRandomCoords`;
            # Landrum recommends it https://github.com/rdkit/rdkit/issues/3764#issuecomment-769367489
            fallback_conformer_generation_status = AllChem.EmbedMultipleConfs(
                rdmol,
                numConfs=n_conformers,
                pruneRmsThresh=rms_cutoff.m_as(unit.angstrom),
                randomSeed=1,
                useRandomCoords=True,
            )

            if not fallback_conformer_generation_status:
                raise ConformerGenerationError("RDKit conformer generation failed.")

        molecule2 = self.from_rdkit(
            rdmol, allow_undefined_stereo=True, _cls=molecule.__class__
        )

        if clear_existing:
            molecule._conformers = list()

        for conformer in molecule2._conformers:
            molecule._add_conformer(conformer)

        if make_carboxylic_acids_cis:
            molecule._make_carboxylic_acids_cis(toolkit_registry=self)

    def assign_partial_charges(
        self,
        molecule: "Molecule",
        partial_charge_method: Optional[str] = None,
        use_conformers: Optional[List[Quantity]] = None,
        strict_n_conformers: bool = False,
        normalize_partial_charges: bool = True,
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
            The charge model to use. One of ['mmff94', 'gasteiger']. If None, 'mmff94' will be used.

            * 'mmff94': Applies partial charges using the Merck Molecular Force Field
                        (MMFF). This method does not make use of conformers, and hence
                        ``use_conformers`` and ``strict_n_conformers`` will not impact
                        the partial charges produced.
        use_conformers : iterable of unit-wrapped numpy arrays, each with
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

        SUPPORTED_CHARGE_METHODS = {"mmff94", "gasteiger"}

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
            charges = [
                mmff_properties.GetMMFFPartialCharge(i) for i in range(molecule.n_atoms)
            ]
        elif partial_charge_method == "gasteiger":
            AllChem.ComputeGasteigerCharges(rdkit_molecule)
            charges = [
                float(rdatom.GetProp("_GasteigerCharge"))
                for rdatom in rdkit_molecule.GetAtoms()
            ]

        molecule.partial_charges = unit.Quantity(
            np.asarray(charges), unit.elementary_charge
        )

        if normalize_partial_charges:
            molecule._normalize_partial_charges()

    @classmethod
    def _elf_is_problematic_conformer(
        cls,
        molecule: "Molecule",
        conformer: unit.Quantity,
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
        cls,
        molecule: "Molecule",
        conformer: unit.Quantity,
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
            molecule.partial_charges.m_as(unit.elementary_charge)
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
        coordinates = conformer.m_as(unit.angstrom)

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
        selected_indices = [0]
        angstrom_tol = rms_tolerance.m_as(unit.angstrom)
        for _ in range(min(limit, molecule.n_conformers) - 1):
            selected_rms = rms_matrix[selected_indices]
            any_too_close = np.any(selected_rms < angstrom_tol, axis=0)
            if np.all(any_too_close):
                # stop if all conformers remaining are within RMS
                # threshold of any selected conformer
                break

            # add the next conformer with the largest summed RMS distance
            # to current selected conformers
            rmsdist = np.where(any_too_close, -np.inf, selected_rms.sum(axis=0))
            selected_indices.append(int(rmsdist.argmax()))

        return [ranked_conformers[i] for i in selected_indices]

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
        allow_undefined_stereo: bool = False,
        hydrogens_are_explicit: bool = False,
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
        # We need the elements of the lanthanides, actinides, and transition
        # metals as we don't want to exclude radicals in these blocks.
        d_and_f_block_elements = {
            *range(21, 31),
            *range(39, 49),
            *range(57, 81),
            *range(89, 113),
        }
        for rda in rdmol.GetAtoms():
            # See issues #1075 for some discussion on radicals
            if (
                rda.GetAtomicNum() not in d_and_f_block_elements
                and rda.GetNumRadicalElectrons() != 0
            ):
                raise RadicalsNotSupportedError(
                    "The OpenFF Toolkit does not currently support parsing molecules with S- and P-block radicals. "
                    f"Found {rda.GetNumRadicalElectrons()} radical electrons on molecule {Chem.MolToSmiles(rdmol)}."
                )

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
            # implicit units of elementary charge
            formal_charge = rda.GetFormalCharge()
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

            res = rda.GetPDBResidueInfo()
            metadata = dict()
            if res is not None:
                metadata["residue_name"] = res.GetResidueName()
                metadata["residue_number"] = res.GetResidueNumber()
                metadata["insertion_code"] = res.GetInsertionCode()
                metadata["chain_id"] = res.GetChainId()

            atom_index = offmol._add_atom(
                atomic_number,
                formal_charge,
                is_aromatic,
                name=name,
                stereochemistry=stereochemistry,
                metadata=metadata,
                invalidate_cache=False,
            )
            map_atoms[rd_idx] = atom_index
            atom_mapping[atom_index] = map_id

        offmol._invalidate_cached_properties()

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
                map_atoms[a1],
                map_atoms[a2],
                order,
                is_aromatic=is_aromatic,
                invalidate_cache=False,
            )
            map_bonds[rdb_idx] = bond_index

        offmol._invalidate_cached_properties()

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
                # Here we assume this always be angstrom
                positions = np.zeros((n_atoms, 3))
                for rd_idx, off_idx in map_atoms.items():
                    atom_coords = conf.GetPositions()[rd_idx, :]
                    positions[off_idx, :] = atom_coords
                offmol._add_conformer(unit.Quantity(positions, unit.angstrom))

        partial_charges = np.zeros(shape=offmol.n_atoms, dtype=np.float64)

        any_atom_has_partial_charge = False
        for rd_idx, rd_atom in enumerate(rdmol.GetAtoms()):
            off_idx = map_atoms[rd_idx]
            if rd_atom.HasProp("PartialCharge"):
                charge = rd_atom.GetDoubleProp("PartialCharge")
                partial_charges[off_idx] = charge
                any_atom_has_partial_charge = True
            else:
                # If some other atoms had partial charges but this one doesn't, raise an Exception
                if any_atom_has_partial_charge:
                    raise ValueError(
                        "Some atoms in rdmol have partial charges, but others do not."
                    )
        if any_atom_has_partial_charge:
            offmol.partial_charges = unit.Quantity(
                partial_charges, unit.elementary_charge
            )
        else:
            offmol.partial_charges = None
        return offmol

    to_rdkit_cache = LRUCache(maxsize=4096)

    @cached(to_rdkit_cache, key=base_wrapper._mol_to_ctab_and_aro_key)
    def _connection_table_to_rdkit(
        self, molecule, aromaticity_model=DEFAULT_AROMATICITY_MODEL
    ):
        from rdkit import Chem

        if aromaticity_model not in ALLOWED_AROMATICITY_MODELS:
            raise InvalidAromaticityModelError(
                f"Given aromaticity model {aromaticity_model} which is not in the set of allowed aromaticity models: "
                f"{ALLOWED_AROMATICITY_MODELS}"
            )

        # Create an editable RDKit molecule
        rdmol = Chem.RWMol()

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
            rdatom.SetFormalCharge(atom.formal_charge.m_as(unit.elementary_charge))
            rdatom.SetIsAromatic(atom.is_aromatic)

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

        if aromaticity_model == "OEAroModel_MDL":
            Chem.SetAromaticity(rdmol, Chem.AromaticityModel.AROMATICITY_MDL)
        else:
            raise InvalidAromaticityModelError(
                f"Given aromaticity model {aromaticity_model} which is not in the set of allowed aromaticity models:"
                f"{ALLOWED_AROMATICITY_MODELS}"
            )

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
        self._assign_rdmol_bonds_stereo(molecule, rdmol)

        # Cleanup the rdmol
        rdmol.UpdatePropertyCache(strict=False)
        Chem.GetSSSR(rdmol)

        # Forcefully assign stereo information on the atoms that RDKit
        # can't figure out. This must be done last as calling AssignStereochemistry
        # again will delete these properties (see #196).
        for rdatom, stereochemistry in undefined_stereo_atoms.items():
            rdatom.SetProp("_CIPCode", stereochemistry)

        return rdmol

    def to_rdkit(
        self, molecule: "Molecule", aromaticity_model: str = DEFAULT_AROMATICITY_MODEL
    ):
        """
        Create an RDKit molecule
        Requires the RDKit to be installed.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------

        aromaticity_model : str, optional, default="OEAroModel_MDL"
            The aromaticity model to use. Only OEAroModel_MDL is supported.

        Returns
        -------

        rdmol : rkit.RDMol
            An RDKit molecule

        Examples
        --------

        Convert a molecule to RDKit
        >>> from openff.toolkit import Molecule
        >>> ethanol = Molecule.from_smiles('CCO')
        >>> rdmol = ethanol.to_rdkit()
        """
        from rdkit import Chem, Geometry

        if aromaticity_model not in ALLOWED_AROMATICITY_MODELS:
            raise InvalidAromaticityModelError(
                f"Given aromaticity model {aromaticity_model} which is not in the set of allowed aromaticity models: "
                f"{ALLOWED_AROMATICITY_MODELS}."
            )

        # Convert the OFF molecule's connectivity table to RDKit, returning a cached rdmol if possible
        rdmol = self._connection_table_to_rdkit(
            molecule, aromaticity_model=aromaticity_model
        )
        # In case a cached rdmol was returned, make a copy of it
        rdmol = Chem.RWMol(rdmol)
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

        for index, atom in enumerate(molecule.atoms):
            rdatom = rdmol.GetAtomWithIdx(index)
            rdatom.SetProp("_Name", atom.name)

            if rdatom.GetPDBResidueInfo() is None:
                res = Chem.AtomPDBResidueInfo()
            else:
                res = rdatom.GetPDBResidueInfo()

            atom_has_any_metadata = False
            # RDKit is very naive about PDB atom names - Needs them to be exactly
            # 4 characters or the columns won't comply with PDB specification
            res.SetName(atom.name.center(4)[:4])
            if "residue_name" in atom.metadata:
                atom_has_any_metadata = True
                res.SetResidueName(atom.metadata["residue_name"])

            if "residue_number" in atom.metadata:
                atom_has_any_metadata = True
                res.SetResidueNumber(int(atom.metadata["residue_number"]))

            if "insertion_code" in atom.metadata:
                atom_has_any_metadata = True
                res.SetInsertionCode(atom.metadata["insertion_code"])

            if "chain_id" in atom.metadata:
                atom_has_any_metadata = True
                res.SetChainId(atom.metadata["chain_id"])

            if atom_has_any_metadata:
                rdatom.SetPDBResidueInfo(res)

        for bond in molecule.bonds:
            atom_indices = (
                bond.atom1.molecule_atom_index,
                bond.atom2.molecule_atom_index,
            )
            rdbond = rdmol.GetBondBetweenAtoms(*atom_indices)
            if not (bond.fractional_bond_order is None):
                rdbond.SetDoubleProp(
                    "fractional_bond_order", bond.fractional_bond_order
                )

        # Set coordinates if we have them
        if molecule._conformers:
            for conformer in molecule._conformers:
                rdmol_conformer = Chem.Conformer()
                for atom_idx in range(molecule.n_atoms):
                    x, y, z = conformer[atom_idx, :].m_as(unit.angstrom)
                    rdmol_conformer.SetAtomPosition(atom_idx, Geometry.Point3D(x, y, z))
                rdmol.AddConformer(rdmol_conformer, assignId=True)

        # Retain charges, if present
        if not (molecule._partial_charges is None):
            rdk_indexed_charges = np.zeros(shape=molecule.n_atoms, dtype=float)
            for atom_idx, charge in enumerate(molecule._partial_charges):
                charge_unitless = charge.m_as(unit.elementary_charge)
                rdk_indexed_charges[atom_idx] = charge_unitless
            for atom_idx, rdk_atom in enumerate(rdmol.GetAtoms()):
                rdk_atom.SetDoubleProp("PartialCharge", rdk_indexed_charges[atom_idx])

            # Note: We could put this outside the "if" statement, which would result in all partial charges in the
            #       resulting file being set to "n/a" if they weren't set in the Open Force Field Toolkit ``Molecule``
            Chem.CreateAtomDoublePropertyList(rdmol, "PartialCharge")

        # Return non-editable version
        return Chem.Mol(rdmol)

    def to_inchi(self, molecule: "Molecule", fixed_hydrogens: bool = False):
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

    def to_inchikey(self, molecule: "Molecule", fixed_hydrogens: bool = False) -> str:
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

    def get_tagged_smarts_connectivity(self, smarts: str):
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

        from openff.toolkit.utils.exceptions import SMIRKSParsingError

        ss = Chem.MolFromSmarts(smarts)

        if ss is None:
            # This is a SMIRKS or SMARTS parsing error? The argument and exception disagree
            raise SMIRKSParsingError(
                f"RDKit was unable to parse SMIRKS/SMARTS {smarts}"
            )

        _unique_tags = set()
        _connections = set()
        for at1 in ss.GetAtoms():
            if at1.GetAtomMapNum() == 0:
                continue
            _unique_tags.add(at1.GetAtomMapNum())
            for at2 in at1.GetNeighbors():
                if at2.GetAtomMapNum() == 0:
                    continue
                cxn_to_add = sorted([at1.GetAtomMapNum(), at2.GetAtomMapNum()])
                _connections.add(tuple(cxn_to_add))
        connections = tuple(sorted(list(_connections)))
        unique_tags = tuple(sorted(list(_unique_tags)))
        return unique_tags, connections

    @staticmethod
    def _find_smarts_matches(
        rdmol,
        smarts: str,
        aromaticity_model: str = "OEAroModel_MDL",
        unique: bool = False,
    ) -> List[Tuple[int, ...]]:
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
        unique : bool, default=False
            If True, only return unique matches. If False, return all matches. This is passed to
            RDKit's ``GetSubstructMatches`` as ``uniquify``.

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

        # Set up query.
        qmol = Chem.MolFromSmarts(smarts)  # cannot catch the error
        if qmol is None:
            raise ValueError(
                'RDKit could not parse the SMIRKS string "{}"'.format(smarts)
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
        molecule: "Molecule",
        smarts: str,
        aromaticity_model: str = "OEAroModel_MDL",
        unique: bool = False,
    ) -> List[Tuple[int, ...]]:
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
        unique : bool, default=False
            If True, only return unique matches. If False, return all matches.

        .. note :: Currently, the only supported ``aromaticity_model`` is ``OEAroModel_MDL``

        """
        rdmol = self._connection_table_to_rdkit(
            molecule, aromaticity_model=aromaticity_model
        )
        return self._find_smarts_matches(
            rdmol,
            smarts,
            aromaticity_model="OEAroModel_MDL",
            unique=unique,
        )

    def atom_is_in_ring(self, atom: "Atom") -> bool:
        """Return whether or not an atom is in a ring.

        It is assumed that this atom is in molecule.

        Parameters
        ----------
        atom : openff.toolkit.topology.molecule.Atom
            The molecule containing the atom of interest

        Returns
        -------
        is_in_ring : bool
            Whether or not the atom is in a ring.

        Raises
        ------
        NotAttachedToMoleculeError
        """
        if atom.molecule is None:
            raise NotAttachedToMoleculeError(
                "This Atom does not belong to a Molecule object"
            )

        molecule = atom.molecule
        atom_index = atom.molecule_atom_index

        rdmol = molecule.to_rdkit()
        rdatom = rdmol.GetAtomWithIdx(atom_index)

        is_in_ring = rdatom.IsInRing()

        return is_in_ring

    def bond_is_in_ring(self, bond: "Bond") -> bool:
        """Return whether or not a bond is in a ring.

        It is assumed that this atom is in molecule.

        Parameters
        ----------
        bond : openff.toolkit.topology.molecule.Bond
            The molecule containing the atom of interest

        Returns
        -------
        is_in_ring : bool
            Whether or not the bond of index `bond_index` is in a ring

        Raises
        ------
        NotAttachedToMoleculeError
        """
        if bond.molecule is None:
            raise NotAttachedToMoleculeError(
                "This Bond does not belong to a Molecule object"
            )

        molecule = bond.molecule
        rdmol = molecule.to_rdkit()

        # Molecule.to_rdkit() is NOT guaranteed to preserve bond ordering,
        # so we must look up the corresponding bond via its constituent atom indices
        rdbond = rdmol.GetBondBetweenAtoms(bond.atom1_index, bond.atom2_index)
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
    def _detect_undefined_stereo(
        cls,
        rdmol,
        err_msg_prefix="",
        raise_warning=False,
    ):
        """Raise UndefinedStereochemistryError if the RDMol has undefined stereochemistry.

        Parameters
        ----------
        rdmol : rdkit.Chem.Mol
            The RDKit molecule.
        err_msg_prefix : str, optional
            A string to prepend to the error message (but not the warning).
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
            msg = "RDMol has unspecified stereochemistry. "
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
                msg = err_msg_prefix + msg
                raise UndefinedStereochemistryError(msg)

    @staticmethod
    def _constrain_end_directions(
        *values: int, bond_indices: List[int], flip_direction: Dict[int, bool]
    ) -> bool:
        """A constraint applied when mapping global E/Z stereochemistry into local RDKit
        bond directions that ensures that the 'left' bonds point in opposite directions
        (i.e. one has to be up and one has to be down) and likewise for the 'right' bonds
        """
        # Account for bond "direction" using flip_directions dict, see more thorough comment
        # in _constrain_rank for details.
        bond_directions = [
            (value if not flip_direction[i] else 1 - value)
            for i, value in zip(bond_indices, values)
        ]
        unique_bond_directions = set(bond_directions)
        return len(unique_bond_directions) == len(values)

    @staticmethod
    def _constrain_rank(
        *values: int,
        bond_indices: List[int],
        flip_direction: Dict[int, bool],
        expected_stereo: str,
    ) -> bool:
        """A constraint applied when mapping global E/Z stereochemistry into local RDKit
        bond directions that ensures that the 'left' bond with the highest CIP rank
        and the 'right' bond with the highest CIP rank point either in the same direction
        if Z stereo or opposite directions if E.
        """
        # The "value" for each bond is ultimately set to either 0 (down) or 1 (up).
        # However, we also need to know the "direction" of the bond to make sense
        # of this - Is it coming FROM, or going TO this double bond while going down or up?
        # This information is in the flipped_values dict. This code assumes that the
        # "normal" case is when the neighboring bonds are coming FROM the double bond,
        # and where that's not true, the following line switches the
        # meaning of "down" and "up" to give us the desired meaning.
        bond_directions = [
            (value if not flip_direction[i] else 1 - value)
            for i, value in zip(bond_indices, values)
        ]

        # Test for equality of items in flipped_values by turning it into a set and
        # counting how many values remain.
        unique_bond_directions = set(bond_directions)
        if expected_stereo == "E":
            return len(unique_bond_directions) == min(2, len(bond_indices))
        elif expected_stereo == "Z":
            return len(unique_bond_directions) == 1
        else:
            raise NotImplementedError()

    @classmethod
    def _assign_rdmol_bonds_stereo(cls, off_molecule: "Molecule", rd_molecule):
        """Copy the info about bonds stereochemistry from the OFF Molecule to RDKit Mol.
        The method proceeds by formulating mapping global E/Z stereo information onto
        local 'bond directions' as a constraint satisfaction problem (CSP).
        In this formalism, the variables correspond to the indices of the bonds
        neighbouring a stereogenic bond, the domain for each variable is 0 (down)
        or 1 (up), and the constraints are designed to ensure the correct E/Z stereo
        is yielded after assignment of the directions.
        """
        from constraint import Problem
        from rdkit import Chem

        _RD_STEREO_TO_STR = {
            Chem.BondStereo.STEREOE: "E",
            Chem.BondStereo.STEREOZ: "Z",
        }

        stereogenic_bonds = [
            bond for bond in off_molecule.bonds if bond.stereochemistry
        ]

        if len(stereogenic_bonds) == 0:
            return

        # Needed to ensure the _CIPRank is present. Note that, despite the kwargs that look like
        # they could mangle existing stereo, it is actually preserved.
        Chem.AssignStereochemistry(
            rd_molecule, cleanIt=True, force=True, flagPossibleStereoCenters=True
        )

        csp_problem = Problem()
        csp_variables = set()

        for bond in stereogenic_bonds:
            # Here we use a notation where atoms 'b' and 'c' are the two atoms involved
            # in the double bond, while 'a' corresponds to a neighbour of 'b' and 'd' a
            # neighbour of 'c'.
            atom_b, atom_c = bond.atom1, bond.atom2
            index_b, index_c = atom_b.molecule_atom_index, atom_c.molecule_atom_index

            indices_a = [
                n.molecule_atom_index for n in atom_b.bonded_atoms if n != atom_c
            ]
            indices_d = [
                n.molecule_atom_index for n in atom_c.bonded_atoms if n != atom_b
            ]
            # A stereogenic double bond should either involve atoms with degree 3
            # (e.g. carbon) or degree 2 (e.g. divalent nitrogen).
            assert 1 <= len(indices_a) <= 2 and 1 <= len(indices_d) <= 2

            # Identify the highest CIP-ranked bond coming off each side of the double
            # bond. This lets us later add a constraint to ensure that we have the
            # correct E/Z stereochemistry
            ranks_a = [
                int(rd_molecule.GetAtomWithIdx(i).GetProp("_CIPRank"))
                for i in indices_a
            ]
            index_a = indices_a[np.argmax(ranks_a)]
            ranks_d = [
                int(rd_molecule.GetAtomWithIdx(i).GetProp("_CIPRank"))
                for i in indices_d
            ]
            index_d = indices_d[np.argmax(ranks_d)]

            index_ab = rd_molecule.GetBondBetweenAtoms(index_a, index_b).GetIdx()
            index_cd = rd_molecule.GetBondBetweenAtoms(index_c, index_d).GetIdx()

            flip_direction = {}

            # Collect lists of the indices of the bonds that appear to the 'left' of
            # and 'right' of the stereogenic bond so we can constrain their directions
            # so that all 'left' bonds do not, for example, point up.
            constraints_ab: List[int] = []
            constraints_cd: List[int] = []

            for index_pair, constraints_list in [
                ((index_a, index_b), constraints_ab) for index_a in indices_a
            ] + [((index_d, index_c), constraints_cd) for index_d in indices_d]:
                # Each single bond neighboring a double bond needs to be defined as a
                # "variable" in the CSP problem. Here, each bond is identified by its
                # bond index in the RDMol (note: this is not guaranteed to be the same
                # as the corresponding bond's index in the OFFMol). This also ensures
                # that we don't define the same bond twice.
                rd_bond = rd_molecule.GetBondBetweenAtoms(*index_pair)
                rd_bond_index = rd_bond.GetIdx()

                if rd_bond_index not in csp_variables:
                    csp_problem.addVariable(rd_bond_index, [1, 0])  # 0 = down, 1 = up
                    csp_variables.add(rd_bond_index)

                # The direction of the bond should point from the double bond to its
                # neighbour. If the bond is pointing from the neighbour to the double
                # bond instead, we need to flip the direction of the bond. See
                # rdkit/Code/GraphMol/Chirality.cpp:findAtomNeighborDirHelper for more
                # details.
                flip_direction[rd_bond_index] = (
                    rd_bond.GetBeginAtomIdx() != index_pair[1]
                )
                constraints_list.append(rd_bond_index)

            # Add one constraint corresponding to the highest-ranked bond on OPPOSITE
            # sides of the double bond, to ensure that they are oriented up/down to
            # achieve the correct E/Z value of the double bond.
            csp_problem.addConstraint(
                functools.partial(
                    cls._constrain_rank,
                    bond_indices=[index_ab, index_cd],
                    flip_direction=flip_direction,
                    expected_stereo=bond.stereochemistry,
                ),
                [index_ab, index_cd],
            )
            # Add constraint(s) corresponding ALL bonds on the SAME side of the double
            # bond, to ensure that they do not all take the same value (if one is "up",
            # the other can not also be "up").
            csp_problem.addConstraint(
                functools.partial(
                    cls._constrain_end_directions,
                    bond_indices=constraints_ab,
                    flip_direction=flip_direction,
                ),
                constraints_ab,
            )
            csp_problem.addConstraint(
                functools.partial(
                    cls._constrain_end_directions,
                    bond_indices=constraints_cd,
                    flip_direction=flip_direction,
                ),
                constraints_cd,
            )

        # Do not assume that every solution found by the solver is valid.
        # Iterate through the solutions and ensure that the
        # desired E/Z marks have really been achieved.
        has_solution = False

        for solution in csp_problem.getSolutionIter():
            for rd_bond_index, direction in solution.items():
                rd_bond = rd_molecule.GetBondWithIdx(rd_bond_index)
                if direction == 0:
                    rd_bond.SetBondDir(Chem.BondDir.ENDDOWNRIGHT)
                elif direction == 1:
                    rd_bond.SetBondDir(Chem.BondDir.ENDUPRIGHT)
                else:
                    raise NotImplementedError()

            Chem.AssignStereochemistry(rd_molecule, cleanIt=True, force=True)

            # Verify that there are no stereo mismatches between the original
            # OFFMol and the newly-assigned RDMol
            stereo_mismatch = False
            for off_bond in stereogenic_bonds:
                rd_bond = rd_molecule.GetBondBetweenAtoms(
                    off_bond.atom1_index, off_bond.atom2_index
                )
                rd_stereo_string = _RD_STEREO_TO_STR.get(rd_bond.GetStereo(), None)
                if off_bond.stereochemistry != rd_stereo_string:
                    stereo_mismatch = True
                    break

            if stereo_mismatch:
                continue

            has_solution = True
            break

        assert has_solution, "E/Z stereo could not be converted to local stereo"
