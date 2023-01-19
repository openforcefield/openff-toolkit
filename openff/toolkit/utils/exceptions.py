from typing import TYPE_CHECKING, DefaultDict, Dict, List, Mapping, Optional, Set, Tuple

if TYPE_CHECKING:
    from openmm.app import Atom as OpenMMAtom
    from openmm.app import Residue as OpenMMResidue
    from openmm.app import Topology as OpenMMTopology


class OpenFFToolkitException(Exception):
    """Base exception for custom exceptions raised by the OpenFF Toolkit"""

    def __init__(self, msg):
        super().__init__(msg)
        self.msg = msg

    def __str__(self):
        return self.msg


class IncompatibleUnitError(OpenFFToolkitException):
    """
    Exception for when a parameter is in the wrong units for a ParameterHandler's unit system
    """


class MissingPackageError(OpenFFToolkitException):
    """This function requires a package that is not installed."""


class ToolkitUnavailableException(OpenFFToolkitException):
    """The requested toolkit is unavailable."""

    # TODO: Allow toolkit to be specified and used in formatting/printing exception.


class LicenseError(ToolkitUnavailableException):
    """This function requires a license that cannot be found."""


class InvalidToolkitError(OpenFFToolkitException):
    """A non-toolkit object was received when a toolkit object was expected"""


class InvalidToolkitRegistryError(OpenFFToolkitException):
    """An object other than a ToolkitRegistry or toolkit wrapper was received"""


class UndefinedStereochemistryError(OpenFFToolkitException):
    """A molecule was attempted to be loaded with undefined stereochemistry"""


class GAFFAtomTypeWarning(RuntimeWarning):
    """A warning raised if a loaded mol2 file possibly uses GAFF atom types."""


class ChargeMethodUnavailableError(OpenFFToolkitException):
    """A toolkit does not support the requested partial_charge_method combination"""


class IncorrectNumConformersError(OpenFFToolkitException):
    """The requested partial_charge_method expects a different number of conformers than was provided"""


class IncorrectNumConformersWarning(Warning):
    """The requested partial_charge_method expects a different number of conformers than was provided"""


class ChargeCalculationError(OpenFFToolkitException):
    """An unhandled error occurred in an external toolkit during charge calculation"""


class ConformerGenerationError(OpenFFToolkitException):
    """Conformer generation via a wrapped toolkit failed."""


class InvalidIUPACNameError(OpenFFToolkitException):
    """Failed to parse IUPAC name"""


class AntechamberNotFoundError(OpenFFToolkitException):
    """The antechamber executable was not found"""


class SMILESParseError(OpenFFToolkitException, ValueError):
    """The record could not be parsed into the given format"""


class RadicalsNotSupportedError(OpenFFToolkitException):
    """The OpenFF Toolkit does not currently support parsing molecules with radicals."""


class InvalidConformerError(OpenFFToolkitException):
    """
    This error is raised when the conformer added to the molecule
    has a different connectivity to that already defined.
    or any other conformer related issues.
    """


# TODO: Remove in favor of SMILESParseError? They are used in different modules
class SmilesParsingError(OpenFFToolkitException):
    """
    This error is raised when parsing a SMILES string results in an error.
    """


class NotAttachedToMoleculeError(OpenFFToolkitException):
    """Exception for when a component does not belong to a Molecule object, but is queried"""


class NotInTopologyError(OpenFFToolkitException):
    """An object was not found in a topology."""


class RemapIndexError(OpenFFToolkitException):
    """An error with indices used to remap a molecule"""


class AtomNotInTopologyError(NotInTopologyError):
    """An atom was not found in a topology."""


class MoleculeNotInTopologyError(NotInTopologyError):
    """A molecule was not found in a topology."""


class InvalidAtomMetadataError(OpenFFToolkitException):
    """The program attempted to set atom metadata to an invalid type"""


class DuplicateUniqueMoleculeError(OpenFFToolkitException):
    """
    Exception for when the user provides indistinguishable unique molecules when trying to identify atoms from a PDB
    """


class NotBondedError(OpenFFToolkitException):
    """
    Exception for when a function requires a bond between two atoms, but none is present
    """


class InvalidBondOrderError(OpenFFToolkitException):
    """
    Exception for passing a non-int to `Molecule.bond_order`
    """


class InvalidBoxVectorsError(OpenFFToolkitException):
    """
    Exception for setting invalid box vectors
    """


class InvalidPeriodicityError(OpenFFToolkitException):
    """
    Exception for setting invalid periodicity
    """


class MissingUniqueMoleculesError(OpenFFToolkitException):
    """
    Exception for when unique_molecules is required but not found
    """


class SMIRKSMismatchError(OpenFFToolkitException):
    """
    Exception for cases where smirks are inappropriate
    for the environment type they are being parsed into
    """


# TODO: This should be renamed to SMARTSParsingError
class SMIRKSParsingError(OpenFFToolkitException):
    """
    Exception for when SMIRKS are not parseable for any environment
    """


class ParameterHandlerRegistrationError(OpenFFToolkitException):
    """
    Exception for errors in ParameterHandler registration
    """


class SMIRNOFFVersionError(OpenFFToolkitException):
    """
    Exception thrown when an incompatible SMIRNOFF version data structure is attempted to be read.
    """


class SMIRNOFFAromaticityError(OpenFFToolkitException):
    """
    Exception thrown when an incompatible SMIRNOFF aromaticity model is checked for compatibility.
    """


class InvalidAromaticityModelError(OpenFFToolkitException, ValueError):
    """
    General exception for errors while setting the aromaticity model of a Topology.
    """


class SMIRNOFFParseError(OpenFFToolkitException):
    """
    Error for when a SMIRNOFF data structure is not parseable by a ForceField
    """

    # TODO: Remove ParseError altogether by v0.11.0


class PartialChargeVirtualSitesError(OpenFFToolkitException):
    """
    Exception thrown when partial charges cannot be computed for a Molecule because the ForceField applies virtual
    sites.
    """


class SMIRNOFFSpecError(OpenFFToolkitException):
    """
    Exception for when data is noncompliant with the SMIRNOFF data specification.
    """


class SMIRNOFFSpecUnimplementedError(OpenFFToolkitException):
    """
    Exception for when a portion of the SMIRNOFF specification is not yet implemented.
    """


class FractionalBondOrderInterpolationMethodUnsupportedError(OpenFFToolkitException):
    """
    Exception for when an unsupported fractional bond order interpolation assignment method is called.
    """


class NotEnoughPointsForInterpolationError(OpenFFToolkitException):
    """Exception for when less than two points are provided for interpolation"""


class IncompatibleParameterError(OpenFFToolkitException):
    """
    Exception for when a set of parameters is scientifically/technically incompatible with another
    """


class UnassignedValenceParameterException(OpenFFToolkitException):
    """Exception raised when there are valence terms for which a ParameterHandler can't find parameters."""


class UnassignedBondParameterException(UnassignedValenceParameterException):
    """Exception raised when there are bond terms for which a ParameterHandler can't find parameters."""


class UnassignedAngleParameterException(UnassignedValenceParameterException):
    """Exception raised when there are angle terms for which a ParameterHandler can't find parameters."""


class UnassignedProperTorsionParameterException(UnassignedValenceParameterException):
    """Exception raised when there are proper torsion terms for which a ParameterHandler can't find parameters."""


class UnassignedMoleculeChargeException(OpenFFToolkitException):
    """Exception raised when no charge method is able to assign charges to a molecule."""


class DuplicateParameterError(OpenFFToolkitException):
    """Exception raised when trying to add a ParameterType that already exists"""


class ParameterLookupError(OpenFFToolkitException):
    """Exception raised when something goes wrong in a parameter lookup in
    ParameterHandler.__getitem__"""


class DuplicateVirtualSiteTypeException(OpenFFToolkitException):
    """Exception raised when trying to register two different virtual site classes with the same 'type'"""


class CallbackRegistrationError(OpenFFToolkitException, TypeError):
    """Error raised when callback registration fails."""


class HierarchySchemeWithIteratorNameAlreadyRegisteredException(OpenFFToolkitException):
    """Exception raised when trying to add a HierarchyScheme to a molecule
    that already has one with the same iterator name"""


# TODO: Should be a subclass of KeyError? Should be replaced by KeyError?
class HierarchySchemeNotFoundException(OpenFFToolkitException):
    """Exception raised when trying to access a HierarchyScheme from a molecule
    that doesn't have one with the given iterator name"""


class MissingIndexedAttributeError(
    OpenFFToolkitException, IndexError, KeyError, AttributeError
):
    """Error raised when an indexed attribute does not exist"""


class MissingPartialChargesError(OpenFFToolkitException, ValueError):
    """Error raised when a molecule is missing partial charges in a context in which it is expected to have them."""


class UnsupportedMoleculeConversionError(OpenFFToolkitException):
    """Error raised when attempting to instantiate a Molecule with insufficient inputs."""


class InconsistentStereochemistryError(OpenFFToolkitException):
    """
    Error raised when stereochemistry is inconsistent before and after conversions between molecule representations.
    """


class UnsupportedFileTypeError(OpenFFToolkitException):
    """Error raised when attempting to parse an unsupported file type."""


class MultipleMoleculesInPDBError(OpenFFToolkitException):
    """Error raised when a multiple molecules are found when one was expected"""


class WrongShapeError(OpenFFToolkitException):
    """Error raised when an array of the wrong shape is found"""


class UnassignedChemistryInPDBError(OpenFFToolkitException, ValueError):
    """
    Error raised when a bond or atom in a polymer could not be assigned chemical information.
    """

    def __init__(
        self,
        msg: Optional[str] = None,
        substructure_library: Optional[Dict[str, Tuple[str, List[str]]]] = None,
        omm_top: Optional["OpenMMTopology"] = None,
        unassigned_bonds: Optional[List[Tuple[int, int]]] = None,
        unassigned_atoms: Optional[List[int]] = None,
        matches: Optional[DefaultDict[int, List[str]]] = None,
    ):
        if omm_top is not None:
            self.omm_top = omm_top
            self._atoms: List["OpenMMAtom"] = list(omm_top.atoms())
            self._bonds: List[Tuple["OpenMMAtom", "OpenMMAtom"]] = list(omm_top.bonds())

        self.substructure_library = substructure_library
        self.unassigned_bonds = [] if unassigned_bonds is None else unassigned_bonds
        self.unassigned_atoms = [] if unassigned_atoms is None else unassigned_atoms
        self.matches = matches

        message = (
            ["Some bonds or atoms in the input could not be identified.", ""]
            if msg is None
            else [msg, ""]
        )

        message += [
            *self.missing_hydrogens_hint(),
            *self.multiple_chains_hint(),
            *self.mismatched_atom_names_hint(),
            *self.unknown_residue_hint(),
            *self.assigned_residue_name_mismatch_hint(),
            *self.unassigned_atoms_err(),
            *self.unassigned_bonds_err(),
        ]

        super().__init__("\n".join(message))

    def residue_of_atom_as_str(self, atom_index: int) -> str:
        res = self._atoms[atom_index].residue
        return self.fmt_residue(res.name, res.id)

    @staticmethod
    def fmt_residue(name: str, num: int) -> str:
        return f"{name}#{num:0>4}"

    def unassigned_atoms_err(self) -> List[str]:
        if self.unassigned_atoms and self.omm_top:
            return [
                (
                    f"Error: The following {len(self.unassigned_atoms)} atoms exist in the input "
                    + "but could not be assigned chemical information from the "
                    + "substructure library:"
                ),
                *(
                    f"    Atom {i: >5} ({self._atoms[i].name}) in residue "
                    + f"{self.residue_of_atom_as_str(i)}"
                    for i in self.unassigned_atoms
                ),
                "",
            ]
        return []

    def unassigned_bonds_err(self) -> List[str]:
        if not (self.unassigned_bonds and self.unassigned_atoms and self.omm_top):
            return []

        unassigned_bonds_with_assigned_atoms = [
            (i_a, i_b)
            for i_a, i_b in self.unassigned_bonds
            if i_a not in self.unassigned_atoms and i_b not in self.unassigned_atoms
        ]

        elided_bonds_hint = []
        if unassigned_bonds_with_assigned_atoms != self.unassigned_bonds:
            elided_bonds_hint.append(
                "Note: Bonds involving atoms that could not be assigned chemical "
                + "information have been elided from this list."
            )

        if unassigned_bonds_with_assigned_atoms:
            return [
                (
                    f"Error: The following {len(unassigned_bonds_with_assigned_atoms)} "
                    + "bonds exist in the input but could not be assigned "
                    + "chemical information from the substructure library:"
                ),
                *(
                    f"    Bond between atom {i_a: >5} ({self._atoms[i_a].name}) "
                    + f"in {self.residue_of_atom_as_str(i_a)} "
                    + f"and atom {i_b: >5} ({self._atoms[i_b].name}) "
                    + f"in {self.residue_of_atom_as_str(i_b)}"
                    for i_a, i_b in unassigned_bonds_with_assigned_atoms
                ),
                *elided_bonds_hint,
                "",
            ]
        return []

    def missing_hydrogens_hint(self) -> List[str]:
        if self.omm_top and all(atom.element.symbol != "H" for atom in self._atoms):
            return [
                "Hint: There are no hydrogens in the input. The OpenFF Toolkit "
                + "requires explicit hydrogens to avoid ambiguities in protonation "
                + "state or bond order. Try generating hydrogens with another package "
                + "and trying again.",
                "",
            ]
        return []

    def unknown_residue_hint(self) -> List[str]:
        if self.substructure_library and self.omm_top:
            unassigned_resnames = [
                self._atoms[i].residue.name for i in self.unassigned_atoms
            ]
            unknown_resnames = set(
                [
                    resname
                    for resname in unassigned_resnames
                    if resname not in self.substructure_library
                ]
            )

            if "HOH" in unknown_resnames:
                solvent_note = [
                    "Note: 'HOH' is a residue code for water. You may have "
                    + "crystallographic waters in your PDB file. Please remove "
                    + "these before proceeding; they can be added back to the "
                    + "topology later."
                ]
            else:
                solvent_note = [""]

            if unknown_resnames:
                return [
                    "Hint: The following residue names with unassigned atoms were not "
                    + "found in the substructure library. While the OpenFF Toolkit "
                    + "identifies residues by matching chemical substructures rather "
                    + "than by residue name, it currently only supports the 20 "
                    + "'canonical' amino acids.",
                    *(f"    {resname}" for resname in unknown_resnames),
                    *solvent_note,
                    "",
                ]
        return []

    def multiple_chains_hint(self) -> List[str]:
        if self.omm_top.getNumChains() > 1:
            return [
                "Hint: The input has multiple chain identifiers. The OpenFF "
                + "Toolkit only supports single-molecule PDB files. Please "
                + "split the file into individual chains and load each "
                + "seperately.",
                "",
            ]

        return []

    def assigned_residue_name_mismatch_hint(self) -> List[str]:
        from collections import defaultdict

        if not self.matches:
            return []

        # Construct a map from input residues to assigned resnames
        residues: Mapping[Tuple[str, str], Set[str]] = defaultdict(set)
        for atom in self.omm_top.atoms():
            input_resname: str = atom.residue.name
            input_resnum: str = atom.residue.id
            matched_resnames = self.matches[atom.index]
            # Only the first match is assigned, so throw out the others
            assigned_resname = next(iter(matched_resnames), "No match")

            residues[(input_resname, input_resnum)].add(assigned_resname)

        # Filter out residues where assigned resname doesn't match the input
        residues = {
            (input_resname, input_resnum): assigned_resnames
            for (input_resname, input_resnum), assigned_resnames in residues.items()
            if set([input_resname]) != assigned_resnames
        }

        if residues:
            return [
                "Hint: The following residues were assigned names that do not "
                + "match the residue name in the input, or could not be assigned "
                + "residue names at all. This may indicate that atoms are "
                + "missing from the input or some other error. The OpenFF "
                + "Toolkit requires all atoms, including hydrogens, to be "
                + "explicit in the input to avoid ambiguities in protonation "
                + "state or bond order:",
                *(
                    (
                        f"    Input residue {self.fmt_residue(resname, int(resnum))} "
                        + f"contains atoms matching substructures {assigned_resnames}"
                    )
                    for (
                        resname,
                        resnum,
                    ), assigned_resnames in residues.items()
                ),
                "",
            ]

        return []

    def mismatched_atom_names_hint(self) -> List[str]:
        from collections import defaultdict

        if not (self.omm_top and self.substructure_library):
            return []

        # Collect all the unassigned atoms by residue
        unassigned_residues: Mapping[OpenMMResidue, List[OpenMMAtom]] = defaultdict(
            list
        )
        for i in self.unassigned_atoms:
            atom = self._atoms[i]
            res: OpenMMResidue = atom.residue
            unassigned_residues[res].append(atom)

        # Mark residues that don't have the right number and elements of atoms
        # by clearing their atoms lists
        for res, atoms in unassigned_residues.items():
            try:
                library_res = self.substructure_library[res.name]
            except KeyError:
                # Residue is not in substructure library at all!
                atoms.clear()
                continue
            for library_elements, names in library_res:
                residue_elements = sorted(atom.element.symbol for atom in res.atoms())
                if library_elements == residue_elements:
                    # Prune the atoms down to just those whose names don't match
                    atoms[:] = [atom for atom in atoms if atom.name not in names]
                    # We're done with this residue
                    break
            else:
                # smiles for loop did not break, so clear the atom list
                atoms.clear()
        misnamed = {res: atoms for res, atoms in unassigned_residues.items() if atoms}

        if misnamed:
            return [
                "Hint: The following residues have the right numbers of the "
                + "right elements to match a substructure with the same name as "
                + "the input residue, but did not match. This most likely "
                + "suggests that their atom names do not match those in the "
                + "substructure library. Try renaming misnamed atoms according "
                + "to the PDB Chemical Component Dictionary.",
                *(
                    (
                        f"    Input residue {self.fmt_residue(res.name, int(res.id))} "
                        + f"has misnamed atoms {', '.join(a.name for a in atoms)}"
                    )
                    for res, atoms in misnamed.items()
                ),
                "",
            ]

        return []
