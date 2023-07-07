"""Interfaces with Amber."""
import textwrap
from collections import defaultdict
from collections.abc import Iterable
from copy import deepcopy
from pathlib import Path
from typing import Optional, Union

import numpy as np
from openff.toolkit import Topology
from openff.units import unit

from openff.interchange import Interchange
from openff.interchange.components.toolkit import _get_num_h_bonds
from openff.interchange.constants import (
    _PME,
    AMBER_COULOMBS_CONSTANT,
    kcal_mol,
    kcal_mol_a2,
    kcal_mol_rad2,
)
from openff.interchange.exceptions import (
    UnsupportedExportError,
    UnsupportedMixingRuleError,
)
from openff.interchange.models import PotentialKey


def _write_text_blob(file, blob):
    def chunker(seq, size):
        return (seq[pos : pos + size] for pos in range(0, len(seq), size))

    if blob == "":
        file.write("\n")
    else:
        for line in chunker(blob, 80):
            file.write(line + "\n")


def _flatten(list_of_lists: Iterable[list[int]]) -> list[int]:
    return [item for sublist in list_of_lists for item in sublist]


def _get_per_atom_exclusion_lists(
    topology: "Topology",
) -> dict[int, defaultdict[int, list[int]]]:
    """
    Get the excluded atoms of each atom in this topology.

    Parameters
    ----------
    topology: Topology
        OpenFF Topology

    Returns
    -------
    per_atom_exclusions: dict[int, list[int]]
        keys: atom indices (OpenFF atoms, zero-indexed)
        values: defaultdict[int, list[int]]
            list of atom indices (OpenFF atoms, zero-indexed) that are excluded from this atom,
            keyed by their separation distance in the graph (1 for bonded, 2 for in an angle, 3
            for in a dihedral)

    """
    per_atom_exclusions: dict[int, defaultdict[int, list[int]]] = dict()

    for atom1 in topology.atoms:
        # Excluded atoms _on this atom_
        this_atom_exclusions = defaultdict(list)
        index1 = topology.atom_index(atom1)

        for atom2 in atom1.bonded_atoms:
            index2 = topology.atom_index(atom2)

            if index2 > index1:
                # These atoms are bonded
                this_atom_exclusions[1].append(index2)

            for atom3 in atom2.bonded_atoms:
                atom3_index = topology.atom_index(atom3)

                if atom3_index > index1:
                    # These atoms are included in an angle
                    this_atom_exclusions[2].append(atom3_index)

                for atom4 in atom3.bonded_atoms:
                    atom4_index = topology.atom_index(atom4)
                    if atom4_index > index1 and atom4_index:
                        # These atoms are included in a torsion
                        this_atom_exclusions[3].append(atom4_index)

        per_atom_exclusions[index1] = this_atom_exclusions

    return per_atom_exclusions


def _get_exclusion_lists(
    per_atom_exclusions: dict[int, defaultdict[int, list[int]]],
) -> tuple[list[int], list[int]]:
    """
    Convert a per-atom exclusion dict to Amer structures.

    Parameters
    ----------
    per_atom_exclusions: dict[int, list[int]]
        keys: atom indices (OpenFF atoms, zero-indexed)
        values: list of atom indices (OpenFF atoms, zero-indexed) that are excluded from this atom

        See _getper_atom_exclusion_lists

    Returns
    -------
    number_excluded_atoms: list[int]
        Number of excluded atoms for each atom (OpenFF atoms, zero-indexed)
    excluded_atoms_list: list[int]
        Flattened list of per-atom exclusions (Amber atoms, one-indexed).
        See EXCLUDED_ATOMS_LIST in https://ambermd.org/prmtop.pdf

    """
    excluded_atoms_list: list[int] = list()
    number_excluded_atoms: list[int] = list()

    for this_atom, this_atom_exclusions in per_atom_exclusions.items():
        if len(this_atom_exclusions) == 0:
            # Amber expects this list to have a 1 in it, pointing to a non-existent
            # atom with index 0, when an atom has no exclusions
            tmp = {0}

        else:
            tmp = set()
            for other_atom in _flatten(this_atom_exclusions.values()):
                if other_atom > this_atom:
                    tmp.add(other_atom + 1)

        for val in sorted(tmp):
            excluded_atoms_list.append(val)

        number_excluded_atoms.append(len(tmp))

    return number_excluded_atoms, excluded_atoms_list


def _get_bond_lists(
    interchange: "Interchange",
    atomic_numbers: tuple,
    potential_key_to_bond_type_mapping: dict[PotentialKey, int],
) -> tuple[list[int], list[int]]:
    bonds_inc_hydrogen: list[int] = list()
    bonds_without_hydrogen: list[int] = list()

    # TODO: Should probably build bond lists and exclusions without assuming bond physics
    for bond, key_ in interchange["Bonds"].key_map.items():
        bond_type_index = potential_key_to_bond_type_mapping[key_]

        these_atomic_numbers = tuple(
            atomic_numbers[index] for index in bond.atom_indices
        )

        bond_indices = tuple(sorted(bond.atom_indices))

        bonds_list = (
            bonds_inc_hydrogen if 1 in these_atomic_numbers else bonds_without_hydrogen
        )

        bonds_list.append(bond_indices[0] * 3)
        bonds_list.append(bond_indices[1] * 3)
        bonds_list.append(bond_type_index + 1)

    return bonds_inc_hydrogen, bonds_without_hydrogen


def _get_angle_lists(
    interchange: "Interchange",
    atomic_numbers: tuple,
    potential_key_to_angle_type_mapping: dict[PotentialKey, int],
) -> tuple[list[int], list[int]]:
    angles_inc_hydrogen: list[int] = list()
    angles_without_hydrogen: list[int] = list()

    for angle, key_ in interchange["Angles"].key_map.items():
        angle_type_index = potential_key_to_angle_type_mapping[key_]

        these_atomic_numbers = tuple(
            atomic_numbers[index] for index in angle.atom_indices
        )

        angle_indices = list(angle.atom_indices)
        angle_indices = (
            angle_indices
            if angle_indices[0] < angle_indices[-1]
            else list(reversed(angle_indices))
        )

        angles_list = (
            angles_inc_hydrogen
            if 1 in these_atomic_numbers
            else angles_without_hydrogen
        )

        angles_list.append(angle_indices[0] * 3)
        angles_list.append(angle_indices[1] * 3)
        angles_list.append(angle_indices[2] * 3)
        angles_list.append(angle_type_index + 1)

    return angles_inc_hydrogen, angles_without_hydrogen


def _get_dihedral_lists(
    interchange: "Interchange",
    atomic_numbers: tuple,
    potential_key_to_dihedral_type_mapping: dict[PotentialKey, int],
    already_counted: set[tuple[int, ...]],
    exclusions: dict[int, defaultdict[int, list[int]]],  # per-atom exclusions
) -> tuple[list[int], list[int]]:
    dihedrals_inc_hydrogen: list[int] = list()
    dihedrals_without_hydrogen: list[int] = list()

    if "ProperTorsions" in interchange.collections:
        for dihedral, proper_key in interchange["ProperTorsions"].key_map.items():
            dihedral_type_index = potential_key_to_dihedral_type_mapping[proper_key]

            atom1_index, atom2_index, atom3_index, atom4_index = dihedral.atom_indices

            # Since 0 can't be negative, attempt to re-arrange this torsion
            # such that the third atom listed is negative.
            # This should only be strictly necessary when _14_tag is -1, but
            # ParmEd likes to always flip it, and always flipping should be harmless.
            # Could put this in an if block if desired.
            if atom3_index == 0:
                atom1_index, atom2_index, atom3_index, atom4_index = reversed(
                    (
                        atom1_index,
                        atom2_index,
                        atom3_index,
                        atom4_index,
                    ),
                )

            these_atomic_numbers = tuple(
                atomic_numbers[index] for index in dihedral.atom_indices
            )

            # See if the non-bonded interactions of the first and fourth atoms need to be exluded

            atom1_neighbors = exclusions[atom1_index][1] + exclusions[atom1_index][2]
            atom4_neighbors = exclusions[atom4_index][1] + exclusions[atom4_index][2]

            _sorted = tuple(sorted([atom1_index, atom4_index]))

            if atom4_index in atom1_neighbors:
                # Exclude because these atoms are separated by only one or two bonds
                exclude_nonbonded = -1

            elif atom1_index in atom4_neighbors:
                # Exclude because these atoms are separated by only one or two bonds
                exclude_nonbonded = -1

            elif _sorted in already_counted:
                # Exclude because these were already counted once in another torsion
                exclude_nonbonded = -1

            else:
                exclude_nonbonded = 1

            dihedrals_list = (
                dihedrals_inc_hydrogen
                if 1 in these_atomic_numbers
                else dihedrals_without_hydrogen
            )

            dihedrals_list.append(atom1_index * 3)
            dihedrals_list.append(atom2_index * 3)
            dihedrals_list.append(atom3_index * 3 * exclude_nonbonded)
            dihedrals_list.append(atom4_index * 3)
            dihedrals_list.append(dihedral_type_index + 1)

            already_counted.add(_sorted)

    if "ImproperTorsions" in interchange.collections:
        for improper, improper_key in interchange["ImproperTorsions"].key_map.items():
            dihedral_type_index = potential_key_to_dihedral_type_mapping[improper_key]

            atom1_index, atom2_index, atom3_index, atom4_index = improper.atom_indices

            # Assume that no improper torsions include 1-4 pairs, so don't check nor track them

            these_atomic_numbers = tuple(
                atomic_numbers[index] for index in improper.atom_indices
            )

            dihedrals_list = (
                dihedrals_inc_hydrogen
                if 1 in these_atomic_numbers
                else dihedrals_without_hydrogen
            )

            dihedrals_list.append(atom1_index * 3)
            dihedrals_list.append(atom2_index * 3)
            dihedrals_list.append(atom3_index * 3 * -1)
            dihedrals_list.append(atom4_index * 3 * -1)
            dihedrals_list.append(dihedral_type_index + 1)

    return dihedrals_inc_hydrogen, dihedrals_without_hydrogen


# TODO: Split this mono-function into smaller functions
def to_prmtop(interchange: "Interchange", file_path: Union[Path, str]):
    """
    Write a .prmtop file. See http://ambermd.org/prmtop.pdf for details.

    """
    if isinstance(file_path, str):
        path = Path(file_path)
    if isinstance(file_path, Path):
        path = file_path

    if interchange["vdW"].mixing_rule != "lorentz-berthelot":
        raise UnsupportedMixingRuleError(interchange["vdW"].mixing_rule)

    if interchange.box is None:
        if interchange["Electrostatics"].periodic_potential != _PME:
            raise UnsupportedExportError(
                f'Electrostatics method PME (`"{_PME}"`) is not valid for a non-periodic system. ',
            )

    with open(path, "w") as prmtop:
        import datetime

        now = datetime.datetime.now()
        prmtop.write(
            "%VERSION  VERSION_STAMP = V0001.000  DATE = "
            f"{now.month:02d}/{now.day:02d}/{(now.year % 100):02}  "
            f"{now.hour:02d}:{now.minute:02d}:{now.second:02d}\n"
            "%FLAG TITLE\n"
            "%FORMAT(20a4)\n"
            "\n",
        )

        from openff.interchange.interop.common import _build_typemap

        typemap = _build_typemap(interchange)  # noqa

        potential_key_to_atom_type_mapping: dict[PotentialKey, int] = {
            key: i for i, key in enumerate(interchange["vdW"].potentials)
        }
        atom_type_indices = [
            potential_key_to_atom_type_mapping[potential_key]
            for potential_key in interchange["vdW"].key_map.values()
        ]

        potential_key_to_bond_type_mapping: dict[PotentialKey, int] = {
            key: i for i, key in enumerate(interchange["Bonds"].potentials)
        }

        potential_key_to_angle_type_mapping: dict[PotentialKey, int] = {
            key: i for i, key in enumerate(interchange["Angles"].potentials)
        }

        dihedral_potentials = dict()
        for key in ["ProperTorsions", "ImproperTorsions"]:
            if key in interchange.collections:
                dihedral_potentials.update(deepcopy(interchange[key].potentials))
                dihedral_potentials.update(deepcopy(interchange[key].potentials))

        potential_key_to_dihedral_type_mapping: dict[PotentialKey, int] = {
            key: i for i, key in enumerate(dihedral_potentials)
        }

        # This is only used to track 1-4 pairs that are already counted. Since we
        # only care about membership and never indexing and the membership check
        # will happen roughly once per dihedral, use a set for speed.
        already_counted: set[tuple[int, ...]] = set()

        per_atom_exclusion_lists = _get_per_atom_exclusion_lists(interchange.topology)

        number_excluded_atoms, excluded_atoms_list = _get_exclusion_lists(
            per_atom_exclusion_lists,
        )
        atomic_numbers = tuple(
            atom.atomic_number for atom in interchange.topology.atoms
        )

        bonds_inc_hydrogen, bonds_without_hydrogen = _get_bond_lists(
            interchange,
            atomic_numbers,
            potential_key_to_bond_type_mapping,
        )

        angles_inc_hydrogen, angles_without_hydrogen = _get_angle_lists(
            interchange,
            atomic_numbers,
            potential_key_to_angle_type_mapping,
        )

        dihedrals_inc_hydrogen, dihedrals_without_hydrogen = _get_dihedral_lists(
            interchange,
            atomic_numbers,
            potential_key_to_dihedral_type_mapping,
            already_counted,
            per_atom_exclusion_lists,
        )

        # total number of atoms
        NATOM = interchange.topology.n_atoms
        # total number of distinct atom types
        NTYPES = len(interchange["vdW"].potentials)
        # number of bonds containing hydrogen
        NBONH = _get_num_h_bonds(interchange.topology)
        # number of bonds not containing hydrogen
        MBONA = interchange.topology.n_bonds - NBONH
        # number of angles containing hydrogen
        NTHETH = int(len(angles_inc_hydrogen) / 4)
        # number of angles not containing hydrogen
        MTHETA = int(len(angles_without_hydrogen) / 4)
        # number of dihedrals containing hydrogen
        NPHIH = int(len(dihedrals_inc_hydrogen) / 5)
        # number of dihedrals not containing hydrogen
        MPHIA = int(len(dihedrals_without_hydrogen) / 5)
        NHPARM = 0  # : currently not used
        NPARM = 0  # : used to determine if addles created prmtop
        # number of excluded atoms
        NNB = len(excluded_atoms_list)
        # number of residues
        NRES = max(1, len([*interchange.topology.hierarchy_iterator("residues")]))
        NBONA = MBONA  # : MBONA + number of constraint bonds
        NTHETA = MTHETA  # : MTHETA + number of constraint angles
        NPHIA = MPHIA  # : MPHIA + number of constraint dihedrals
        # number of unique bond types
        NUMBND = len(potential_key_to_bond_type_mapping)
        # number of unique angle types
        NUMANG = len(potential_key_to_angle_type_mapping)
        # number of unique dihedral types
        NPTRA = len(potential_key_to_dihedral_type_mapping)
        # number of atom types in parameter file, see SOLTY below
        # this appears to be unused, but ParmEd writes a 1 here (?)
        NATYP = 1
        NPHB = 0  # : number of distinct 10-12 hydrogen bond pair types
        IFPERT = 0  # : set to 1 if perturbation info is to be read in
        NBPER = 0  # : number of bonds to be perturbed
        NGPER = 0  # : number of angles to be perturbed
        NDPER = 0  # : number of dihedrals to be perturbed
        MBPER = 0  # : number of bonds with atoms completely in perturbed group
        MGPER = 0  # : number of angles with atoms completely in perturbed group
        MDPER = 0  # : number of dihedrals with atoms completely in perturbed groups
        # set to 1 if standard periodic box, 2 when truncated octahedral
        IFBOX = 0 if interchange.box is None else 1
        # number of atoms in the largest residue
        NMXRS = interchange.topology.n_atoms
        IFCAP = 0  # : set to 1 if the CAP option from edit was specified
        NUMEXTRA = 0  # : number of extra points found in topology
        # number of PIMD slices / number of beads
        # ParmEd does not seem to write this _at all_ in most/all cases
        NCOPY = 0

        pointers = [
            NATOM,
            NTYPES,
            NBONH,
            MBONA,
            NTHETH,
            MTHETA,
            NPHIH,
            MPHIA,
            NHPARM,
            NPARM,
            NNB,
            NRES,
            NBONA,
            NTHETA,
            NPHIA,
            NUMBND,
            NUMANG,
            NPTRA,
            NATYP,
            NPHB,
            IFPERT,
            NBPER,
            NGPER,
            NDPER,
            MBPER,
            MGPER,
            MDPER,
            IFBOX,
            NMXRS,
            IFCAP,
            NUMEXTRA,
            NCOPY,
        ]

        prmtop.write("%FLAG POINTERS\n" "%FORMAT(10I8)\n")
        text_blob = "".join([str(val).rjust(8) for val in pointers])
        _write_text_blob(prmtop, text_blob)

        prmtop.write("%FLAG ATOM_NAME\n" "%FORMAT(20a4)\n")

        atom_names: list[str] = list()

        for atom in interchange.topology.atoms:
            if hasattr(atom, "name"):
                if atom.name:
                    atom_names.append(atom.name)
                else:
                    atom_names.append(atom.symbol)
            else:
                atom_names.append(atom.symbol)

        text_blob = "".join([str(val).rjust(4) for val in atom_names])
        _write_text_blob(prmtop, text_blob)

        prmtop.write("%FLAG CHARGE\n" "%FORMAT(5E16.8)\n")
        charges = [
            charge.m_as(unit.e) * AMBER_COULOMBS_CONSTANT
            for charge in interchange["Electrostatics"].charges.values()
        ]
        text_blob = "".join([f"{val:16.8E}" for val in charges])
        _write_text_blob(prmtop, text_blob)

        prmtop.write("%FLAG ATOMIC_NUMBER\n" "%FORMAT(10I8)\n")
        text_blob = "".join([str(val).rjust(8) for val in atomic_numbers])
        _write_text_blob(prmtop, text_blob)

        prmtop.write("%FLAG MASS\n" "%FORMAT(5E16.8)\n")
        masses = [a.mass.m for a in interchange.topology.atoms]
        text_blob = "".join([f"{val:16.8E}" for val in masses])
        _write_text_blob(prmtop, text_blob)

        prmtop.write("%FLAG ATOM_TYPE_INDEX\n" "%FORMAT(10I8)\n")
        text_blob = "".join([str(val + 1).rjust(8) for val in atom_type_indices])
        _write_text_blob(prmtop, text_blob)

        prmtop.write("%FLAG NUMBER_EXCLUDED_ATOMS\n" "%FORMAT(10I8)\n")
        # https://ambermd.org/prmtop.pdf says this section is ignored (!?)
        text_blob = "".join([str(val).rjust(8) for val in number_excluded_atoms])
        _write_text_blob(prmtop, text_blob)

        acoefs = [None] * int((NTYPES + 1) * NTYPES / 2)
        bcoefs = [None] * int((NTYPES + 1) * NTYPES / 2)

        nonbonded_parm_indices: list[Optional[int]] = [None] * (NTYPES * NTYPES)

        for key_i, i in potential_key_to_atom_type_mapping.items():
            for key_j, j in potential_key_to_atom_type_mapping.items():
                if j < i:
                    # Only need to handle the lower triangle as everything symmetric.
                    continue

                atom_type_index_i = i + 1
                atom_type_index_j = j + 1

                coeff_index = int(i + (j + 1) * j / 2) + 1  # FORTRAN IDX

                # index = NONBONDED PARM INDEX [NTYPES × (ATOM TYPE INDEX(i) − 1) + ATOM TYPE INDEX(j)]
                parm_index_fwd = (
                    NTYPES * (atom_type_index_i - 1) + atom_type_index_j  # FORTRAN IDX
                )
                parm_index_rev = (
                    NTYPES * (atom_type_index_j - 1) + atom_type_index_i  # FORTRAN IDX
                )

                # TODO: Figure out the right way to map cross-interactions, using the
                #       key_i and key_j objects as lookups to parameters
                sigma_i = interchange["vdW"].potentials[key_i].parameters["sigma"]
                sigma_j = interchange["vdW"].potentials[key_j].parameters["sigma"]
                epsilon_i = interchange["vdW"].potentials[key_i].parameters["epsilon"]
                epsilon_j = interchange["vdW"].potentials[key_j].parameters["epsilon"]

                sigma = (sigma_i + sigma_j) * 0.5
                epsilon = (epsilon_i * epsilon_j) ** 0.5

                acoef = (4 * epsilon * sigma**12).m_as(kcal_mol * unit.angstrom**12)
                bcoef = (4 * epsilon * sigma**6).m_as(kcal_mol * unit.angstrom**6)

                acoefs[coeff_index - 1] = acoef
                bcoefs[coeff_index - 1] = bcoef

                nonbonded_parm_indices[parm_index_fwd - 1] = coeff_index
                nonbonded_parm_indices[parm_index_rev - 1] = coeff_index

        assert all(
            value is not None
            for values in [acoefs, bcoefs, nonbonded_parm_indices]
            for value in values  # type: ignore
        ), "an internal error occurred"

        prmtop.write("%FLAG NONBONDED_PARM_INDEX\n" "%FORMAT(10I8)\n")
        text_blob = "".join([str(val).rjust(8) for val in nonbonded_parm_indices])
        _write_text_blob(prmtop, text_blob)

        residue_names = [
            getattr(residue, "residue_name", "RES")
            for residue in interchange.topology.hierarchy_iterator("residues")
        ]

        # Avoid blank residue labels, though this is a band-aid. See issue #652
        if residue_names == list():
            residue_names = ["RES"]

        prmtop.write("%FLAG RESIDUE_LABEL\n" "%FORMAT(20a4)\n")
        text_blob = "".join([val.ljust(4) for val in residue_names])
        _write_text_blob(prmtop, text_blob)

        residue_pointers = (
            [
                interchange.topology.atom_index([*residue.atoms][0])
                for residue in interchange.topology.hierarchy_iterator("residues")
            ]
            if NRES > 1
            else [0]
        )
        prmtop.write("%FLAG RESIDUE_POINTER\n" "%FORMAT(10I8)\n")
        text_blob = "".join([str(val + 1).rjust(8) for val in residue_pointers])
        _write_text_blob(prmtop, text_blob)

        # TODO: Exclude (?) bonds containing hydrogens
        prmtop.write("%FLAG BOND_FORCE_CONSTANT\n" "%FORMAT(5E16.8)\n")
        bond_k = [
            interchange["Bonds"].potentials[key].parameters["k"].m_as(kcal_mol_a2) / 2
            for key in potential_key_to_bond_type_mapping
        ]
        text_blob = "".join([f"{val:16.8E}" for val in bond_k])
        _write_text_blob(prmtop, text_blob)

        prmtop.write("%FLAG BOND_EQUIL_VALUE\n" "%FORMAT(5E16.8)\n")
        bond_length = [
            interchange["Bonds"]
            .potentials[key]
            .parameters["length"]
            .m_as(unit.angstrom)
            for key in potential_key_to_bond_type_mapping
        ]
        text_blob = "".join([f"{val:16.8E}" for val in bond_length])
        _write_text_blob(prmtop, text_blob)

        prmtop.write("%FLAG ANGLE_FORCE_CONSTANT\n" "%FORMAT(5E16.8)\n")
        angle_k = [
            interchange["Angles"].potentials[key].parameters["k"].m_as(kcal_mol_rad2)
            / 2  # noqa
            for key in potential_key_to_angle_type_mapping
        ]
        text_blob = "".join([f"{val:16.8E}" for val in angle_k])
        _write_text_blob(prmtop, text_blob)

        prmtop.write("%FLAG ANGLE_EQUIL_VALUE\n" "%FORMAT(5E16.8)\n")
        angle_theta = [
            interchange["Angles"].potentials[key].parameters["angle"].m_as(unit.radian)
            for key in potential_key_to_angle_type_mapping
        ]
        text_blob = "".join([f"{val:16.8E}" for val in angle_theta])
        _write_text_blob(prmtop, text_blob)

        dihedral_k: list[int] = list()
        dihedral_periodicity: list[int] = list()
        dihedral_phase: list[int] = list()

        for key_ in potential_key_to_dihedral_type_mapping:
            params = interchange[key_.associated_handler].potentials[key_].parameters  # type: ignore
            idivf = int(params["idivf"]) if "idivf" in params else 1
            dihedral_k.append((params["k"] / idivf).m_as(kcal_mol))
            dihedral_periodicity.append(params["periodicity"].m_as(unit.dimensionless))
            dihedral_phase.append(params["phase"].m_as(unit.radian))

        prmtop.write("%FLAG DIHEDRAL_FORCE_CONSTANT\n" "%FORMAT(5E16.8)\n")
        text_blob = "".join([f"{val:16.8E}" for val in dihedral_k])
        _write_text_blob(prmtop, text_blob)

        prmtop.write("%FLAG DIHEDRAL_PERIODICITY\n" "%FORMAT(5E16.8)\n")
        text_blob = "".join([f"{val:16.8E}" for val in dihedral_periodicity])
        _write_text_blob(prmtop, text_blob)

        prmtop.write("%FLAG DIHEDRAL_PHASE\n" "%FORMAT(5E16.8)\n")
        text_blob = "".join([f"{val:16.8E}" for val in dihedral_phase])
        _write_text_blob(prmtop, text_blob)

        prmtop.write("%FLAG SCEE_SCALE_FACTOR\n" "%FORMAT(5E16.8)\n")
        scee = NPTRA * [1.2]
        text_blob = "".join([f"{val:16.8E}" for val in scee])
        _write_text_blob(prmtop, text_blob)

        prmtop.write("%FLAG SCNB_SCALE_FACTOR\n" "%FORMAT(5E16.8)\n")
        scnb = NPTRA * [2.0]
        text_blob = "".join([f"{val:16.8E}" for val in scnb])
        _write_text_blob(prmtop, text_blob)

        prmtop.write("%FLAG SOLTY\n" "%FORMAT(5E16.8)\n")
        prmtop.write(f"{0:16.8E}\n")

        prmtop.write("%FLAG LENNARD_JONES_ACOEF\n" "%FORMAT(5E16.8)\n")
        text_blob = "".join([f"{val:16.8E}" for val in acoefs])
        _write_text_blob(prmtop, text_blob)

        prmtop.write("%FLAG LENNARD_JONES_BCOEF\n" "%FORMAT(5E16.8)\n")
        text_blob = "".join([f"{val:16.8E}" for val in bcoefs])
        _write_text_blob(prmtop, text_blob)

        prmtop.write("%FLAG BONDS_INC_HYDROGEN\n" "%FORMAT(10I8)\n")
        text_blob = "".join([str(val).rjust(8) for val in bonds_inc_hydrogen])
        _write_text_blob(prmtop, text_blob)

        prmtop.write("%FLAG BONDS_WITHOUT_HYDROGEN\n" "%FORMAT(10I8)\n")
        text_blob = "".join([str(val).rjust(8) for val in bonds_without_hydrogen])
        _write_text_blob(prmtop, text_blob)

        prmtop.write("%FLAG ANGLES_INC_HYDROGEN\n" "%FORMAT(10I8)\n")
        text_blob = "".join([str(val).rjust(8) for val in angles_inc_hydrogen])
        _write_text_blob(prmtop, text_blob)

        prmtop.write("%FLAG ANGLES_WITHOUT_HYDROGEN\n" "%FORMAT(10I8)\n")
        text_blob = "".join([str(val).rjust(8) for val in angles_without_hydrogen])
        _write_text_blob(prmtop, text_blob)

        prmtop.write("%FLAG DIHEDRALS_INC_HYDROGEN\n" "%FORMAT(10I8)\n")
        text_blob = "".join([str(val).rjust(8) for val in dihedrals_inc_hydrogen])
        _write_text_blob(prmtop, text_blob)

        prmtop.write("%FLAG DIHEDRALS_WITHOUT_HYDROGEN\n" "%FORMAT(10I8)\n")
        text_blob = "".join([str(val).rjust(8) for val in dihedrals_without_hydrogen])
        _write_text_blob(prmtop, text_blob)

        prmtop.write("%FLAG EXCLUDED_ATOMS_LIST\n" "%FORMAT(10I8)\n")
        text_blob = "".join([str(val).rjust(8) for val in excluded_atoms_list])
        _write_text_blob(prmtop, text_blob)

        prmtop.write("%FLAG HBOND_ACOEF\n" "%FORMAT(5E16.8)\n")
        _write_text_blob(prmtop, "")

        prmtop.write("%FLAG HBOND_BCOEF\n" "%FORMAT(5E16.8)\n")
        _write_text_blob(prmtop, "")

        prmtop.write("%FLAG HBCUT\n" "%FORMAT(5E16.8)\n")
        _write_text_blob(prmtop, "")

        prmtop.write("%FLAG AMBER_ATOM_TYPE\n" "%FORMAT(20a4)\n")
        text_blob = "".join([val.ljust(4) for val in typemap.values()])
        _write_text_blob(prmtop, text_blob)

        prmtop.write("%FLAG TREE_CHAIN_CLASSIFICATION\n" "%FORMAT(20a4)\n")
        blahs = NATOM * ["BLA"]
        text_blob = "".join([val.ljust(4) for val in blahs])
        _write_text_blob(prmtop, text_blob)

        prmtop.write("%FLAG JOIN_ARRAY\n" "%FORMAT(10I8)\n")
        _ = NATOM * [0]
        text_blob = "".join([str(val).rjust(8) for val in _])
        _write_text_blob(prmtop, text_blob)

        prmtop.write("%FLAG IROTAT\n" "%FORMAT(10I8)\n")
        _ = NATOM * [0]
        text_blob = "".join([str(val).rjust(8) for val in _])
        _write_text_blob(prmtop, text_blob)

        if IFBOX == 1:
            prmtop.write("%FLAG SOLVENT_POINTERS\n" "%FORMAT(3I8)\n")
            prmtop.write("       1       1       2\n")

            # TODO: No easy way to accurately export this section while
            #       using an MDTraj topology
            prmtop.write("%FLAG ATOMS_PER_MOLECULE\n" "%FORMAT(10I8)\n")
            prmtop.write(str(interchange.topology.n_atoms).rjust(8))
            prmtop.write("\n")

            prmtop.write("%FLAG BOX_DIMENSIONS\n" "%FORMAT(5E16.8)\n")
            box = [90.0]
            for i in range(3):
                box.append(interchange.box[i, i].m_as(unit.angstrom))  # type: ignore
            text_blob = "".join([f"{val:16.8E}" for val in box])
            _write_text_blob(prmtop, text_blob)

        prmtop.write("%FLAG RADIUS_SET\n" "%FORMAT(1a80)\n")
        prmtop.write("0\n")

        prmtop.write("%FLAG RADII\n" "%FORMAT(5E16.8)\n")
        radii = NATOM * [0]
        text_blob = "".join([f"{val:16.8E}" for val in radii])
        _write_text_blob(prmtop, text_blob)

        prmtop.write("%FLAG SCREEN\n" "%FORMAT(5E16.8)\n")
        screen = NATOM * [0]
        text_blob = "".join([f"{val:16.8E}" for val in screen])
        _write_text_blob(prmtop, text_blob)

        prmtop.write("%FLAG IPOL\n" "%FORMAT(1I8)\n")
        prmtop.write("       0\n")


def to_inpcrd(interchange: "Interchange", file_path: Union[Path, str]):
    """
    Write a .prmtop file. See https://ambermd.org/FileFormats.php#restart for details.

    """
    if isinstance(file_path, str):
        path = Path(file_path)
    if isinstance(file_path, Path):
        path = file_path

    n_atoms = interchange.topology.n_atoms
    time = 0.0

    with open(path, "w") as inpcrd:
        inpcrd.write(f"\n{n_atoms:5d}{time:15.7e}\n")

        coords = interchange.positions.m_as(unit.angstrom)
        blob = "".join([f"{val:12.7f}".rjust(12) for val in coords.flatten()])

        for line in textwrap.wrap(blob, width=72, drop_whitespace=False):
            inpcrd.write(line + "\n")

        if interchange.box is not None:
            box = interchange.box.to(unit.angstrom).magnitude
            if (box == np.diag(np.diagonal(box))).all():
                for i in range(3):
                    inpcrd.write(f"{box[i, i]:12.7f}")
                for _ in range(3):
                    inpcrd.write("  90.0000000")
            else:
                # TODO: Handle non-rectangular
                raise NotImplementedError

        inpcrd.write("\n")
