import pathlib

import numpy
from openff.units import unit

from openff.interchange._experimental import experimental
from openff.interchange.interop.gromacs.models.models import (
    GROMACSAngle,
    GROMACSAtom,
    GROMACSAtomType,
    GROMACSBond,
    GROMACSDihedral,
    GROMACSExclusion,
    GROMACSMolecule,
    GROMACSPair,
    GROMACSSettles,
    GROMACSSystem,
    LennardJonesAtomType,
    PeriodicImproperDihedral,
    PeriodicProperDihedral,
    RyckaertBellemansDihedral,
)


@experimental
def from_files(top_file, gro_file, cls=GROMACSSystem) -> GROMACSSystem:
    """
    Parse a GROMACS topology file. Adapted from Intermol.

    https://github.com/shirtsgroup/InterMol/blob/v0.1.2/intermol/gromacs/gromacs_parser.py
    """
    with open(top_file) as f:
        for line in f:
            stripped = line.split(";")[0].strip()

            if len(stripped) == 0:
                continue

            if stripped.startswith(";"):
                continue

            if stripped.startswith("["):
                if not len(stripped.split()) == 3 and stripped.endswith("]"):
                    raise ValueError("Invalid GROMACS topology file")

                current_directive = stripped[1:-1].strip()

                continue

            if current_directive == "defaults":
                (
                    nonbonded_function,
                    combination_rule,
                    gen_pairs,
                    vdw_14,
                    coul_14,
                ) = _process_defaults(line)

                system = cls(
                    nonbonded_function=nonbonded_function,
                    combination_rule=combination_rule,
                    gen_pairs=gen_pairs,
                    vdw_14=vdw_14,
                    coul_14=coul_14,
                )

            elif current_directive == "atomtypes":
                atom_type = _process_atomtype(line)
                system.atom_types[atom_type.name] = atom_type

            elif current_directive == "moleculetype":
                molecule_type = _process_moleculetype(line)
                system.molecule_types[molecule_type.name] = molecule_type

                current_molecule = molecule_type.name

            elif current_directive == "atoms":
                system.molecule_types[current_molecule].atoms.append(
                    _process_atom(line),
                )

            elif current_directive == "pairs":
                system.molecule_types[current_molecule].pairs.append(
                    _process_pair(line),
                )

            elif current_directive == "settles":
                system.molecule_types[current_molecule].settles.append(
                    _process_settles(line),
                )

            elif current_directive == "bonds":
                system.molecule_types[current_molecule].bonds.append(
                    _process_bond(line),
                )

            elif current_directive == "angles":
                system.molecule_types[current_molecule].angles.append(
                    _process_angle(line),
                )

            elif current_directive == "dihedrals":
                system.molecule_types[current_molecule].dihedrals.append(
                    _process_dihedral(line),
                )

            elif current_directive in ["exclusions"]:
                system.molecule_types[current_molecule].exclusions.append(
                    _process_exclusion(line),
                )

            elif current_directive == "system":
                system.name = _process_system(line)

            elif current_directive == "molecules":
                molecule_name, number_of_copies = _process_molecule(line)

                system.molecules[molecule_name] = number_of_copies

            else:
                raise ValueError(f"Invalid directive {current_directive}")

    for molecule_type in system.molecule_types.values():
        this_molecule_atom_type_names = tuple(
            atom.atom_type for atom in molecule_type.atoms
        )

        molecule_type._contained_atom_types = {
            atom_type_name: system.atom_types[atom_type_name]
            for atom_type_name in this_molecule_atom_type_names
        }

    system.positions = _read_coordinates(gro_file)
    system.box = _read_box(gro_file)

    return system


def _process_defaults(line: str) -> tuple[int, int, str, float, float]:
    split = line.split()

    nonbonded_function = int(split[0])

    if nonbonded_function != 1:
        raise ValueError("Only LJ nonbonded functions are supported.")

    combination_rule = int(split[1])

    if combination_rule not in (2, 3):
        raise ValueError("combination rule 1 not supported.")

    gen_pairs = split[2]
    lj_14 = float(split[3])
    coul_14 = float(split[4])

    return nonbonded_function, combination_rule, gen_pairs, lj_14, coul_14


def _process_atomtype(
    line: str,
) -> GROMACSAtomType:
    split = line.split()

    # The second (bonding type) and third (atomic number) columns are optional. Handle these cases
    # explicitly by assuming that the ptype column is always a 1-length string.
    # See https://github.com/shirtsgroup/InterMol/blob/v0.1.2/intermol/gromacs/gromacs_parser.py#L1357-L1368

    if len(split[3]) == 1:
        # Bonded type and atomic number are both missing.
        split.insert(1, None)  # type: ignore[arg-type]
        split.insert(1, None)  # type: ignore[arg-type]

    elif len(split[4]) == 1 and len(split[5]) >= 1:
        if split[1][0].isalpha():
            # Atomic number is missing.
            split.insert(2, None)  # type: ignore[arg-type]
        else:
            # Bonded type is missing.
            split.insert(1, None)  # type: ignore[arg-type]

    atom_type = split[0]
    bonding_type = split[1]

    atomic_number = int(split[2]) if split[2] is not None else None
    mass = unit.Quantity(float(split[3]), unit.dalton)

    charge = unit.Quantity(float(split[4]), unit.elementary_charge)

    particle_type = split[5]

    if particle_type == "A":
        sigma = unit.Quantity(float(split[6]), unit.nanometer)
        epsilon = unit.Quantity(float(split[7]), unit.kilojoule_per_mole)
    else:
        raise ValueError(f"Particle type must be A, parsed {particle_type}.")

    return LennardJonesAtomType(
        name=atom_type,
        bonding_type=bonding_type,
        atomic_number=atomic_number,
        mass=mass,
        charge=charge,
        particle_type=particle_type,
        sigma=sigma,
        epsilon=epsilon,
    )


def _process_moleculetype(line: str) -> GROMACSMolecule:
    split = line.split()

    molecule_type = split[0]
    nrexcl = int(split[1])

    return GROMACSMolecule(name=molecule_type, nrexcl=nrexcl)


def _process_atom(
    line: str,
) -> GROMACSAtom:
    split = line.split()

    atom_number = int(split[0])
    atom_type = split[1]
    residue_number = int(split[2])
    residue_name = split[3]
    atom_name = split[4]
    charge_group_number = int(split[5])
    charge = unit.Quantity(float(split[6]), unit.elementary_charge)
    mass = unit.Quantity(float(split[7]), unit.amu)

    return GROMACSAtom(
        index=atom_number,
        atom_type=atom_type,
        name=atom_name,
        residue_index=residue_number,
        residue_name=residue_name,
        charge_group_number=charge_group_number,
        charge=charge,
        mass=mass,
    )


def _process_pair(line: str) -> GROMACSPair:
    split = line.split()

    atom1 = int(split[0])
    atom2 = int(split[1])

    nonbonded_function = int(split[2])

    if nonbonded_function == 1:
        return GROMACSPair(
            atom1=atom1,
            atom2=atom2,
        )
    else:
        raise ValueError(
            f"Nonbonded function must be 1, parsed a pair with {nonbonded_function}.",
        )


def _process_settles(line: str) -> GROMACSSettles:
    split = line.split()

    first_atom = int(split[0])

    oxygen_hydrogen_distance = unit.Quantity(float(split[2]), unit.nanometer)
    hydrogen_hydrogen_distance = unit.Quantity(float(split[3]), unit.nanometer)

    return GROMACSSettles(
        first_atom=first_atom,
        oxygen_hydrogen_distance=oxygen_hydrogen_distance,
        hydrogen_hydrogen_distance=hydrogen_hydrogen_distance,
    )


def _process_bond(line: str) -> GROMACSBond:
    split = line.split()

    atom1 = int(split[0])
    atom2 = int(split[1])

    bond_function = int(split[2])

    if bond_function == 1:
        bond_length = unit.Quantity(float(split[3]), unit.nanometer)
        bond_k = unit.Quantity(
            float(split[4]),
            unit.kilojoule_per_mole / unit.nanometer**2,
        )

        return GROMACSBond(
            atom1=atom1,
            atom2=atom2,
            function=bond_function,
            length=bond_length,
            k=bond_k,
        )

    else:
        raise ValueError(f"Bond function must be 1, parsed {bond_function}.")


def _process_angle(
    line: str,
) -> GROMACSAngle:
    split = line.split()

    atom1 = int(split[0])
    atom2 = int(split[1])
    atom3 = int(split[2])

    angle_function = int(split[3])

    if angle_function == 1:
        angle = unit.Quantity(float(split[4]), unit.degrees)
        k = unit.Quantity(float(split[5]), unit.kilojoule_per_mole)
    else:
        raise ValueError(f"Angle function must be 1, parsed {angle_function}.")

    return GROMACSAngle(
        atom1=atom1,
        atom2=atom2,
        atom3=atom3,
        angle=angle,
        k=k,
    )


def _process_dihedral(
    line: str,
) -> GROMACSDihedral:
    split = line.split()

    atom1 = int(split[0])
    atom2 = int(split[1])
    atom3 = int(split[2])
    atom4 = int(split[3])

    dihedral_function = int(split[4])

    if dihedral_function == 1:
        return PeriodicProperDihedral(
            atom1=atom1,
            atom2=atom2,
            atom3=atom3,
            atom4=atom4,
            phi=unit.Quantity(float(split[5]), unit.degrees),
            k=unit.Quantity(float(split[6]), unit.kilojoule_per_mole),
            multiplicity=int(float(split[7])),
        )

    elif dihedral_function == 3:
        return RyckaertBellemansDihedral(
            atom1=atom1,
            atom2=atom2,
            atom3=atom3,
            atom4=atom4,
            c0=unit.Quantity(float(split[5]), unit.kilojoule_per_mole),
            c1=unit.Quantity(float(split[6]), unit.kilojoule_per_mole),
            c2=unit.Quantity(float(split[7]), unit.kilojoule_per_mole),
            c3=unit.Quantity(float(split[8]), unit.kilojoule_per_mole),
            c4=unit.Quantity(float(split[9]), unit.kilojoule_per_mole),
            c5=unit.Quantity(float(split[10]), unit.kilojoule_per_mole),
        )

    elif dihedral_function == 4:
        return PeriodicImproperDihedral(
            atom1=atom1,
            atom2=atom2,
            atom3=atom3,
            atom4=atom4,
            phi=unit.Quantity(float(split[5]), unit.degrees),
            k=unit.Quantity(float(split[6]), unit.kilojoule_per_mole),
            multiplicity=int(float(split[7])),
        )

    else:
        raise ValueError(
            f"Dihedral functions 1, 3, and 4 supported, parsed {dihedral_function}.",
        )


def _process_exclusion(line: str) -> GROMACSExclusion:
    split = line.split()

    return GROMACSExclusion(
        first_atom=int(split[0]),
        other_atoms=[int(atom) for atom in split[1:]],
    )


def _process_molecule(line: str) -> tuple[str, int]:
    split = line.split()

    molecule_name = split[0]
    number_of_molecules = int(split[1])

    return molecule_name, number_of_molecules


def _process_system(line: str) -> str:
    split = line.split()

    system_name = split[0]

    return system_name


def _read_coordinates(file_path: pathlib.Path) -> unit.Quantity:
    def _infer_coord_precision(file_path: pathlib.Path) -> int:
        """
        Infer decimal precision of coordinates by parsing periods in atoms lines.
        """
        with open(file_path) as file_in:
            file_in.readline()
            file_in.readline()
            atom_line = file_in.readline()
            period_indices = [i for i, x in enumerate(atom_line) if x == "."]
            spacing_between_periods = period_indices[-1] - period_indices[-2]
            precision = spacing_between_periods - 5
            return precision

    precision = _infer_coord_precision(file_path)
    coordinate_width = precision + 5
    # Column numbers in file separating x, y, z coords of each atom.
    # Default (3 decimals of precision -> 8 columns) are 20, 28, 36, 44
    coordinate_columns = [
        20,
        20 + coordinate_width,
        20 + 2 * coordinate_width,
        20 + 3 * coordinate_width,
    ]

    with open(file_path) as gro_file:
        # Throw away comment / name line
        gro_file.readline()
        n_atoms = int(gro_file.readline())

        unitless_coordinates = numpy.zeros((n_atoms, 3))
        for coordinate_index in range(n_atoms):
            line = gro_file.readline()
            _ = int(line[:5])  # residue_index
            _ = line[5:10]  # residue_name
            _ = line[10:15]  # atom_name
            _ = int(line[15:20])  # atom_index
            x = float(line[coordinate_columns[0] : coordinate_columns[1]])
            y = float(line[coordinate_columns[1] : coordinate_columns[2]])
            z = float(line[coordinate_columns[2] : coordinate_columns[3]])
            unitless_coordinates[coordinate_index] = numpy.array([x, y, z])

        coordinates = unitless_coordinates * unit.nanometer

    return coordinates


def _read_box(file_path: pathlib.Path) -> unit.Quantity:
    with open(file_path) as gro_file:
        # Throw away comment / name line
        gro_file.readline()
        n_atoms = int(gro_file.readline())

        box_line = gro_file.readlines()[n_atoms]

    parsed_box = [float(val) for val in box_line.split()]

    if len(parsed_box) == 3:
        box = parsed_box * numpy.eye(3) * unit.nanometer

    return box
