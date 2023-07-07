"""Interfaces with LAMMPS."""
from pathlib import Path
from typing import IO, Union

import numpy as np
from openff.units import unit

from openff.interchange import Interchange
from openff.interchange.exceptions import UnsupportedExportError
from openff.interchange.models import AngleKey, BondKey


def to_lammps(interchange: Interchange, file_path: Union[Path, str]):
    """Write an Interchange object to a LAMMPS data file."""
    if isinstance(file_path, str):
        path = Path(file_path)
    if isinstance(file_path, Path):
        path = file_path

    n_atoms = interchange.topology.n_atoms
    if "Bonds" in interchange.collections:
        n_bonds = len(interchange["Bonds"].key_map.keys())
    else:
        n_bonds = 0
    if "Angles" in interchange.collections:
        n_angles = len(interchange["Angles"].key_map.keys())
    else:
        n_angles = 0
    if "ProperTorsions" in interchange.collections:
        n_propers = len(interchange["ProperTorsions"].key_map.keys())
    else:
        n_propers = 0
    if "ImproperTorsions" in interchange.collections:
        n_impropers = len(interchange["ImproperTorsions"].key_map.keys())
    else:
        n_impropers = 0

    with open(path, "w") as lmp_file:
        lmp_file.write("Title\n\n")

        lmp_file.write(f"{n_atoms} atoms\n")
        lmp_file.write(f"{n_bonds} bonds\n")
        lmp_file.write(f"{n_angles} angles\n")
        lmp_file.write(f"{n_propers} dihedrals\n")
        lmp_file.write(f"{n_impropers} impropers\n")

        lmp_file.write(f"\n{len(interchange['vdW'].potentials)} atom types")
        if n_bonds > 0:
            lmp_file.write(f"\n{len(interchange['Bonds'].potentials)} bond types")
        if n_angles > 0:
            lmp_file.write(f"\n{len(interchange['Angles'].potentials)} angle types")
        if n_propers > 0:
            lmp_file.write(
                f"\n{len(interchange['ProperTorsions'].potentials)} dihedral types",
            )
        if n_impropers > 0:
            lmp_file.write(
                f"\n{len(interchange['ImproperTorsions'].potentials)} improper types",
            )

        lmp_file.write("\n")

        # write types section

        x_min, y_min, z_min = np.min(
            interchange.positions.to(unit.angstrom),
            axis=0,
        ).magnitude
        if interchange.box is None:
            L_x, L_y, L_z = 100, 100, 100
        else:
            L_x, L_y, L_z = np.diag(interchange.box.to(unit.angstrom).magnitude)

        lmp_file.write(
            "{:.10g} {:.10g} xlo xhi\n"
            "{:.10g} {:.10g} ylo yhi\n"
            "{:.10g} {:.10g} zlo zhi\n".format(
                x_min,
                x_min + L_x,
                y_min,
                y_min + L_y,
                z_min,
                z_min + L_z,
            ),
        )

        lmp_file.write("0.0 0.0 0.0 xy xz yz\n")

        lmp_file.write("\nMasses\n\n")

        vdw_handler = interchange["vdW"]
        atom_type_map = dict(enumerate(vdw_handler.potentials))
        key_map_inv = dict({v: k for k, v in vdw_handler.key_map.items()})

        for atom_type_idx, smirks in atom_type_map.items():
            # Find just one topology atom matching this SMIRKS by vdW
            matched_atom_idx = key_map_inv[smirks].atom_indices[0]
            matched_atom = interchange.topology.atom(matched_atom_idx)
            mass = matched_atom.mass.m

            lmp_file.write(f"{atom_type_idx + 1:d}\t{mass:.8g}\n")

        lmp_file.write("\n\n")

        _write_pair_coeffs(
            lmp_file=lmp_file,
            interchange=interchange,
            atom_type_map=atom_type_map,
        )

        if n_bonds > 0:
            _write_bond_coeffs(lmp_file=lmp_file, interchange=interchange)
        if n_angles > 0:
            _write_angle_coeffs(lmp_file=lmp_file, interchange=interchange)
        if n_propers > 0:
            _write_proper_coeffs(lmp_file=lmp_file, interchange=interchange)
        if n_impropers > 0:
            _write_improper_coeffs(lmp_file=lmp_file, interchange=interchange)

        _write_atoms(
            lmp_file=lmp_file,
            interchange=interchange,
            atom_type_map=atom_type_map,
        )
        if n_bonds > 0:
            _write_bonds(lmp_file=lmp_file, interchange=interchange)
        if n_angles > 0:
            _write_angles(lmp_file=lmp_file, interchange=interchange)
        if n_propers > 0:
            _write_propers(lmp_file=lmp_file, interchange=interchange)
        if n_impropers > 0:
            _write_impropers(lmp_file=lmp_file, interchange=interchange)


def _write_pair_coeffs(lmp_file: IO, interchange: Interchange, atom_type_map: dict):
    """Write the Pair Coeffs section of a LAMMPS data file."""
    lmp_file.write("Pair Coeffs\n\n")

    vdw_handler = interchange["vdW"]

    for atom_type_idx, smirks in atom_type_map.items():
        params = vdw_handler.potentials[smirks].parameters

        sigma = params["sigma"].to(unit.angstrom).magnitude
        epsilon = params["epsilon"].to(unit.Unit("kilocalorie / mole")).magnitude

        lmp_file.write(f"{atom_type_idx + 1:d}\t{epsilon:.8g}\t{sigma:.8g}\n")

    lmp_file.write("\n")


def _write_bond_coeffs(lmp_file: IO, interchange: Interchange):
    """Write the Bond Coeffs section of a LAMMPS data file."""
    lmp_file.write("Bond Coeffs\n\n")

    bond_handler = interchange["Bonds"]
    bond_type_map = dict(enumerate(bond_handler.potentials))

    for bond_type_idx, smirks in bond_type_map.items():
        params = bond_handler.potentials[smirks].parameters

        k = params["k"].to(unit.Unit("kilocalorie / mole / angstrom ** 2")).magnitude
        k = k * 0.5  # Account for LAMMPS wrapping 1/2 into k
        length = params["length"].to(unit.angstrom).magnitude

        lmp_file.write(f"{bond_type_idx+1:d} harmonic\t{k:.16g}\t{length:.16g}\n")

    lmp_file.write("\n")


def _write_angle_coeffs(lmp_file: IO, interchange: Interchange):
    """Write the Angle Coeffs section of a LAMMPS data file."""
    lmp_file.write("\nAngle Coeffs\n\n")

    angle_handler = interchange["Angles"]
    angle_type_map = dict(enumerate(angle_handler.potentials))

    for angle_type_idx, smirks in angle_type_map.items():
        params = angle_handler.potentials[smirks].parameters

        k = params["k"].to(unit.Unit("kilocalorie / mole / radian ** 2")).magnitude
        k = k * 0.5  # Account for LAMMPS wrapping 1/2 into k
        theta = params["angle"].to(unit.degree).magnitude

        lmp_file.write(f"{angle_type_idx+1:d} harmonic\t{k:.16g}\t{theta:.16g}\n")

    lmp_file.write("\n")


def _write_proper_coeffs(lmp_file: IO, interchange: Interchange):
    """Write the Dihedral Coeffs section of a LAMMPS data file."""
    lmp_file.write("\nDihedral Coeffs\n\n")

    proper_handler = interchange["ProperTorsions"]
    proper_type_map = dict(enumerate(proper_handler.potentials))

    for proper_type_idx, smirks in proper_type_map.items():
        params = proper_handler.potentials[smirks].parameters

        k = params["k"].to(unit.Unit("kilocalorie / mole")).magnitude
        n = int(params["periodicity"])
        phase = params["phase"].to(unit.degree).magnitude
        idivf = int(params["idivf"])
        k = k / idivf

        lmp_file.write(
            f"{proper_type_idx+1:d} fourier 1\t{k:.16g}\t{n:d}\t{phase:.16g}\n",
        )

    lmp_file.write("\n")


def _write_improper_coeffs(lmp_file: IO, interchange: Interchange):
    """Write the Improper Coeffs section of a LAMMPS data file."""
    lmp_file.write("\nImproper Coeffs\n\n")

    improper_handler = interchange["ImproperTorsions"]
    improper_type_map = dict(enumerate(improper_handler.potentials))

    for improper_type_idx, smirks in improper_type_map.items():
        params = improper_handler.potentials[smirks].parameters

        k = params["k"].to(unit.Unit("kilocalorie / mole")).magnitude
        n = int(params["periodicity"])
        phase = params["phase"].to(unit.degree).magnitude
        idivf = int(params["idivf"])
        k = k / idivf

        # See https://lammps.sandia.gov/doc/improper_cvff.html
        # E_periodic = k * (1 + cos(n * theta - phase))
        # E_cvff = k * (1 + d * cos(n * theta))
        # k_periodic = k_cvff
        # if phase = 0,
        #   d_cvff = 1
        #   n_periodic = n_cvff
        # if phase = 180,
        #   cos(n * x - pi) == - cos(n * x)
        #   d_cvff = -1
        #   n_periodic = n_cvff
        # k * (1 + cos(n * phi - pi / 2)) == k * (1 - cos(n * phi))

        if phase == 0:
            k_cvff = k
            d_cvff = 1
            n_cvff = n
        elif phase == 180:
            k_cvff = k
            d_cvff = -1
            n_cvff = n
        else:
            raise UnsupportedExportError(
                "Improper exports to LAMMPS are funky and not well-supported, the only compatibility"
                "found between periodidic impropers is with improper_style cvff when phase = 0 or 180 degrees",
            )

        lmp_file.write(
            f"{improper_type_idx+1:d} {k_cvff:.16g}\t{d_cvff:d}\t{n_cvff:.16g}\n",
        )

    lmp_file.write("\n")


def _write_atoms(lmp_file: IO, interchange: Interchange, atom_type_map: dict):
    """Write the Atoms section of a LAMMPS data file."""
    lmp_file.write("\nAtoms\n\n")

    atom_type_map_inv = dict({v: k for k, v in atom_type_map.items()})

    electrostatics_handler = interchange["Electrostatics"]
    vdw_hander = interchange["vdW"]

    charges = electrostatics_handler.charges

    for atom in interchange.topology.atoms:
        atom_index = interchange.topology.atom_index(atom)
        try:
            molecule_index = int(atom.metadata["residue_number"])
        except KeyError:
            # TODO: Is there a mapping between molecules and their
            #       "molecule index" somewhere?
            molecule_index = 0

        top_key = AngleKey(atom_indices=(atom_index,))
        pot_key = vdw_hander.key_map[top_key]
        atom_type = atom_type_map_inv[pot_key]

        charge = charges[top_key].m_as(unit.e)
        pos = interchange.positions[atom_index].to(unit.angstrom).magnitude
        lmp_file.write(
            "{:d}\t{:d}\t{:d}\t{:.8g}\t{:.8g}\t{:.8g}\t{:.8g}\n".format(
                atom_index + 1,
                molecule_index + 1,
                atom_type + 1,
                charge,
                pos[0],
                pos[1],
                pos[2],
            ),
        )


def _write_bonds(lmp_file: IO, interchange: Interchange):
    """Write the Bonds section of a LAMMPS data file."""
    lmp_file.write("\nBonds\n\n")

    bond_handler = interchange["Bonds"]
    bond_type_map = dict(enumerate(bond_handler.potentials))

    bond_type_map_inv = dict({v: k for k, v in bond_type_map.items()})

    for bond_idx, bond in enumerate(interchange.topology.bonds):
        indices = (
            interchange.topology.atom_index(bond.atom1),
            interchange.topology.atom_index(bond.atom2),
        )
        top_key = BondKey(atom_indices=indices)
        if top_key in bond_handler.key_map:
            pot_key = bond_handler.key_map[top_key]
        else:
            top_key = BondKey(atom_indices=indices[::-1])
            pot_key = bond_handler.key_map[top_key]

        bond_type = bond_type_map_inv[pot_key]

        lmp_file.write(
            "{:d}\t{:d}\t{:d}\t{:d}\n".format(
                bond_idx + 1,
                bond_type + 1,
                indices[0] + 1,
                indices[1] + 1,
            ),
        )


def _write_angles(lmp_file: IO, interchange: Interchange):
    """Write the Angles section of a LAMMPS data file."""
    lmp_file.write("\nAngles\n\n")

    angle_handler = interchange["Angles"]
    angle_type_map = dict(enumerate(angle_handler.potentials))

    angle_type_map_inv = dict({v: k for k, v in angle_type_map.items()})

    for angle_idx, angle in enumerate(interchange.topology.angles):
        # These are "topology indices"
        indices = tuple(interchange.topology.atom_index(a) for a in angle)
        top_key = AngleKey(atom_indices=indices)
        pot_key = angle_handler.key_map[top_key]
        angle_type = angle_type_map_inv[pot_key]

        lmp_file.write(
            "{:d}\t{:d}\t{:d}\t{:d}\t{:d}\n".format(
                angle_idx + 1,
                angle_type + 1,
                indices[0] + 1,
                indices[1] + 1,
                indices[2] + 1,
            ),
        )


def _write_propers(lmp_file: IO, interchange: Interchange):
    """Write the Dihedrals section of a LAMMPS data file."""
    lmp_file.write("\nDihedrals\n\n")

    proper_handler = interchange["ProperTorsions"]
    proper_type_map = dict(enumerate(proper_handler.potentials))

    proper_type_map_inv = dict({v: k for k, v in proper_type_map.items()})

    for proper_idx, proper in enumerate(interchange.topology.propers):
        indices = tuple(interchange.topology.atom_index(a) for a in proper)

        for top_key, pot_key in proper_handler.key_map.items():
            if indices == top_key.atom_indices:
                proper_type_idx = proper_type_map_inv[pot_key]

                lmp_file.write(
                    "{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\n".format(
                        proper_idx + 1,
                        proper_type_idx + 1,
                        indices[0] + 1,
                        indices[1] + 1,
                        indices[2] + 1,
                        indices[3] + 1,
                    ),
                )


def _write_impropers(lmp_file: IO, interchange: Interchange):
    """Write the Impropers section of a LAMMPS data file."""
    lmp_file.write("\nImpropers\n\n")

    improper_handler = interchange["ImproperTorsions"]
    improper_type_map = dict(enumerate(improper_handler.potentials))

    improper_type_map_inv = dict({v: k for k, v in improper_type_map.items()})

    # Molecule/Topology.impropers lists the central atom **second** ...
    for improper_idx, improper in enumerate(interchange.topology.impropers):
        indices = tuple(interchange.topology.atom_index(a) for a in improper)

        # ... so the tuple must be modified to list the central atom **first**,
        # which is how the improper handler's slot map is built up
        _indices = tuple((indices[1], indices[0], indices[2], indices[3]))

        for top_key, pot_key in improper_handler.key_map.items():
            if _indices == top_key.atom_indices:
                improper_type_idx = improper_type_map_inv[pot_key]

                # https://github.com/openforcefield/openff-interchange/issues/544
                # LAMMPS, at least with `improper_style cvff`, lists the
                # central atom FIRST, whereas `indices` lists it SECOND
                lmp_file.write(
                    "{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\n".format(
                        improper_idx + 1,
                        improper_type_idx + 1,
                        indices[1] + 1,
                        indices[0] + 1,
                        indices[2] + 1,
                        indices[3] + 1,
                    ),
                )
