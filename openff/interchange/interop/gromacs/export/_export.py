import pathlib
import warnings
from typing import Optional, Union

import numpy
from openff.models.models import DefaultModel
from openff.units import unit

from openff.interchange.exceptions import MissingPositionsError
from openff.interchange.interop.gromacs.models.models import (
    GROMACSSystem,
    LennardJonesAtomType,
    PeriodicImproperDihedral,
    PeriodicProperDihedral,
    RyckaertBellemansDihedral,
)


class GROMACSWriter(DefaultModel):
    """Thin wrapper for writing GROMACS systems."""

    system: GROMACSSystem
    top_file: Optional[Union[pathlib.Path, str]] = None
    gro_file: Optional[Union[pathlib.Path, str]] = None

    def to_top(self):
        """Write a GROMACS topology file."""
        if self.top_file is None:
            raise ValueError("No TOP file specified.")

        with open(self.top_file, "w") as top:
            self._write_defaults(top)
            self._write_atomtypes(top)

            self._write_moleculetypes(top)

            self._write_system(top)
            self._write_molecules(top)

    def to_gro(self, decimal: int = 3):
        """Write a GROMACS coordinate file."""
        if self.gro_file is None:
            raise ValueError("No GRO file specified.")

        with open(self.gro_file, "w") as gro:
            self._write_gro(gro, decimal)

    def _write_defaults(self, top):
        top.write("[ defaults ]\n")
        top.write("; nbfunc\tcomb-rule\tgen-pairs\tfudgeLJ\tfudgeQQ\n")

        gen_pairs = "yes" if self.system.gen_pairs else "no"
        top.write(
            f"{self.system.nonbonded_function:6d}\t"
            f"{self.system.combination_rule:6d}\t"
            f"{gen_pairs:6s}\t"
            f"{self.system.vdw_14:8.6f}\t"
            f"{self.system.coul_14:8.6f}\n\n",
        ),

    def _write_atomtypes(self, top):
        top.write("[ atomtypes ]\n")
        top.write(
            ";type, bondingtype, atomic_number, mass, charge, ptype, sigma, epsilon\n",
        )

        for atom_type in self.system.atom_types.values():
            if not isinstance(atom_type, LennardJonesAtomType):
                raise NotImplementedError(
                    "Only Lennard-Jones atom types are currently supported.",
                )

            top.write(
                f"{atom_type.name :<11s}\t"
                f"{atom_type.atomic_number :6d}\t"
                f"{atom_type.mass.m :.16g}\t"
                f"{atom_type.charge.m :.16f}\t"
                f"{atom_type.particle_type :5s}\t"
                f"{atom_type.sigma.m :.16g}\t"
                f"{atom_type.epsilon.m :.16g}\n",
            )

        top.write("\n")

    def _write_moleculetypes(self, top):
        for molecule_name, molecule_type in self.system.molecule_types.items():
            top.write("[ moleculetype ]\n")

            top.write(
                f"{molecule_name.replace(' ', '_')}\t"
                f"{molecule_type.nrexcl:10d}\n\n",
            )

            self._write_atoms(top, molecule_type)
            self._write_pairs(top, molecule_type)
            self._write_bonds(top, molecule_type)
            self._write_angles(top, molecule_type)
            self._write_dihedrals(top, molecule_type)
            self._write_settles(top, molecule_type)
            self._write_exclusions(top, molecule_type)

        top.write("\n")

    def _write_atoms(self, top, molecule_type):
        top.write("[ atoms ]\n")
        top.write(";index, atom type, resnum, resname, name, cgnr, charge, mass\n")

        for atom in molecule_type.atoms:
            top.write(
                f"{atom.index :6d} "
                f"{atom.atom_type :6s}"
                f"{atom.residue_index :8d} "
                f"{atom.residue_name :8s} "
                f"{atom.name :6s}"
                f"{atom.charge_group_number :6d}"
                f"{atom.charge.m :20.12f}"
                f"{atom.mass.m :20.12f}\n",
            )

        top.write("\n")

    def _write_pairs(self, top, molecule_type):
        top.write("[ pairs ]\n")
        top.write(";ai    aj   funct\n")

        function = 1

        for pair in molecule_type.pairs:
            top.write(
                f"{pair.atom1 :6d}\t{pair.atom2 :6d}\t{function :6d}\n",
            )

        top.write("\n")

    def _write_bonds(self, top, molecule_type):
        top.write("[ bonds ]\n")
        top.write(";ai    aj   funct r k\n")

        function = 1

        for bond in molecule_type.bonds:
            top.write(
                f"{bond.atom1 :6d} "
                f"{bond.atom2 :6d} "
                f"{function :6d}"
                f"{bond.length.m :20.12f} "
                f"{bond.k.m :20.12f} ",
            )

            top.write("\n")

        top.write("\n")

    def _write_angles(self, top, molecule_type):
        top.write("[ angles ]\n")
        top.write(";ai    aj   ak   funct theta  k\n")

        function = 1

        for angle in molecule_type.angles:
            top.write(
                f"{angle.atom1 :6d} "
                f"{angle.atom2 :6d} "
                f"{angle.atom3 :6d} "
                f"{function :6d} "
                f"{angle.angle.m :20.12f} "
                f"{angle.k.m :20.12f} ",
            )

            top.write("\n")

        top.write("\n")

    def _write_dihedrals(self, top, molecule_type):
        top.write("[ dihedrals ]\n")
        top.write(";ai    aj   ak   al   funct phi  k\n")

        functions = {
            PeriodicProperDihedral: 1,
            RyckaertBellemansDihedral: 3,
            PeriodicImproperDihedral: 4,
        }

        for dihedral in molecule_type.dihedrals:
            function = functions[type(dihedral)]

            top.write(
                f"{dihedral.atom1 :6d}"
                f"{dihedral.atom2 :6d}"
                f"{dihedral.atom3 :6d}"
                f"{dihedral.atom4 :6d}"
                f"{functions[type(dihedral)] :6d}",
            )

            if function in [1, 4]:
                top.write(
                    f"{dihedral.phi.m :20.12f}"
                    f"{dihedral.k.m :20.12f}"
                    f"{dihedral.multiplicity :18d}",
                )

            elif function == 3:
                top.write(
                    f"{dihedral.c0.m :20.12f}"
                    f"{dihedral.c1.m :20.12f}"
                    f"{dihedral.c2.m :20.12f}"
                    f"{dihedral.c3.m :20.12f}"
                    f"{dihedral.c4.m :20.12f}"
                    f"{dihedral.c5.m :20.12f}",
                )

            else:
                raise ValueError(f"Invalid dihedral function {function}.")

            top.write("\n")

        top.write("\n")

    def _write_exclusions(self, top, molecule_type):
        top.write("[ exclusions ]\n")
        top.write(";ai    aj\n")

        for exclusion in molecule_type.exclusions:
            top.write(
                f"{exclusion.first_atom :6d}",
            )
            for other_atom in exclusion.other_atoms:
                top.write(
                    f"{other_atom :6d}",
                )

            top.write("\n")

        top.write("\n")

    def _write_settles(self, top, molecule_type):
        top.write("[ settles ]\n")
        top.write(";i  funct   dOH  dHH\n")

        function = 1

        for settle in molecule_type.settles:
            top.write(
                f"{settle.first_atom :6d}\t"
                f"{function :6d}\t"
                f"{settle.oxygen_hydrogen_distance.m_as(unit.nanometer) :20.12f}\t"
                f"{settle.hydrogen_hydrogen_distance.m_as(unit.nanometer):20.12f}\n",
            )

        top.write("\n")

    def _write_system(self, top):
        top.write("[ system ]\n")
        top.write(";name\n")

        top.write(f"{self.system.name}\n")

        top.write("\n")

    def _write_molecules(self, top):
        top.write("[ molecules ]\n")
        top.write(";name\tnumber\n")

        for molecule_name, n_molecules in self.system.molecules.items():
            top.write(f"{molecule_name.replace(' ', '_')}\t{n_molecules}\n")

        top.write("\n")

    def _write_gro(self, gro, decimal: int):
        if self.system.positions is None:
            raise MissingPositionsError(
                "Positions are required to write a `.gro` file but found None.",
            )
        elif numpy.allclose(self.system.positions, 0):
            warnings.warn(
                "Positions seem to all be zero. Result coordinate file may be non-physical.",
                UserWarning,
            )

        n_particles = sum(
            len(molecule_type.atoms) * self.system.molecules[molecule_name]
            for molecule_name, molecule_type in self.system.molecule_types.items()
        )

        assert n_particles == self.system.positions.shape[0], (
            n_particles,
            self.system.positions.shape[0],
        )

        # Explicitly round here to avoid ambiguous things in string formatting
        positions = numpy.round(self.system.positions, decimal).m_as(unit.nanometer)

        gro.write("Generated by Interchange\n")
        gro.write(f"{n_particles}\n")

        count = 0
        for molecule_name, molecule in self.system.molecule_types.items():
            n_copies = self.system.molecules[molecule_name]

            for _ in range(n_copies):
                for atom in molecule.atoms:
                    gro.write(
                        f"%5d%-5s%5s%5d"
                        f"%{decimal+5}.{decimal}f"
                        f"%{decimal+5}.{decimal}f"
                        f"%{decimal+5}.{decimal}f\n"
                        % (
                            atom.residue_index,  # This needs to be looked up from a different data structure
                            atom.residue_name,
                            atom.name,
                            (count + 1) % 100000,
                            positions[count, 0],
                            positions[count, 1],
                            positions[count, 2],
                        ),
                    )

                    count += 1

        if self.system.box is None:
            warnings.warn(
                "WARNING: System defined with no box vectors, which GROMACS does not offically "
                "support in versions 2020 or newer (see "
                "https://gitlab.com/gromacs/gromacs/-/issues/3526). Setting box vectors to a 5 "
                " nm cube.",
            )
            box = 5 * numpy.eye(3)
        else:
            box = self.system.box.m_as(unit.nanometer)

        # Check for rectangular
        if (box == numpy.diag(numpy.diagonal(box))).all():
            for i in range(3):
                gro.write(f"{box[i, i]:11.7f}")
        else:
            raise NotImplementedError(
                "Non-rectangular boxes are not yet tested. Please open an issue on GitHub if you "
                "need this feature.",
            )

        gro.write("\n")
