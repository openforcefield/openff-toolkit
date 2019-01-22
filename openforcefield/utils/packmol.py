#!/usr/bin/env python

# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Packmol API.

Authors
-------
* Levi N. Naden <levi.naden@choderalab.org>
* John D. Chodera <john.chodera@choderalab.org>
* Simon Boothroyd <simon.boothroyd@choderalab.org>

Based on SolvationToolkit (https://github.com/MobleyLab/SolvationToolkit)

"""


# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

import copy
import logging
import os
import random
import string
import subprocess
import tempfile
import shutil

from distutils.spawn import find_executable

from simtk import openmm
from simtk import unit
from simtk.openmm import app

# =============================================================================================
# Module Constants
# =============================================================================================

PACKMOL_PATH = find_executable("packmol") or shutil.which("packmol")

_HEADER_TEMPLATE = """
# Mixture
tolerance {0:f}
filetype pdb
output {1:s}
"""
# add_amber_ter

_BOX_TEMPLATE = """
structure {0:s}
  number {1:d}
  inside box 0. 0. 0. {2:f} {3:f} {4:f}
end structure
"""


# =============================================================================================
# Packmol
# =============================================================================================

def pack_box(molecules,
             n_copies,
             tolerance=2.0,
             box_size=None,
             mass_density=None,
             verbose=False):

    """Run packmol to generate a box containing a mixture of molecules.

    Parameters
    ----------
    molecules : list of OEMol
        Molecules in the system (with 3D geometries)
    n_copies : list of int (same length as 'molecules')
        Number of copies of the molecules
    tolerance : float, optional, default=2.0
        The mininum spacing between molecules during packing.  In ANGSTROMS!
    box_size : simtk.unit.Quantity in units compatible with angstroms
        The size of the box to generate.
        Default generates boxes that are very large for increased stability.
        May require extra time for energy minimization and equilibration.
    mass_density : simtk.unit.Quantity with units compatible with grams/milliliters, optional, default = 1.0*grams/milliliters
        Target mass density for final system, if available.
    verbose : bool, optional, default=False
        If True, verbose output is written.

    Returns
    -------
    topology : simtk.openmm.Topology
        Topology of the resulting system
    positions : simtk.unit.Quantity wrapped [natoms,3] numpy array with units compatible with angstroms
        Single frame trajectory with mixture box.

    """
    from openeye import oechem

    if len(molecules) != len(n_copies):

        logging.error("Length of 'molecules' and 'n_copies' must be identical")
        return None, None

    # Create PDB files for all components
    pdb_filenames = list()

    pdb_flavor = (oechem.OEOFlavor_PDB_CurrentResidues
                  | oechem.OEOFlavor_PDB_ELEMENT
                  | oechem.OEOFlavor_PDB_BONDS
                  | oechem.OEOFlavor_PDB_HETBONDS
                  | oechem.OEOFlavor_PDB_BOTH
                  | oechem.OEIFlavor_PDB_Connect)

    for molecule in molecules:

        tmp_filename = tempfile.mktemp(suffix=".pdb")
        pdb_filenames.append(tmp_filename)

        # Write PDB file
        mol_copy = copy.deepcopy(molecule)
        ofs = oechem.oemolostream(tmp_filename)
        ofs.SetFlavor(oechem.OEFormat_PDB, pdb_flavor)

        # Fix residue names
        residue_name = "".join([random.choice(string.ascii_uppercase) for i in range(3)])

        for atom in mol_copy.GetAtoms():

            residue = oechem.OEAtomGetResidue(atom)
            residue.SetName(residue_name)
            oechem.OEAtomSetResidue(atom, residue)

        oechem.OEWriteMolecule(ofs, mol_copy)
        ofs.close()

    # Run packmol
    if PACKMOL_PATH is None:
        raise IOError("Packmol not found, cannot run pack_box()")

    output_filename = tempfile.mktemp(suffix=".pdb")

    # Approximate volume to initialize box
    if box_size is None:

        if mass_density is not None:
            # Estimate box_size from mass density.
            box_size = approximate_volume_by_density(molecules, n_copies, mass_density=mass_density)
        else:
            # Use vdW radii to estimate box_size
            box_size = approximate_volume(molecules, n_copies)

    unitless_box_angstrom = box_size / box_size.unit

    header = _HEADER_TEMPLATE.format(tolerance, output_filename)

    for (pdb_filename, molecule, count) in zip(pdb_filenames, molecules, n_copies):

        header += _BOX_TEMPLATE.format(pdb_filename,
                                       count,
                                       unitless_box_angstrom,
                                       unitless_box_angstrom,
                                       unitless_box_angstrom)

    if verbose:
        print(header)

    # Write packmol input
    packmol_filename = tempfile.mktemp(suffix=".txt")

    with open(packmol_filename, 'w') as file_handle:
        file_handle.write(header)

    packmol_succeeded = False

    with open(packmol_filename) as file_handle:

        result = subprocess.check_output(PACKMOL_PATH,
                                         stdin=file_handle,
                                         stderr=subprocess.STDOUT).decode("utf-8")

        if verbose:
            print(result)

        packmol_succeeded = result.find('Success!') > 0

    os.unlink(packmol_filename)

    for filename in pdb_filenames:
        os.unlink(filename)

    if not packmol_succeeded:

        logging.error("Packmol failed to converge")
        os.unlink(output_filename)

        return None, None

    # Append missing connect statements to the end of the
    # output file.
    _append_connect_statements(output_filename, molecules, n_copies)

    # Read the resulting PDB file.
    pdbfile = app.PDBFile(output_filename)
    os.unlink(output_filename)

    # Extract topology and positions
    topology = pdbfile.getTopology()
    positions = pdbfile.getPositions()

    unitless_box_nm = box_size / unit.nanometers
    # import numpy as np
    # box_vectors = np.diag([box_size]*3)

    box_vector_x = openmm.Vec3(unitless_box_nm, 0, 0)
    box_vector_y = openmm.Vec3(0, unitless_box_nm, 0)
    box_vector_z = openmm.Vec3(0, 0, unitless_box_nm)

    # Set the periodic box vectors.
    topology.setPeriodicBoxVectors([box_vector_x, box_vector_y, box_vector_z] * unit.nanometers)
    # topology.setPeriodicBoxVectors(box_vectors)

    return topology, positions


def approximate_volume(molecules,
                       n_copies,
                       box_scaleup_factor=2.0):
    """Approximate the appropriate box size based on the number and types of atoms present.

    Parameters
    ----------
    molecules : list of OEMol
        Molecules in the system (with 3D geometries)
    n_copies : list of int (same length as 'molecules')
        Number of copies of the molecules
    box_scaleup_factor : float, optional, default = 2.0
        Factor by which the estimated box size is increased

    Returns
    -------
    box_size : simtk.unit.Quantity with units compatible with angstroms
        The size of the box to generate.

    Notes
    -----
    By default, boxes are very large for increased stability, and therefore may
    require extra time for energy minimization and equilibration.

    """
    from openeye import oechem

    volume = 0.0 * unit.angstrom**3

    for (molecule, number) in zip(molecules, n_copies):

        molecule_volume = 0.0 * unit.angstrom**3

        for atom in molecule.GetAtoms():
            molecule_volume += oechem.OEGetBondiVdWRadius(atom.GetAtomicNum()) * unit.angstrom**3

        volume += molecule_volume * number

    # Add 2 angs to help ease PBC issues.
    box_edge = volume**(1.0/3.0) * box_scaleup_factor + 2.0 * unit.angstrom

    return box_edge


def approximate_volume_by_density(molecules,
                                  n_copies,
                                  mass_density=1.0*unit.grams/unit.milliliters,
                                  box_scaleup_factor=1.1):
    """Generate an approximate box size based on the number and molecular weight of molecules present, and a target
    density for the final solvated mixture. If no density is specified, the target density is assumed to be 1 g/ml.

    Parameters
    ----------
    molecules : list of OEMol
        Molecules in the system (with 3D geometries)
    n_copies : list of int (same length as 'molecules')
        Number of copies of the molecules.
    box_scaleup_factor : float, optional, default = 1.1
        Factor by which the estimated box size is increased
    mass_density : simtk.unit.Quantity with units compatible with grams/milliliters, optional, default = 1.0*grams/milliliters
        Target mass density for final system, if available.

    Returns
    -------
    box_edge : simtk.unit.Quantity with units compatible with angstroms
        The size (edge length) of the box to generate.

    Notes
    -----
    By default, boxes are only modestly large. This approach has not been extensively tested for stability but has been
    used in the Mobley lab for perhaps ~100 different systems without substantial problems.

    """
    from openeye import oechem

    # Load molecules to get molecular weights
    volume = 0.0 * unit.angstrom**3

    for (molecule, number) in zip(molecules, n_copies):

        molecule_mass = oechem.OECalculateMolecularWeight(molecule) * \
                        unit.grams / unit.mole / unit.AVOGADRO_CONSTANT_NA

        molecule_volume = molecule_mass / mass_density

        volume += molecule_volume * number

    # Add 2 angs to help ease PBC issues.
    box_edge = volume**(1.0/3.0) * box_scaleup_factor + 2.0 * unit.angstrom

    return box_edge


def _append_connect_statements(file_name, molecules, n_copies):

    lines = []

    with open(file_name, 'r') as file:
        lines = file.readlines()

    if lines[len(lines) - 1].find('END') == 0:
        lines.pop()

    atom_counter = 0

    # TODO: Does packmol always give the exact number of mols asked for?
    # In future may be better way to figure this out.
    for (molecule, count) in zip(molecules, n_copies):

        bonds = {}

        for bond in molecule.GetBonds():

            index_A = bond.GetBgnIdx()
            index_B = bond.GetEndIdx()

            if index_A not in bonds:
                bonds[index_A] = []
            if index_B not in bonds:
                bonds[index_B] = []

            bonds[index_A].append(index_B)
            bonds[index_B].append(index_A)

        for i in range(count):

            for index_A in bonds:

                if len(bonds[index_A]) == 0:
                    continue

                connect_string = 'CONECT' + "%5d" % (index_A + atom_counter + 1)

                for j in range(len(bonds[index_A])):
                    connect_string += "%5d" % (bonds[index_A][j] + atom_counter + 1)

                connect_string += '\n'

                lines.append(connect_string)

            atom_counter += molecule.NumAtoms()

    lines.append('END')

    with open(file_name, 'w') as file:
        file.writelines(lines)


