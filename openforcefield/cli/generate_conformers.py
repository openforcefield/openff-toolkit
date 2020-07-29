import argparse
from copy import deepcopy

from openforcefield.typing.engines.smirnoff import ForceField
from openforcefield.topology import Molecule, Topology
from openforcefield.topology.molecule import UndefinedStereochemistryError
from openforcefield.utils.toolkits import ToolkitRegistry

# from rdkit.Chem import rdMolAlign
import numpy as np
from simtk import openmm, unit


def generate_conformers(
    molecule, toolkit, forcefield, rms_cutoff=None, constrained=False,
):
    if toolkit.lower() == "openeye":
        from openforcefield.utils.toolkits import OpenEyeToolkitWrapper

        toolkit_registry = ToolkitRegistry(toolkit_precedence=[OpenEyeToolkitWrapper])
    elif toolkit.lower() == "rdkit":
        from openforcefield.utils.toolkits import RDKitToolkitWrapper

        toolkit_registry = ToolkitRegistry(toolkit_precedence=[RDKitToolkitWrapper])

    ff_name = forcefield
    if constrained:
        tmp = ff_name.split("-")
        tmp[0] = "openff_unconstrained"
        ff_name = "-".join(tmp)
    if not ff_name.endswith(".offxml"):
        ff_name += ".offxml"
    try:
        ff = ForceField(ff_name)
    except Exception as e:
        print("catch speifics here")
        raise e

    ambiguous_stereochemistry = False
    try:
        mols = toolkit_registry.call(
            "from_file", molecule, file_format=molecule.split(".")[-1],
        )
    except UndefinedStereochemistryError:
        ambiguous_stereochemistry = True
        raw_mols = toolkit_registry.call(
            "from_file", molecule, file_format=molecule.split(".")[-1], allow_undefined_stereo=True,
        )
        # TODO: This is a brute-force approach, it would be better to check stereo
        #  without needing to call enumerate_stereoisomers
        mols = []
        for mol in raw_mols:
            stereoisomers = mol.enumerate_stereoisomers()
            if stereoisomers:
                [mols.append(val) for val in stereoisomers]
            else:
                mols.append(mol)

    if type(mols) != list:
        mols = [mols]

    mols = _collapse_conformers(mols)

    for mol in mols:
        toolkit_registry.call('generate_conformers', molecule=mol, n_conformers=5, rms_cutoff=None)
        top = Topology.from_molecules(mol)
        omm_sys, new_top = ff.create_openmm_system(top, return_topology=True)

    for mol in mols:
        simulation = _build_simulation(molecule=mol, forcefield=ff)

        for conformer in mol.conformers:
            energy, positions = _get_minimized_data(conformer, simulation)

            _save_minimized_conformer(mol, energy, positions, 'a.sdf')


def _collapse_conformers(molecules):
    """
    Collapse conformers of isomorphic molecules into single Molecule objects.

    This is useful when reading in multi-conformer SDF files because the SDF
    reader does not automatically collapse conformers. This function should
    not modify lists of Molecule objects whose conformers are already collapsed.

    Parameters
    ----------
    molecules : list of Molecule
        List of Molecule objects

    Returns
    -------
    collapsed_molecules : list of Molecule
        List of Molecule objects with only one object per isomorphic molecule
        and conformers of isomorphic molecules stored in each

    """
    collapsed_molecules = [molecules[0]]
    for molecule in molecules[1:]:
        if molecule == collapsed_molecules[-1]:
            for conformer in molecule.conformers:
                collapsed_molecules[-1].add_conformer(conformer)
        else:
            collapsed_molecules.append(molecule)

    return collapsed_molecules


def _build_simulation(molecule, forcefield):
    """Given a Molecule and ForceField, initialize a barebones OpenMM Simulation."""
    mol_copy = deepcopy(molecule)
    off_top = molecule.to_topology()
    system = forcefield.create_openmm_system(off_top)

    # Use OpenMM to compute initial and minimized energy for all conformers
    integrator = openmm.VerletIntegrator(1 * unit.femtoseconds)
    platform = openmm.Platform.getPlatformByName('Reference')
    omm_top = off_top.to_openmm()
    simulation = openmm.app.Simulation(omm_top, system, integrator, platform)

    return simulation


def _get_minimized_data(conformer, simulation):
    """Given an OpenMM simulation and conformer, minimze and return an energy"""
    simulation.context.setPositions(conformer)
    simulation.minimizeEnergy()

    min_state = simulation.context.getState(getEnergy=True, getPositions=True)
    min_energy = min_state.getPotentialEnergy()
    min_coords = min_state.getPositions()

    return min_energy, min_coords

def _save_minimized_conformer(mol, energy, positions, filename):
    """Save minimized structures to SDF files"""
    # TODO: Add SD data pairs in output SDF file for
    #   * absolute_energy in kcal/mol
    #   * conformer generation method and version (openeye version or rdkit version)
    #   * partial charges (so, use ForceField.to_openmm_system's return_topology kwarg, and attach the final coordinates to ret_topology._reference_molecules[0], then save that to file)
    mol = deepcopy(mol)
    mol._conformers = None
    min_coords = np.array([[atom.x, atom.y, atom.z] for atom in positions]) * unit.nanometer
    mol.add_conformer(min_coords)
    # Add `energy` to mol
    mol.to_file(filename, file_format='sdf')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate conformers with cheminformatics toolkits"
    )
    parser.add_argument(
        "-t", "--toolkit", type=str, help="Name of the underlying cheminformatics toolkit to use",
    )
    parser.add_argument(
        "-f", "--forcefield", type=str, help="Name of the force field to use",
    )
    parser.add_argument(
        "-m",
        "--molecule",
        type=str,
        help="Path to an input file containing a molecule(s)",
    )
    parser.add_argument(
        "-r",
        "--rms-cutoff",
        type=float,
        default=0.25,
        help="The redundancy cutoff between pre-minimized conformers",
    )
    args = parser.parse_args()  # ['-t', '-f', '-m'])

    generate_conformers(
        molecule=args.molecule, toolkit=args.toolkit, forcefield=args.forcefield
    )
