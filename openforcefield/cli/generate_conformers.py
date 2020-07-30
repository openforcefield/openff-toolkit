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
        print("catch specifics here")
        raise e

    # TODO: Parse prefix of molecule file
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

    # If no conformers were found (i.e. SMI files), generate some
    # TODO: How many should be generated?
    # TODO: If 1 or a few conformers are found, should more be generated?
    for mol in mols:
        if not mol.conformers:
            mol.generate_conformers(toolkit_registry=toolkit_registry, n_conformers=1)

    if type(mols) != list:
        mols = [mols]

    # TODO: What happens if some molecules in a multi-molecule file have charges, others don't?
    mols_with_charges = []
    for mol in mols:
        if mol.partial_charges is not None:
            mols_with_charges.append(mol)

    mols = _collapse_conformers(mols)

    mols_out = []
    for mol in mols:
        simulation, partial_charges = _build_simulation(molecule=mol, forcefield=ff, mols_with_charge=mols_with_charges)
        mol._partial_charges = partial_charges

        for conformer in mol.conformers:
            energy, positions = _get_minimized_data(conformer, simulation)

            mols_out.append(_reconstruct_mol_from_conformer(mol, energy, positions))
            # _save_minimized_conformer(mol, energy, positions, 'a.sdf')

    return mols_out

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


def _build_simulation(molecule, forcefield, mols_with_charge):
    """Given a Molecule and ForceField, initialize a barebones OpenMM Simulation."""
    mol_copy = deepcopy(molecule)
    off_top = molecule.to_topology()
    system, top_with_charges = forcefield.create_openmm_system(
        off_top,
        charge_from_molecules=mols_with_charge,
        allow_nonintegral_charges=True,
        return_topology=True
    )

    # Use OpenMM to compute initial and minimized energy for all conformers
    integrator = openmm.VerletIntegrator(1 * unit.femtoseconds)
    platform = openmm.Platform.getPlatformByName('Reference')
    omm_top = off_top.to_openmm()
    simulation = openmm.app.Simulation(omm_top, system, integrator, platform)

    # TODO: return_topology does not produce a top with charges, either make it do so
    # or reconstruct partial charges from OpenMM system
    # partial_charges = [*top_with_charges.reference_molecules][0].partial_charges
    partial_charges = [system.getForces()[0].getParticleParameters(i)[0] for i in range(mol_copy.n_atoms)]
    # Unwrap list of Quantity objects into a single Quantity that contains a list
    # Surely there's a simpler way to to this?
    partial_charges = unit.Quantity(
        np.asarray([val.value_in_unit(unit.elementary_charge) for val in partial_charges]),
        unit=unit.elementary_charge,
    )

    return simulation, partial_charges


def _get_minimized_data(conformer, simulation):
    """Given an OpenMM simulation and conformer, minimze and return an energy"""
    simulation.context.setPositions(conformer)
    simulation.minimizeEnergy()

    min_state = simulation.context.getState(getEnergy=True, getPositions=True)
    min_energy = min_state.getPotentialEnergy()
    min_coords = min_state.getPositions()

    return min_energy, min_coords

def _reconstruct_mol_from_conformer(mol, energy, positions):
    # TODO: Add SD data pairs in output SDF file for
    #   * absolute_energy in kcal/mol
    #   * conformer generation method and version (openeye version or rdkit version)
    #   * partial charges (so, use ForceField.to_openmm_system's return_topology kwarg, and attach the final coordinates to ret_topology._reference_molecules[0], then save that to file)
    mol = deepcopy(mol)
    mol._conformers = None
    min_coords = np.array([[atom.x, atom.y, atom.z] for atom in positions]) * unit.nanometer
    mol.add_conformer(min_coords)
    return mol

def _save_minimized_conformer(mol, energy, positions, filename):
    """Save minimized structures to SDF files"""
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
