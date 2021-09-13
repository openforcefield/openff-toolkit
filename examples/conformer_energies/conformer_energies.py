import argparse
from openff.toolkit.topology import Molecule, Topology
from openff.toolkit.utils import RDKitToolkitWrapper
from rdkit.Chem import rdMolAlign
import numpy as np
from simtk import openmm, unit


def compute_conformer_energies_from_file(filename):
    # Load in the molecule and its conformers.
    # Note that all conformers of the same molecule are loaded as separate Molecule objects
    # If using a OFF Toolkit version before 0.7.0, loading SDFs through RDKit and OpenEye may provide
    # different behavior in some cases. So, here we force loading through RDKit to ensure the correct behavior
    rdktkw = RDKitToolkitWrapper()
    loaded_molecules = Molecule.from_file(filename, toolkit_registry=rdktkw)
    # The logic below only works for lists of molecules, so if a
    # single molecule was loaded, cast it to list
    if type(loaded_molecules) is not list:
        loaded_molecules = [loaded_molecules]
    # Collatate all conformers of the same molecule
    # NOTE: This isn't necessary if you have already loaded or created multi-conformer molecules;
    # it is just needed because our SDF reader does not automatically collapse conformers.
    molecules = [loaded_molecules[0]]
    for molecule in loaded_molecules[1:]:
        if molecule == molecules[-1]:
            for conformer in molecule.conformers:
                molecules[-1].add_conformer(conformer)
        else:
            molecules.append(molecule)

    n_molecules = len(molecules)
    n_conformers = sum([mol.n_conformers for mol in molecules])
    print(
        f"{n_molecules} unique molecule(s) loaded, with {n_conformers} total conformers"
    )

    # Load the openff-2.0.0 force field appropriate for vacuum calculations (without constraints)
    from openff.toolkit.typing.engines.smirnoff import ForceField

    forcefield_str = "openff_unconstrained-2.0.0.offxml"
    forcefield = ForceField(forcefield_str)
    # Loop over molecules and minimize each conformer
    for molecule in molecules:
        # If the molecule doesn't have a name, set mol.name to be the hill formula
        if molecule.name == "":
            molecule.name = Topology._networkx_to_hill_formula(molecule.to_networkx())
            print("%s : %d conformers" % (molecule.name, molecule.n_conformers))
            # Make a temporary copy of the molecule that we can update for each minimization
        mol_copy = Molecule(molecule)
        # Make an OpenFF Topology so we can parameterize the system
        off_top = molecule.to_topology()
        print(f"Parametrizing {molecule.name} (may take a moment to calculate charges)")
        system = forcefield.create_openmm_system(off_top)
        # Use OpenMM to compute initial and minimized energy for all conformers
        integrator = openmm.VerletIntegrator(1 * unit.femtoseconds)
        platform = openmm.Platform.getPlatformByName("Reference")
        omm_top = off_top.to_openmm()
        simulation = openmm.app.Simulation(omm_top, system, integrator, platform)

        # Print text header
        print(
            "Using force field",
            forcefield_str,
            "\nConformer         Initial PE         Minimized PE       RMS between initial and minimized conformer",
        )
        output = [
            [
                "Conformer",
                "Initial PE (kcal/mol)",
                "Minimized PE (kcal/mol)",
                "RMS between initial and minimized conformer (Angstrom)",
            ]
        ]
        for conformer_index, conformer in enumerate(molecule.conformers):
            simulation.context.setPositions(conformer)
            orig_potential = simulation.context.getState(
                getEnergy=True
            ).getPotentialEnergy()
            simulation.minimizeEnergy()
            min_state = simulation.context.getState(getEnergy=True, getPositions=True)
            min_potential = min_state.getPotentialEnergy()

            # Calculate the RMSD between the initial and minimized conformer
            min_coords = min_state.getPositions()
            min_coords = (
                np.array([[atom.x, atom.y, atom.z] for atom in min_coords])
                * unit.nanometer
            )
            mol_copy._conformers = None
            mol_copy.add_conformer(conformer)
            mol_copy.add_conformer(min_coords)
            rdmol = mol_copy.to_rdkit()
            rmslist = []
            rdMolAlign.AlignMolConformers(rdmol, RMSlist=rmslist)
            minimization_rms = rmslist[0]

            # Save the minimized conformer to file
            mol_copy._conformers = None
            mol_copy.add_conformer(min_coords)
            mol_copy.to_file(
                f"{molecule.name}_conf{conformer_index+1}_minimized.sdf",
                file_format="sdf",
            )
            print(
                "%5d / %5d : %8.3f kcal/mol %8.3f kcal/mol  %8.3f Angstroms"
                % (
                    conformer_index + 1,
                    molecule.n_conformers,
                    orig_potential / unit.kilocalories_per_mole,
                    min_potential / unit.kilocalories_per_mole,
                    minimization_rms,
                )
            )
            output.append(
                [
                    str(conformer_index + 1),
                    f"{orig_potential/unit.kilocalories_per_mole:.3f}",
                    f"{min_potential/unit.kilocalories_per_mole:.3f}",
                    f"{minimization_rms:.3f}",
                ]
            )
            # Write the results out to CSV
        with open(f"{molecule.name}.csv", "w") as of:
            for line in output:
                of.write(",".join(line) + "\n")
                # Clean up OpenMM Simulation
        del simulation, integrator


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Perform energy minimization on a molecule, potentially with many conformers. For each conformer, this script will provide the initial energy, minimized energy, and RMSD between the initial and minimized structure (both as STDOUT and a csv file). The minimized conformers will be written out to SDF."
    )
    parser.add_argument(
        "-f", "--filename", help="Name of an input file containing conformers"
    )
    args = parser.parse_args()

    compute_conformer_energies_from_file(args.filename)
