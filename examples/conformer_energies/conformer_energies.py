import argparse

import openmm
from openff.interchange import Interchange
from openff.units.openmm import from_openmm
from rdkit.Chem import rdMolAlign

from openff.toolkit import ForceField, Molecule, RDKitToolkitWrapper


def compute_conformer_energies_from_file(filename):
    # Load in the molecule and its conformers.
    # Note that all conformers of the same molecule are loaded as separate Molecule objects
    # If using a OFF Toolkit version before 0.7.0, loading SDFs through RDKit and OpenEye may provide
    # different behavior in some cases. So, here we force loading through RDKit to ensure the correct behavior
    rdktkw = RDKitToolkitWrapper()
    loaded_molecules = Molecule.from_file(filename, toolkit_registry=rdktkw)
    # The logic below only works for lists of molecules, so if a
    # single molecule was loaded, cast it to list
    try:
        loaded_molecules = [*loaded_molecules]
    except TypeError:
        loaded_molecules = [loaded_molecules]

    # Collatate all conformers of the same molecule
    # NOTE: This isn't necessary if you have already loaded or created multi-conformer molecules;
    # it is just needed because our SDF reader does not automatically collapse conformers.
    molecule = loaded_molecules.pop(0)
    for next_molecule in loaded_molecules:
        if next_molecule == molecule:
            for conformer in next_molecule.conformers:
                molecule.add_conformer(conformer)
        else:
            # We're assuming the SDF just has multiple conformers of the
            # same molecule, so raise an error if that's not the case
            raise ValueError("Multiple chemical species loaded")

    # Make sure the molecule has a name
    if not molecule.name:
        molecule.name = molecule.to_hill_formula()

    print(
        f"Loaded {molecule.n_conformers} conformers"
        + f" of {molecule.to_smiles(explicit_hydrogens=False)!r}"
        + f" ({molecule.name})"
    )

    # Load the openff-2.1.0 force field appropriate for vacuum calculations (without constraints)
    forcefield = ForceField("openff_unconstrained-2.1.0.offxml")
    print(f"Parametrizing {molecule.name} (may take a moment to calculate charges)...")
    interchange = Interchange.from_smirnoff(forcefield, [molecule])
    print("Done.")
    integrator = openmm.VerletIntegrator(1 * openmm.unit.femtoseconds)
    simulation = interchange.to_openmm_simulation(integrator)

    # We'll store energies in two lists
    initial_energies = []
    minimized_energies = []

    # And minimized conformers in a second molecule
    minimized_molecule = Molecule(molecule)
    minimized_molecule.conformers.clear()

    for conformer in molecule.conformers:
        # Tell the OpenMM Simulation the positions of this conformer
        simulation.context.setPositions(conformer.to_openmm())

        # Keep a record of the initial energy
        initial_energies.append(
            simulation.context.getState(getEnergy=True).getPotentialEnergy()
        )

        # Perform the minimization
        simulation.minimizeEnergy()

        # Record minimized energy and positions
        min_state = simulation.context.getState(getEnergy=True, getPositions=True)

        minimized_energies.append(min_state.getPotentialEnergy())
        minimized_molecule.add_conformer(from_openmm(min_state.getPositions()))

    n_confs = molecule.n_conformers
    print(f"{molecule.name}: {n_confs} conformers")

    # Create a copy of the molecule so we can work on it
    working_mol = Molecule(molecule)

    # Print text header
    print("Conformer         Initial PE        Minimized PE        RMSD")
    output = [
        [
            "Conformer",
            "Initial PE (kcal/mol)",
            "Minimized PE (kcal/mol)",
            "RMSD between initial and minimized conformer (Angstrom)",
        ]
    ]

    for i, (init_energy, init_coords, min_energy, min_coords) in enumerate(
        zip(
            initial_energies,
            molecule.conformers,
            minimized_energies,
            minimized_molecule.conformers,
        )
    ):
        # Clear the conformers from the working molecule
        working_mol.conformers.clear()
        # Save the minimized conformer to file
        working_mol.add_conformer(min_coords)
        working_mol.to_file(
            f"{molecule.name}_conf{i+1}_minimized.sdf",
            file_format="sdf",
        )

        # Calculate the RMSD between the initial and minimized conformer
        working_mol.add_conformer(init_coords)
        rdmol = working_mol.to_rdkit()
        rmslist = []
        rdMolAlign.AlignMolConformers(rdmol, RMSlist=rmslist)
        minimization_rms = rmslist[0]

        # Record the results
        output.append(
            [
                i + 1,
                init_energy.value_in_unit(openmm.unit.kilocalories_per_mole),
                min_energy.value_in_unit(openmm.unit.kilocalories_per_mole),
                minimization_rms,
            ]
        )
        print(
            f"{{:5d}} / {n_confs:5d} :  {{:8.3f}} kcal/mol {{:8.3f}} kcal/mol {{:8.3f}} Å".format(
                *output[-1]
            )
        )

    # Write the results out to CSV
    with open(f"{molecule.name}.csv", "w") as of:
        of.write(", ".join(output.pop(0)) + "\n")
        for line in output:
            of.write("{}, {:.3f}, {:.3f}, {:.3f}".format(*line) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Perform energy minimization on a molecule, potentially with many conformers. For each conformer, this script will provide the initial energy, minimized energy, and RMSD between the initial and minimized structure (both as STDOUT and a csv file). The minimized conformers will be written out to SDF."
    )
    parser.add_argument(
        "-f", "--filename", help="Name of an input file containing conformers"
    )
    args = parser.parse_args()

    compute_conformer_energies_from_file(args.filename)
