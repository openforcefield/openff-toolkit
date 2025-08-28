import argparse

from openff.interchange.drivers.openmm import get_openmm_energies
from rdkit.Chem import rdMolAlign

from openff.toolkit import ForceField, Molecule
from openff.toolkit.utils import get_data_file_path


def compute_conformer_energies_from_file(filename):
    # First, load conformers from an SDF file.
    loaded_molecules = Molecule.from_file(
        get_data_file_path("molecules/ruxolitinib_conformers.sdf"),
    )

    # Normalize to list
    try:
        loaded_molecules = [*loaded_molecules]
    except TypeError:
        loaded_molecules = [loaded_molecules]

    # from_file loads each entry in the SDF into its own molecule,
    # so collapse conformers into the same molecule
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

    # Load the openff-2.2.1 force field appropriate for vacuum calculations (without constraints)
    forcefield = ForceField("openff_unconstrained-2.2.1.offxml")
    # Create an Interchange object, which stores the result of parametrizing with this force field
    print(f"Parametrizing {molecule.name} (may take a moment to calculate charges)...")
    interchange = forcefield.create_interchange(molecule.to_topology())
    print("Done.")

    # We'll store energies in two lists
    initial_energies = []
    minimized_energies = []

    # And minimized conformers in a second molecule
    minimized_molecule = Molecule(molecule)
    minimized_molecule.conformers.clear()

    for conformer in molecule.conformers:
        # Use this conformer to update the positions of the Interchange object
        interchange.positions = conformer

        # Get the (total) initial energy from this conformer and store it
        initial_energies.append(get_openmm_energies(interchange).total_energy)

        # Minimize using Interchange.minimize, which wraps OpenMM
        interchange.minimize(engine="openmm")

        # Record minimized energy and positions
        minimized_energies.append(get_openmm_energies(interchange).total_energy)
        minimized_molecule.add_conformer(interchange.positions.to("angstrom"))

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
            f"{molecule.name}_conf{i + 1}_minimized.sdf",
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
                init_energy.m_as("kilocalories_per_mole"),
                min_energy.m_as("kilocalories_per_mole"),
                minimization_rms,
            ]
        )
        print(f"{{:5d}} / {n_confs:5d} :  {{:8.3f}} kcal/mol {{:8.3f}} kcal/mol {{:8.3f}} Å".format(*output[-1]))

    # Write the results out to CSV
    with open(f"{molecule.name}.csv", "w") as of:
        of.write(", ".join(output.pop(0)) + "\n")
        for line in output:
            of.write("{}, {:.3f}, {:.3f}, {:.3f}".format(*line) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Perform energy minimization on a molecule, potentially with many conformers. For each conformer, this "
            "script will provide the initial energy, minimized energy, and RMSD between the initial and minimized "
            "structure (both as STDOUT and a csv file). The minimized conformers will be written out to SDF."
        ),
    )
    parser.add_argument("-f", "--filename", help="Name of an input file containing conformers")
    args = parser.parse_args()

    compute_conformer_energies_from_file(args.filename)
