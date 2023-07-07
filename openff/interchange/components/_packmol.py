"""
A wrapper around PACKMOL. Taken from OpenFF Evaluator v0.4.3.
"""
import os
import random
import shutil
import string
import subprocess
import tempfile
import warnings
from collections import defaultdict
from distutils.spawn import find_executable
from typing import Optional

import mdtraj
import numpy
from openff.toolkit import Molecule, Topology
from openff.toolkit.utils.rdkit_wrapper import RDKitToolkitWrapper
from openff.units import Quantity, unit
from openff.utilities.utilities import requires_package, temporary_cd

from openff.interchange.exceptions import PACKMOLRuntimeError, PACKMOLValueError


def _find_packmol() -> Optional[str]:
    """
    Attempt to find the path to the `packmol` binary.

    Returns
    -------
    str, optional
        The path to the packmol binary if it could be found, otherwise
        `None`.

    """
    return (
        find_executable("packmol") or shutil.which("packmol") or None
        if "PACKMOL" not in os.environ
        else os.environ["PACKMOL"]
    )


def _validate_inputs(
    molecules: list[Molecule],
    number_of_copies: list[int],
    structure_to_solvate: Optional[Topology],
    box_aspect_ratio: Optional[Quantity],
    box_size: Optional[Quantity],
    mass_density: Optional[Quantity],
):
    """
    Validate the inputs which were passed to the main pack method.

    Parameters
    ----------
    molecules : list of openff.toolkit.topology.Molecule
        The molecules in the system.
    number_of_copies : list of int
        A list of the number of copies of each molecule type, of length
        equal to the length of `molecules`.
    structure_to_solvate: Topology, optional
        A topology to be solvated.
    box_size : openff.units.Quantity, optional
        The size of the box to generate in units compatible with angstroms.
        If `None`, `mass_density` must be provided.
    mass_density : openff.units.Quantity, optional
        Target mass density for final system with units compatible with g / mL.
         If `None`, `box_size` must be provided.
    box_aspect_ratio: list of float, optional
        The aspect ratio of the simulation box, used in conjunction with
        the `mass_density` parameter.

    """
    if box_size is None and mass_density is None:
        raise PACKMOLValueError(
            "Either a `box_size` or `mass_density` must be specified.",
        )

    if box_size is not None and len(box_size) != 3:
        raise PACKMOLValueError(
            "`box_size` must be a openff.units.unit.Quantity wrapped list of length 3",
        )

    if box_aspect_ratio is not None:
        assert len(box_aspect_ratio) == 3
        assert all(x > 0.0 for x in box_aspect_ratio)

    if len(molecules) != len(number_of_copies):
        raise PACKMOLValueError(
            "The length of `molecules` and `number_of_copies` must be identical.",
        )

    if structure_to_solvate is not None:
        if not isinstance(structure_to_solvate, Topology):
            raise PACKMOLValueError(
                "`structure_to_solvate` must be a openff.toolkit.topology.Topology",
            )

        positions = structure_to_solvate.get_positions()

        try:
            assert positions is not None
            assert positions.shape[0] == structure_to_solvate.n_atoms
        except AssertionError:
            raise PACKMOLValueError(
                "`structure_to_solvate` missing some atomic positions.",
            )


def _approximate_box_size_by_density(
    molecules: list[Molecule],
    n_copies: list[int],
    mass_density: Quantity,
    box_aspect_ratio: list[float],
    box_scaleup_factor: float = 1.1,
) -> Quantity:
    """
    Approximate box size.

    Generate an approximate box size based on the number and molecular
    weight of the molecules present, and a target density for the final
    solvated mixture.

    Parameters
    ----------
    molecules : list of openff.toolkit.topology.Molecule
        The molecules in the system.
    n_copies : list of int
        The number of copies of each molecule.
    mass_density : openff.units.Quantity
        The target mass density for final system. It should have units
        compatible with g / mL.
    box_aspect_ratio: List of float
        The aspect ratio of the simulation box, used in conjunction with
        the `mass_density` parameter.
    box_scaleup_factor : float
        The factor by which the estimated box size should be
        increased.

    Returns
    -------
    openff.units.Quantity
        A list of the three box lengths in units compatible with angstroms.

    """
    total_mass = 0.0 * unit.dalton
    for molecule, number in zip(molecules, n_copies):
        for atom in molecule.atoms:
            total_mass += number * atom.mass
    volume = total_mass / mass_density

    box_length = volume ** (1.0 / 3.0) * box_scaleup_factor
    box_length_angstrom = box_length.to(unit.angstrom).magnitude

    aspect_ratio_normalizer = (
        box_aspect_ratio[0] * box_aspect_ratio[1] * box_aspect_ratio[2]
    ) ** (1.0 / 3.0)

    box_size = [
        box_length_angstrom * box_aspect_ratio[0],
        box_length_angstrom * box_aspect_ratio[1],
        box_length_angstrom * box_aspect_ratio[2],
    ] * unit.angstrom

    box_size /= aspect_ratio_normalizer

    return box_size


def _standardize_smiles(smiles: str) -> str:
    """
    Standardize a SMILES pattern using the OpenFF Toolkit and RDKit.

    Taken from OpenFF Evaluator v0.4.3.

    See openff.evaluator.substances.components::Component._standardize_smiles

    """
    from openff.toolkit.utils.rdkit_wrapper import RDKitToolkitWrapper
    from openff.toolkit.utils.toolkit_registry import ToolkitRegistry

    # This parsing was previously done with `cmiles.utils.load_molecule`, which
    # * did NOT enforce stereochemistry while parsing SMILES and
    # * implicitly used the same toolkit to write the SMILES back from an object
    # This is hard-coded to keep test results consistent across OpenEye status
    # and compared to older versions; if desired this could be relaxed
    rdkit_registry = ToolkitRegistry(toolkit_precedence=[RDKitToolkitWrapper()])

    molecule = Molecule.from_smiles(
        smiles,
        toolkit_registry=rdkit_registry,
        allow_undefined_stereo=True,
    )

    try:
        # Try to make the smiles isomeric.
        return molecule.to_smiles(
            isomeric=True,
            explicit_hydrogens=False,
            mapped=False,
            toolkit_registry=rdkit_registry,
        )
    except ValueError:
        # Fall-back to non-isomeric.
        return molecule.to_smiles(
            isomeric=False,
            explicit_hydrogens=False,
            mapped=False,
            toolkit_registry=rdkit_registry,
        )


def _generate_residue_name(residue, smiles):
    """
    Generate residue names corresponding to a particular smiles pattern.

    Where possible (i.e for amino acids and ions) a standard residue
    name will be returned, otherwise a random name will be used.

    Parameters
    ----------
    residue: mdtraj.core.topology.Residue
        The residue to assign the name to.
    smiles: str
        The SMILES pattern to generate a resiude name for.

    """
    from mdtraj.core import residue_names
    from openff.toolkit.topology import Molecule

    # Define the set of residue names which should be discarded
    # if randomly generated as they have a reserved meaning.
    forbidden_residue_names = [
        *residue_names._AMINO_ACID_CODES,
        *residue_names._SOLVENT_TYPES,
        *residue_names._WATER_RESIDUES,
        "ADE",
        "CYT",
        "CYX",
        "DAD",
        "DGU",
        "FOR",
        "GUA",
        "HID",
        "HIE",
        "HIH",
        "HSD",
        "HSH",
        "HSP",
        "NMA",
        "THY",
        "URA",
    ]

    amino_residue_mappings = {
        "C[C@H](N)C(=O)O": "ALA",
        "N=C(N)NCCC[C@H](N)C(=O)O": "ARG",
        "NC(=O)C[C@H](N)C(=O)O": "ASN",
        "N[C@@H](CC(=O)O)C(=O)O": "ASP",
        "N[C@@H](CS)C(=O)O": "CYS",
        "N[C@@H](CCC(=O)O)C(=O)O": "GLU",
        "NC(=O)CC[C@H](N)C(=O)O": "GLN",
        "NCC(=O)O": "GLY",
        "N[C@@H](Cc1c[nH]cn1)C(=O)O": "HIS",
        "CC[C@H](C)[C@H](N)C(=O)O": "ILE",
        "CC(C)C[C@H](N)C(=O)O": "LEU",
        "NCCCC[C@H](N)C(=O)O": "LYS",
        "CSCC[C@H](N)C(=O)O": "MET",
        "N[C@@H](Cc1ccccc1)C(=O)O": "PHE",
        "O=C(O)[C@@H]1CCCN1": "PRO",
        "N[C@@H](CO)C(=O)O": "SER",
        "C[C@@H](O)[C@H](N)C(=O)O": "THR",
        "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O": "TRP",
        "N[C@@H](Cc1ccc(O)cc1)C(=O)O": "TYR",
        "CC(C)[C@H](N)C(=O)O": "VAL",
    }

    standardized_smiles = _standardize_smiles(smiles)

    # Check for amino acids.
    if standardized_smiles in amino_residue_mappings:
        residue.name = amino_residue_mappings[standardized_smiles]
        return

    # Check for water
    if standardized_smiles == "O":
        residue.name = "HOH"

        # Re-assign the water atom names. These need to be set to get
        # correct CONECT statements.
        h_counter = 1

        for atom in residue.atoms:
            if atom.element.symbol == "O":
                atom.name = "O1"
            else:
                atom.name = f"H{h_counter}"
                h_counter += 1

        return

    # Check for ions
    openff_molecule = Molecule.from_smiles(smiles, allow_undefined_stereo=True)

    if openff_molecule.n_atoms == 1:
        residue.name = _ion_residue_name(openff_molecule)
        residue.atom(0).name = residue.name

        return

    # Randomly generate a name
    random_residue_name = "".join(
        [random.choice(string.ascii_uppercase) for _ in range(3)],
    )

    while random_residue_name in forbidden_residue_names:
        # Re-choose the residue name until we find a safe one.
        random_residue_name = "".join(
            [random.choice(string.ascii_uppercase) for _ in range(3)],
        )

    residue.name = random_residue_name

    # Assign unique atom names.
    element_counter = defaultdict(int)

    for atom in residue.atoms:
        atom.name = f"{atom.element.symbol}{element_counter[atom.element.symbol] + 1}"
        element_counter[atom.element.symbol] += 1


def _ion_residue_name(molecule: Molecule) -> str:
    """
    Generate a residue name for a monatomic ion.

    Parameters
    ----------
    molecule: openff.toolkit.topology.Molecule
        The monoatomic ion to generate a resiude name for.

    Returns
    -------
    residue_name: str
        The residue name of the ion

    """
    element_symbol = molecule.atoms[0].symbol
    charge_symbol = ""

    formal_charge = molecule.atoms[0].formal_charge

    if isinstance(formal_charge, unit.Quantity):
        formal_charge = formal_charge.m_as(unit.elementary_charge)

    formal_charge = int(formal_charge)

    if formal_charge != 0:
        charge_symbol = "-" if formal_charge < 0 else "+"
        formal_charge = abs(formal_charge)

        if formal_charge > 1:
            charge_symbol = f"{formal_charge}{charge_symbol}"

    residue_name = f"{element_symbol}{charge_symbol}"

    # This truncates i.e. Mg2+ to Mg. What's supposed to happen here?
    residue_name = residue_name[:3]

    return residue_name


def _create_trajectory(molecule: Molecule) -> mdtraj.Trajectory:
    """
    Create an `mdtraj` topology from a molecule object.

    Parameters
    ----------
    molecule: openff.toolkit.topology.Molecule
        The SMILES pattern.

    Returns
    -------
    mdtraj.Trajectory
        The created trajectory.

    """
    import mdtraj

    # Check whether the molecule has a configuration defined, and if not,
    # define one.
    if molecule.n_conformers <= 0:
        molecule.generate_conformers(n_conformers=1)

    # We need to save out the molecule and then reload it as the toolkit
    # will not always save the atoms in the same order that they are
    # present in the molecule object.
    with tempfile.NamedTemporaryFile(suffix=".pdb") as file:
        # Toolkit's PDB writer is untrustworthy (for non-biopolymers) with OpenEyeToolktiWrapper
        # See https://github.com/openforcefield/openff-toolkit/issues/1307 and linked issues
        molecule.to_file(
            file.name,
            file_format="PDB",
            toolkit_registry=RDKitToolkitWrapper(),
        )
        # Load the pdb into an mdtraj object.
        mdtraj_trajectory = mdtraj.load_pdb(file.name)

    # Change the assigned residue name (sometimes molecules are assigned
    # an amino acid residue name even if that molecule is not an amino acid,
    # e.g. C(CO)N is not Gly) and save the altered object as a pdb.
    for residue in mdtraj_trajectory.topology.residues:
        _generate_residue_name(residue, molecule.to_smiles())

    return mdtraj_trajectory


def _build_input_file(
    molecule_file_names: list[str],
    molecule_counts: list[int],
    structure_to_solvate: Optional[Topology],
    center_solute: bool,
    box_size: Quantity,
    tolerance: Quantity,
    output_file_name: str,
) -> str:
    """
    Construct the packmol input file.

    Parameters
    ----------
    molecule_file_names: list of str
        The paths to the molecule pdb files.
    molecule_counts: list of int
        The number of each molecule to add.
    structure_to_solvate: Topology, optional
        The topology to solvate.
    center_solute: bool
        If `True`, the structure to solvate will be centered in the
        simulation box.
    box_size: openff.units.Quantity
        The lengths of each box vector.
    tolerance: openff.units.Quantity
        The packmol convergence tolerance.
    output_file_name: str
        The path to save the packed pdb to.

    Returns
    -------
    str
        The path to the input file.

    """
    box_size = box_size.m_as(unit.angstrom)
    tolerance = tolerance.m_as(unit.angstrom)

    # Add the global header options.
    input_lines = [
        f"tolerance {tolerance:f}",
        "filetype pdb",
        f"output {output_file_name}",
        "",
    ]

    # Add the section of the molecule to solvate if provided.
    if structure_to_solvate is not None:
        structure_to_solvate.to_file("_structure.pdb")
        solute_position = [0.0] * 3

        if center_solute:
            solute_position = [box_size[i] / 2.0 for i in range(3)]

        input_lines.extend(
            [
                "structure _structure.pdb",
                "  number 1",
                "  fixed "
                f"{solute_position[0]} "
                f"{solute_position[1]} "
                f"{solute_position[2]} "
                "0. 0. 0.",
                "centerofmass" if center_solute else "",
                "end structure",
                "",
            ],
        )

    # Add a section for each type of molecule to add.
    for file_name, count in zip(molecule_file_names, molecule_counts):
        input_lines.extend(
            [
                f"structure {file_name}",
                f"  number {count}",
                f"  inside box 0. 0. 0. {box_size[0]} {box_size[1]} {box_size[2]}",
                "end structure",
                "",
            ],
        )

    packmol_input = "\n".join(input_lines)

    # Write packmol input
    packmol_file_name = "packmol_input.txt"

    with open(packmol_file_name, "w") as file_handle:
        file_handle.write(packmol_input)

    return packmol_file_name


def _correct_packmol_output(
    file_path: str,
    molecule_topologies: list[mdtraj.Topology],
    number_of_copies: list[int],
    structure_to_solvate: Optional[Topology] = None,
) -> mdtraj.Trajectory:
    """
    Corrects the PDB file output by packmol, namely by adding full connectivity information.

    Parameters
    ----------
    file_path: str
        The file path to the packmol output file.
    molecule_topologies: list of mdtraj.Topology
        A list of topologies for the molecules which packmol has
        added.
    number_of_copies: list of int
        The total number of each molecule which packmol should have
        created.
    structure_to_solvate: Topology, optional
        The topology that was solvated.

    Returns
    -------
    mdtraj.Trajectory
        A trajectory containing the packed system with full connectivity.

    """
    import mdtraj

    with warnings.catch_warnings():
        if structure_to_solvate is not None:
            # Catch the known warning which is fixed in the next section.
            warnings.filterwarnings(
                "ignore",
                message="WARNING: two consecutive residues with same number",
            )

        trajectory = mdtraj.load(file_path)

    all_topologies = []
    all_n_copies = []

    if structure_to_solvate is not None:
        structure_to_solvate.to_file("structure.pdb")
        solvated_trajectory = mdtraj.load("structure.pdb")

        all_topologies.append(solvated_trajectory.topology)
        all_n_copies.append(1)

        # We have to split the topology to ensure the structure to solvate
        # ends up in its own chain.
        n_solvent_atoms = trajectory.n_atoms - solvated_trajectory.n_atoms
        solvent_indices = numpy.arange(n_solvent_atoms) + solvated_trajectory.n_atoms

        solvent_topology = trajectory.topology.subset(solvent_indices)

        full_topology = solvated_trajectory.topology.join(solvent_topology)
        trajectory.topology = full_topology

    all_topologies.extend(molecule_topologies)
    all_n_copies.extend(number_of_copies)

    all_bonds = []
    offset = 0

    for molecule_topology, count in zip(all_topologies, all_n_copies):
        _, molecule_bonds = molecule_topology.to_dataframe()

        for _ in range(count):
            for bond in molecule_bonds:
                all_bonds.append(
                    [int(bond[0].item()) + offset, int(bond[1].item()) + offset],
                )

            offset += molecule_topology.n_atoms

    if len(all_bonds) > 0:
        all_bonds = numpy.unique(all_bonds, axis=0).tolist()

    # We have to check whether there are any existing bonds, because mdtraj
    # will sometimes automatically detect some based on residue names (e.g HOH),
    # and this behaviour cannot be disabled.
    existing_bonds = []

    for bond in trajectory.topology.bonds:
        existing_bonds.append(bond)

    for bond in all_bonds:
        atom_a = trajectory.topology.atom(bond[0])
        atom_b = trajectory.topology.atom(bond[1])

        bond_exists = False

        for existing_bond in existing_bonds:
            if (existing_bond.atom1 == atom_a and existing_bond.atom2 == atom_b) or (
                existing_bond.atom2 == atom_a and existing_bond.atom1 == atom_b
            ):
                bond_exists = True
                break

        if bond_exists:
            continue

        trajectory.topology.add_bond(atom_a, atom_b)

    return trajectory


@requires_package("mdtraj")
def pack_box(
    molecules: list[Molecule],
    number_of_copies: list[int],
    structure_to_solvate: Optional[Topology] = None,
    center_solute: bool = True,
    tolerance: Quantity = 2.0 * unit.angstrom,
    box_size: Optional[Quantity] = None,
    mass_density: Optional[Quantity] = None,
    box_aspect_ratio: Optional[list[float]] = None,
    working_directory: Optional[str] = None,
    retain_working_files: bool = False,
    add_box_buffers: bool = True,
) -> Topology:
    """
    Run packmol to generate a box containing a mixture of molecules.

    Parameters
    ----------
    molecules : list of openff.toolkit.topology.Molecule
        The molecules in the system.
    number_of_copies : list of int
        A list of the number of copies of each molecule type, of length
        equal to the length of ``molecules``.
    structure_to_solvate: Topology, optional
        A topology to be solvated.
    center_solute: bool
        If ``True``, the structure to solvate will be placed in the center of
        the simulation box. This option is only applied when
        ``structure_to_solvate`` is set.
    tolerance : openff.units.Quantity
        The minimum spacing between molecules during packing in units of
        distance.
    box_size : openff.units.Quantity, optional
        The size of the box to generate in units compatible with angstroms.
        If ``None``, ``mass_density`` must be provided.
    mass_density : openff.units.Quantity, optional
        Target mass density for final system with units compatible with g / mL.
        If ``None``, ``box_size`` must be provided.
    box_aspect_ratio: list of float, optional
        The aspect ratio of the simulation box, used in conjunction with
        the ``mass_density`` parameter. If none, an isotropic ratio (i.e.
        ``[1.0, 1.0, 1.0]``) is used.
    verbose : bool
        If ``True``, verbose output is written.
    working_directory: str, optional
        The directory in which to generate the temporary working files. If
        ``None``, a temporary one will be created.
    retain_working_files: bool
        If ``True`` all of the working files, such as individual molecule
        coordinate files, will be retained.
    add_box_buffers: bool
        If ``True`` (the default), buffers will be added to the box size to
        help reduce clashes.

    Returns
    -------
    mdtraj.Trajectory
        The packed box encoded in an mdtraj trajectory.
    list of str
        The residue names which were assigned to each of the
        molecules in the `molecules` list.

    Raises
    ------
    PACKMOLRuntimeError
        When packmol fails to execute / converge.

    """
    if mass_density is not None and box_aspect_ratio is None:
        box_aspect_ratio = [1.0, 1.0, 1.0]

    # Make sure packmol can be found.
    packmol_path = _find_packmol()

    if packmol_path is None:
        raise OSError("Packmol not found, cannot run pack_box()")

    # Validate the inputs.
    _validate_inputs(
        molecules,
        number_of_copies,
        structure_to_solvate,
        box_aspect_ratio,
        box_size,
        mass_density,
    )

    unique_molecules = set(molecules)

    if structure_to_solvate is not None:
        for unique_molecule in structure_to_solvate.unique_molecules:
            unique_molecules.add(unique_molecule)

    # Estimate the box_size from mass density if one is not provided.
    if box_size is None:
        box_size = _approximate_box_size_by_density(
            molecules,
            number_of_copies,
            mass_density,
            box_aspect_ratio,  # type: ignore[arg-type]
            box_scaleup_factor=1.1 if add_box_buffers else 1.0,
        )

    # Set up the directory to create the working files in.
    temporary_directory = False

    if working_directory is None:
        working_directory = tempfile.mkdtemp()
        temporary_directory = True

    if len(working_directory) > 0:
        os.makedirs(working_directory, exist_ok=True)

    assigned_residue_names = []

    with temporary_cd(working_directory):
        # Create PDB files for all of the molecules.
        pdb_file_names = []
        mdtraj_topologies = []

        for index, molecule in enumerate(molecules):
            mdtraj_trajectory = _create_trajectory(molecule)

            pdb_file_name = f"{index}.pdb"
            pdb_file_names.append(pdb_file_name)

            mdtraj_trajectory.save_pdb(pdb_file_name)
            mdtraj_topologies.append(mdtraj_trajectory.topology)

            residue_name = mdtraj_trajectory.topology.residue(0).name
            assigned_residue_names.append(residue_name)

        # Generate the input file.
        output_file_name = "packmol_output.pdb"

        input_file_path = _build_input_file(
            pdb_file_names,
            number_of_copies,
            structure_to_solvate,
            center_solute,
            box_size,
            tolerance,
            output_file_name,
        )

        with open(input_file_path) as file_handle:
            result = subprocess.check_output(
                packmol_path,
                stdin=file_handle,
                stderr=subprocess.STDOUT,
            ).decode("utf-8")

            packmol_succeeded = result.find("Success!") > 0

        if not retain_working_files:
            os.unlink(input_file_path)

            for file_path in pdb_file_names:
                os.unlink(file_path)

        if not packmol_succeeded:
            if os.path.isfile(output_file_name):
                os.unlink(output_file_name)

            if temporary_directory and not retain_working_files:
                shutil.rmtree(working_directory)

            raise PACKMOLRuntimeError(result)

        # Add a 2 angstrom buffer to help alleviate PBC issues.
        if add_box_buffers:
            box_size = [
                (x + 2.0 * unit.angstrom).to(unit.nanometer).magnitude for x in box_size
            ]

        # Append missing connect statements to the end of the
        # output file.
        trajectory = _correct_packmol_output(
            output_file_name,
            mdtraj_topologies,
            number_of_copies,
            structure_to_solvate,
        )
        trajectory.unitcell_lengths = box_size
        trajectory.unitcell_angles = numpy.array([90.0] * 3)

        if not retain_working_files:
            os.unlink(output_file_name)

    if temporary_directory and not retain_working_files:
        shutil.rmtree(working_directory)

    topology = Topology.from_openmm(
        openmm_topology=trajectory.topology.to_openmm(),
        unique_molecules=unique_molecules,
        positions=unit.Quantity(trajectory.xyz[0], unit.nanometer),
    )

    # MDTraj's Topology.to_openmm() doesn't set the box vectors
    topology.box_vectors = unit.Quantity(trajectory.unitcell_vectors[0], unit.nanometer)

    return topology
