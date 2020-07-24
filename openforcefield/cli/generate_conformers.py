import argparse

from openforcefield.typing.engines.smirnoff import ForceField
from openforcefield.topology import Molecule, Topology
from openforcefield.utils.toolkits import ToolkitRegistry

# from rdkit.Chem import rdMolAlign
# import numpy as np
# from simtk import openmm, unit


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
    except:
        print("catch speifics here")

    mols = toolkit_registry.call(
        "from_file", molecule, file_format=molecule.split(".")[-1], allow_undefined_stereo=True,
    )

    if type(mols) != list:
        mols = [mols]
    """
    The conformers should be output as separate SDFs with SD data pairs for
    absolute_energy in kcal/mol
    conformer generation method (openeye version or rdkit version)
    partial charges (so, use ForceField.to_openmm_system's return_topology kwarg, and attach the final coordinates to ret_topology._reference_molecules[0], then save that to file)
    """
    for mol in mols:
        toolkit_registry.call('generate_conformers', molecule=mol, n_conformers=5, rms_cutoff=None)
        top = Topology.from_molecules(mol)
        omm_sys, new_top = ff.create_openmm_system(top, return_topology=True)

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
