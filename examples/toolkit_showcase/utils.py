from tempfile import NamedTemporaryFile
from simtk import openmm, unit
from simtk.openmm import app
import nglview
import mdtraj as mdt


def find_clashing_water(pmd_struct, lig_resname, distance=0.15):
    """
    Find waters that are sterically clashing with a ligand.

    Parameters
    ----------
    pmd_struct : parmed.Structure
        The structure to analyze.
    lig_resname : str
        The up-to-three character residue name.
    distance : float
        The distance cutoff (in nanometers) for clash detection.

    Returns
    -------
    water_resnums : Iterable[int]
        The residue numbers of waters that are clashing with the ligand.

    """
    with NamedTemporaryFile(suffix=".pdb") as tf:
        app.PDBFile.writeFile(
            pmd_struct.topology, pmd_struct.positions, open(tf.name, "w")
        )
        traj = mdt.load(tf.name)
    top = traj.topology
    lig_atom_idxs = top.select(f"resname {lig_resname}")
    lig_res_idx = top.atom(lig_atom_idxs[1]).residue.index
    wat_atom_idxs = top.select("resname HOH and name O")
    wat_res_idxs = [top.atom(i).residue.index for i in wat_atom_idxs]
    potential_contacts = [(lig_res_idx, wat_res_idx) for wat_res_idx in wat_res_idxs]
    contacts = mdt.compute_contacts(
        traj, contacts=potential_contacts, scheme="closest", ignore_nonprotein=False
    )

    # Note that this is 0-indexed, while the parmed structure is
    # 1-indexed, therefore we add 1 before returning
    clash_res_idx = [i[1] + 1 for i in contacts[1][(contacts[0] < distance)[0, :]]]
    return clash_res_idx


def minimize_and_visualize(molecule, forcefield):
    """
    Minimize the molecule with the force field and then visualize it with nglview

    Parameters
    ----------
    molecule : openff.toolkit.topology.molecule.Molecule
        The molecule to minimize
    forcefield : openff.toolkit.typing.engines.smirnoff.forcefield.ForceField
        The force field to minimize against

    Returns
    -------
    view : nglview.widget.NGLWidget
        A view of the molecule
    """
    mol_topology = molecule.to_topology()
    mol_system = forcefield.create_openmm_system(
        mol_topology, charge_from_molecules=[molecule]
    )

    integrator = openmm.LangevinIntegrator(
        300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds
    )
    simulation = openmm.app.Simulation(mol_topology.to_openmm(), mol_system, integrator)
    simulation.context.setPositions(molecule.conformers[0])
    simulation.minimizeEnergy()
    with NamedTemporaryFile(suffix=".pdb") as tf:
        openmm.app.PDBFile.writeModel(
            simulation.topology,
            simulation.context.getState(getPositions=True).getPositions(),
            open(tf.name, "w"),
        )
        view = nglview.show_file(tf.name)
    return view
