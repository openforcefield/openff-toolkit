from functools import partial
from smarty import ForceField
import smarty
import openeye
import os
from smarty.utils import get_data_filename
from simtk.openmm.app import element as elem
from simtk.openmm.app import Topology
from simtk import unit, openmm
import numpy as np

def positions_from_oemol(mol):
    """
    Extract OpenMM positions from OEMol.

    Parameters
    ----------
    mol : oechem.openeye.OEMol
        OpenEye molecule from which to extract coordinates.

    Returns
    -------
    positions : simtk.unit.Quantity of dimension (nparticles,3)

    """
    if mol.GetDimension() != 3:
        # Assign coordinates
        omega = oeomega.OEOmega()
        omega.SetMaxConfs(1)
        omega.SetIncludeInput(False)
        omega.SetCanonOrder(False)
        omega.SetSampleHydrogens(True)  # Word to the wise: skipping this step can lead to significantly different charges!
        status = omega(mol)  # generate conformation

    coordinates = mol.GetCoords()
    natoms = len(coordinates)
    positions = np.zeros([natoms,3], np.float32)
    for index in range(natoms):
        (x,y,z) = coordinates[index]
        positions[index,0] = x
        positions[index,1] = y
        positions[index,2] = z
    positions = unit.Quantity(positions, unit.angstroms)
    return positions

def check_energy_is_finite(system, positions):
    """
    Check the potential energy is not NaN.

    Parameters
    ----------
    system : simtk.openmm.System
        The system to check
    positions : simtk.unit.Quantity of dimension (natoms,3) with units of length
        The positions to use

    """
    integrator = openmm.VerletIntegrator(1.0 * unit.femtoseconds)
    context = openmm.Context(system, integrator)
    context.setPositions(positions)
    state = context.getState(getEnergy=True)
    energy = state.getPotentialEnergy() / unit.kilocalories_per_mole
    if np.isnan(energy):
        raise Exception('Potential energy is NaN')

def test_read_ffxml():
    """Test reading of ffxml files.
    """
    forcefield = ForceField(get_data_filename('forcefield/Frosst_AlkEtOH.ffxml'))

def check_system_creation_from_molecule(forcefield, mol, verbose=False):
    """
    Generate a System from the given OEMol and SMIRFF forcefield and check that its energy is finite.

    Parameters
    ----------
    forcefield : smarty.ForceField
        SMIRFF forcefield
    mol : oechem.OEMol
        Molecule to test (need not have coordinates)

    """

    from smarty.forcefield import generateTopologyFromOEMol
    topology = generateTopologyFromOEMol(mol)
    system = forcefield.createSystem(topology, [mol], verbose=verbose)
    # Test energy computation.
    positions = positions_from_oemol(mol)
    check_energy_is_finite(system, positions)

def check_system_creation_from_topology(forcefield, topology, mols, positions, verbose=False):
    """
    Generate a System from the given OEMol and SMIRFF forcefield and check that its energy is finite.

    Parameters
    ----------
    forcefield : smarty.ForceField
        SMIRFF forcefield
    topology : simtk.openmm.app.Topology
        Topology to construct system from
    mols : list of oechem.OEMol
        Reference molecules
    positions : simtk.unit.Quantity with dimension (natoms,3) with units of length

    """
    from smarty.forcefield import CutoffPeriodic
    system = forcefield.createSystem(topology, mols, verbose=verbose, nonbondedMethod=CutoffPeriodic)
    # Test energy computation.
    check_energy_is_finite(system, positions)

def test_create_system_molecules(verbose=False):
    """Test creation of a System object from small molecules.
    """
    forcefield = ForceField(get_data_filename('forcefield/Frosst_AlkEtOH.ffxml'))

    from openeye import oechem
    ifs = oechem.oemolistream(get_data_filename('molecules/AlkEtOH-tripos.mol2.gz'))
    mol = oechem.OEGraphMol()
    while oechem.OEReadMolecule(ifs, mol):
        args = { 'verbose' : verbose }
        f = partial(check_system_creation_from_molecule, forcefield, mol, **args)
        f.description ='Testing creation of system object from small molecules (%s)' % mol.GetTitle()
        yield f

def test_create_system_boxes(verbose=False):
    """Test creation of a System object from some boxes of mixed solvents.
    """
    forcefield = ForceField(get_data_filename('forcefield/Frosst_AlkEtOH.ffxml'))

    # Read monomers
    mols = list()
    monomers = ['cyclohexane', 'ethanol', 'propane', 'methane', 'butanol']
    from openeye import oechem
    mol = oechem.OEGraphMol()
    for monomer in monomers:
        filename = get_data_filename(os.path.join('systems', 'monomers', monomer + '.sdf'))
        ifs = oechem.oemolistream(filename)
        while oechem.OEReadMolecule(ifs, mol):
            oechem.OETriposAtomNames(mol)
            mols.append( oechem.OEGraphMol(mol) )
    print('%d reference molecules loaded' % len(mols))

    # Read systems.
    boxes = ['cyclohexane_ethanol_0.4_0.6.pdb', 'propane_methane_butanol_0.2_0.3_0.5.pdb']
    from simtk.openmm.app import PDBFile
    for box in boxes:
        filename = get_data_filename(os.path.join('systems', 'packmol_boxes', box))
        pdbfile = PDBFile(filename)
        f = partial(check_system_creation_from_topology, forcefield, pdbfile.topology, mols, pdbfile.positions, verbose=verbose)
        f.description = 'Test creation of System object from %s' % box
        yield f

if __name__ == '__main__':
    #test_smirks()
    test_create_system_boxes(verbose=True)
