from functools import partial
from smarty import ForceField
import smarty
import openeye
import os
from smarty.utils import get_data_filename
from simtk.openmm import app
from simtk.openmm.app import element as elem
from simtk.openmm.app import Topology
from simtk import unit, openmm
import numpy as np
from openmoltools import system_checker
import copy

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
    Generate a System from the given topology, OEMols matching the contents of the topology, and SMIRFF forcefield and check that its energy is finite.

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

def create_system_from_amber( prmtop_filename, crd_filename ):
    """Utility function. Create and return an OpenMM System given a prmtop and
       crd, AMBER format files.
    
    Parameters
    ----------
    prmtop_filename : str (filename)
        Filename of input AMBER format prmtop file
    crd_filename : str (filename)
        Filename of input AMBER format crd file

    Returns
    _______
    topology : OpenMM Topology
    system : OpenMM System
    positions : initial atomic positions (OpenMM)
"""
    
    # Create System object
    prmtop = app.AmberPrmtopFile(prmtop_filename)
    topology = prmtop.topology
    system = prmtop.createSystem(nonbondedMethod = app.NoCutoff, constraints = None, implicitSolvent = None )
    
    # Read coordinates
    crd = app.AmberInpcrdFile( crd_filename )
    positions = pcrd.getPositions()
     
    return (topology, system, positions)

def create_system_from_molecule(forcefield, mol, verbose=False):
    """
    Generate a System from the given OEMol and SMIRFF forcefield, return the resulting System.

    Parameters
    ----------
    forcefield : smarty.ForceField
        SMIRFF forcefield
    mol : oechem.OEMol
        Molecule to test (must have coordinates)


    Returns
    ----------
    topology : OpenMM Topology
    system : OpenMM System
    positions : initial atomic positions (OpenMM)
    """
    # Create system
    from smarty.forcefield import generateTopologyFromOEMol
    topology = generateTopologyFromOEMol(mol)
    system = forcefield.createSystem(topology, [mol], verbose=verbose)
    
    # Get positions
    coordinates = mol.GetCoords()
    natoms = len(coordinates)
    positions = np.zeros([natoms,3], np.float32)
    for index in range(natoms):
        (x,y,z) = coordinates[index]
        positions[index,0] = x
        positions[index,1] = y
        positions[index,2] = z
    positions = unit.Quantity(positions, unit.angstroms)
    
    return topology, system, positions

def compare_system_energies( system0, system1, positions0, positions1=None, 
                                label0="AMBER system", label1 = "SMIRFF system", verbose = True )
    """
    Given two OpenMM systems, check that their energies and component-wise 
    energies are consistent, and return these. The same positions will be used
    for both systems unless a second set of positions is provided.

    Parameters
    ----------
    topology0 : OpenMM Topology
        Topology of first system
    topology1 : OpenMM Topology
        Topology of second system
    system0 : OpenMM System 
        First system for comparison (usually from AMBER)
    system1 : OpenMM System
        Second system for comparison (usually from SMIRFF)
    positions0 : simtk.unit.Quantity wrapped
        Positions to use for energy evaluation comparison
    positions1 (optional) : simtk.unit.Quantity wrapped (optional)
        Positions to use for second OpenMM system; original positions are used 
        if this is not provided
    label0 (optional) : str
        String labeling system0 for output. Default, "AMBER system"
    label1 (optional) : str
        String labeling system1 for output. Default, "SMIRFF system"
    verbose (optional) : bool
        Print out info on energies, True/False (default True)

    Returns
    ----------
        groups0 : dict
            As returned by openmoltools.system_checker.check_energy_groups, 
            a dictionary with keys "bond", "angle", "nb", "torsion" and values
            corresponding to the energies of these components for the first simulation object
        groups1 : dict
            As returned by openmoltools.system_checker.check_energy_groups, 
            a dictionary with keys "bond", "angle", "nb", "torsion" and values
            corresponding to the energies of these components for the second simulation object
        energy0 : simtk.unit.Quantity 
            Energy of first system
        energy1 : simtk.unit.Quantity
            Energy of second system

    TO DO:
        Allow energy extraction/comparison of terms specified by particular
        SMARTS queries i.e. for specific bond, angle, or torsional terms.
    """

    # Create integrator
    timestep = 1.0 * units.femtoseconds
    integrator = simtk.openmm.VerletIntegrator( timestep )     

    # Grab second positions
    if positions1 == None:
        positions1 = copy.deepcopy( positions0 ) 

    # Create simulations
    platform = simtk.openmm.Platform.getPlatformByName("Reference")
    simulation0 = app.Simulation( topology0, system0, integrator, platform = platform ) 
    simulation.context.setPositions(positions0)
    simulation1 = app.Simulation( topology1, system1, integrator, platform = platform ) 
    simulation.context.setPositions(positions1)

   
    # Do energy comparison, print info if desired
    syscheck = system_checker.SystemChecker( simulation0, simulation1 )
    syscheck.check_force_parameters()
    groups0, groups1 = syscheck.check_energy_groups()
    energy0, energy1 = syscheck.check_energies()   
    if verbose:
        print("Energy of %s: " % label0, energy0 )
        print("Energy of %s: " % label1, energy1 )
        print("Components of %s: \n" % label0 )
        for key in groups0.keys():
            print("%s: ", groups0[key] ) 
        print("Components of %s: \n" % label1 )
        for key in groups1.keys():
            print("%s: ", groups1[key] ) 

    # Return
    return groups0, groups1, energy0, energy1 


def compare_molecule_energies( prmtop, crd, forcefield, mol ):
    """
    Compare energies for OpenMM Systems/topologies created from an AMBER prmtop
    and crd versus from a SMIRFF forcefield file and OEMol which should
    parameterize the same system with same parameters.

    
    Parameters
    ----------
    prmtop_filename : str (filename)
        Filename of input AMBER format prmtop file
    crd_filename : str (filename)
        Filename of input AMBER format crd file
    forcefield : smarty.ForceField
        SMIRFF forcefield
    mol : oechem.OEMol
        Molecule to test
    verbose (optional): Bool
        Print out info. Default: True

    Returns
    --------
        groups0 : dict
            As returned by openmoltools.system_checker.check_energy_groups, 
            a dictionary with keys "bond", "angle", "nb", "torsion" and values
            corresponding to the energies of these components for the first simulation object
        groups1 : dict
            As returned by openmoltools.system_checker.check_energy_groups, 
            a dictionary with keys "bond", "angle", "nb", "torsion" and values
            corresponding to the energies of these components for the second simulation object
        energy0 : simtk.unit.Quantity 
            Energy of first system
        energy1 : simtk.unit.Quantity
            Energy of second system
    """ 

    ambertop, ambersys, amberpos = create_system_from_amber( prmtop, crd )
    smirfftop, smirffsys, smirffpos = create_system_from_molecule(forcefield, mol)

    groups0, groups1, energy0, energy1 = compare_system_energies( system0, 
                    system1, amberpos, verbose = verbose )

    return groups0, groups1, energy0, energy1

if __name__ == '__main__':
    #test_smirks()
    test_create_system_boxes(verbose=True)
