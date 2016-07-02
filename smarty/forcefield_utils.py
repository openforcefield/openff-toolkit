#!/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
forcefield_utils.py

OpenMM ForceField replacement using SMIRKS-based matching.

AUTHORS

David L. Mobley <dmobley@mobleylab.org>

Based loosely on code from github.com/choderalab/openmoltools, and especially
parts from John Chodera and Kyle Beauchamp.
"""
#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import os
import smarty
from smarty import ForceField
from smarty.utils import get_data_filename
import simtk.openmm
from simtk.openmm import app
from simtk.openmm.app import element as elem
from simtk.openmm.app import Topology
import numpy as np
from openmoltools import system_checker
import copy

import openeye.oechem
import openeye.oeomega
import openeye.oequacpac
from openeye import oechem

from simtk import openmm, unit


#=============================================================================
# UTILITY FUNCTIONS
#=============================================================================


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
    positions = crd.getPositions()
     
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

def compare_system_energies( topology0, topology1, system0, system1, positions0, positions1=None, label0="AMBER system", label1 = "SMIRFF system", verbose = True ):
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
    timestep = 1.0 * unit.femtoseconds
    integrator0 = simtk.openmm.VerletIntegrator( timestep )     
    integrator1 = simtk.openmm.VerletIntegrator( timestep )     

    # Grab second positions
    if positions1 == None:
        positions1 = copy.deepcopy( positions0 ) 

    # Create simulations
    platform = simtk.openmm.Platform.getPlatformByName("Reference")
    simulation0 = app.Simulation( topology0, system0, integrator0, platform = platform ) 
    simulation0.context.setPositions(positions0)
    simulation1 = app.Simulation( topology1, system1, integrator1, platform = platform ) 
    simulation1.context.setPositions(positions1)

   
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
            print("%s: " % key, groups0[key] ) 
        print("Components of %s: \n" % label1 )
        for key in groups1.keys():
            print("%s: " % key, groups1[key] ) 

    # Return
    return groups0, groups1, energy0, energy1 


def compare_molecule_energies( prmtop, crd, forcefield, mol, verbose = True ):
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

    groups0, groups1, energy0, energy1 = compare_system_energies( ambertop, 
               smirfftop, ambersys, smirffsys, amberpos, verbose = verbose )

    return groups0, groups1, energy0, energy1

