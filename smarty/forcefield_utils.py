#!/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
forcefield_utils.py

Utilities relating to OpenMM ForceField replacement using SMIRKS-based matching.

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
import simtk.openmm as mm
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
import parmed

#=============================================================================
# UTILITY FUNCTIONS
#=============================================================================


def create_system_from_amber( prmtop_filename, crd_filename, verbose = False ):
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

def compare_system_energies( topology0, topology1, system0, system1, positions0, positions1=None, label0="AMBER system", label1 = "SMIRFF system", verbose = True, skip_assert = False, skip_improper = False ):
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
    skip_assert (optional) : bool
        Skip assertion that energies must be equal within specified tolerance. Default False.
    skip_improper (optional) : bool
        Skip detailed checking of force terms on impropers (helpful here if comparing with AMBER force fields using different definitions of impropers.) Default False.

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

    # Print what torsions were found if verbose
    if verbose:
        # Build list of atoms for debugging info
        atoms0 = [ atom for atom in simulation0.topology.atoms() ]
        atoms1 = [ atom for atom in simulation1.topology.atoms() ]
        # Loop over first system and print torsion info
        for force in simulation0.system.getForces():
            if type(force) == mm.PeriodicTorsionForce:
                print("Num (type) \t Num (type) \t Num (type) \t Num (type) \t per \t phase \t k0")
                for k in range(force.getNumTorsions()):
                    i0, i1, i2, i3, per, phase, k0 = force.getTorsionParameters(k)
                    print("%3s (%3s)- %3s (%3s)- \t %s (%3s)- \t %3s (%3s)- \t %f \t %f \t %f " % (i0, atoms0[i0].name, i1, atoms0[i1].name, i2, atoms0[i2].name, i3, atoms0[i3].name, per, phase/unit.degree, k0/unit.kilojoule_per_mole) )
        for force in simulation1.system.getForces():
            if type(force) == mm.PeriodicTorsionForce:
                print("Num (type) \t Num (type) \t Num (type) \t Num (type) \t per \t phase \t k0")
                for k in range(force.getNumTorsions()):
                    i0, i1, i2, i3, per, phase, k0 = force.getTorsionParameters(k)
                    print("%3s (%3s)- %3s (%3s)- %3s (%3s)- %3s (%3s) - %f \t %f \t %f " % (i0, atoms1[i0].name, i1, atoms1[i1].name, i2, atoms1[i2].name, i3, atoms1[i3].name, per, phase/unit.degree, k0/unit.kilojoule_per_mole) )




    # Do energy comparison, print info if desired
    syscheck = system_checker.SystemChecker( simulation0, simulation1 )
    if not skip_assert:
        # Only check force terms if we want to make sure energies are identical
        syscheck.check_force_parameters(skipImpropers = skip_improper)
    groups0, groups1 = syscheck.check_energy_groups(skip_assert = skip_assert)
    energy0, energy1 = syscheck.check_energies(skip_assert = skip_assert)
    if verbose:
        print("Energy of %s: " % label0, energy0 )
        print("Energy of %s: " % label1, energy1 )
        print("\nComponents of %s:" % label0 )
        for key in groups0.keys():
            print("%s: " % key, groups0[key] )
        print("\nComponents of %s:" % label1 )
        for key in groups1.keys():
            print("%s: " % key, groups1[key] )

    # Return
    return groups0, groups1, energy0, energy1


def compare_molecule_energies( prmtop, crd, forcefield, mol, verbose = True, skip_assert=False, skip_improper = False):
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
    skip_assert : bool
        Skip assertion that energies must be equal within tolerance. Default, False.
    skip_improper (optional) : bool
        Skip detailed checking of force terms on impropers (helpful here if comparing with AMBER force fields using different definitions of impropers.) Default False.


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
    smirfftop, smirffsys, smirffpos = create_system_from_molecule(forcefield, mol, verbose = verbose)

    groups0, groups1, energy0, energy1 = compare_system_energies( ambertop,
               smirfftop, ambersys, smirffsys, amberpos, verbose = verbose, skip_assert = skip_assert, skip_improper = skip_improper )

    return groups0, groups1, energy0, energy1


def get_molecule_parameterIDs( oemols, ffxml):
    """Process a list of oemols with a specified SMIRFF ffxml file and determine which parameters are used by which molecules, returning collated results.


    Parameters
    ----------
    oemols : list
        List of OpenEye OEChem molecules to parse; must have explicit hydrogens.

    Returns
    -------
    parameters_by_molecule : dict
        Parameter IDs used in each molecule, keyed by isomeric SMILES
        generated from provided OEMols. Each entry in the dict is a list
        which does not necessarily have unique entries; i.e. parameter IDs
        which are used more than once will occur multiple times.

    parameters_by_ID : dict
        Molecules in which each parameter ID occur, keyed by parameter ID.
        Each entry in the dict is a set of isomeric SMILES for molecules
        in which that parameter occurs. No frequency information is stored.

    """

    # Create storage
    parameters_by_molecule = {}
    parameters_by_ID = {}

    # Generate isomeric SMILES
    isosmiles = list()
    for mol in oemols:
        smi = oechem.OECreateIsoSmiString(mol)
        if not smi in isosmiles:
            isosmiles.append(smi)
        # If the molecule is already here, raise exception
        else:
            raise ValueError("Error: get_molecule_parameterIDs has been provided a list of oemols which contains the same molecule, having isomeric smiles %s, more than once." % smi )
    # Label molecules
    ff = ForceField( ffxml )
    labels = ff.labelMolecules( oemols )

    # Organize labels into output dictionary by looping over all molecules/smiles
    for idx in range(len(isosmiles)):
        # Pull smiles, initialize storage
        smi = isosmiles[idx]
        parameters_by_molecule[smi] = []

        # Organize data for this molecule
        data = labels[idx]
        for force_type in data.keys():
            for (atom_indices, pid, smirks) in data[force_type]:
                # Store pid to molecule
                parameters_by_molecule[smi].append(pid)

                # Store which molecule this pid occurred in
                if pid not in parameters_by_ID:
                    parameters_by_ID[pid] = set()
                    parameters_by_ID[pid].add(smi)
                else:
                    parameters_by_ID[pid].add(smi)

    return parameters_by_molecule, parameters_by_ID

def getMolParamIDToAtomIndex( oemol, ff):
    """Take an OEMol and a SMIRFF forcefield object and return a dictionary, keyed by parameter ID, where each entry is a tuple of ( smirks, [[atom1, ... atomN], [atom1, ... atomN]) giving the SMIRKS corresponding to that parameter ID and a list of the atom groups in that molecule that parameter is applied to.

    Parameters
    ----------
    oemol : OEMol
        OpenEye OEMol with the molecule to investigate.
    ff : ForceField
        SMIRFF ForceField object (obtained from an ffxml via ForceField(ffxml)) containing FF of interest.

    Returns
    -------
    param_usage : dictionary
        Dictionary, keyed by parameter ID, where each entry is a tuple of ( smirks, [[atom1, ... atomN], [atom1, ... atomN]) giving the SMIRKS corresponding to that parameter ID and a list of the atom groups in that molecule that parameter is applied to.

    """

    labels = ff.labelMolecules([oemol])
    param_usage = {}
    for mol_entry in range(len(labels)):
        for force in labels[mol_entry].keys():
            for (atom_indices, pid, smirks) in labels[mol_entry][force]:
                if not pid in param_usage:
                    param_usage[pid] = (smirks, [atom_indices])
                else:
                    param_usage[pid][1].append( atom_indices )

    return param_usage


def save_system_to_amber( topology, system, positions, prmtop, crd ):
    """Save an OpenMM System, with provided topology and positions, to AMBER prmtop and coordinate files.

    Parameters
    ----------
    topology : OpenMM Topology
        Topology of the system to be saved, perhaps as loaded from a PDB file or similar.
    system : OpenMM System
        Parameterized System to be saved, containing components represented by Topology
    positions : unit.Quantity position array
        Position array containing positions of atoms in topology/system
    prmtop : filename
        AMBER parameter file name to write
    crd : filename
        AMBER coordinate file name (ASCII crd format) to write

    """

    structure = parmed.openmm.topsystem.load_topology( topology, system, positions )
    structure.save( prmtop, overwrite = True, format="amber" )
    structure.save( crd, format='rst7', overwrite = True)


def save_system_to_gromacs( topology, system, positions, top, gro ):
    """Save an OpenMM System, with provided topology and positions, to AMBER prmtop and coordinate files.

    Parameters
    ----------
    topology : OpenMM Topology
        Topology of the system to be saved, perhaps as loaded from a PDB file or similar.
    system : OpenMM System
        Parameterized System to be saved, containing components represented by Topology
    positions : unit.Quantity position array
        Position array containing positions of atoms in topology/system
    top : filename
        GROMACS topology file name to write
    gro : filename
        GROMACS coordinate file name (.gro format) to write

    """

    structure = parmed.openmm.topsystem.load_topology( topology, system, positions )
    structure.save( prmtop, overwrite = True, format="gromacs")
    structure.save( gro, overwrite = True, format="gro")
