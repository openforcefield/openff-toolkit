#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Utilities for testing

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

from openforcefield.utils import get_data_filename #, generateTopologyFromOEMol, read_molecules
#from openforcefield.utils import check_energy_is_finite, get_energy

from simtk import unit, openmm
from simtk.openmm import app

#=============================================================================================
# UTILITIES
#=============================================================================================


def get_amber_system(prefix='cyclohexane_ethanol_0.4_0.6'):
    """Get AMBER prmtop and inpcrd test data filenames

    Parameters
    ----------
    prefix : str, optional, default='cyclohexane_ethanol_0.4_0.6'
        The prefix of .prmtop and .inpcrd files to retrieve from testdata/systems/amber

    Returns
    -------
    prmtop_filename : str
        Absolute path to the AMBER prmtop filename in testdata/systems/amber
    inpcrd_filename : str
        Absolute path to the AMBER inpcrd filename in testdata/systems/amber
    """
    prefix = os.path.join('systems', 'amber', prefix)
    prmtop_filename = get_data_filename(prefix+'.prmtop')
    inpcrd_filename = get_data_filename(prefix+'.inpcrd')
    return prmtop_filename, inpcrd_filename

def get_packmol_pdbfile(prefix='cyclohexane_ethanol_0.4_0.6'):
    """Get PDB filename for a packmol-generated box

    Parameters
    ----------
    prefix : str, optional, default='cyclohexane_ethanol_0.4_0.6'
        The prefix of .pdb file to retrieve from testdata/systems/packmol_boxes

    Returns
    -------
    pdb_filename : str
        Absolute path to the PDB file
    """
    prefix = os.path.join('systems', 'packmol_boxes', prefix)
    pdb_filename = get_data_filename(prefix+'.pdb')
    return pdb_filename

def get_monomer_mol2file(prefix='ethanol'):
    """Get absolute filepath for a mol2 file denoting a small molecule monomer in testdata

    Parameters
    ----------
    prefix : str, optional, default='ethanol'
        The prefix of .mol2 file to retrieve from systems/monomers/

    Returns
    -------
    mol2_filename : str
        Absolute path to the mol2 file
    """
    prefix = os.path.join('systems', 'monomers', prefix)
    mol2_filename = get_data_filename(prefix+'.pdb')
    return mol2_filename

def compare_system_energies(self, topology0, topology1, system0, system1, positions0, positions1=None, label0="AMBER system", label1="SMIRNOFF system", verbose=True, skip_assert=False, skip_improper=False):
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
        Second system for comparison (usually from SMIRNOFF)
    positions0 : simtk.unit.Quantity wrapped
        Positions to use for energy evaluation comparison
    positions1 (optional) : simtk.unit.Quantity wrapped (optional)
        Positions to use for second OpenMM system; original positions are used
        if this is not provided
    label0 (optional) : str
        String labeling system0 for output. Default, "AMBER system"
    label1 (optional) : str
        String labeling system1 for output. Default, "SMIRNOFF system"
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
    and crd versus from a SMIRNOFF forcefield file and OEMol which should
    parameterize the same system with same parameters.
    Parameters
    ----------
    prmtop_filename : str (filename)
        Filename of input AMBER format prmtop file
    crd_filename : str (filename)
        Filename of input AMBER format crd file
    forcefield : ForceField
        SMIRNOFF forcefield
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
