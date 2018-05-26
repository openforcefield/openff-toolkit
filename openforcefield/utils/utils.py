#!/usr/bin/env python

"""
Utility subroutines for open forcefield tools

"""
#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import sys
import string

from optparse import OptionParser # For parsing of command line arguments

import os
import math
import copy
import re
import numpy
import random
import parmed

from simtk import unit, openmm
from simtk.openmm import app
from simtk.openmm.app import element as elem
from simtk.openmm.app import Topology

from openmoltools import system_checker

import time
from simtk import unit

#=============================================================================================
# UTILITY ROUTINES
#=============================================================================================

def all_subclasses(cls):
    return cls.__subclasses__() + [ g for s in cls.__subclasses__() for g in all_subclasses(s) ]

def checkCharges(molecule):
    # Check that molecule is charged.
    #is_charged = False
    for atom in molecule.GetAtoms():
        if atom.GetPartialCharge() != 0.0:
            return True
        else:
            print('WARNING: Molecule %s has no charges; input molecules must be charged.' % molecule.GetTitle())
            return False

def generateSMIRNOFFStructure(molecule):
    """
    Given an OpenEye molecule (oechem.OEMol), create an OpenMM System and use to
    generate a ParmEd structure using the SMIRNOFF forcefield parameters.
    """
    from openforcefield.typing.engines.smirnoff import ForceField
    from openforcefield.typing.engines.smirnoff.forcefield_utils import create_system_from_molecule

    ff = get_data_filename('forcefield/smirnoff99Frosst.offxml')
    with open(ff) as ffxml:
        mol_ff = ForceField(ffxml)

    if not checkCharges(molecule):
        from openmoltools.openeye import get_charges
        print("Assigning charges to molecule.")
        charged_molecule = get_charges(molecule)
    else:
        charged_molecule = molecule
    mol_top, mol_sys, mol_pos = create_system_from_molecule(mol_ff, charged_molecule)
    molecule_structure = parmed.openmm.load_topology(mol_top, mol_sys, xyz=mol_pos)

    return molecule_structure

def generateProteinStructure(proteinpdb, protein_forcefield='amber99sbildn.xml', solvent_forcefield='tip3p.xml'):
    """
    Given an OpenMM PDBFile, create the OpenMM System of the protein and
    then generate the parametrized ParmEd Structure of the protein.
    Parameters
    ----------
    proteinpdb : openmm.app.PDBFile object,
        Loaded PDBFile object of the protein.
    protein_forcefield : xml file, default='amber99sbildn.xml'
        Forcefield parameters for protein
    solvent_forcefield : xml file, default='tip3p.xml'
        Forcefield parameters for solvent
    Returns
    -------
    solv_structure : parmed.structure.Structure
        The parameterized Structure of the protein with solvent molecules. (No ligand).
    """
    #Generate protein Structure object
    forcefield = app.ForceField(protein_forcefield, solvent_forcefield)
    protein_system = forcefield.createSystem( proteinpdb.topology )
    protein_structure = parmed.openmm.load_topology(proteinpdb.topology,
                                                    protein_system,
                                                    xyz=proteinpdb.positions)
    return protein_structure

def combinePostions(proteinPositions, molPositions):
    """
    Loops through the positions from the ParmEd structures of the protein and ligand,
    divides by unit.angstroms which will ensure both positions arrays are in the same units.
    Parameters
    ----------
    proteinPositions : list of 3-element Quantity tuples.
        Positions list taken directly from the protein Structure.
    molPositions : list of 3-element Quantity tuples.
        Positions list taken directly from the molecule Structure.
    Returns
    -------
    positions : list of 3-element Quantity tuples.
        ex. unit.Quantity(positions, positions_unit)
        Combined positions of the protein and molecule Structures.
    """
    positions_unit = unit.angstroms
    positions0_dimensionless = numpy.array(proteinPositions / positions_unit)
    positions1_dimensionless = numpy.array(molPositions / positions_unit)
    coordinates = numpy.vstack(
        (positions0_dimensionless, positions1_dimensionless))
    natoms = len(coordinates)
    positions = numpy.zeros([natoms, 3], numpy.float32)
    for index in range(natoms):
            (x, y, z) = coordinates[index]
            positions[index, 0] = x
            positions[index, 1] = y
            positions[index, 2] = z
    positions = unit.Quantity(positions, positions_unit)
    return positions

def mergeStructure(proteinStructure, molStructure):
    """
    Combines the parametrized ParmEd structures of the protein and ligand to
    create the Structure for the protein:ligand complex, while retaining the SMIRNOFF
    parameters on the ligand. Preserves positions and box vectors.
    (Not as easily achieved using native OpenMM tools).
    Parameters
    ----------
    proteinStructure : parmed.structure.Structure
        The parametrized structure of the protein.
    moleculeStructure : parmed.structure.Structure
        The parametrized structure of the ligand.
    Returns
    -------
    structure : parmed.structure.Structure
        The parametrized structure of the protein:ligand complex.
    """
    structure = proteinStructure + molStructure
    positions = combinePostions(proteinStructure.positions, molStructure.positions)
    # Concatenate positions arrays (ensures same units)
    structure.positions = positions
    # Restore original box vectors
    structure.box = proteinStructure.box
    return structure


def generateTopologyFromOEMol(molecule):
    """
    Generate an OpenMM Topology object from an OEMol molecule.

    Parameters
    ----------
    molecule : openeye.oechem.OEMol
        The molecule from which a Topology object is to be generated.

    Returns
    -------
    topology : simtk.openmm.app.Topology
        The Topology object generated from `molecule`.

    """
    from openeye import oechem

    # Avoid manipulating the molecule
    mol = oechem.OEMol(molecule)

    # Create a Topology object with one Chain and one Residue.
    from simtk.openmm.app import Topology
    topology = Topology()
    chain = topology.addChain()
    resname = mol.GetTitle()
    residue = topology.addResidue(resname, chain)

    # Make sure the atoms have names, otherwise bonds won't be created properly below
    if any([atom.GetName() =='' for atom in mol.GetAtoms()]):
        oechem.OETriposAtomNames(mol)
    # Check names are unique; non-unique names will also cause a problem
    atomnames = [ atom.GetName() for atom in mol.GetAtoms() ]
    if any( atomnames.count(atom.GetName())>1 for atom in mol.GetAtoms()):
        raise Exception("Error: Reference molecule must have unique atom names in order to create a Topology.")

    # Create atoms in the residue.
    for atom in mol.GetAtoms():
        name = atom.GetName()
        element = elem.Element.getByAtomicNumber(atom.GetAtomicNum())
        openmm_atom = topology.addAtom(name, element, residue)

    # Create bonds.
    atoms = { atom.name : atom for atom in topology.atoms() }
    for bond in mol.GetBonds():
        aromatic = None
        if bond.IsAromatic(): aromatic = 'Aromatic'
        # Add bond, preserving order assessed by OEChem
        topology.addBond(atoms[bond.GetBgn().GetName()], atoms[bond.GetEnd().GetName()], type=aromatic, order=bond.GetOrder())

    return topology

def get_data_filename(relative_path):
    """Get the full path to one of the reference files in testsystems.

    In the source distribution, these files are in ``openforcefield/data/``,
    but on installation, they're moved to somewhere in the user's python
    site-packages directory.

    Parameters
    ----------
    name : str
        Name of the file to load (with respect to the repex folder).

    """

    from pkg_resources import resource_filename
    fn = resource_filename('openforcefield', os.path.join('data', relative_path))

    if not os.path.exists(fn):
        raise ValueError("Sorry! %s does not exist. If you just added it, you'll have to re-install" % fn)

    return fn

def normalize_molecules(molecules):
    """
    Normalize all molecules in specified set.

    ARGUMENTS

    molecules (list of OEMol) - molecules to be normalized (in place)

    """
    from openeye import oechem

    # Add explicit hydrogens.
    for molecule in molecules:
        oechem.OEAddExplicitHydrogens(molecule)

    # Build a conformation for all molecules with Omega.
    print("Building conformations for all molecules...")
    from openeye import oeomega
    omega = oeomega.OEOmega()
    omega.SetMaxConfs(1)
    omega.SetFromCT(True)
    for molecule in molecules:
        #omega.SetFixMol(molecule)
        omega(molecule)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print("%.3f s elapsed" % elapsed_time)

    # Regularize all molecules through writing as mol2.
    print("Regularizing all molecules...")
    ligand_mol2_dirname  = os.path.dirname(mcmcDbName) + '/mol2'
    if( not os.path.exists( ligand_mol2_dirname ) ):
        os.makedirs(ligand_mol2_dirname)
    ligand_mol2_filename = ligand_mol2_dirname + '/temp' + os.path.basename(mcmcDbName) + '.mol2'
    start_time = time.time()
    omolstream = oechem.oemolostream(ligand_mol2_filename)
    for molecule in molecules:
        # Write molecule as mol2, changing molecule through normalization.
        oechem.OEWriteMolecule(omolstream, molecule)
    omolstream.close()
    end_time = time.time()
    elapsed_time = end_time - start_time
    print("%.3f s elapsed" % elapsed_time)

    # Assign AM1-BCC charges.
    print("Assigning AM1-BCC charges...")
    start_time = time.time()
    for molecule in molecules:
        # Assign AM1-BCC charges.
        if molecule.NumAtoms() == 1:
            # Use formal charges for ions.
            OEFormalPartialCharges(molecule)
        else:
            # Assign AM1-BCC charges for multiatom molecules.
            OEAssignPartialCharges(molecule, OECharges_AM1BCC, False) # use explicit hydrogens
        # Check to make sure we ended up with partial charges.
        if OEHasPartialCharges(molecule) == False:
            print("No charges on molecule: '%s'" % molecule.GetTitle())
            print("IUPAC name: %s" % OECreateIUPACName(molecule))
            # TODO: Write molecule out
            # Delete themolecule.
            molecules.remove(molecule)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print("%.3f s elapsed" % elapsed_time)
    print("%d molecules remaining" % len(molecules))

    return

def read_molecules(filename, verbose=True):
    """
    Read molecules from an OpenEye-supported file.

    Parameters
    ----------
    filename : str
        Filename from which molecules are to be read (e.g. mol2, sdf)

    Returns
    -------
    molecules : list of OEMol
        List of molecules read from file

    """
    from openeye import oechem

    if not os.path.exists(filename):
        built_in = get_data_filename('molecules/%s' % filename)
        if not os.path.exists(built_in):
            raise Exception("File '%s' not found." % filename)
        filename = built_in

    if verbose: print("Loading molecules from '%s'..." % filename)
    start_time = time.time()
    molecules = list()
    input_molstream = oechem.oemolistream(filename)

    flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield
    input_molstream.SetFlavor(oechem.OEFormat_MOL2, flavor)

    molecule = oechem.OECreateOEGraphMol()
    while oechem.OEReadMolecule(input_molstream, molecule):
        # If molecule has no title, try getting SD 'name' tag
        if molecule.GetTitle() == '':
            name = oechem.OEGetSDData(molecule, 'name').strip()
            molecule.SetTitle(name)
        # Append to list.
        molecule_copy = oechem.OEMol(molecule)
        molecules.append(molecule_copy)
    input_molstream.close()
    if verbose: print("%d molecules read" % len(molecules))
    end_time = time.time()
    elapsed_time = end_time - start_time
    if verbose: print("%.3f s elapsed" % elapsed_time)

    return molecules

def setPositionsInOEMol(molecule, positions):
    """Set the positions in an OEMol using a position array with units from simtk.unit, i.e. from OpenMM. Atoms must have same order.

    Arguments:
    ---------
    molecule : OEMol
        OpenEye molecule
    positions : Nx3 array
        Unit-bearing via simtk.unit Nx3 array of coordinates
    """
    from openeye import oechem

    if molecule.NumAtoms() != len(positions): raise ValueError("Number of atoms in molecule does not match length of position array.")
    pos_unitless = positions/unit.angstroms

    coordlist = []
    for idx in range(len(pos_unitless)):
        for j in range(3):
            coordlist.append( pos_unitless[idx][j])
    molecule.SetCoords(oechem.OEFloatArray(coordlist))

def extractPositionsFromOEMol(molecule):
    """Get the positions from an OEMol and return in a position array with units via simtk.unit, i.e. foramtted for OpenMM.
    Adapted from choderalab/openmoltools test function extractPositionsFromOEMOL

    Arguments:
    ----------
    molecule : OEMol
        OpenEye molecule

    Returns:
    --------
    positions : Nx3 array
        Unit-bearing via simtk.unit Nx3 array of coordinates
    """

    positions = unit.Quantity(numpy.zeros([molecule.NumAtoms(), 3], numpy.float32), unit.angstroms)
    coords = molecule.GetCoords()
    for index in range(molecule.NumAtoms()):
        positions[index,:] = unit.Quantity(coords[index], unit.angstroms)
    return positions

def read_typelist(filename):
    """
    Read a parameter type or decorator list from a file.
    Lines in these files have the format
    "SMARTS/SMIRKS  shorthand"
    lines beginning with '%' are ignored

    Parameters
    ----------
    filename : str
        Path and name of file to be read
        Could be file in openforcefield/data/

    Returns
    -------
    typelist : list of tuples
        Typelist[i] is element i of the typelist in format (smarts, shorthand)
    """
    if filename is None:
        return None

    if not os.path.exists(filename):
        built_in = get_data_filename(filename)
        if not os.path.exists(built_in):
            raise Exception("File '%s' not found." % filename)
        filename = built_in

    typelist = list()
    ifs = open(filename)
    lines = ifs.readlines()
    used_typenames = list()

    for line in lines:
        # Strip trailing comments
        index = line.find('%')
        if index != -1:
            line = line[0:index]

        # Split into tokens.
        tokens = line.split()
        # Process if we have enough tokens
        if len(tokens) >= 2:
            smarts = tokens[0]
            typename = ' '.join(tokens[1:])
            if typename not in used_typenames:
                typelist.append([smarts,typename])
                used_typenames.append(typename)
            else:
                raise Exception("Error in file '%s' -- each entry must "
                        "have a unique name." % filename )

    ifs.close()

    return typelist

# TODO: Fix me
def is_openeye_installed(oetools=('oechem', 'oequacpac', 'oeiupac', 'oeomega')):
    """
    Check if a given OpenEye tool is installed and Licensed

    If the OpenEye toolkit is not installed, returns False

    Parameters
    ----------
    oetools : str or iterable of strings, Optional, Default: ('oechem', 'oequacpac', 'oeiupac', 'oeomega')
        Set of tools to check by their string name. Defaults to the complete set that YANK *could* use, depending on
        feature requested.

        Only checks the subset of tools if passed. Also accepts a single tool to check as a string instead of an
        iterable of length 1.

    Returns
    -------
    all_installed : bool
        True if all tools in ``oetools`` are installed and licensed, False otherwise
    """
    # Complete list of module: License check
    tools_license = {
        'oechem': 'OEChemIsLicensed',
        'oequacpac': 'OEQuacPacIsLicensed',
        'oeiupac': 'OEIUPACIsLicensed',
        'oeomega': 'OEOmegaIsLicensed'
        }
    tool_keys = tools_license.keys()
    # Cast oetools to tuple if its a single string
    if type(oetools) is str:
        oetools = (oetools,)
    tool_set = set(oetools)
    valid_tool_set = set(tool_keys)
    if tool_set & valid_tool_set == set():
        # Check for empty set intersection
        raise ValueError("Expected OpenEye tools to have at least of the following {}, "
                         "but instead got {}".format(tool_keys, oetools))
    try:
        for tool in oetools:
            if tool in tool_keys:
                # Try loading the module
                try:
                    module = importlib.import_module('openeye', tool)
                except SystemError: # Python 3.4 relative import fix
                    module = importlib.import_module('openeye.' + tool)
                # Check that we have the license
                if not getattr(module, tools_license[tool])():
                    raise ImportError
    except ImportError:
        return False
    return True

# TODO: Fix me
class LicenseError(Exception):
    pass

# TODO: Fix me
def requires_openeye_licenses(f, *args):
    """
    Decorator to check that OpenEye licenses are found, raising LicenseError if valid license not found

    """
    if not is_openeye_installed(oetools=args):
        # TODO: Include more informative error messages
        raise LicenseError()


# TODO: This only works with the OpenEye toolkit installed; replace with Molecule API
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
    from openeye import oechem, oeomega
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

def get_energy(system, positions):
    """
    Return the potential energy.

    Parameters
    ----------
    system : simtk.openmm.System
        The system to check
    positions : simtk.unit.Quantity of dimension (natoms,3) with units of length
        The positions to use
    Returns
    ---------
    energy
    """

    integrator = openmm.VerletIntegrator(1.0 * unit.femtoseconds)
    context = openmm.Context(system, integrator)
    context.setPositions(positions)
    state = context.getState(getEnergy=True)
    energy = state.getPotentialEnergy() / unit.kilocalories_per_mole
    return energy

#=============================================================================
# OPENMM MERGING AND EXPORTING UTILITY FUNCTIONS
#=============================================================================

# TODO: Reorganize this file, moving exporters to openforcefield.exporters

def create_system_from_amber(prmtop_filename, crd_filename, verbose = False):
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
    Generate a System from the given OEMol and SMIRNOFF forcefield, return the resulting System.

    Parameters
    ----------
    forcefield : ForceField
        SMIRNOFF forcefield
    mol : oechem.OEMol
        Molecule to test (must have coordinates)


    Returns
    ----------
    topology : OpenMM Topology
    system : OpenMM System
    positions : initial atomic positions (OpenMM)
    """
    # Create system
    from openforcefield.utils import generateTopologyFromOEMol
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

def compare_system_energies( topology0, topology1, system0, system1, positions0, positions1=None, label0="AMBER system", label1 = "SMIRNOFF system", verbose = True, skip_assert = False, skip_improper = False ):
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

def get_molecule_parameterIDs( oemols, ffxml):
    """Process a list of oemols with a specified SMIRNOFF ffxml file and determine which parameters are used by which molecules, returning collated results.


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
    """Take an OEMol and a SMIRNOFF forcefield object and return a dictionary, keyed by parameter ID, where each entry is a tuple of ( smirks, [[atom1, ... atomN], [atom1, ... atomN]) giving the SMIRKS corresponding to that parameter ID and a list of the atom groups in that molecule that parameter is applied to.

    Parameters
    ----------
    oemol : OEMol
        OpenEye OEMol with the molecule to investigate.
    ff : ForceField
        SMIRNOFF ForceField object (obtained from an ffxml via ForceField(ffxml)) containing FF of interest.

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

def merge_system( topology0, topology1, system0, system1, positions0, positions1, label0="AMBER system", label1 = "SMIRNOFF system", verbose = True):
    """Merge two given OpenMM systems. Returns the merged OpenMM System.

    Parameters
    ----------
    topology0 : OpenMM Topology
        Topology of first system (i.e. a protein)
    topology1 : OpenMM Topology
        Topology of second system (i.e. a ligand)
    system0 : OpenMM System
        First system for merging (usually from AMBER)
    system1 : OpenMM System
        Second system for merging (usually from SMIRNOFF)
    positions0 : simtk.unit.Quantity wrapped
        Positions to use for energy evaluation comparison
    positions1 (optional) : simtk.unit.Quantity wrapped (optional)
        Positions to use for second OpenMM system
    label0 (optional) : str
        String labeling system0 for output. Default, "AMBER system"
    label1 (optional) : str
        String labeling system1 for output. Default, "SMIRNOFF system"
    verbose (optional) : bool
        Print out info on topologies, True/False (default True)

    Returns
    ----------
    topology : OpenMM Topology
    system : OpenMM System
    positions: unit.Quantity position array
    """

    #Load OpenMM Systems to ParmEd Structures
    structure0 = parmed.openmm.load_topology( topology0, system0 )
    structure1 = parmed.openmm.load_topology( topology1, system1 )

    #Merge parameterized Structure
    structure = structure0 + structure1
    topology = structure.topology

    #Concatenate positions arrays
    positions_unit = unit.angstroms
    positions0_dimensionless = np.array( positions0 / positions_unit )
    positions1_dimensionless = np.array( positions1 / positions_unit )

    coordinates = np.vstack((positions0_dimensionless,positions1_dimensionless))
    natoms = len(coordinates)
    positions = np.zeros([natoms,3], np.float32)
    for index in range(natoms):
        (x,y,z) = coordinates[index]
        positions[index,0] = x
        positions[index,1] = y
        positions[index,2] = z
    positions = unit.Quantity(positions, positions_unit)

    #Generate merged OpenMM system
    system = structure.createSystem()

    if verbose:
        print("Generating ParmEd Structures...\n \t{}: {}\n \t{}: {}\n".format(label0, structure0, label1, structure1))
        print("Merged ParmEd Structure: {}".format( structure ))

    return topology, system, positions

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
    structure.save( top, overwrite = True, format="gromacs")
    structure.save( gro, overwrite = True, format="gro")
