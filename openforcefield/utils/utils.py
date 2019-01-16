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

from simtk.openmm import app
from simtk.openmm.app import element as elem
from simtk.openmm.app import Topology

import time
from simtk import unit

import logging
import functools

#=============================================================================================
# UTILITY DECORATORS
#=============================================================================================


def deprecated(func):
    """A useful decorator to mark functions as deprecanted.

    TODO: append sphinx markup to functions docstring?
    """
    @functools.wraps(func)
    def new_func(*args, **kwargs):

        logging.warning("A deprecated function has been called {}.".format(func.__name__))
        return func(*args, **kwargs)

    return new_func

#=============================================================================================
# UTILITY ROUTINES
#=============================================================================================
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
    from openeye.oechem import OEMol

    # Avoid manipulating the molecule
    mol = OEMol(molecule)

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

    import openeye.oechem
    import openeye.oeomega
    import openeye.oequacpac

    from openeye.oechem import OEFormalPartialCharges, OEHasPartialCharges
    from openeye.oequacpac import OEAssignPartialCharges, OECharges_AM1BCC
    from openeye.oeiupac import OECreateIUPACName

    # Add explicit hydrogens.
    for molecule in molecules:
        openeye.oechem.OEAddExplicitHydrogens(molecule)

    # Build a conformation for all molecules with Omega.
    print("Building conformations for all molecules...")
    start_time = time.time()
    import openeye.oeomega
    omega = openeye.oeomega.OEOmega()
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
    mcmcDbName = '' # TODO: Fix line.
    ligand_mol2_dirname  = os.path.dirname(mcmcDbName) + '/mol2'
    if( not os.path.exists( ligand_mol2_dirname ) ):
        os.makedirs(ligand_mol2_dirname)
    ligand_mol2_filename = ligand_mol2_dirname + '/temp' + os.path.basename(mcmcDbName) + '.mol2'
    start_time = time.time()
    omolstream = openeye.oechem.oemolostream(ligand_mol2_filename)
    for molecule in molecules:
        # Write molecule as mol2, changing molecule through normalization.
        openeye.oechem.OEWriteMolecule(omolstream, molecule)
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
    from openeye.oechem import OECreateOEGraphMol, oemolistream, OEReadMolecule, OEGetSDData, OEMol

    if not os.path.exists(filename):
        built_in = get_data_filename('molecules/%s' % filename)
        if not os.path.exists(built_in):
            raise Exception("File '%s' not found." % filename)
        filename = built_in

    if verbose: print("Loading molecules from '%s'..." % filename)
    start_time = time.time()
    molecules = list()
    input_molstream = oemolistream(filename)

    flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield
    input_molstream.SetFlavor(oechem.OEFormat_MOL2, flavor)

    molecule = OECreateOEGraphMol()

    while OEReadMolecule(input_molstream, molecule):
        # If molecule has no title, try getting SD 'name' tag
        if molecule.GetTitle() == '':
            name = OEGetSDData(molecule, 'name').strip()
            molecule.SetTitle(name)
        # Append to list.
        molecule_copy = OEMol(molecule)
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
    from openeye.oechem import OEFloatArray

    if molecule.NumAtoms() != len(positions): raise ValueError("Number of atoms in molecule does not match length of position array.")
    pos_unitless = positions/unit.angstroms

    coordlist = []
    for idx in range(len(pos_unitless)):
        for j in range(3):
            coordlist.append( pos_unitless[idx][j])

    molecule.SetCoords(OEFloatArray(coordlist))

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
