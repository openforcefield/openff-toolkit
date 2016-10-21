#!/usr/bin/env python

"""
Utility subroutines for SMARTY atom type sampling

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

import openeye.oechem
import openeye.oeomega
import openeye.oequacpac

from openeye.oechem import *
from openeye.oeomega import *
from openeye.oequacpac import *

import time
from simtk import unit

#=============================================================================================
# UTILITY ROUTINES
#=============================================================================================

def get_data_filename(relative_path):
    """Get the full path to one of the reference files in testsystems.

    In the source distribution, these files are in ``smarty/data/``,
    but on installation, they're moved to somewhere in the user's python
    site-packages directory.

    Parameters
    ----------
    name : str
        Name of the file to load (with respect to the repex folder).

    """

    from pkg_resources import resource_filename
    fn = resource_filename('smarty', os.path.join('data', relative_path))

    if not os.path.exists(fn):
        raise ValueError("Sorry! %s does not exist. If you just added it, you'll have to re-install" % fn)

    return fn

def normalize_molecules(molecules):
    """
    Normalize all molecules in specified set.

    ARGUMENTS

    molecules (list of OEMol) - molecules to be normalized (in place)

    """

    # Add explicit hydrogens.
    for molecule in molecules:
        openeye.oechem.OEAddExplicitHydrogens(molecule)

    # Build a conformation for all molecules with Omega.
    print("Building conformations for all molecules...")
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

    if not os.path.exists(filename):
        built_in = get_data_filename('molecules/%s' % filename)
        if not os.path.exists(built_in):
            raise Exception("File '%s' not found." % filename)
        filename = built_in

    if verbose: print("Loading molecules from '%s'..." % filename)
    start_time = time.time()
    molecules = list()
    input_molstream = oemolistream(filename)

    from openeye import oechem
    flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield
    input_molstream.SetFlavor(oechem.OEFormat_MOL2, flavor)

    molecule = OECreateOEGraphMol()
    while OEReadMolecule(input_molstream, molecule):
        # Get molecule name.
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

def parse_odds_file(filename, verbose = False):
    """
    parses files that have the form
    decorator       odds
    if only one column odds will be assumed equally probable

    Parameters
    -----------
    filename: string or file object
    may be an absolute file path, a path relative to the current working directory, a path relative to this module's data subdirectory (for built in decorator files), or an opten file-like object with a readlines() method.

    Returns
    --------
    choices: 2-tuple of the form ( [decorators], [odds] )
    """
    if verbose:
        if isinstance(filename, file):
            print("Attempting to parse file '%s'" % filename.name)
        else:
            print("Attempting to parse file '%s'" % filename)

    # if no file return None
    if filename is None:
        return None

    # if input is a file object
    try:
        input_lines = filename.readlines()
        if verbose: print("Attempting to parse file '%s'" % filename.name)
    except AttributeError:
        if verbose: print("Attempting to parse file '%s'" % filename)
        try:
            ifs = open(filename, 'r')
            input_lines = ifs.readlines()
        except IOError:
            ifs = get_data_filename(filename)
            ifs = open(ifs, 'r')
            input_lines = ifs.readlines()
        except Exception as e:
            raise Exception("%s\nProvided file (%s) could not be parsed" % (str(e), filename))
    except Exception as e:
        msg = str(e) + '\n'
        msg += "Could not read data from file %s" % filename
        raise Exception(msg)

    # close file
    ifs.close()

    decorators = []
    odds = []
    noOdds = False
    for l in input_lines:
        # skip empty lines
        if len(l) == 0:
            continue
        # check for and remove comments
        comment = l.find('%')
        if comment == -1: # no comment
            entry = l.split()
        elif comment > 0: # remove trailing comment
            entry = l[:comment].split()
        else: # whole line is a comment skip
            continue

        # add decorator
        if entry[0] == "''" or entry[0] == '""':
            decorators.append('')
        else:
            decorators.append(entry[0])

        if len(entry) == 2:
            odds.append(float(entry[1]))
        elif len(entry) == 1:
            noOdds = True
        else:
            raise Exception("Error entry (%s) in decorator file '%s' is invalid" % (l, filename))

    if (odds.count(0) == len(odds)) or noOdds:
        odds = None
        #TODO: handle case where 1 line is missing odds entry

    return (decorators, odds)


def setPositionsInOEMol(molecule, positions):
    """Set the positions in an OEMol using a position array with units from simtk.unit, i.e. from OpenMM. Atoms must have same order.

    Arguments:
    ---------
    molecule : OEMol
        OpenEye molecule
    positions : Nx3 array
        Unit-bearing via simtk.unit Nx3 array of coordinates
    """
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
