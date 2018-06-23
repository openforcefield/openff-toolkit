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

from openforcefield.utils import get_testdata_filename, generateTopologyFromOEMol, read_molecules
from openforcefield.utils import check_energy_is_finite, get_energy

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
    prmtop_filename = get_testdata_filename(prefix+'.prmtop')
    inpcrd_filename = get_testdata_filename(prefix+'.inpcrd')
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
    pdb_filename = get_testdata_filename(prefix+'.pdb')
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
    mol2_filename = get_testdata_filename(prefix+'.pdb')
    return mol2_filename
