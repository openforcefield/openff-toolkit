from functools import partial
import smarty
import openeye
from openeye.oechem import *
import os
from smarty.utils import get_data_filename
import numpy as np
from smarty.forcefield_labeler import *


def test_read_ffxml():
    """Test reading of ffxml files.
    """
    labeler = ForceField_labeler(get_data_filename('forcefield/Frosst_AlkEtOH.ffxml'))

def test_molecule_labeling(verbose = False):
    """Test using ForceField_labeler to provide force terms applied to an oemol."""
    mol = OEMol()
    OEParseSmiles(mol, 'CCC')
    OEAddExplicitHydrogens(mol)
    labeler = ForceField_labeler(get_data_filename('forcefield/Frosst_AlkEtOH.ffxml'))
    labels = labeler.labelMolecules( [mol], verbose = verbose)

    # Check that force terms aren't empty
    if not 'HarmonicBondForce' in labels[0].keys(): 
        raise Exception("No force term assigned for harmonic bonds.")
    if not 'HarmonicAngleForce' in labels[0].keys(): 
        raise Exception("No force term assigned for harmonic angles.")
    if not 'PeriodicTorsionForce' in labels[0].keys(): 
        raise Exception("No force term assigned for periodic torsions.")
    if not 'NonbondedForce' in labels[0].keys(): 
        raise Exception("No nonbonded force term assigned.")


