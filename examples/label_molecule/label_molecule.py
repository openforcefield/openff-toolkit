#!/bin/env python

from openeye.oechem import *
from smarty.utils import get_data_filename

# Create an oemol
mol = OEMol()
OEParseSmiles(mol, 'CCC')
OEAddExplicitHydrogens(mol)

from smarty.forcefield_labeler import *


labeler = ForceField_labeler( get_data_filename('forcefield/Frosst_AlkEtOH.ffxml') )

print labeler.labelMolecules( [mol], verbose = True )
