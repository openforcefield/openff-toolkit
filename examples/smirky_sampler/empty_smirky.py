from smarty.environment import *
from smarty.sampler_smirky import *
from smarty.utils import *

# General things
ORs = ['X1', 'X2', 'X3', 'X4']
ANDs = ['+0']
mol2file = 'test_filt1_tripos.mol2'
SMIRFF = "forcefield/Frosst_AlkEtOH.ffxml"
verbose = True
#typetag = 'Bond'  
typetag = 'Angle'
#typetag = 'VdW'
#typetag = 'Torsion'

# load molecules:
molecules = read_molecules(mol2file, verbose = verbose)

# Try to initialize a Bond sampler:
print "Attempting to sample %sS" % typetag.upper()
sampler = TypeSampler(molecules, ORs, ANDs, typetag, None, SMIRFF, 0.0, True)


