from smarty.environment import *
from smarty.sampler_smirky import *
from smarty.utils import *
import time

# General things
temperature = 1.0*10**(-6)
Decors = ['X1', 'X2', 'X3', 'X4', 'H0', 'H1', 'H2', 'H3','+0', '+1', '-1', 'a', 'A']
mol2file = 'test_filt1_tripos.mol2'
SMIRFF = "forcefield/Frosst_AlkEtOH.ffxml"
verbose = True
replacements = None
elementList = ['#%i' % i for i in range(1,119)]
typetag = 'Angle'
traj = "startEmpty_Angle.csv"
niterations = 100

# load molecules:
molecules = read_molecules(mol2file, verbose = verbose)

# Try to initialize a Bond sampler:
print "Attempting to sample %sS" % typetag.upper()
init_time = time.time()
sampler = TypeSampler(molecules, typetag, elementList, Decors, Decors, replacements, None, SMIRFF, temperature, True)
sampler.run(niterations, traj)
final_time = time.time()
elapsed = (final_time-init_time) / 60.0
print "%i iterations took %.2f minutes" % (niterations, elapsed)

