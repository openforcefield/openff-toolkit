from smarty.environment import *
from smarty.sampler_smirky import *
from smarty.utils import *
import copy
import sys # used to exit while testing

# General things
ORs = ['X1', 'X2', 'X3', 'X4']
ANDs = ['+0']
mol2file = 'test_filt1_tripos.mol2'
SMIRFF = "forcefield/Frosst_AlkEtOH.ffxml"
verbose = True
typetag = 'VdW'

# load molecules:
molecules = read_molecules(mol2file, verbose = verbose)

# default atomlists
H = [['#1'], None]
bond = [None, None]
C = [['#6'], None]
O = [['#8'], None]

# make atom types using environments
# This is going to be tedius to type out, but we can copy from here in the future
H_default = AtomChemicalEnvironment(H)

HC = copy.deepcopy(H_default)
HC.addAtom(HC.atom1, None, None, C[0])

H1 = copy.deepcopy(H_default)
c = H1.addAtom(H1.atom1, None, None, C[0])
H1.addAtom(c, None, None, O[0])


H2 = copy.deepcopy(H_default)
c = H2.addAtom(H2.atom1, None, None, C[0])
H2.addAtom(c, None, None, O[0])
H2.addAtom(c, None, None, O[0])

H3 = copy.deepcopy(H_default)
c = H3.addAtom(H3.atom1, None, None, C[0])
H3.addAtom(c, None, None, O[0])
H3.addAtom(c, None, None, O[0])
H3.addAtom(c, None, None, O[0])

HO = copy.deepcopy(H_default)
HO.addAtom(HO.atom1, None, None, O[0])

CT = AtomChemicalEnvironment([['#6X4'], None])

OS = AtomChemicalEnvironment([['#8X2'], None])

OH = AtomChemicalEnvironment([['#8X2+0'], None])
OH.addAtom(OH.atom1, None, None, ['#1'])

initialList = [HC, H1, H2, H3, HO, CT, OS, OH]

# Begin sampling
print "Attempting to sample %sS" % typetag.upper()
sampler = TypeSampler(molecules, ORs, ANDs, typetag, initialList, SMIRFF, 0.0, True)


