from smarty.environment import *
from smarty.sampler_smirky import *
from smarty.utils import *
from operator import itemgetter, attrgetter
from openeye.oechem import *
import copy
import sys # used to exit while testing

# General things
ORs = ['X1', 'X2', 'X3', 'X4']
ANDs = ['+0']
mol2file = 'test_filt1_tripos.mol2'
SMIRFF = "forcefield/Frosst_AlkEtOH.ffxml"
verbose = True
typetag = 'Bond'
replacements = None
forcetype = "NonbondedForce"
elements = ["[#%i]" %i for i in range(1,119)]

# load molecules:
molecules = read_molecules(mol2file, verbose = verbose)

# default atomlists
H = [['#1'], None]
empty = [None, None]
bond = [['-'], None]
CT = [['#6X4'], None]
O = [['#8'], None]

# make bond types using environments
# This is going to be tedius to type out, but we can copy from here in the future
CTCT = BondChemicalEnvironment(CT,bond, CT)
CTCT.label = 'CTCT'

CTH = BondChemicalEnvironment(CT, bond, H)
CTH.label = 'CTH_'

OH = BondChemicalEnvironment(O, empty, H)
OH.label = 'OHHO'

CTOH = BondChemicalEnvironment(CT, bond, O)
CTOH.atom2.setANDtypes(['X2','H1'])
CTOH.label = 'CTOH'

CTOS = BondChemicalEnvironment(CT, bond, O)
CTOS.atom2.setANDtypes(['X2', 'H0'])
CTOS.label = 'CTOS'

initialList = [CTCT, CTH, OH, CTOH, CTOS]
for env in initialList:
    print env.asSMIRKS()

initialList = [[e.asSMIRKS(), e.label] for e in initialList]
# Begin sampling
print "Attempting to sample %sS" % typetag.upper()
sampler = TypeSampler(molecules, typetag, elements, ORs, ANDs, replacements, initialList, SMIRFF, 0.0, True)
frac = sampler.run(2)
print frac

