from smarty.environment import *
from smarty.forcefield_labeler import *
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
typetag = 'Torsion'
replacements = None
elements = ["[#%i]" %i for i in range(1,119)]

# load molecules:
molecules = read_molecules(mol2file, verbose = verbose)

# Make molecule for depict notebook
smile = "C[C@@H]1[C@@H]([C@@H](O[C@@H](OCCO1)O)O)O"
mol = OEMol()
if not OEParseSmiles(mol, smile):
    raise Exception("Error couldn't parse smiles: %s" % smile)
OEAddExplicitHydrogens(mol)
#molecules = [mol]

# default atomlists
H = [['#1'], None]
empty = [None, None]
bond = [['-'], None]
CT = [['#6X4'], None]
O = [['#8'], None]
OX2 = [['#8X2'], None]
nH = [['!#1'], None]
OH0 = [['#8'], ['X2','H0']]

# This is going to be tedius to type out, but we can copy from here in the future
t1 = TorsionChemicalEnvironment(empty, bond, CT, bond, CT, bond, empty)

t2 = TorsionChemicalEnvironment(empty, bond, CT, bond, OX2, bond, H)

t3 = TorsionChemicalEnvironment(empty, bond, CT, bond, OX2, bond, nH)

t4 = TorsionChemicalEnvironment( H, bond, CT, bond, CT, bond, H)

t5 = TorsionChemicalEnvironment(H, bond, CT, bond, CT, bond, CT)

t6 = TorsionChemicalEnvironment(CT, bond, CT, bond, OX2, bond, H)

t7 = TorsionChemicalEnvironment(CT, bond, CT, bond, CT, bond, CT)

t8 = TorsionChemicalEnvironment(CT, bond, CT, bond, OX2, bond, CT)

t9 = TorsionChemicalEnvironment(CT, bond, OX2, bond, CT,bond, OH0)

t10 = TorsionChemicalEnvironment(OX2, bond, CT, bond, CT, bond, OX2)

t11 = TorsionChemicalEnvironment(OX2, bond, CT, bond, CT, bond, H)

t12 = TorsionChemicalEnvironment(H, bond, CT, bond, CT, bond,OX2)

initialList = [t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12]

# Begin sampling
print "Attempting to sample %sS" % typetag.upper()
sampler = TypeSampler(molecules, typetag, elements, ORs, ANDs, replacements, initialList, SMIRFF, 0.0, True)


