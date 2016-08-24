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
typetag = 'Angle'
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

# make bond types using environments
# This is going to be tedius to type out, but we can copy from here in the future

aCTa = AngleChemicalEnvironment(empty, bond, CT, bond, empty)
aCTa.label = 'aCTa'

HCTH = AngleChemicalEnvironment(H, bond, CT, bond, H)
HCTH.label = 'HCTH'

CTCTCT = AngleChemicalEnvironment(CT, bond, CT, bond, CT)
CTCTCT.label = 'CTCTCT'

OCTO = AngleChemicalEnvironment(O, bond, CT, bond, O)
OCTO.label = 'OCTO'

CTOHHO = AngleChemicalEnvironment(CT, bond, O, bond, H)
CTOHHO.atom2.setORtypes = ['#8X2']
CTOHHO.label = 'CTOHHO'

CTOSCT = AngleChemicalEnvironment(CT, bond, O, bond, CT)
CTOSCT.atom2.setORtypes = ['#8X2']
CTOSCT.label = 'CTOSCT'

initialList = [aCTa, HCTH, CTCTCT, OCTO, CTOHHO, CTOSCT]

# Begin sampling
print "Attempting to sample %sS" % typetag.upper()
sampler = TypeSampler(molecules, typetag, elements, ORs, ANDs, replacements, initialList, SMIRFF, 0.0, True)



