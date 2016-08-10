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
typetag = 'VdW'
forcetype = "NonbondedForce"
elements = ["[#%i]" %i for i in range(1,119)]
replacements = None

# load molecules:
# molecules = read_molecules(mol2file, verbose = verbose)

# Make molecule for depict notebook
smile = "C[C@@H]1[C@@H]([C@@H](O[C@@H](OCCO1)O)O)O"
mol = OEMol()
if not OEParseSmiles(mol, smile):
    raise Exception("Error couldn't parse smiles: %s" % smile)
OEAddExplicitHydrogens(mol)
molecules = [mol]

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
HC.label = 'HC'

H1 = copy.deepcopy(H_default)
c = H1.addAtom(H1.atom1, None, None, C[0])
H1.addAtom(c, None, None, O[0])
H1.label = 'H1'


H2 = copy.deepcopy(H_default)
c = H2.addAtom(H2.atom1, None, None, C[0])
H2.addAtom(c, None, None, O[0])
H2.addAtom(c, None, None, O[0])
H2.label = 'H2'

H3 = copy.deepcopy(H_default)
c = H3.addAtom(H3.atom1, None, None, C[0])
H3.addAtom(c, None, None, O[0])
H3.addAtom(c, None, None, O[0])
H3.addAtom(c, None, None, O[0])
H3.label = 'H3'

HO = copy.deepcopy(H_default)
HO.addAtom(HO.atom1, None, None, O[0])
HO.label = 'HO'

CT = AtomChemicalEnvironment([['#6X4'], None])
CT.label = 'CT'

OS = AtomChemicalEnvironment([['#8X2'], None])
OS.label = 'OS'

OH = AtomChemicalEnvironment([['#8X2+0'], None])
OH.addAtom(OH.atom1, None, None, ['#1'])
OH.label = 'OH'

initialList = [HC, H1, H2, H3, HO, CT, OS, OH]

# Begin sampling
print "Attempting to sample %sS" % typetag.upper()
sampler = TypeSampler(molecules, typetag, elements, ORs, ANDs, replacements, initialList, SMIRFF, 0.0, True)

