# Chemical Environment usage examples

This directory shows how to use Chemical Environments. First we will outline how they are structured with some example code snippets. 
For a more hands on tutorial you can look at the jupyter notebook `using_environments.py`. 

### Initiating ChemicalEnvironments

All Chemical Environments can be initated using SMIRKS strings. 
The general ChemicalEnvironment, with no input string is completely empty, however there
are 5 subtypes of ChemicalEnvironments that match the types of parameters found in the SMIRNOFF format:

* AtomChemicalEnvironment
    - expects 1 indexed atom
    - default/generic SMIRKS `"[*:1]"`
* BondChemicalEnvironment
    - expects 2 indexed atoms
    - default/generic SMIRKS: `"[*:1]~[*:2]"`
* Angle
    - expects 3 indexed atoms
    - default/generic SMIRKS: `"[*:1]~[*:2]~[*:3]"`
* Torsion
    - expects 4 indexed atoms in a proper dihedral angle
    - default/generic SMIRKS: `"[*:1]~[*:2]~[*:3]~[*:4]"`
* Improper 
    - expects 4 indexed atoms in an improper dihedral angle
    - default/generic SMIRKS: `"[*:1]~[*:2](~[*:3])~[*:4]"` 

Here is an example for initiating a chemical environment:
```python
from openforcefield.typing.chemistry import environment as env
smirks = "[#6X3,#7;+0:1]~;@[#8;r:2]~;@[#6X3,#7;+0:3]"
angle = env.AngleChemicalEnvironment(smirks = smirks)
```

### ANDtypes and ORtypes

ChemicalEnvironments store information about the SMIRKS string in the same way a 
chemist might think about a fragment, that is there are atoms connected by bonds. 
The descriptors for both atoms and bonds are broken down into ORtypes, decorators OR'd together with a ','
and ANDtypes, decorators AND'd to the end of an atom or bond with a ';'. 
ORtypes are further broken down into bases and decorators. A base is typically an atomic number or symbol.  

As an example, consider atom 1 above ('[#6X3,#7X2,#7X3;+0:1]'):
* ORtypes (base, [list of decorators] )
    - `('#6', ['X3'])`
    - `('#7', [])`
* ANDtypes
    - `['+0']`

You can set and get ORtypes and ANDtypes for atoms and bonds. You can also addORtypes and addANDtypes for atoms and bonds. 
For example:
```python
new_ORtypes = [ ('#6', ['X3']), ('#7', ['X2']) ]
atom1.setORtypes(new_ORtypes)
print("New Atom 1: %s " % atom1.asSMIRKS())
# New Atom 1: [#6X3,#7X2;+0:1]

# Change atom2's AND and OR types with the add*type methods
atom2 = angle.selectAtom(2)
atom2.addANDtype('R1')
atom2.addORtype('#7', ['X3', '+0'])
print("New Atom 2: %s" % atom2.asSMIRKS())
# New Atom 2: [#8,#7+0X3;R1:2]

print("\nNew SMIRKS: %s" % angle.asSMIRKS())
# New SMIRKS: [#6X3,#7X2:1]~;@[#8,#7+0X3;R1:2]~;@[#6X3,#7:3]
```

### Adding and removing atoms

You can add new atoms to an environment by specifying an atom to connect
to and the OR and AND types for the new bond and new atom. 

You can also remove existing atoms, this move returns True if the atom was removed and False if it wasn't. 
Removing an atom removes its associated bonds. Calling removeAtom on an indexed atom or an atom connecting two 
other atoms will always return False. 

```python
atom3 = angle.selectAtom(3)

alpha_ORtypes = [('#8', ['X2', 'H1'])]
alpha_ANDtypes = ['R0']
alpha_bondANDtypes = ['!@']
alpha = angle.addAtom(atom3, bondANDtypes = alpha_bondANDtypes, 
                      newORtypes = alpha_ORtypes, newANDtypes = alpha_ANDtypes)
print("Alpha Atom SMIRKS: %s" % alpha.asSMIRKS())
# Alpha Atom SMIRKS: [#8X2H1;R0]

beta_ORtypes = [('#1', [])]
beta = angle.addAtom(alpha, newORtypes = beta_ORtypes)
print("Beta Atom SMIRKS: %s" % beta.asSMIRKS())
# Beta Atom SMIRKS: [#1]

print("New overall SMIRKS: %s" % angle.asSMIRKS())
# New overall SMIRKS: [#6X3,#7X2:1]~;@[#8,#7+0X3;R1:2]~;@[#6X3,#7:3]~;!@[#8X2H1;R0]~[#1]

removed = angle.removeAtom(beta)
print("The hydrogen beta to atom3 was remove: ", removed)
# The hydrogen beta to atom3 was remove:  True
print("Updated SMIRKS string: %s" % angle.asSMIRKS())
# print("Updated SMIRKS string: %s" % angle.asSMIRKS())
```

### Other ChemicalEnvironment Methods

0. Selecting atoms and bonds
```python
# when given an integer selectAtom or selectBond returns that atom or bond with that index 
# otherwise it returns a random atom or bond that meets the specified requirement 
# If no atom or bond with that type is found then None is returned
alpha = angle.selectAtom('alpha')
atom4 = angle.selectAtom(4) # None, angle has 3 indexed atoms
beta_bond = angle.selectBond('beta') # bond between an alpha and beta atom
```

1. Getting information about an atom or bond in an environment (i.e. isAlpha returns a boolean)
```python
# Check if the alpha atom above is any of the following
angle.isIndexed(alpha) # False
angle.isUnindexed(alpha) # True 
angle.isAlpha(alpha) # True
angle.isBeta(alpha) # False
# NOTE - These methods can take an atom or a bond as an argument
```
2. Get atoms or bonds in each type of position, for example getIndexedAtoms or getAlphaBonds
```python
# We will print the SMIRKS for each indexed atom:
all_atoms = angle.getAtoms()
indexed_list = angle.getIndexedAtoms() # list of indexed atoms
unindexed_list = angle.getUnindexedAtoms() # list of atoms with no index
alpha_list = angle.getAlphaAtoms()
beta_list = angle.getBetaAtoms
# Note - Atoms can be replaced with Bonds in each of the above examples
```
3. Report the minimum order of a bond with Bond.getOrder
```python
bond1 = angle.selectBond(1)
bond1.getOrder() # returns 1
```
4. Report the valence and bond order around an atom can be reported with getValence and getBondORder
```python
atom3 = angle.selectAtom(3)
angle.getValence(atom3) # number of neighbors atom 3 has
angle.getBondOrder(atom3) # total bond order around atom 3
```
5. Get a bond between two atoms (or determine if the atoms are bonded) with getBond(atom1, atom2)
```python
# Check for bonds between each pair of indexed atoms
A = 1
B = 2
atomA = angle.selectAtom(A)
atomB = angle.selectAtom(B)
bondAB = angle.getBond(atomA, atomB) # returns None if no bond found
```
6. Get atoms bound to a specified atom with getNieghbors
```python
# get the neighbors for each indexed atom
atomA = angle.selectAtom(A)
neighbor_list = angle.getNeighbors(atomA) 
```

### `using_environments.ipynb`

This notebook should be your first stop for interactively understanding and using ChemicalEnvironemnts. 
It includes all of the examples shown above and more in a jupyter notebook so you can play with the different method options. 

### Making moves with Chemical Environments

Here our goal is to show how to use Chemical Environments to make "moves" in chemical perception space. First a list of weighted moves were generated and saved to output files. Then, a seconds notebook is used to make these moves in chemical space. 

These notebooks were generated in the early planning stages for smirky, the first attempt at sampling chemical perception for SMIRNOFF parameter types.  

* `create_move_types_and_weights.ipynb` - used to create lists of moves and their probabilities based on the type of parameter (vdw, bond, angle, torsion, improper). Created the text files `moveTrees.uniq.*.txt` where `*` refers to the parameter type.  

* `moves_using_environments.ipynb` - uses the `moveTrees.uniq.*.txt` files to illustrate how to make moves in chemical perception space. 

For a more thorough and automated tool for using ChemicalEnvironments to make changes to chemical perception see [smirky](https://github.com/open-forcefield-group/smarty), a tool we developed to rediscover SMIRNOFF parameters. 
