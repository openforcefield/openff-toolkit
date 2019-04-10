#!/usr/bin/env python

#==============================================================================
# MODULE DOCSTRING
#==============================================================================

"""
environment.py

.. warning :: This file is  will be updated to comply with PEP8.

Classes defining a chemical environment for atoms and how they are connected
using networkx graph objects to organize and make changes to the structure.
Output will be in the form of SMARTS and SMIRKS.

AUTHORS

Caitlin Bannan <bannanc@uci.edu>, Mobley Lab, University of California Irvine,
with contributions from John Chodera, Memorial Sloan Kettering Cancer Center
and David Mobley, UC Irvine.

"""

__all__ = [
    'SMIRKSMismatchError',
    'SMIRKSParsingError',
    'ChemicalEnvironment',
    'AtomChemicalEnvironment',
    'BondChemicalEnvironment',
    'AngleChemicalEnvironment',
    'TorsionChemicalEnvironment',
    'ImproperChemicalEnvironment'
]


#==============================================================================
# GLOBAL IMPORTS
#==============================================================================

import re
import copy

import networkx as nx
from numpy import random

import openforcefield.utils


#==============================================================================
# Functions
#==============================================================================

def _find_embedded_brackets(string, in_char, out_char):
    """
    Finds the substring surrounded by the in_char and out_char
    intended use is to identify embedded bracketed sequences

    string - a string you want separated
    in_char - regular expression for the character you're looking for '\(' for '('
    out_char - regular expression for the closing character such as '\)' for ')'

    string = "[#1$(*-C(-[#7,#8,F,#16,Cl,Br])-[#7,#8,F,#16,Cl,Br]):1]"
    sub_string, in_idx, out_idx = _find_embedded_brackets(string, '\(','\)')
    # sub_string = (*-C(-[#7,#8,F,#16,Cl,Br])-[#7,#8,F,#16,Cl,Br])  in_idx = 4, out_idx = 50
    """
    in_list = [m.start() for m in re.finditer(in_char, string)]
    out_list = [m.start() for m in re.finditer(out_char, string)]
    # If no occurance of the in_char return an empty string
    if len(in_list) == 0:
        return "", -1, -1
    # If no out_char returns the first in_char to the end
    if len(out_list) == 0:
        return string[in_list[0]:], in_list[0], -1

    # Otherwise find closure from the first in_char
    list_idx = 0
    while list_idx < len(in_list) - 1:
        if in_list[list_idx+1] > out_list[list_idx]:
            break
        list_idx+=1
    in_idx = in_list[0]
    out_idx = out_list[list_idx]
    return string[in_idx:out_idx+1], in_idx, out_idx

def _convert_embedded_SMIRKS(smirks):
    """
    Converts a SMIRKS string with the $(...) in an atom to the
    form expected by the environment parser

    smirks = any smirks string, if no $(...) then the original smirks is returned

    initial_smirks = "[#1$(*~[#6]):1]"
    new_smirks = _convert_embedded_SMIRKS(initial_smirks)
    # new_smirks = [#1:1]~[#6]
    """
    a_out = 0
    while smirks.find('$(') != -1:
        # Find first atom
        atom, a_in, a_out = _find_embedded_brackets(smirks, r'\[', r'\]')
        d = atom.find('$(')
        # Find atom with the $ string embedded
        while d == -1:
            atom, temp_in, temp_out = _find_embedded_brackets(smirks[a_out+1:], r'\[', r'\]')
            a_in = a_out + temp_in + 1
            a_out += temp_out + 1
            d = atom.find('$(')

        # Store the smirks pattern before and after the relevant atom
        pre_smirks = smirks[:a_in]
        post_smirks = smirks[a_out+1:]

        # Check for ring index, i.e. the 1s in "[#6:1]1-CCCCC1"
        match = re.match(r'(\d+)',post_smirks)
        if match is not None: # leftover starts with int
            ring_out = re.findall(r'(\d+)',post_smirks)[0]
            # update post_smirks
            post_smirks = post_smirks[match.end():]
        else:
            ring_out = ''

        embedded, p_in, p_out = _find_embedded_brackets(atom, r'\(', r'\)')
        # two forms of embedded strings $(*~stuff) or $([..]~stuff)
        # in the latter case the first atom refers the current atom
        if embedded[1] == '[':
            first, f_in, f_out = _find_embedded_brackets(embedded, r'\[', r'\]')
            first = _convert_embedded_SMIRKS(first)
            new_atom = atom[:d]+first[1:-1]+atom[p_out+1:]
            embedded = embedded[f_out+1:]
            # if embedded is empty between brackets, remove it
            if embedded.replace('(','').replace(')','') == '':
                embedded = ''

        elif embedded[1] == '*': # embedded[1] = *
            new_atom = atom[:d]+atom[p_out+1:]
            embedded = embedded[2:]

        else: # embedded starts with a "no bracket" atom such as 'C'
            embedded = embedded[1:] # remove leading '('
            # atoms by symbol don't need brackets, this covers atomic symbols and aromatic atoms
            no_bracket = r'(!?[A-Z][a-z]?|!?[cnops])'
            match = re.match(no_bracket, embedded)
            if match is not None:
                new_atom = atom[:d]+embedded[:match.end()]+atom[p_out+1:]
                embedded = embedded[match.end():]
            else:
                new_atom = atom[:d]+atom[p_out+1]

        # Look for ring insided embedded SMIRKS "[#6$(*1CCC1)]"
        match = re.match(r'(\d+)', embedded)
        if match is not None: # embedded starts with an int
            ring_in = re.findall(r'(\d+)', embedded)[0]
            embedded = '(' + embedded[match.end():]
        else:
            ring_in = ''
            if embedded != '':
                embedded = '(' + embedded

        # Make new smirks
        smirks = pre_smirks+new_atom+ring_out+ring_in+embedded+post_smirks

    return smirks

def _remove_blanks_repeats(init_list, remove_list = ['']):
    """
    Returns the input list 'init_list'
    without any repeating entries or blank strings ''
    """
    final_list = [item for item in init_list if item not in remove_list]
    return list( set(final_list) )


class SMIRKSMismatchError(openforcefield.utils.MessageException):
    """
    Exception for cases where smirks are inappropriate
    for the environment type they are being parsed into
    """
    pass


class SMIRKSParsingError(openforcefield.utils.MessageException):
    """
    Exception for when SMIRKS are not parseable for any environment
    """
    pass


class ChemicalEnvironment:
    """Chemical environment abstract base class that matches an atom, bond, angle, etc.

    .. warning :: This class is largely redundant with the same one in the Chemper package, and will likely be removed.

    """
    class Atom:
        """Atom representation, which may have some ORtypes and ANDtypes properties.

        Attributes
        ----------
        ORtypes : list of tuples in the form (base, [list of decorators])
            where bases and decorators are both strings
            The descriptor types that will be combined with logical OR
        ANDtypes : list of string
            The descriptor types  that will be combined with logical AND
        """
        def __init__(self, ORtypes = None, ANDtypes = None, index = None, ring = None):
            """Initialize an Atom object with optional descriptors.

            Parameters
            -----------
            ORtypes: list of tuples for ORbases and ORdecorators,
                in the form (base, [list of decorators])
                optional, default = []
            ANDtypes: list of str,
                strings that will be AND'd together in a SMARTS
                optional, default = None
            index : int, optional, default=None
                If not None, the specified index will be attached as a SMIRKS index (e.g. '[#6:1]')
            ring : int, optional, default = None
                If not None, the specified ring index will be attached at the end of the atom i.e. '[#6:1]1'
            """
            # List of 2 tuples in the form [ (ORbase, ORdecorator), ...]
            if ORtypes is not None:
                self.ORtypes = copy.deepcopy(ORtypes)
            else:
                self.ORtypes = list()

            # Set of strings that will be AND'd to the the end
            if ANDtypes is not None:
                self.ANDtypes = list(copy.deepcopy(ANDtypes))
            else:
                self.ANDtypes = list()

            self.index = index
            self.ring = ring
            self._atom = True

        def asSMARTS(self):
            """Return the atom representation as SMARTS.

            Returns
            --------
            smarts : str
            The SMARTS string for this atom
            """

            smarts = '['

            # Add the OR'd features
            if self.ORtypes:
                ORList = list()
                for (base, ORdecorators) in self.ORtypes:
                    if len(base) > 0 and base[0] == '$':
                        # after a $base an explicit '&' is necessary
                        if ORdecorators:
                            OR = base+'&'+''.join(ORdecorators)
                        else:
                            OR = base
                    else: # base doesn't start with $
                        OR = base+''.join(ORdecorators)
                    ORList.append(OR)
                smarts += ','.join(ORList)
            else:
                smarts += '*'

            if len(self.ANDtypes) > 0:
                smarts += ';' + ';'.join(self.ANDtypes)

            if self.ring is not None:
                return smarts + ']' + str(self.ring)
            else:
                return smarts + ']'

        def asSMIRKS(self):
            """Return the atom representation as SMIRKS.

            Returns
            --------
            smirks : str
            The SMIRKS string for this atom
            """
            smirks = self.asSMARTS()

            # No index specified so SMIRKS = SMARTS
            if self.index is None:
                return smirks

            # Add label to the end of SMARTS
            else:
                sub_string, start, end = _find_embedded_brackets(smirks, r'\[', r'\]')
                if self.ring is not None:
                    return sub_string[:-1] + ':' + str(self.index) + ']'+str(self.ring)
                else:
                    return sub_string[:-1] + ':' + str(self.index) + ']'

        def addORtype(self, ORbase, ORdecorators):
            """
            Adds ORtype to the set for this atom.

            Parameters
            ----------
            ORbase: string, such as '#6'
            ORdecorators: list of strings, such as ['X4','+0']
            """
            ORdecorators = _remove_blanks_repeats(ORdecorators, ['',ORbase])
            self.ORtypes.append((ORbase, ORdecorators))

        def addANDtype(self, ANDtype):
            """
            Adds ANDtype to the set for this atom.

            Parameters
            ----------
            ANDtype: string
                added to the list of ANDtypes for this atom
            """
            self.ANDtypes.append(ANDtype)
            self.ANDtypes = _remove_blanks_repeats(self.ANDtypes)

        def getORtypes(self):
            """
            returns a copy of the dictionary of ORtypes for this atom
            """
            return copy.deepcopy(self.ORtypes)

        def setORtypes(self, newORtypes):
            """
            sets new ORtypes for this atom

            Parameters
            ----------
            newORtypes: list of tuples in the form (base, [ORdecorators])
                for example: ('#6', ['X4','H0','+0']) --> '#6X4H0+0'
            """
            self.ORtypes = list()
            if newORtypes is not None:
                for (base, decs) in newORtypes:
                    adjusted_decs = _remove_blanks_repeats(decs, ['', base])
                    self.ORtypes.append( (base, adjusted_decs) )

        def getANDtypes(self):
            """
            returns a copy of the list of ANDtypes for this atom
            """
            return list(copy.deepcopy(self.ANDtypes))

        def setANDtypes(self, newANDtypes):
            """
            sets new ANDtypes for this atom

            Parameters
            ----------
            newANDtypes: list of strings
                strings that will be AND'd together in a SMARTS
            """
            if newANDtypes is None:
                self.ANDtypes = list()
            else:
                self.ANDtypes = _remove_blanks_repeats(newANDtypes)

    class Bond(Atom):
        """Bond representation, which may have ORtype and ANDtype descriptors.

        Attributes
        ----------
        ORtypes : list of tuples of ORbases and ORdecorators
            in form (base: [list of decorators])
            The ORtype types that will be combined with logical OR
        ANDtypes : list of string
            The ANDtypes that will be combined with logical AND

        """
        # Implementation identical to atoms apart from what is put in the asSMARTS/asSMIRKS strings

        def __init__(self, ORtypes = None, ANDtypes = None):
            """
            Parameters
            -----------
            ORtypes: list of tuples, optional, default = None
                tuples have form (base, [ORdecorators])
                bond descriptors that will be OR'd together in a SMARTS
            ANDtypes: list of str, optional, default = None
                strings that will be AND'd together in a SMARTS
            index: integer, default = None
                This is for book keeping inside environments and will not be shown in SMARTS or SMIRKS
                example: bond1 in a Bond is the bond between atom1 and atom2
            """
            super(ChemicalEnvironment.Bond,self).__init__(ORtypes, ANDtypes, None, None)
            self._atom = False
            return

        def asSMARTS(self):
            """Return the atom representation as SMARTS.

            Returns
            -------
            smarts : str
                The SMARTS string for just this atom
            """
            if self.ORtypes:
                ORcombos = list()
                for (ORbase, ORdecorators) in self.ORtypes:
                    ORcombos.append(ORbase+''.join(ORdecorators))
                smarts = ','.join(ORcombos)
            else:
                smarts = '~'

            if len(self.ANDtypes) > 0:
                smarts += ';' + ';'.join(self.ANDtypes)

            return smarts

        def asSMIRKS(self):
            """
            Returns
            -------
            smarts : str
                The SMIRKS string for just this bond
            """
            #the same as asSMARTS()
            #    for consistency asSMARTS() or asSMIRKS() can be called
            #    for all environment objects
            return self.asSMARTS()

        def getOrder(self):
            """
            Returns a float for the order of this bond
            for multiple ORtypes or ~ it returns the minimum possible order
            the intended application is for checking valence around a given atom
            """
            # Minimum order for empty ORtypes is 1:
            if not self.ORtypes:
                return 1

            orderDict = {'~':1.,
                    '-':1., ':': 1.5, '=':2., '#':3.,
                    '!-':1.5, '!:':1., '!=':1., '!#':1.}
            orderList = [orderDict[base] for (base, decor) in self.ORtypes]
            return min(orderList)

    @staticmethod
    def validate(smirks, ensure_valence_type=None, toolkit='openeye'):
        """Validate the provided SMIRKS string is valid, and if requested, tags atoms appropriate to the specified valence type.

        Parameters
        ----------
        smirks : str
            The SMIRKS expression to validate
        ensure_valence_type : str, optional, default=None
            If specified, ensure the tagged atoms are appropriate to the specified valence type

        This method will raise a :class:`SMIRKSParsingError` if the provided SMIRKS string is not valid.

        """
        chemenv = ChemicalEnvironment(smirks, toolkit=toolkit)

        if ensure_valence_type:
            valence_type = chemenv.getType()
            if valence_type != ensure_valence_type:
                raise SMIRKSParsingError("Tagged atoms in SMARTS string '%s' specifies valence type '%s', expected '%s'." % (smirks, valence_type, ensure_valence_type))

    def __init__(self, smirks = None, label = None, replacements = None, toolkit='openeye'):
        """Initialize a chemical environment abstract base class.

        smirks = string, optional
            if smirks is not None, a chemical environment is built
            from the provided SMIRKS string
        label = anything, optional
            intended to be used to label this chemical environment
            could be a string, int, or float, or anything
        replacements = list of lists, optional,
            [substitution, smarts] form for parsing SMIRKS
        """
        # TODO: Refactor all this class to use the ToolkitRegistry API.
        if toolkit.lower() == 'openeye' and openforcefield.utils.OpenEyeToolkitWrapper.is_available():
            self.toolkit = 'openeye'
        elif toolkit.lower() == 'rdkit' and openforcefield.utils.RDKitToolkitWrapper.is_available():
            self.toolkit = 'rdkit'
        else:
            raise ValueError("Could not find toolkit {}, please use/install "
                             "openeye or rdkit.".format(toolkit))

        # Define the regular expressions used for all SMIRKS decorators
        # There are a limited number of descriptors for smirks string they are:
        # That is a # followed by one or more ints w/or w/o at ! in front '!#16'
        element_num = r"!?[#]\d+"
        # covers element symbols, i.e. N,C,O,Br not followed by a number
        element_sym = "!?[A-Z][a-z]?"
        # covers element symbols that are aromatic:
        aro_sym = "!?[cnops]"
        # replacement strings
        replace_str = r"\$\w+"
        # a or A w/ or w/o a ! in front 'A'
        aro_ali = "!?[aA]"
        # the decorators (D,H,j,r,V,X,^) followed by one or more integers
        needs_int = r"!?[DHjrVX^]\d+"
        # R(x), +, - do not need to be followed by a integer w/ or w/o a ! 'R2'
        optional_int = r"!?[Rx+-]\d*"
        # chirality options, "@", "@@", "@int" w/ or w/o a ! in front
        chirality = r"!?[@]\d+|!?[@]@?"

        # Generate RegEx string for decorators:
        self.no_bracket_atom_reg = r'('+'|'.join([element_sym, aro_sym, replace_str])+')'
        self.atom_reg = '|'.join([element_num, aro_ali, needs_int,
            optional_int, chirality, replace_str, element_sym, aro_sym])
        self.atom_reg = r'('+self.atom_reg+')'

        # Define bond regular expression options below in order:
        # single, double, triple, aromatic, directional up bond, directional down bond
        # Each can have ! in from and directional can have ? after
        # up and down bonds have lots of \ to fit the python requirements
        self.bond_regs = ['!?[-]', '!?[=]', '!?[#]', '!?[:]', '!?[@]', '!?[\\\\]\\??', '!?[\\/]\\??']
        self.bond_regs = r'('+'|'.join(self.bond_regs)+')'
        # Note, not looking for ~ because that is used for empty bonds

        # Create an empty graph which will store Atom objects.
        self._graph = nx.Graph()
        self.label = label
        self.replacements = replacements

        if smirks is not None:
            # Check that it is a valid SMIRKS
            if not self.isValid(smirks):
                raise SMIRKSParsingError("Error Provided SMIRKS ('%s') was \
not parseable with %s tools" % (smirks, self.toolkit))

            # Check for SMIRKS not supported by Chemical Environments
            if smirks.find('.') != -1:
                raise SMIRKSParsingError("Error: Provided SMIRKS ('%s') \
contains a '.' indicating multiple molecules in the same pattern. This type \
of pattern is not parseable into ChemicalEnvironments" % smirks)
            if smirks.find('>') != -1:
                raise SMIRKSParsingError("Error: Provided SMIRKS ('%s') \
contains a '>' indicating a reaction. This type of pattern is not parseable \
into ChemicalEnvironments." % smirks)

            # try parsing into environment object
            try:
                self._parse_smirks(smirks)
            except:
                raise SMIRKSParsingError("Error SMIRKS (%s) was not parseable\
                        into a ChemicalEnvironment" % smirks)

        # Check that the created Environment is valid
        if not self.isValid():
            raise SMIRKSParsingError("Input SMIRKS (%s), converted to %s \
                    is now invalid" % (smirks, self.asSMIRKS()))

        return

    def _graph_remove_node(self, node):
        self._graph.remove_node(node)
        return True

    def _graph_nodes(self, data=False):
        """
        When data is False returns a list of nodes in graph
        otherwise returns a dictionary in the form {node: data}
        """
        if data:
            return dict(self._graph.nodes(data=True))
        return list(self._graph.nodes())

    def _graph_edges(self, data=False, node=None):
        """
        returns a list of tuples,
        If data is False it has the form [(node1, node2)]
        Otherwise it includes the data [(node1, node2, data_dictionary)]
        If node is not None then the list includes only edges connected to that node
        """
        if node is None:
            return list(self._graph.edges(data=data))
        return list(self._graph.edges(node, data=data))

    def _graph_neighbors(self, node):
        """
        returns a list of neighbors for the given node
        """
        return list(self._graph.neighbors(node))

    def _graph_get_edge_data(self, node1, node2):
        """
        returns a dictionary for the data at the edged connecting
        node1 and node2 in graph
        """
        return self._graph.get_edge_data(node1, node2)

    def isValid(self, smirks=None):
        """
        Returns if the environment is valid, that is if it
        creates a parseable SMIRKS string.
        """
        if smirks is None:
            smirks = self._asSMIRKS()
        if self.toolkit == 'openeye':
            return self._oe_isValid(smirks)
        elif self.toolkit == 'rdkit':
            return self._rdk_isValid(smirks)
        else:
            raise Exception("Could not import openeye.oechem or rdkit.Chem")

    def _rdk_isValid(self, smirks):
        from rdkit import Chem
        if self.replacements is not None:
            for substring, replace_with in self.replacements:
                smirks = smirks.replace(substring, '('+replace_with+')')
        ss = Chem.MolFromSmarts(smirks)
        if ss is None:
            print(smirks, 'not parsed')
        return ss is not None

    def _oe_isValid(self, smirks):
        """
        Returns if the atom is valid, that is if it
        creates a parseable SMIRKS string.
        """
        from openeye import oechem
        qmol = oechem.OEQMol()
        if self.replacements is not None:
            smirks = oechem.OESmartsLexReplace(smirks, self.replacements)
        return oechem.OEParseSmarts(qmol, smirks)

    def _parse_smirks(self,input_smirks):
        """
        This function converts a smirks string to a Chemical Environment
        """
        smirks = _convert_embedded_SMIRKS(input_smirks)
        atoms = dict() # store created atom
        idx = 1 # current atom being created
        store = list() # to store indices while branching
        bondingTo = idx # which atom are we going to bond to

        atom_string, start, end = _find_embedded_brackets(smirks, r'\[', r'\]')

        if start != 0: # first atom is not in square brackets
            if start != -1:
                start_string = smirks[:start]
            else:
                start_string = smirks

            # Check for atoms not between square brackets
            split = re.split(self.no_bracket_atom_reg, start_string)
            atom_string = split[1]

            # update leftover for this condition
            if start != -1: # there is at least 1 more bracketed atom
                leftover = ''.join(split[2:])+smirks[start:]
            else:
                leftover = ''.join(split[2:])

        else: # First atom is in square brackets
            leftover = smirks[end+1:]
            # remove square brackets for parsing
            atom_string = atom_string[1:-1]

        # Check for ring index, i.e. the 1s in "[#6:1]1-CCCCC1"
        match = re.match(r'(\d+)',leftover)
        if match is not None: # leftover starts with int
            ring = re.findall(r'(\d+)',leftover)[0]
            leftover = leftover[match.end():]
        else:
            ring = None

        # Get atom information and create first atom
        OR, AND, index = self._getAtomInfo(atom_string)
        new_atom = self.addAtom(None, newORtypes = OR, newANDtypes = AND,
                newAtomIndex = index, newAtomRing = ring, beyondBeta = True)
        atoms[idx] = new_atom

        while len(leftover) > 0:
            idx += 1

            # Check for branching
            if leftover[0] == ')':
                bondingTo = store.pop()
                leftover = leftover[1:]
                continue

            if leftover[0] == '(':
                store.append(bondingTo)
                leftover = leftover[1:]
                continue

            # find beginning and end of next [atom]
            atom_string, start, end = _find_embedded_brackets(leftover, r'\[', r'\]')

            if start != -1: # no more square brackets
                bond_string = leftover[:start]
            else:
                bond_string = leftover

            # Check for atoms not between square brackets
            bond_split = re.split(self.no_bracket_atom_reg, bond_string)
            # Next atom is not in brackets for example C in "[#7:1]-C"
            if len(bond_split) > 1:
                bond_string = bond_split[0]
                atom_string = '['+bond_split[1]+']'
                # update leftover for this condition
                if start != -1: # ther is at least 1 more bracketed atom
                    leftover = ''.join(bond_split[2:])+leftover[start:]
                else:
                    leftover = ''.join(bond_split[2:])

            else: # next atom is in the brackets [atom]
                # bond and atom string stay the same, update leftover
                leftover = leftover[end+1:]

            # Get bond and atom info
            bOR, bAND = self._getBondInfo(bond_string)
            aOR, aAND, index = self._getAtomInfo(atom_string[1:-1])

            # Check for ring index, i.e. the 1s in "[#6:1]1-CCCCC1"
            match = re.match(r'(\d+)',leftover)
            if match is not None: # leftover starts with int
                ring = re.findall(r'(\d+)',leftover)[0]
                leftover = leftover[match.end():]
            else:
                ring = None

            # create new atom
            new_atom = self.addAtom(atoms[bondingTo], bondORtypes=bOR,
                    bondANDtypes=bAND, newORtypes=aOR, newANDtypes=aAND,
                    newAtomIndex=index, newAtomRing=ring, beyondBeta=True)

            # update state
            atoms[idx] = new_atom
            bondingTo = idx
        return

    def _getAtomInfo(self, atom):
        """
        given atom string, returns ORtypes, ANDtypes, and index
        """
        # Find atom index
        colon = atom.find(':')
        if colon == -1:
            index = None
        else:
            index = int(atom[colon+1:])
            atom = atom[:colon]

        split = atom.split(';')

        # Get ANDtypes (and split them if they don't use ;)
        ANDtypes = list()
        for a in split[1:]:
            ANDtypes += re.findall(self.atom_reg, a)

        # Get ORtypes
        ORList = split[0].split(',')
        ORtypes = list()
        # Separate ORtypes into bases and decorators
        for OR in ORList:
            ORbase, ORdecors = self._separateORtypes(OR)
            if ORbase != None:
                ORtypes.append( (ORbase, ORdecors) )

        return ORtypes, ANDtypes, index

    def _separateORtypes(self, ORtype):
        """
        Separates ORtype (i.e. "#6X4R+0") into
        a base and decorators (i.e. '#6', ['X4','R','+0'] )
        """
        # special case 1: wild card
        if ORtype == '*':
            return None, []

        # if ORbase is a wildcard
        if ORtype[0] == '*':
            return '*', re.findall(self.atom_reg, ORtype[1:])

        # Split up decorators by RegEx strings for atoms
        split = re.findall(self.atom_reg, ORtype)
        if len(split) == 0:
            return None, []

        base = split[0]
        decs = _remove_blanks_repeats(split[1:], ['',base])
        return base, decs

    def _getBondInfo(self, bond):
        """
        given bond strings returns ORtypes and ANDtypes
        """
        # blank bond string is single or aromatic
        # empty ORtypes in Chemical Environments are treated as ~ bonds
        if bond == "":
            ANDtypes = list()
            ORtypes = [ ('-', []), (':', []) ]
            return ORtypes, ANDtypes

        # AND types indicated by ; at the end
        split = bond.split(';')
        ANDtypes = list()
        for a in split[1:]:
            ANDtypes += re.findall(self.bond_regs, a)

        # ORtypes are divided by ,
        ORList = split[0].split(',')
        ORtypes = list()
        for OR in ORList:
            if OR == '~':
                continue
            or_divide = re.findall(self.bond_regs, OR)
            if len(or_divide) > 0:
                ORtypes.append( (or_divide[0], or_divide[1:]))

        return ORtypes, ANDtypes

    def asSMIRKS(self, smarts = False):
        """
        Returns a SMIRKS representation of the chemical environment

        Parameters
        -----------
        smarts: optional, boolean
            if True, returns a SMARTS instead of SMIRKS without index labels
        """
        return self._asSMIRKS(None, None, smarts)

    def _asSMIRKS(self, initialAtom = None, neighbors = None, smarts = False):
        """Return a SMIRKS representation of the chemical environment.

        Parameters
        ----------
        initalAtom = optional, atom object
            This is randomly selected if not chosen.
        neighbors = optional, list of atom objects
            This is all of the initalAtom neighbors if not specified
            generally this is used only for the recursive calls
            so initial atoms are not reprinted
        smarts = optional, boolean
            if True, returns a SMARTS string instead of SMIRKS
        """
        # If empty chemical environment
        if len(self._graph_nodes()) == 0:
            return ""

        if initialAtom is None:
            initialAtom = self.getAtoms()[0]

        if neighbors is None:
            neighbors = self._graph_neighbors(initialAtom)

        # sort neighbors to guarantee order is constant
        neighbors = sorted(neighbors, key=lambda atom: atom.asSMIRKS())

        # initialize smirks for starting atom
        if smarts:
            smirks = initialAtom.asSMARTS()
        else:
            smirks = initialAtom.asSMIRKS()

        # loop through neighbors
        for idx, neighbor in enumerate(neighbors):
            # get the SMIRKS for the bond between these atoms
            # bonds are the same if smarts or smirks
            bond_edge = self._graph_get_edge_data(initialAtom, neighbor)
            bondSMIRKS = bond_edge['bond'].asSMIRKS()

            # Get the neighbors for this neighbor
            new_neighbors = self._graph_neighbors(neighbor)
            # Remove initialAtom so it doesn't get reprinted
            new_neighbors.remove(initialAtom)

            # Call asSMIRKS again to get the details for that atom
            atomSMIRKS = self._asSMIRKS(neighbor, new_neighbors, smarts)

            # Use ( ) for branch atoms (all but last)
            if idx < len(neighbors) - 1:
                smirks += '(' + bondSMIRKS + atomSMIRKS + ')'
            # This is for the atoms that are a part of the main chain
            else:
                smirks += bondSMIRKS + atomSMIRKS

        return smirks

    def selectAtom(self, descriptor = None):
        """Select a random atom fitting the descriptor.

        Parameters
        ----------
        descriptor: optional, None
            None - returns any atom with equal probability
            int - will return an atom with that index
            'Indexed' - returns a random indexed atom
            'Unindexed' - returns a random unindexed atom
            'Alpha' - returns a random alpha atom
            'Beta' - returns a random beta atom

        Returns
        --------
        a single Atom object fitting the description
        or None if no such atom exists
        """
        if descriptor is None:
            return random.choice(self._graph_nodes())

        # TODO: Is there a better way to do this?
        try:
            descriptor = int(descriptor)
        except:
            pass

        if type(descriptor) is int:
            for atom in self.getAtoms():
                if atom.index == descriptor:
                    return atom
            return None

        atoms = self.getComponentList('atom',descriptor)
        if len(atoms) == 0:
            return None

        return random.choice(atoms)

    def getComponentList(self, component_type, descriptor = None):
        """
        Returns a list of atoms or bonds matching the descriptor

        Parameters
        ----------
        component_type: string: 'atom' or 'bond'
        descriptor: string, optional
            'all', 'Indexed', 'Unindexed', 'Alpha', 'Beta'

        Returns
        -------
        """
        if descriptor is not None:
            d = descriptor.lower()
        else:
            d = None

        if not component_type.lower() in ['atom', 'bond']:
            raise Exception("Error: 'getComponentList()' component_type must be 'atom' or 'bond'")

        if component_type.lower() == 'atom':
            if d == 'indexed':
                return self.getIndexedAtoms()
            elif d == 'unindexed':
                return self.getUnindexedAtoms()
            elif d == 'alpha':
                return self.getAlphaAtoms()
            elif d == 'beta':
                return self.getBetaAtoms()
            else:
                return self.getAtoms()

        elif component_type.lower() == 'bond':
            if d == 'indexed':
                return self.getIndexedBonds()
            elif d == 'unindexed':
                return self.getUnindexedBonds()
            elif d == 'alpha':
                return self.getAlphaBonds()
            elif d == 'beta':
                return self.getBetaBonds()

            return self.getBonds()

        return None

    def selectBond(self, descriptor = None):
        """Select a random bond fitting the descriptor.

        Parameters
        ----------
        descriptor: optional, None
            None - returns any bond with equal probability
            int - will return an bond with that index
            'Indexed' - returns a random indexed bond
            'Unindexed' - returns a random unindexed bond
            'Alpha' - returns a random alpha bond
            'Beta' - returns a random beta bond

        Returns
        -------
        a single Bond object fitting the description
        or None if no such atom exists
        """
        # TODO: Is there a better way to do this?
        try:
            descriptor = int(descriptor)
        except:
            pass

        if type(descriptor) is int:
            for bond in self.getBonds():
                if bond._bond_type == descriptor:
                    return bond
            return None

        bonds = self.getComponentList('bond', descriptor)
        if len(bonds) == 0:
            return None

        return random.choice(bonds)

    def addAtom(self, bondToAtom, bondORtypes= None, bondANDtypes = None,
            newORtypes = None, newANDtypes = None, newAtomIndex = None,
            newAtomRing = None, beyondBeta = False):
        """Add an atom to the specified target atom.

        Parameters
        -----------
        bondToAtom: atom object, required
            atom the new atom will be bound to
        bondORtypes: list of tuples, optional
            strings that will be used for the ORtypes for the new bond
        bondANDtypes: list of strings, optional
            strings that will be used for the ANDtypes for the new bond
        newORtypes: list of strings, optional
            strings that will be used for the ORtypes for the new atom
        newANDtypes: list of strings, optional
            strings that will be used for the ANDtypes for the new atom
        newAtomIndex: int, optional
            integer label that could be used to index the atom in a SMIRKS string
        beyondBeta: boolean, optional
            if True, allows bonding beyond beta position

        Returns
        --------
        newAtom: atom object for the newly created atom
        """
        if bondToAtom is None:
            if len(self._graph_nodes()) > 0:
                return None
            newType = newAtomIndex
            if newType is None:
                newType = 0

            newAtom = self.Atom(newORtypes, newANDtypes, newAtomIndex, newAtomRing)
            self._graph.add_node(newAtom, atom_type = newType)
            return newAtom

        # Check if we can get past beta position
        bondToType = self._graph_nodes(data=True)[bondToAtom]['atom_type']
        if bondToType < 0 and not beyondBeta:
            return None

        # determine the type integer for the new atom and bond
        if newAtomIndex != None:
            newType = newAtomIndex
            bondType = max(newType, bondToType) - 1
        else:
            if bondToType > 0:
                newType = 0
            else:
                newType = bondToType - 1
            bondType = newType

        # create new bond
        newBond = self.Bond(bondORtypes, bondANDtypes)

        # create new atom
        newAtom = self.Atom(newORtypes, newANDtypes, newAtomIndex, newAtomRing)

        # Add node for newAtom
        self._graph.add_node(newAtom, atom_type = newType)

        # Connect original atom and new atom
        self._graph.add_edge(bondToAtom, newAtom, bond = newBond, bond_type = bondType)
        newBond._bond_type = bondType

        return newAtom

    def removeAtom(self, atom, onlyEmpty = False):
        """Remove the specified atom from the chemical environment.
        if the atom is not indexed for the SMIRKS string or
        used to connect two other atoms.

        Parameters
        ----------
        atom: atom object, required
            atom to be removed if it meets the conditions.
        onlyEmpty: boolean, optional
            True only an atom with no ANDtypes and 1 ORtype can be removed

        Returns
        --------
        Boolean True: atom was removed, False: atom was not removed
        """
        # labeled atoms can't be removed
        if atom.index is not None:
            return False

        # Atom connected to more than one other atom cannot be removed
        if len(self._graph_neighbors(atom)) > 1:
            return False

        # if you can remove "decorated atoms" remove it
        if not onlyEmpty:
            self._graph_remove_node(atom)
            return True

        if len(atom.ANDtypes) > 0:
            return False
        elif len(atom.ORtypes) > 1:
            return False

        self._graph_remove_node(atom)
        return True

    def getAtoms(self):
        """
        Returns
        -------
        list of atoms in the environment
        """
        return self._graph_nodes()

    def getBonds(self, atom = None):
        """
        Parameters
        ----------
        atom: Atom object, optional, returns bonds connected to atom
        returns all bonds in fragment if atom is None

        Returns
        --------
        a complete list of bonds in the fragment
        """
        if atom is None:
            edge_list = self._graph_edges(data=True)
            bonds = [data['bond'] for a1, a2, data in edge_list]
        else:
            bonds = []
            for (a1, a2, info) in self._graph_edges(data=True, node=atom):
                bonds.append(info['bond'])

        return bonds

    def getBond(self, atom1, atom2):
        """
        Get bond betwen two atoms

        Parameters
        -----------
        atom1 and atom2: atom objects

        Returns
        --------
        bond object between the atoms or None if no bond there
        """
        if atom2 in self._graph_neighbors(atom1):
            return self._graph_get_edge_data(atom1, atom2)['bond']
        else:
            return None

    def getIndexedAtoms(self):
        """
        returns the list of Atom objects with an index
        """
        index_atoms = []
        for atom, info in self._graph_nodes(data=True).items():
            if info['atom_type'] > 0:
                index_atoms.append(atom)
        return index_atoms

    def getUnindexedAtoms(self):
        """
        returns a list of Atom objects that are not indexed
        """
        unindexed_atoms = []
        for atom, info in self._graph_nodes(data=True).items():
            if info['atom_type'] < 1:
                unindexed_atoms.append(atom)
        return unindexed_atoms

    def getAlphaAtoms(self):
        """
        Returns a list of atoms alpha to any indexed atom
            that are not also indexed
        """
        alpha_atoms = []
        for atom, info in self._graph_nodes(data=True).items():
            if info['atom_type'] == 0:
                alpha_atoms.append(atom)

        return alpha_atoms

    def getBetaAtoms(self):
        """
        Returns a list of atoms beta to any indexed atom
            that are not alpha or indexed atoms
        """
        beta_atoms = []
        for atom, info in self._graph_nodes(data=True).items():
            if info['atom_type'] == -1:
                beta_atoms.append(atom)
        return beta_atoms

    def getIndexedBonds(self):
        """
        Returns a list of Bond objects that connect two indexed atoms
        """
        indexedBonds = []
        for bond in self.getBonds():
            if bond._bond_type > 0:
                indexedBonds.append(bond)
        return indexedBonds

    def getUnindexedBonds(self):
        """
        Returns a list of Bond objects that connect
            an indexed atom to an unindexed atom
            two unindexed atoms
        """
        unindexedBonds = []
        for bond in self.getBonds():
            if bond._bond_type < 1:
                unindexedBonds.append(bond)
        return unindexedBonds

    def getAlphaBonds(self):
        """
        Returns a list of Bond objects that connect
            an indexed atom to alpha atoms
        """
        alphaBonds = []
        for bond in self.getBonds():
            if bond._bond_type == 0:
                alphaBonds.append(bond)
        return alphaBonds

    def getBetaBonds(self):
        """
        Returns a list of Bond objects that connect
            alpha atoms to beta atoms
        """
        betaBonds = []
        for bond in self.getBonds():
            if bond._bond_type == -1:
                betaBonds.append(bond)
        return betaBonds

    def isAlpha(self, component):
        """
        Takes an atom or bond are returns True if it is alpha to an indexed atom
        """
        if component._atom:
            return self._graph_nodes(data=True)[component]['atom_type'] == 0
        else:
            return component._bond_type == 0

    def isUnindexed(self, component):
        """
        returns True if the atom or bond is not indexed
        """
        if component._atom:
            return component.index is None
        else:
            return component._bond_type < 1

    def isIndexed(self, component):
        """
        returns True if the atom or bond is indexed
        """
        if component._atom:
            return component.index != None
        else:
            return component._bond_type > 0

    def isBeta(self, component):
        """
        Takes an atom or bond are returns True if it is beta to an indexed atom
        """
        if component._atom:
            return self._graph_nodes(data=True)[component]['atom_type'] == -1
        else:
            return component._bond_type == -1

    # TODO: We may want to overhaul ChemicalEnvironment.getType() to return one of ['atom', 'bond', 'angle', 'proper', 'improper']
    # and check to make sure the expected connectivity is represented in the SMIRKS expression.
    def getType(self):
        """
        Uses number of indexed atoms and bond connectivity
        to determine the type of chemical environment

        Returns
        -------
        chemical environemnt type:
            'Atom', 'Bond', 'Angle', 'ProperTorsion', 'ImproperTorsion'
            None if number of indexed atoms is 0 or > 4
        """
        index_atoms = self.getIndexedAtoms()
        natoms = len(index_atoms)

        if natoms == 1:
            return "Atom"
        if natoms == 2:
            return "Bond"
        if natoms == 3:
            return "Angle"
        if natoms == 4:
            atom2 = self.selectAtom(2)
            atom4 = self.selectAtom(4)
            bond24 = self.getBond(atom2, atom4)
            if bond24 != None:
                return "ImproperTorsion"
            return "ProperTorsion"
        else:
            return None

    def getNeighbors(self, atom):
        """
        Returns atoms that are bound to the given atom
        in the form of a list of Atom objects
        """
        return self._graph_neighbors(atom)

    def getValence(self, atom):
        """
        Returns the valence (number of neighboring atoms)
        around the given atom
        """
        return len(self._graph_neighbors(atom))

    def getBondOrder(self, atom):
        """
        Returns minimum bond order around a given atom
        0 if atom has no neighbors
        aromatic bonds count as 1.5
        any bond counts as 1.0
        """
        order = 0.
        for a1, a2, info in self._graph_edges(data=True, node=atom):
            order += info['bond'].getOrder()
        return order

class AtomChemicalEnvironment(ChemicalEnvironment):
    """Chemical environment matching one labeled atom.

    """
    def __init__(self, smirks = "[*:1]", label = None, replacements = None, toolkit='openeye'):
        """Initialize a chemical environment corresponding to matching a single atom.

        Parameters
        -----------
        smirks: string, optional
            the default is an empty Atom corresponding to "[*:1]"
        label = anything, optional
            intended to be used to label this chemical environment
            could be a string, int, or float, or anything
        replacements = list of lists, optional,
            [substitution, smarts] form for parsing SMIRKS

        For example:
            # create an atom that is carbon, nitrogen, or oxygen with no formal charge
            atom = AtomChemicalEnvironment([['#6', '#7', '#8'], ['+0']])
            print atom.asSMIRKS()
            # prints: "[#6,#7,#8;+0:1]"
        """
        # Initialize base class
        super(AtomChemicalEnvironment,self).__init__(smirks, label, replacements, toolkit)
        correct, expected = self._checkType()
        if not correct:
            assigned = self.getType()
            raise SMIRKSMismatchError("The SMIRKS (%s) was assigned the type %s when %s was expected" % (smirks, assigned, expected))
        self.atom1 = self.selectAtom(1)

    def _checkType(self):
        return (self.getType() == 'Atom'), 'Atom'

    def asSMIRKS(self, smarts = False):
        """
        Returns a SMIRKS representation of the chemical environment

        Parameters
        -----------
        smarts: optional, boolean
            if True, returns a SMARTS instead of SMIRKS without index labels
        """
        return self._asSMIRKS(self.atom1, None, smarts)

    def asAtomtypeSMARTS(self):
        """
        Makes SMARTS string for one atom

        Returns
        --------
        string for the SMARTS string for the first atom (labeled with index :1)

        This is a single atom with neighbors as decorators in the form:
        [atom1$(*~neighbor1)$(*~neighbor2)...]
        """
        # smarts for atom1 without last ']'
        smarts = self.atom1.asSMARTS()[:-1]

        for idx, neighbor in enumerate(self._graph_neighbors(self.atom1)):
            new_neighbors = self._graph_neighbors(neighbor)
            new_neighbors.remove(self.atom1)

            bondSMARTS = self._graph_get_edge_data(self.atom1, neighbor)['bond'].asSMARTS()
            neighborSMARTS = self._asSMIRKS(neighbor, new_neighbors, True)

            smarts += '$(*' + bondSMARTS + neighborSMARTS + ')'

        return smarts + ']'

class BondChemicalEnvironment(AtomChemicalEnvironment):
    """Chemical environment matching two labeled atoms (or a bond).
    """
    def __init__(self, smirks = "[*:1]~[*:2]", label = None, replacements = None, toolkit='openeye'):
        """Initialize a chemical environment corresponding to matching two atoms (bond).

        Parameters
        -----------
        smirks: string, optional
            the default is an empty Bond corresponding to "[*:1]~[*:2]"
        label = anything, optional
            intended to be used to label this chemical environment
            could be a string, int, or float, or anything
        replacements = list of lists, optional,
            [substitution, smarts] form for parsing SMIRKS

        """
        # Initialize base class
        super(BondChemicalEnvironment,self).__init__(smirks, label, replacements, toolkit)

        # Add initial atom
        self.atom2 = self.selectAtom(2)
        if self.atom2 is None:
            raise Exception("Error: Bonds need 2 indexed atoms, there were not enough in %s" % smirks)

        self.bond2 = self._graph_get_edge_data(self.atom1, self.atom2)['bond']

    def _checkType(self):
        return (self.getType() == 'Bond'), 'Bond'

class AngleChemicalEnvironment(BondChemicalEnvironment):
    """Chemical environment matching three marked atoms (angle).
    """
    def __init__(self, smirks = "[*:1]~[*:2]~[*:3]", label = None, replacements = None, toolkit='openeye'):

        """Initialize a chemical environment corresponding to matching three atoms.

        Parameters
        -----------
        smirks: string, optional
            the default is an empty Angle corresponding to "[*:1]~[*:2]~[*:3]"
        label = anything, optional
            intended to be used to label this chemical environment
            could be a string, int, or float, or anything
        replacements = list of lists, optional,
            [substitution, smarts] form for parsing SMIRKS
        """
        # Initialize base class
        super(AngleChemicalEnvironment,self).__init__(smirks, label, replacements, toolkit)

        # Add initial atom
        self.atom3 = self.selectAtom(3)
        self.bond2 = self._graph_get_edge_data(self.atom2, self.atom3)['bond']

    def _checkType(self):
        return (self.getType() == 'Angle'), 'Angle'

class TorsionChemicalEnvironment(AngleChemicalEnvironment):
    """Chemical environment matching four marked atoms (torsion).
    """
    def __init__(self, smirks = "[*:1]~[*:2]~[*:3]~[*:4]", label = None, replacements = None, toolkit='openeye'):
        """Initialize a chemical environment corresponding to matching four atoms (torsion).

        Parameters
        -----------
        smirks: string, optional
            the default is an empty Torsion corresponding to
            "[*:1]~[*:2]~[*:3]~[*:4]"
        label = anything, optional
            intended to be used to label this chemical environment
            could be a string, int, or float, or anything
        replacements = list of lists, optional,
            [substitution, smarts] form for parsing SMIRKS
        """
        # Initialize base class
        super(TorsionChemicalEnvironment,self).__init__(smirks, label, replacements, toolkit)

        # Add initial atom
        self.atom4 = self.selectAtom(4)
        self.bond3 = self._graph_get_edge_data(self.atom3, self.atom4)['bond']

    def _checkType(self):
        return (self.getType() == 'ProperTorsion'), 'ProperTorsion'

class ImproperChemicalEnvironment(AngleChemicalEnvironment):
    """Chemical environment matching four marked atoms (improper).
    """
    def __init__(self, smirks = "[*:1]~[*:2](~[*:3])~[*:4]", label = None, replacements = None, toolkit='openeye'):
        """Initialize a chemical environment corresponding four atoms (improper).

        Parameters
        -----------
        smirks: string, optional
            the default is an empty Improper corresponding to
            "[*:1]~[*:2](~[*:3])~[*:4]"
        label = anything, optional
            intended to be used to label this chemical environment
            could be a string, int, or float, or anything
        """
        # Initialize base class
        super(ImproperChemicalEnvironment,self).__init__(smirks, label, replacements, toolkit)

        # Add initial atom
        self.atom4 = self.selectAtom(4)
        self.bond3 = self._graph_get_edge_data(self.atom2, self.atom4)['bond']

    def _checkType(self):
        return (self.getType() == 'Improper'), 'Improper'
