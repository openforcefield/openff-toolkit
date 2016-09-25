#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
smarty.py

Example illustrating a scheme to create and destroy atom types automatically using SMARTS.

AUTHORS

John Chodera <john.chodera@choderalab.org>, Memorial Sloan Kettering Cancer Center.
Additional contributions from the Mobley lab, UC Irvine, including David Mobley, Caitlin Bannan, and Camila Zanette.

The AtomTyper class is based on 'patty' by Pat Walters, Vertex Pharmaceuticals.

"""
#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import sys
import string

from optparse import OptionParser # For parsing of command line arguments

import os
import math
import copy
import re
import numpy
import random

import openeye.oechem
import openeye.oeomega
import openeye.oequacpac

from openeye.oechem import *
from openeye.oeomega import *
from openeye.oequacpac import *

import networkx

import time

from . import AtomTyper
from score_utils import load_trajectory
from score_utils import scores_vs_time

#=============================================================================================
# ATOMTYPE SAMPLER
#=============================================================================================

class AtomTypeSamplerElemental(object):
    """
    Atom type sampler.

    """
    def __init__(self, molecules, basetypes_filename, initialtypes_filename, decorators_filename, replacements_filename=None, reference_typed_molecules=None, temperature=0.1, verbose=False, decorator_behavior='combinatorial-decorators', element=None):
        """
        Initialize an atom type sampler.

        ARGUMENTS

        molecules : list of molecules for typing
            List of molecules for typing
        basetypes_filename : str
            File defining base/generic atom types (which cannot be destroyed); often these are elemental types
        initialtypes_filename : 
            File defining initial atom types (which CAN be destroyed, except for those which occur in basetypes_filename
        decorators_filename : str
            File containing decorators that can be added to existing types to generate subtypes
        replacements_filename : str, optional, default=None
            If specified, SMARTS replacement definitions will be read from this file
        reference_typed_molecules : list of OEMol, optional, default=None
            List of molecules with reference types for use in Monte Carlo acceptance.
            If specified, the likelihood function will utilize the maximal number of matched atom types with these molecules.
            If not specified, no likelihood function will be employed.
        temperature : float, optional, default=0.1
            Temperature for Monte Carlo acceptance/rejection
        verbose : bool, optional, default=False
            If True, verbose output will be printed.

        Notes
        -----
        This is just a proof of concept.  No scoring of molecular properties is performed.

        """
        
        
        # Choose one element per run to find its atomtypes (string format)
        self.initial_element = element

        self.verbose = verbose
        
        self.decorator_behavior = decorator_behavior

        # Define internal typing tag.
        self.typetag = 'atomtype'

        # Save bond list to use throughout
        self.bondset = [("-","simply"), ("=", "doubly"), ("#","triply"), (":", "aromatic"), ("~","any")]

        # Read atomtypes (initial and base) and decorators.
        self.atomtypes = AtomTyper.read_typelist(initialtypes_filename)
        self.basetypes = AtomTyper.read_typelist(basetypes_filename)
        self.decorators = AtomTyper.read_typelist(decorators_filename)
        self.replacements = AtomTyper.read_typelist(replacements_filename)
        # Try to ensure base/initial types have unique names as name 
        # clashes between initial and target types will cause problems
        for idx, [smarts, typename] in enumerate(self.atomtypes):
            self.atomtypes[idx] = [smarts, 'c_'+typename]
        for idx, [smarts, typename] in enumerate(self.basetypes):
            self.basetypes[idx] = [smarts, 'c_'+typename]
        self.unmatched_atomtypes = copy.deepcopy(self.atomtypes)

        # Store smarts for basetypes
        self.basetypes_smarts = [ smarts for (smarts, name) in self.basetypes ]

        # Ensure all base types are in initial types (and add if not) as 
        # base types are generics (such as elemental) and need to be present 
        # at the start

        initial_smarts = [ smarts for (smarts, name) in self.atomtypes ]
        for [smarts, typename] in self.basetypes:
            if smarts not in initial_smarts:
                self.atomtypes = [[smarts, typename]] + self.atomtypes
                if self.verbose: print("Added base (generic) type `%s`, name %s, to initial types." % (smarts, typename) )
        # Store initially populated base types, as these will be retained even 
        # if they have zero occupancy (whereas unpopulated base types
        # need never be used ever and can be deleted- i.e. if we have no 
        # phosphorous in the set we don't need a phosphorous base type)
        self.used_basetypes = []
        self.used_atomtypes = []
        
        # Creat dictionary to store children of initial atom types
        #self.parents = dict()
        #for [smarts, typename] in self.atomtypes:
        #    #store empty list of chlidren for each atomtype
        #    self.parents[smarts] = []

        # Store a deep copy of the molecules since they will be annotated
        self.molecules = copy.deepcopy(molecules)

        # Type all molecules with current typelist to ensure that starting types are sufficient.
        self.type_molecules(self.atomtypes, self.molecules)

        # Compute atomtype statistics on molecules.
        [atom_typecounts, molecule_typecounts] = self.compute_type_statistics(self.atomtypes, self.molecules)
        if self.verbose: self.show_type_statistics(self.atomtypes, atom_typecounts, molecule_typecounts)
        # For use later, also see which base types are used (get those stats) - which means I need to type a copy of molecules then recompute stats
        tmpmolecules = copy.deepcopy(molecules)
        self.type_molecules(self.basetypes, tmpmolecules)
        [ basetype_typecounts, molecule_basetype_typecounts] = self.compute_type_statistics( self.basetypes, tmpmolecules )

        # Compute total atoms
        self.total_atoms = 0.0
        for molecule in self.molecules:
            for atom in molecule.GetAtoms():
                self.total_atoms += 1.0

        # Compute total atoms of the specif element
        self.total_element = self.calculate_number_element(self.atomtypes, self.molecules)

        # Store reference molecules
        self.reference_typed_molecules = None
        self.reference_atomtypes = set()
        self.current_atom_matches = None
        self.temperature = temperature
        if reference_typed_molecules is not None:
            self.reference_typed_molecules = copy.deepcopy(reference_typed_molecules)
            # Extract list of reference atom types
            for molecule in reference_typed_molecules:
                for atom in molecule.GetAtoms():
                    atomtype = atom.GetType()
                    self.reference_atomtypes.add(atomtype)
            self.reference_atomtypes = list(self.reference_atomtypes)
            # Compute current atom matches
            [self.atom_type_matches, self.total_atom_type_matches] = self.best_match_reference_types(self.atomtypes, self.molecules)
            # Count atom types.
            self.reference_atomtypes_atomcount = { atomtype : 0 for atomtype in self.reference_atomtypes }
            for molecule in reference_typed_molecules:
                for atom in molecule.GetAtoms():
                    atomtype = atom.GetType()
                    self.reference_atomtypes_atomcount[atomtype] += 1
        
        # Maintain a list of SMARTS matches without any atom type matches in the dataset
        # This is used for efficiency.
        self.atomtypes_with_no_matches = set()

        # Track used vs unused base types - unused base types are not retained
        for (smarts, atom_type) in self.basetypes:
            # If this type is used, then track it
            if basetype_typecounts[atom_type] > 0:
                self.used_basetypes.append( [ smarts, atom_type] )
                if self.verbose: print("Storing used base type `%s`, name `%s` with count %s..." % (smarts, atom_type, atom_typecounts[atom_type] )) 
            # If unused, it matches nothing in the set
            else:  
                self.atomtypes_with_no_matches.add( smarts )
                if self.verbose: print("Storing base atom type `%s`, which is unused, so that it will not be proposed further." % smarts )
        # Atom basetypes to create new smart strings
        self.atom_basetype = copy.deepcopy(self.used_basetypes)
        
        # Track used vs unused atomtypes - unused atomtypes are not retained
        for (smarts, atom_type) in self.atomtypes:
            # If this type is used, then track it
            if atom_typecounts[atom_type] > 0:
                self.used_atomtypes.append( [ smarts, atom_type] )
                if self.verbose: print("Storing used atomtype `%s`, name `%s` with count %s..." % (smarts, atom_type, atom_typecounts[atom_type] )) 
            # If unused, it matches nothing in the set
            elif atom_typecounts[atom_type] == 0 and (smarts not in self.basetypes_smarts): 
                self.atomtypes_with_no_matches.add( smarts )
                if self.verbose: print("Storing base atom type `%s`, which is unused, so that it will not be proposed further." % smarts )
        
        # Atom basetypes to create new smart strings
        self.atom_basetype = copy.deepcopy(self.used_basetypes)
        self.atomtypes = copy.deepcopy(self.used_atomtypes)

        # Track unused initial types that are not base types as we also don't 
        # need to retain those
        #for (smarts, atom_type) in self.atomtypes:
        #    if atom_typecounts[atom_type] == 0 and (smarts not in self.basetypes_smarts):
        #        self.atomtypes_with_no_matches.add( smarts )
        #        if self.verbose: print("Storing initial atom type `%s`, which is unused, so that it will not be proposed further." % smarts )   

        if self.atomtypes != self.atom_basetype:
            self.atomtypes = self.atomtypes + self.atom_basetype

        self.newatomtypes = []
        for [smarts, typename] in self.atomtypes:
            element = re.findall('\d+', smarts)[0]
            if element == self.initial_element:
                self.newatomtypes += [[smarts, typename]]
        self.atomtypes = copy.deepcopy(self.newatomtypes)

        # Creat dictionary to store children of initial atom types
        self.parents = dict()
        for [smarts, typename] in self.newatomtypes:
            #store empty list of chlidren for each atomtype
            self.parents[smarts] = []

        return

    def best_match_reference_types(self, atomtypes, molecules):
        """
        Determine best match for each parameter with reference atom types

        Parameters
        ----------
        atomtypes :
            Current atom types
        molecules : list of OEMol
            Typed molecules, where types are stored in self.atomtypetag string data.

        Returns
        -------
        atom_type_matches : list of tuples (current_atomtype, reference_atomtype, counts)
            Best correspondence between current and reference atomtypes, along with number of atoms equivalently typed in reference molecule set.
        total_atom_type_matches : int
            The total number of correspondingly typed atoms in the reference molecule set.

        * Currently, types for reference typed molecules are accessed via atom.GetType(), while types for current typed molecules are accessed via atom.GetStringData(self.typetag).
          This should be homogenized.

        Contributor:
        * Josh Fass <josh.fass@choderalab.org> contributed this algorithm.

        """
        if self.reference_typed_molecules is None:
            if self.verbose: print('No reference molecules specified, so skipping likelihood calculation.')
            return None

        # Create bipartite graph (U,V,E) matching current atom types U with reference atom types V via edges E with weights equal to number of atoms typed in common.
        if self.verbose: print('Creating graph matching current atom types with reference atom types...')
        initial_time = time.time()
        import networkx as nx
        graph = nx.Graph()

        # Get current atomtypes and reference atom types
        current_atomtypes = [ typename for (smarts, typename) in atomtypes ]
        reference_atomtypes = [ typename for typename in self.reference_atomtypes ]
        # check that current atom types are not in reference atom types
        if set(current_atomtypes) & set(reference_atomtypes):
            raise Exception("Current and reference atom types must be unique")
        # Add current atom types
        for atomtype in current_atomtypes:
            graph.add_node(atomtype, bipartite=0)
        # Add reference atom types
        for atomtype in reference_atomtypes:
            graph.add_node(atomtype, bipartite=1)
        # Add edges.
        atoms_in_common = dict()
        for current_atomtype in current_atomtypes:
            for reference_atomtype in reference_atomtypes:
                atoms_in_common[(current_atomtype,reference_atomtype)] = 0
        for (current_typed_molecule, reference_typed_molecule) in zip(molecules, self.reference_typed_molecules):
            for (current_typed_atom, reference_typed_atom) in zip(current_typed_molecule.GetAtoms(), reference_typed_molecule.GetAtoms()):
                current_atomtype = current_typed_atom.GetStringData(self.typetag)
                reference_atomtype = reference_typed_atom.GetType()
                if (current_atomtype,reference_atomtype) in atoms_in_common.keys():
                    atoms_in_common[(current_atomtype,reference_atomtype)] += 1
        for current_atomtype in current_atomtypes:
            for reference_atomtype in reference_atomtypes:
                weight = atoms_in_common[(current_atomtype,reference_atomtype)]
                graph.add_edge(current_atomtype, reference_atomtype, weight=weight)
        elapsed_time = time.time() - initial_time
        if self.verbose: print('Graph creation took %.3f s' % elapsed_time)

        # Compute maximum match
        if self.verbose: print('Computing maximum weight match...')
        initial_time = time.time()
        mate = nx.algorithms.max_weight_matching(graph, maxcardinality=False)
        elapsed_time = time.time() - initial_time
        if self.verbose: print('Maximum weight match took %.3f s' % elapsed_time)

        # Compute match dictionary and total number of matches.
        atom_type_matches = list()
        total_atom_type_matches = 0
        for current_atomtype in current_atomtypes:
            if current_atomtype in mate:
                reference_atomtype = mate[current_atomtype]
                counts = graph[current_atomtype][reference_atomtype]['weight']
                total_atom_type_matches += counts
                atom_type_matches.append( (current_atomtype, reference_atomtype, counts) )
            else:
                atom_type_matches.append( (current_atomtype, None, None) )

        # Report on matches
        if self.verbose:
            print("PROPOSED:")
            self.show_type_matches(atom_type_matches)

        return (atom_type_matches, total_atom_type_matches)

    def show_type_matches(self, atom_type_matches):
        """
        Show pairing of current to reference atom types.

        atom_type_matches : list of (current_atomtype, reference_atomtype, counts)
            List of atom type matches.

        Returns fraction_matched_atoms, the fractional count of matched atoms

        """
        print('Atom type matches:')
        total_atom_type_matches = 0
        for (current_atomtype, reference_atomtype, counts) in atom_type_matches:
            if reference_atomtype is not None:
                print('%-64s matches %8s : %8d atoms matched' % (current_atomtype, reference_atomtype, counts))
                total_atom_type_matches += counts
            else:
                print('%-64s         no match' % (current_atomtype))

        #fraction_matched_atoms = float(total_atom_type_matches) / float(self.total_atoms)
        fraction_matched_atoms = float(total_atom_type_matches) / float(self.total_element)
        #print('%d / %d total atoms match (%.3f %%)' % (total_atom_type_matches, self.total_atoms, fraction_matched_atoms * 100))
        print('%d / %d total atoms match (%.3f %%)' % (total_atom_type_matches, self.total_element, fraction_matched_atoms * 100))

        return fraction_matched_atoms
    
    
    def AtomDecorator(self, atom1type, decorator):
        """
        Given an atom and a decorator ammend the SMARTS string with that decorator 
        Returns at "atom" which is a tuple of [SMARTS, typeName] for the decorated atom
        """
        if self.HasAlpha(atom1type):
            # decorators should go before the $ sign on the atom
            dollar = atom1type[0].find('$')
            proposed_atomtype = atom1type[0][:dollar] + decorator[0] + atom3[0][dollar:]
            proposed_typename = atom1type[1] + ' ' + decorator[1]
        else:
            # No alpha atom so the decorator goes before the ']'
            proposed_atomtype = atom1type[0][:-1] + decorator[0] + ']'
            proposed_typename = atom1type[1] + ' '  + decorator[1]
        return proposed_atomtype, proposed_typename

    def PickAnAtom(self, atomList):
        """
        takes a desired current set of atomtypes (current, base, populated base, or some combination of the 
        above or something else) and pick one at random. (Here, we will typically but perhaps not always want 
        to use currently populated atomtypes, i.e. the atomtypes list, not basetypes (generic types) or 
        used_basetypes (base types which type anything)). A little silly, but still helpful. 
        Returns a tuple corresponding to the selected (smarts, name).
        """
        atomtype_index = random.randint(0, len(atomList)-1)
        return atomList[atomtype_index]

    def HasAlpha(self, atom1type):
        """
        Takes a specified atomtype tuple (smarts, name) and determines whether or not it already has an alpha 
        substituent, returning True or False.
        """
        if atom1type[0].find("$") != -1:
            return True
        else:
            return False

    def AddAlphaSubstituentAtom(self, atom1type, bondset, atom2type, first_alpha):
        """
        Takes specified atomtypes atom1type and atom2type (where atom1type is a type without an alpha substituent), 
        and a specified bond set, and introduces an alpha substituent involving atom2type (which can be a decorated 
        type such as output by AtomDecorator); returns a tuple of (smarts, name) from the result. Example output 
        for input of atom1type=("[#1]", "hydrogen"), bondset = (['-'], "single"), and atom2type = ("[#6]" is 
        ("[#1$(*-[#6])]", "hydrogen singly bonded to carbon") or something along those lines. This should basically 
        just be adding $(*zw) where z is the selected bond type and w is atom2type. It should raise an exception 
        if an alpha substituent is attempted to be added to an atom1type which already has an alpha substituent.
        """
        if first_alpha:
            result = re.match('\[(.+)\]', atom1type[0])
            proposed_atomtype = '[' + result.groups(1)[0] + '$(*' + bondset[0] + atom2type[0] + ')' + ']'
        else:
            # Add the new alpha at the end
            proposed_atomtype = atom1type[0][:len(atom1type[0])-1] + '$(*' + bondset[0] + atom2type[0] + ')' + ']'
        proposed_typename = atom1type[1] + ' ' + bondset[1] + ' ' + atom2type[1] + ' '
        return proposed_atomtype, proposed_typename

    def AddBetaSubstituentAtom(self, atom1type, bondset, atom2type):
        """
        Takes specified atomtypes atom1type and atom2type (where atom1type is a type WITH an alpha substituent), 
        and a specified bond set, and introduces an alpha substituent involving atom2type (which can be a 
        decorated type such as output by AtomDecorator); returns a tuple of (smarts, name) from the result. 
        Example output for input of atom1type=("[#1$(*-[#6])]", "hydrogen singly bonded to carbon"), 
        bondset = (['-'], "single"), and atom2type = ("[#8]" is ("[#1$(*-[#6]-[#8])]", "hydrogen singly bonded 
        to (carbon singly bonded to oxygen)") or something along those lines. This will have to handle two cases 
        separately -- addition of the first beta substituent (where it is inserted between the characters (] and 
        is inserted without parens), and addition of subsequent beta substituents (where it is inserted after 
        the first set of [] in the substituent and is inserted enclosed in parens, i.e. (-[w]).) It should raise 
        an exception if a beta substituent is attempted to be added to an atom1type which does not have an alpha 
        substituent.
        """

        # counting '[' tells us how many atoms are in the mix
        count = atom1type[0].count('[')
        proposed_atomtype = ""
        number_brackets = 0
        # find closed alpha atom
        closeAlpha = atom1type[0].find(']')
        # This has two atoms (already has an alpha atom)
        if count == 2: 
            proposed_atomtype = atom1type[0][:closeAlpha+1]
            proposed_atomtype += bondset[0] + atom2type[0] + ')]'
            proposed_typename = atom1type[1] + bondset[1] + ' ' + atom2type[1]
            if self.verbose: print("ADD FIRST BETA SUB: proposed --- %s %s" % ( str(proposed_atomtype), str(proposed_typename)))
        elif count > 2:
            # Has an alpha atom with at least 1 beta atom
            proposed_atomtype = atom1type[0][:closeAlpha+1]
            proposed_atomtype += '(' + bondset[0] + atom2type[0] + ')'
            proposed_atomtype += atom1type[0][closeAlpha+1:]
            proposed_typename = atom1type[1] + ' (' + bondset[1] + ' ' + atom2type[1] + ')'
            if self.verbose: print("ADD MORE BETA SUB: proposed --- %s %s" % ( str(proposed_atomtype), str(proposed_typename)))
        else:
            # Has only 1 atom which means there isn't an alpha atom yet, add an alpha atom instead
            proposed_atomtype, proposed_typename = self.AddAlphaSubstituentAtom(atom1type, bondset, atom2type) 
        return proposed_atomtype, proposed_typename


    def sample_atomtypes(self): 
        """
        Perform one step of atom type sampling.

        """
        # Copy current atomtypes for proposal.
        proposed_atomtypes = copy.deepcopy(self.atomtypes)
        proposed_molecules = copy.deepcopy(self.molecules)
        proposed_parents = copy.deepcopy(self.parents)
        natomtypes = len(proposed_atomtypes)
        ndecorators = len(self.decorators)
        natombasetypes = len(self.atom_basetype)

        valid_proposal = True
        
        #creating = False
        
        if random.random() < 0.5:
            # Pick an atom type to destroy.
            atomtype_index = random.randint(0, natomtypes-1)
            (atomtype, typename) = proposed_atomtypes[atomtype_index]
            if self.verbose: print("Attempting to destroy atom type %s : %s..." % (atomtype, typename))
            # Reject deletion of (populated) base types as we want to retain 
            # generics even if empty
            if [atomtype, typename] in self.used_basetypes: 
                if self.verbose: print("Destruction rejected for atom type %s because this is a generic type which was initially populated." % atomtype )
                return False

            # Delete the atomtype.
            proposed_atomtypes.remove([atomtype, typename])

            # update proposed parent dictionary
            for parent, children in proposed_parents.items():
                if atomtype in [at for [at, tn] in children]:
                    children += proposed_parents[atomtype]
                    children.remove([atomtype, typename])

            del proposed_parents[atomtype]

            # Try to type all molecules.
            try:
                self.type_molecules(proposed_atomtypes, proposed_molecules)
            except AtomTyper.TypingException as e:
                # Reject since typing failed.
                if self.verbose: print("Typing failed; rejecting.")
                valid_proposal = False
        else:
            #creating = True
            if self.decorator_behavior == 'simple-decorators':
                # Pick an atomtype to subtype.
                atomtype_index = random.randint(0, natomtypes-1)
                # Pick a decorator to add.
                decorator_index = random.randint(0, ndecorators-1)
                # Create new atomtype to insert by appending decorator with 'and' operator.
                (atomtype, atomtype_typename) = self.atomtypes[atomtype_index]
                (decorator, decorator_typename) = self.decorators[decorator_index]
                result = re.match('\[(.+)\]', atomtype)
                proposed_atomtype = '[' + result.groups(1)[0] + '&' + decorator + ']'
                proposed_typename = atomtype_typename + ' ' + decorator_typename
                if self.verbose: print("Attempting to create new subtype: '%s' (%s) + '%s' (%s) -> '%s' (%s)" % (atomtype, atomtype_typename, decorator, decorator_typename, proposed_atomtype, proposed_typename))
            
                # Update proposed parent dictionary
                proposed_parents[atomtype].append([proposed_atomtype, proposed_typename])
                # Hack to make naming consistent with below
                atom1smarts, atom1typename = atomtype, atomtype_typename

            else:
                # combinatorial-decorators
                nbondset = len(self.bondset)
                # Pick an atomtype
                atom1type = self.PickAnAtom(self.unmatched_atomtypes)
                # Check if it is the element you want, if not, try to get only that atomtype element
                element = re.findall('\d+', atom1type[0])[0]
                if self.verbose: print("****** Got this element '%s'." % element)
                while element != self.initial_element:
                    atom1type = self.PickAnAtom(self.unmatched_atomtypes)
                    element = re.findall('\d+', atom1type[0])[0]
                    if self.verbose: print("****** Got this element '%s'." % element)
                atom1smarts, atom1typename = atom1type
                # Check if we need to add an alfa or beta substituent
                if self.HasAlpha(atom1type):
                    # Has alpha
                    bondset_index = random.randint(0, nbondset-1)
                    atom2type = self.PickAnAtom(self.used_basetypes)
                    if random.random() < 0.5 or atom1type[0][2] == '1': # Add Beta Substituent Atom randomly or when it is Hydrogen
                        proposed_atomtype, proposed_typename = self.AddBetaSubstituentAtom(atom1type, self.bondset[bondset_index], atom2type)
                    else: # Add another Alpha Substituent if it is not a Hydrogen
                        proposed_atomtype, proposed_typename = self.AddAlphaSubstituentAtom(atom1type, self.bondset[bondset_index], atom2type, first_alpha = False)
                    if self.verbose: print("Attempting to create new subtype: -> '%s' (%s)" % (proposed_atomtype, proposed_typename))
                else:
                    # Has no alpha
                    if random.random() < 0.5:
                        # Add a no-bond decorator
                        decorator_index = random.randint(0, ndecorators-1)
                        decorator = self.decorators[decorator_index]
                        proposed_atomtype, proposed_typename = self.AtomDecorator(atom1type, decorator)
                        if self.verbose: print("Attempting to create new subtype: '%s' (%s) + '%s' (%s) -> '%s' (%s)" % (atom1type[0], atom1type[1], decorator[0], decorator[1], proposed_atomtype, proposed_typename))
                    else:
                        bondset_index = random.randint(0, nbondset-1)
                        atom2type = self.PickAnAtom(self.used_basetypes)
                        proposed_atomtype, proposed_typename = self.AddAlphaSubstituentAtom(atom1type, self.bondset[bondset_index], atom2type, first_alpha = True)
                        if self.verbose: print("Attempting to create new subtype: '%s' (%s) -> '%s' (%s)" % (atom1type[0], atom1type[1], proposed_atomtype, proposed_typename))

                # Update proposed parent dictionary
                proposed_parents[atom1type[0]].append([proposed_atomtype, proposed_typename])

            proposed_parents[proposed_atomtype] = []

            # Check that we haven't already determined this atom type isn't matched in the dataset.
            if proposed_atomtype in self.atomtypes_with_no_matches:
                if self.verbose: print("Atom type '%s' (%s) unused in dataset; rejecting." % (proposed_atomtype, proposed_typename))
                return False

            # Check if proposed atomtype is already in set.
            existing_atomtypes = set()
            for (a, b) in self.atomtypes:
                existing_atomtypes.add(a)
            if proposed_atomtype in existing_atomtypes:
                if self.verbose: print("Atom type already exists; rejecting to avoid duplication.")
                valid_proposal = False
            
            # Check for valid proposal before proceeding.
            if not valid_proposal:
                return False

            # Insert atomtype immediately after.
            proposed_atomtypes.insert(natomtypes, [proposed_atomtype, proposed_typename]) # Insert in the end (hierarchy issue)
            # Try to type all molecules.
            try:
                # Type molecules.
                self.type_molecules(proposed_atomtypes, proposed_molecules)
                # Compute updated statistics.
                [proposed_atom_typecounts, proposed_molecule_typecounts] = self.compute_type_statistics(proposed_atomtypes, proposed_molecules)
                # Reject if new type is unused.
                if (proposed_atom_typecounts[proposed_typename] == 0):
                    # Reject because new type is unused in dataset.
                    if self.verbose: print("Atom type '%s' (%s) unused in dataset; rejecting." % (proposed_atomtype, proposed_typename))
                    valid_proposal = False
                    # Store this atomtype to speed up future rejections
                    self.atomtypes_with_no_matches.add(proposed_atomtype)
                # Reject if parent type is now unused, UNLESS it is a base type
                if (proposed_atom_typecounts[atom1typename] == 0) and (atom1smarts not in self.basetypes_smarts):
                    # Reject because new type is unused in dataset.
                    if self.verbose: print("Parent type '%s' (%s) now unused in dataset; rejecting." % (atom1smarts, atom1typename))
                    valid_proposal = False
            except AtomTyper.TypingException as e:
                print("Exception: %s" % str(e))
                # Reject since typing failed.
                if self.verbose: print("Typing failed for one or more molecules using proposed atomtypes; rejecting.")
                valid_proposal = False

        # Check for valid proposal
        if not valid_proposal:
            return False

        if self.verbose: print('Proposal is valid...')


        # IF is creating a new atom type
        #if creating:
        #
        #    # Get the Element we are creating
        #    element = re.findall('\d+', proposed_atomtype)[0]
        #    if self.verbose: print("****** Calculating only for element '%s'." % element)
        #    # Create proposed atomtypes by element
        #    proposed_atomtype_by_element = []
        #    for atomt in proposed_atomtypes:
        #        if re.findall('\d+', atomt[0])[0] == element:
        #            proposed_atomtype_by_element += [atomt]
        #    print proposed_atomtype_by_element
        #    print proposed_atomtypes
        #
        #    # Get only the element atomtypes form the Total Atoms
        #    total_atom_type_element = []
        #    for total_atomt in self.atomtypes:
        #        if re.findall('\d+', total_atomt[0])[0] == element:
        #            total_atom_type_element += [total_atomt]
        #    print total_atom_type_element
        #                #self.type_molecules(proposed_atomtype_by_element, self.molecules)
        #    if self.verbose: print("****** GRAPH for atual atomtypes that are related with element %s." % (element))
        #    [self.atom_type_matches, self.total_atom_type_matches] = self.best_match_reference_types(total_atom_type_element, self.molecules)
        #    if self.verbose: print("****** GRAPH for proposed atomtypes that are related with element %s" % (element))
        #    (proposed_atom_type_matches, proposed_total_atom_type_matches) = self.best_match_reference_types(proposed_atomtype_by_element, proposed_molecules)
        #else:
        #    (proposed_atom_type_matches, proposed_total_atom_type_matches) = self.best_match_reference_types(proposed_atomtypes, proposed_molecules)


        # Accept automatically if no reference molecules
        accept = False
        if self.reference_typed_molecules is None:
            accept = True
        else:
            # Compute effective temperature
            if self.temperature == 0.0:
                effective_temperature = 1
            else:
                effective_temperature = (self.total_atoms * self.temperature)

            # Compute likelihood for accept/reject
            (proposed_atom_type_matches, proposed_total_atom_type_matches) = self.best_match_reference_types(proposed_atomtypes, proposed_molecules)
            log_P_accept = (proposed_total_atom_type_matches - self.total_atom_type_matches) / effective_temperature
            print('Proposal score: %d >> %d : log_P_accept = %.5e' % (self.total_atom_type_matches, proposed_total_atom_type_matches, log_P_accept))
            if (log_P_accept > 0.0) or (numpy.random.uniform() < numpy.exp(log_P_accept)):
                accept = True

        # Accept or reject
        if accept:
            self.atomtypes = proposed_atomtypes
            self.molecules = proposed_molecules
            self.parents = proposed_parents
            (proposed_atom_type_matches, proposed_total_atom_type_matches) = self.best_match_reference_types(proposed_atomtypes, proposed_molecules)
            self.atom_type_matches = proposed_atom_type_matches
            self.total_atom_type_matches = proposed_total_atom_type_matches
            return True
        else:
            [self.atom_type_matches, self.total_atom_type_matches] = self.best_match_reference_types(self.atomtypes, self.molecules)
            return False

    def type_molecules(self, typelist, molecules):
        """
        Type all molecules with the specified typelist.

        """
        # Create an atom typer.
        atomtyper = AtomTyper(typelist, self.typetag, replacements=self.replacements)

        # Type molecules.
        for molecule in molecules:
            atomtyper.assignTypes(molecule, "elemental")

        return


    def calculate_number_element(self, typelist, molecules):
        """
        Compute statistics for numnber of molecules assigned each type.
                    
        ARGUMENTS
                    
        typelist
        molecules
                    
        RETURNS
        #
        atom_typecounts (dict) - counts of number of atoms containing each atomtype
        molecule_typecounds (dict) - counts of number of molecules containing each atom type
        
        """
        # Zero type counts by atom and molecule.
        atom_typecounts = dict()
        number_element = 0
        for [smarts, typename] in typelist:
            # Get statistics for only the element we are looking for
            element = re.findall('\d+', smarts)[0]
            if element == self.initial_element:
                atom_typecounts[typename] = 0
        # Count number of atoms with each type.
        for molecule in molecules:
            types_in_this_molecule = set()
            for atom in molecule.GetAtoms():
                atomtype = atom.GetStringData(self.typetag)
                # Check if the atom key is in the atom_typecounts dictionary
                if (atomtype) in atom_typecounts.keys():
                    types_in_this_molecule.add(atomtype)
                    number_element += 1

        return number_element

    
    def compute_type_statistics(self, typelist, molecules, only_element=False):
        """
        Compute statistics for numnber of molecules assigned each type.

        ARGUMENTS

        typelist
        molecules

        RETURNS
#
        atom_typecounts (dict) - counts of number of atoms containing each atomtype
        molecule_typecounds (dict) - counts of number of molecules containing each atom type

        """
        # Zero type counts by atom and molecule.
        atom_typecounts = dict()
        molecule_typecounts = dict()
        for [smarts, typename] in typelist:
            # Get statistics for only the element we are looking for
            element = re.findall('\d+', smarts)[0]
            if only_element:
                if element == self.initial_element:
                    atom_typecounts[typename] = 0
                    molecule_typecounts[typename] = 0
            else:
                atom_typecounts[typename] = 0
                molecule_typecounts[typename] = 0

        # Count number of atoms with each type.
        for molecule in molecules:
            types_in_this_molecule = set()
            for atom in molecule.GetAtoms():
                atomtype = atom.GetStringData(self.typetag)
                # Check if the atom key is in the atom_typecounts dictionary
                if (atomtype) in atom_typecounts.keys():
                    types_in_this_molecule.add(atomtype)
                    atom_typecounts[atomtype] += 1
            for atomtype in types_in_this_molecule:
                molecule_typecounts[atomtype] += 1

        return (atom_typecounts, molecule_typecounts)

    def show_type_statistics(self, typelist, atom_typecounts, molecule_typecounts, atomtype_matches=None):
        """
        Print atom type statistics.

        """
        index = 1
        natoms = 0

        if atomtype_matches is not None:
            reference_type_info = dict()
            for (typename, reference_atomtype, count) in atomtype_matches:
                reference_type_info[typename] = (reference_atomtype, count)

        # Print header
        if atomtype_matches is not None:
            print "%5s   %10s %10s   %64s %32s %8s %46s" % ('INDEX', 'ATOMS', 'MOLECULES', 'TYPE NAME', 'SMARTS', 'REF TYPE', 'FRACTION OF REF TYPED MOLECULES MATCHED')
        else:
            print "%5s   %10s %10s   %64s %32s" % ('INDEX', 'ATOMS', 'MOLECULES', 'TYPE NAME', 'SMARTS')

        # Print counts
        for [smarts, typename] in typelist:
            element = re.findall('\d+', smarts)[0]
            if element == self.initial_element:
                if atomtype_matches is not None:
                    (reference_atomtype, reference_count) = reference_type_info[typename]
                    if reference_atomtype is not None:
                        reference_total = self.reference_atomtypes_atomcount[reference_atomtype]
                        reference_fraction = float(reference_count) / float(reference_total)
                        print "%5d : %10d %10d | %64s %32s %8s %16d / %16d (%7.3f%%)" % (index, atom_typecounts[typename], molecule_typecounts[typename], typename, smarts, reference_atomtype, reference_count, reference_total, reference_fraction*100)
                    else:
                        print "%5d : %10d %10d | %64s %32s" % (index, atom_typecounts[typename], molecule_typecounts[typename], typename, smarts)
                else:
                    print "%5d : %10d %10d | %64s %32s" % (index, atom_typecounts[typename], molecule_typecounts[typename], typename, smarts)

                natoms += atom_typecounts[typename]
                index += 1

        nmolecules = len(self.molecules)

        if atomtype_matches is not None:
            print "%5s : %10d %10d |  %64s %32s %8d / %8d match (%.3f %%)" % ('TOTAL', natoms, nmolecules, '', '', self.total_atom_type_matches, self.total_element, (float(self.total_atom_type_matches) / float(self.total_element)) * 100)
        else:
            print "%5s : %10d %10d" % ('TOTAL', natoms, nmolecules)

        return

    def save_type_statistics(self, typelist, atom_typecounts, molecule_typecounts, atomtype_matches=None):
        """
        Save "atom type" matches to be output to trajectory
        This isn't the most elegant solution, but it will make an output file we can read back in

        """
        if atomtype_matches is not None:
            reference_type_info = dict()
            for (typename, reference_atomtype, count) in atomtype_matches:
                reference_type_info[typename] = (reference_atomtype, count)

        index = 1
        output = []
        # Print counts
        # INDEX, SMARTS, PARENT INDEX, REF TYPE, MATCHES, MOLECULES, FRACTION, OUT of, PERCENTAGE
        for [smarts, typename] in typelist:
            if atomtype_matches is not None:
                (reference_atomtype, reference_count) = reference_type_info[typename]
                if reference_atomtype is not None:
                    reference_total = self.reference_atomtypes_atomcount[reference_atomtype]
                    reference_fraction = float(reference_count) / float(reference_total)
                    # Save output
                    output.append("%i,'%s',%i,%i,'%s',%i,%i,%i,%i" % (index, smarts, 0, 0, reference_atomtype, atom_typecounts[typename], molecule_typecounts[typename], reference_count, reference_total))
                else:
                    output.append("%i,'%s',%i,%i,'%s',%i,%i,%i,%i" % (index, smarts, 0, 0, 'NONE', atom_typecounts[typename], molecule_typecounts[typename], 0, 0))

            else:
                output.append("%i,'%s',%i,%i,'%s',%i,%i,%i,%i" % (index, smarts, 0, 0, 'NONE', atom_typecounts[typename], molecule_typecounts[typename], 0, 0))
            index += 1
        return output

    def get_unfinishedAtomList(self, atom_typecounts, molecule_typecounts, atomtype_matches = None):
        """
        This method prunes the set of current atomtypes so that if all branches 
        of a base type have been found it no longer tries extending any atom of that base type.  
        """
        # Reset unmatched atom types incase something was destroyed
        self.unmatched_atomtypes = copy.deepcopy(self.atomtypes)

        # If we don't have reference matches, unmatched_atomtypes should be all current atomtypes 
        if atomtype_matches is None:
            return
        else: # store counts for each atom type
            reference_counts = dict()
            for (typename, reference_atomtype, count) in atomtype_matches:
                if reference_atomtype is None:
                    reference_counts[typename] = 0
                else:
                    reference_counts[typename] = count

        # If all of a basetype and it's children match found atoms and reference remove from list
        for [base_smarts, base_typename] in self.used_basetypes:
            includeBase = True
            element = re.findall('\d+', base_smarts)[0]
            if element == self.initial_element:
                # If the number of atoms matches the references are the same for basetypes and their children
                # then we have found all reference types for that element and should stop searching that branch
                if atom_typecounts[base_typename] == reference_counts[base_typename]:
                    includeBase = False
                    for [child_smarts, child_name] in self.parents[base_smarts]:
                        # If any of the children's atom count and reference count don't agree then these should stay in the unmatched_atomtypes
                        if not atom_typecounts[child_name] == reference_counts[child_name]:
                            includeBase = True
                            break
            # Remove atomtypes from completed element branches
            #if not includeBase:
            #    print "INclude Base False"
            #    self.unmatched_atomtypes.remove([base_smarts, base_typename])
            #    for child in self.parents[base_smarts]:
            #        self.unmatched_atomtypes.remove(child)

        return

    def print_parent_tree(self, roots, start=''):
        """
        Recursively prints the parent tree. 

        Parameters
        ----------
        roots = list of smarts strings to print
        """
        for r in roots:
            print("%s%s" % (start, r))
            if r in self.parents.keys():
                new_roots = [smart for [smart, name] in self.parents[r]]
                self.print_parent_tree(new_roots, start+'\t')


    def run(self, niterations, trajFile=None, plotFile=None):
        """
        Run atomtype sampler for the specified number of iterations.

        Parameters
        ----------
        niterations : int
            The specified number of iterations
        trajFile : str, optional, default=None
            Output trajectory filename
        plotFile : str, optional, default=None
            Filename for output of plot of score versus time

        Returns
        ----------
        fraction_matched_atoms : float
            fraction of total atoms matched successfully at end of run

        """
        self.traj = []
        for iteration in range(niterations):
            if self.verbose:
                print("Iteration %d / %d" % (iteration, niterations))

            accepted = self.sample_atomtypes()
            [atom_typecounts, molecule_typecounts] = self.compute_type_statistics(self.atomtypes, self.molecules)
            self.get_unfinishedAtomList(atom_typecounts, molecule_typecounts, atomtype_matches = self.atom_type_matches)

            if trajFile is not None:
                # Get data as list of csv strings
                lines = self.save_type_statistics(self.atomtypes, atom_typecounts, molecule_typecounts, atomtype_matches=self.atom_type_matches)
                # Add lines to trajectory with iteration number:
                for l in lines:
                    self.traj.append('%i,%s \n' % (iteration, l))

            if self.verbose:
                if accepted:
                    print('Accepted.')
                else:
                    print('Rejected.')

                # Compute atomtype statistics on molecules.
                self.show_type_statistics(self.atomtypes, atom_typecounts, molecule_typecounts, atomtype_matches=self.atom_type_matches)
                print('')

                # Print parent tree as it is now.
                roots = self.parents.keys()
                # Remove keys from roots if they are children
                for parent, children in self.parents.items():
                    child_smarts = [smarts for [smarts, name] in children]
                    for child in child_smarts:
                        if child in roots:
                            roots.remove(child)

                print("Atom type hierarchy:")
                self.print_parent_tree(roots, '\t')

        if trajFile is not None:
            # make "trajectory" file
            if os.path.isfile(trajFile):
                print("trajectory file already exists, it was overwritten")
            f = open(trajFile, 'w')
            start = ['Iteration,Index,Smarts,ParNum,ParentParNum,RefType,Matches,Molecules,FractionMatched,Denominator\n']
            f.writelines(start + self.traj)
            f.close()
 
            # Get/print some stats on trajectory
            # Load timeseries
            timeseries = load_trajectory( trajFile )
            time_fractions = scores_vs_time( timeseries )
            print("Maximum score achieved: %.2f" % max(time_fractions['all']))

        # If desired, make plot
        if plotFile:
            import pylab as pl
            if not trajFile:
                raise Exception("Cannot construct plot of trajectory without a trajectory file.")
            # Load timeseries
            timeseries = load_trajectory( trajFile )
            time_fractions = scores_vs_time( timeseries )

            # Plot overall score
            pl.plot( time_fractions['all'], 'k-', linewidth=2.0)

            # Grab reference types other than 'all'
            plot_others = False
            if plot_others:
                reftypes = time_fractions.keys()
                reftypes.remove('all')

                # Plot scores for individual types
                for reftype in reftypes:
                    pl.plot( time_fractions[reftype] )
            
            # Axis labels and legend
            pl.xlabel('Iteration')
            pl.ylabel('Fraction of reference type found')
            if plot_others:
                pl.legend(['all']+reftypes, loc="lower right")
            pl.ylim(-0.1, 1.1)

            # Save
            pl.savefig( plotFile )


        #Compute final type stats
        [atom_typecounts, molecule_typecounts] = self.compute_type_statistics(self.atomtypes, self.molecules)
        fraction_matched_atoms = self.show_type_matches(self.atom_type_matches)

        # If verbose print parent tree:
        if self.verbose: 
            roots = self.parents.keys()
            # Remove keys from roots if they are children
            for parent, children in self.parents.items():
                child_smarts = [smarts for [smarts, name] in children]
                for child in child_smarts:
                    if child in roots:
                        roots.remove(child)

            print("Atom type hierarchy:")
            self.print_parent_tree(roots, '\t')
        return fraction_matched_atoms
