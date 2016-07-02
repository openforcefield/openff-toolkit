#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
smarty.py

Example illustrating a scheme to create and destroy atom types automatically using SMARTS.

AUTHORS

John Chodera <jchodera@berkeley.edu>, University of California, Berkeley

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

#=============================================================================================
# ATOMTYPE SAMPLER
#=============================================================================================

class AtomTypeSampler(object):
    """
    Atom type sampler.

    """
    def __init__(self, molecules, basetypes_filename, decorators_filename, replacements_filename=None, reference_typed_molecules=None, temperature=0.1, verbose=False):
        """
        Initialize an atom type sampler.

        ARGUMENTS

        molecules : list of molecules for typing
            List of molecules for typing
        basetypes_filename : str
            File defining base atom types (which cannot be destroyed)
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

        self.verbose = verbose

        # Define internal typing tag.
        self.typetag = 'atomtype'

        # Read atomtypes and decorators.
        self.atomtypes = AtomTyper.read_typelist(basetypes_filename)
        for idx, [smarts, typename] in enumerate(self.atomtypes):
            self.atomtypes[idx] = [smarts, 'c_'+typename]
        self.decorators = AtomTyper.read_typelist(decorators_filename)
        self.replacements = AtomTyper.read_typelist(replacements_filename)

        # Store a copy of the basetypes, as these (and only these) are allowed 
        # to end up with zero occupancy
        self.basetypes = copy.deepcopy(self.atomtypes)
        # Store smarts for basetypes
        self.basetypes_smarts = [ smarts for (smarts, name) in self.basetypes ]

        # Store a deep copy of the molecules since they will be annotated
        self.molecules = copy.deepcopy(molecules)

        # Type all molecules with current typelist to ensure that basetypes are sufficient.
        self.type_molecules(self.atomtypes, self.molecules)

        # Compute atomtype statistics on molecules.
        [atom_typecounts, molecule_typecounts] = self.compute_type_statistics(self.atomtypes, self.molecules)
        if self.verbose: self.show_type_statistics(self.atomtypes, atom_typecounts, molecule_typecounts)

        # Compute total atoms
        self.total_atoms = 0.0
        for molecule in self.molecules:
            for atom in molecule.GetAtoms():
                self.total_atoms += 1.0

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

        """
        print('Atom type matches:')
        total_atom_type_matches = 0
        for (current_atomtype, reference_atomtype, counts) in atom_type_matches:
            if reference_atomtype is not None:
                print('%-64s matches %8s : %8d atoms matched' % (current_atomtype, reference_atomtype, counts))
                total_atom_type_matches += counts
            else:
                print('%-64s         no match' % (current_atomtype))

        fraction_matched_atoms = float(total_atom_type_matches) / float(self.total_atoms)
        print('%d / %d total atoms match (%.3f %%)' % (total_atom_type_matches, self.total_atoms, fraction_matched_atoms * 100))


    def sample_atomtypes(self):
        """
        Perform one step of atom type sampling.

        """
        # Copy current atomtypes for proposal.
        proposed_atomtypes = copy.deepcopy(self.atomtypes)
        proposed_molecules = copy.deepcopy(self.molecules)
        natomtypes = len(proposed_atomtypes)
        ndecorators = len(self.decorators)

        valid_proposal = True

        if random.random() < 0.5:
            # Pick an atom type to destroy.
            atomtype_index = random.randint(0, natomtypes-1)
            (atomtype, typename) = proposed_atomtypes[atomtype_index]
            if self.verbose: print("Attempting to destroy atom type %s : %s..." % (atomtype, typename))
            # Delete the atomtype.
            proposed_atomtypes.remove([atomtype, typename])
            # Try to type all molecules.
            try:
                self.type_molecules(proposed_atomtypes, proposed_molecules)
            except AtomTyper.TypingException as e:
                #print e
                # Reject since typing failed.
                if self.verbose: print("Typing failed; rejecting.")
                valid_proposal = False
        else:
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
            proposed_atomtypes.insert(atomtype_index+1, [proposed_atomtype, proposed_typename])
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
                if (proposed_atom_typecounts[atomtype_typename] == 0) and (atomtype not in self.basetypes_smarts):
                    # Reject because new type is unused in dataset.
                    if self.verbose: print("Parent type '%s' (%s) now unused in dataset; rejecting." % (atomtype, atomtype_typename))
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
            self.atom_type_matches = proposed_atom_type_matches
            self.total_atom_type_matches = proposed_total_atom_type_matches
            return True
        else:
            return False

    def type_molecules(self, typelist, molecules):
        """
        Type all molecules with the specified typelist.

        """
        # Create an atom typer.
        atomtyper = AtomTyper(typelist, self.typetag, replacements=self.replacements)

        # Type molecules.
        for molecule in molecules:
            atomtyper.assignTypes(molecule)

        return

    def compute_type_statistics(self, typelist, molecules):
        """
        Compute statistics for numnber of molecules assigned each type.

        ARGUMENTS

        typelist
        molecules

        RETURNS

        atom_typecounts (dict) - counts of number of atoms containing each atomtype
        molecule_typecounds (dict) - counts of number of molecules containing each atom type

        """
        # Zero type counts by atom and molecule.
        atom_typecounts = dict()
        molecule_typecounts = dict()
        for [smarts, typename] in typelist:
            atom_typecounts[typename] = 0
            molecule_typecounts[typename] = 0

        # Count number of atoms with each type.
        for molecule in molecules:
            types_in_this_molecule = set()
            for atom in molecule.GetAtoms():
                atomtype = atom.GetStringData(self.typetag)
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
            print "%5s : %10d %10d |  %64s %32s %8d / %8d match (%.3f %%)" % ('TOTAL', natoms, nmolecules, '', '', self.total_atom_type_matches, self.total_atoms, (float(self.total_atom_type_matches) / float(self.total_atoms)) * 100)
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

    def run(self, niterations, trajFile):
        """
        Run atomtype sampler for the specified number of iterations.

        Parameters
        ----------
        niterations : int
            The specified number of iterations

        """
        self.traj = []
        for iteration in range(niterations):
            if self.verbose:
                print("Iteration %d / %d" % (iteration, niterations))

            accepted = self.sample_atomtypes()

            if self.verbose:
                if accepted:
                    print('Accepted.')
                else:
                    print('Rejected.')

                # Compute atomtype statistics on molecules.
                [atom_typecounts, molecule_typecounts] = self.compute_type_statistics(self.atomtypes, self.molecules)
                self.show_type_statistics(self.atomtypes, atom_typecounts, molecule_typecounts, atomtype_matches=self.atom_type_matches)

                # Get data as list of csv strings
                lines = self.save_type_statistics(self.atomtypes, atom_typecounts, molecule_typecounts, atomtype_matches=self.atom_type_matches)
                # Add lines to trajectory with iteration number:
                for l in lines:
                    self.traj.append('%i,%s \n' % (iteration, l))
                print('')
        
        # make "trajectory" file
        if os.path.isfile(trajFile):
            print "trajectory file already exists, it was overwritten"
        f = open(trajFile, 'w')
        start = ['Iteration,Index,Smarts,ParNum,ParentParNum,RefType,Matches,Molecules,FractionMatched,Denominator\n']
        f.writelines(start + self.traj)
        f.close()

