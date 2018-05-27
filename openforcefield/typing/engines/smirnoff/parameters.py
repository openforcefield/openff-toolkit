#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Parameter handlers for the SMIRNOFF force field engine

This file contains standard parameter handlers for the SMIRNOFF force field engine.
These classes implement the object model for self-contained parameter assignment.
New pluggable handlers can be created by creating subclasses of :class:`ParameterHandler`.

.. codeauthor:: John D. Chodera <john.chodera@choderalab.org>
.. codeauthor:: David L. Mobley <dmobley@mobleylab.org>
.. codeauthor:: Peter K. Eastman <peastman@stanford.edu>

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import os
import re
import sys
import math
import copy
import time
import logger
import itertools

from collections import OrderedDict

import numpy as np

import lxml.etree as etree

from simtk import openmm, unit
from simtk.openmm.app import element as elem

from openforcefield.utils import get_data_filename, all_subclasses
from openforcefield.topology import Topology, ValenceDict, ImproperDict
from openforcefield.topology import DEFAULT_AROMATICITY_MODEL
from openforcefield.typing.chemistry import ChemicalEnvironment, SMIRKSParsingError

#=============================================================================================
# CONFIGURE LOGGER
#=============================================================================================

logger = logging.getLogger(__name__)

#=============================================================================================
# PARAMETER HANDLERS
#
# The following classes are Handlers that know how to create Force subclasses and add them to a System that is being
# created.  Each Handler class must define three methods: 1) a static method that takes an etree Element and a ForceField,
# and returns the corresponding Handler object; 2) a createForce() method that constructs the Force object and adds it
# to the System; and 3) a labelForce() method that provides access to which
# terms are applied to which atoms in specified oemols.
# The static method should be added to the parsers map.
#=============================================================================================

class ParameterList(list):
    """Parameter list that also supports accessing items by SMARTS string.
    """
    # TODO: Make this faster by caching SMARTS -> index lookup?

    def __getitem__(self, item):
        """Retrieve item by index or SMIRKS
        """
        if type(item) == str:
            # Try to retrieve by SMARTS
            for result in self:
                if result.smirks == item:
                    return result
        # Try traditional access
        result = list.__getitem__(self, item)
        try:
            return ParameterList(result)
        except TypeError:
            return result
    def __contains__(self, item):
        """Check to see if either Parameter or SMIRKS is contained in parameter list.
        """
        if type(item) == str:
            # Special case for SMIRKS strings
            if item in [result.smirks for result in self]:
                return True
        # Fall back to traditional access
        return list.__contains__(self, item)

# TODO: Rename to better reflect role as parameter base class?
class ParameterType(object):
    """
    Base class for SMIRNOFF parameter types.

    """

    _VALENCE_TYPE = None # ChemicalEnvironment valence type for checking SMIRKS is conformant

    def __init__(self, smirks):
        """
        Create a ParameterType

        Parameters
        ----------
        smirks : str
            The SMIRKS match for the provided parameter type.

        """
        self.smirks = smirks

    @property
    def smirks(self):
        return self.__smirks

    @smirks.setter
    def set_smirks(self, smirks):
        # Validate the SMIRKS string to ensure it matches the expected parameter type,
        # raising an exception if it is invalid or doesn't tag a valid set of atoms
        ChemicalEnvironment.validate(smirks, ensure_valence_type=self.__valence_type)
        self.__smirks = smirks

    # TODO: Can we automatically check unit compatibilities for other parameters we create?
    # For example, if we have a parameter with units energy/distance**2, can we check to make
    # sure the dimensionality is preserved when the parameter is modified?

class ParameterHandler(object):
    """Virtual base class for parameter handlers.

    Parameter handlers are configured with some global parameters for a given section, and

    .. warning

       Parameter handler objects can only belong to a single :class:`ForceField` object.
       If you need to create a copy to attach to a different :class:`ForceField` object, use ``create_copy()``.

    """

    # TODO: Keep track of preferred units for parameter handlers.

    # TODO: Remove these?
    _TAGNAME = None # str of section type handled by this ParameterHandler (XML element name for SMIRNOFF XML representation)
    _VALENCE_TYPE = None # ChemicalEnvironment valence type string expected by SMARTS string for this Handler
    _INFOTYPE = None # container class with type information that will be stored in self._types
    _OPENMMTYPE = None # OpenMM Force class (or None if no equivalent)
    _DEPENDENCIES = None # list of ParameterHandler classes that must precede this, or None
    _DEFAULTS = {} # dict of attributes and their default values at tag-level
    _KWARGS = [] # list of keyword arguments accepted by the force Handler on initialization
    _SMIRNOFF_VERSION_INTRODUCED = 0.0 # the earliest version of SMIRNOFF spec that supports this ParameterHandler
    _SMIRNOFF_VERSION_DEPRECATED = None # if deprecated, the first SMIRNOFF version number it is no longer used
    _REQUIRE_UNITS = None # list of parameters that require units to be defined

    # TODO: Do we need to store the parent forcefield object?
    def __init__(self, forcefield):
        self._forcefield = forcefield # the ForceField object that this ParameterHandler is registered with
        self.parameters = list() # list of ParameterType objects # TODO: Change to method accessor so we can access as list or dict

    # TODO: Do we need to return these, or can we handle this internally
    @property
    def known_kwargs(self):
        """List of kwargs that can be parsed by the function.
        """
        # TODO: Should we use introspection to inspect the function signature instead?
        return set(self._KWARGS)

    def check_compatibility(self):
        """Check to make sure the global parameters for this parameter type are compatible with how this class has currently been configured.

        Raises a ValueError if the parameters are incompatible.
        """
        pass

    # TODO: Can we ensure SMIRKS and other parameters remain valid after manipulation?
    def add_parameter(self):
        """Add a parameter to the forcefield, ensuring all parameters are valid.
        """
        pass

    def get_matches(self, entity):
        """Retrieve all force terms for a chemical entity, which could be a Molecule, group of Molecules, or Topology.

        Parameters
        ----------
        entity : openforcefield.topology.ChemicalEntity
            Chemical entity for which constraints are to be enumerated

        Returns
        ---------
        matches : ValenceDict
            matches[atoms] is the ParameterType object corresponding to the tuple of Atom objects ``Atoms``

        """
        logger.info(self.__class__.__name__) # TODO: Overhaul logging
        matches = ValenceDict()
        for force_type in self._types:
            matches_for_this_type = { atoms : force_type for atoms in topology.chemical_environment_matches(force_type.smirks) }
            matches.update(matches_for_this_type)
            logger.info('{:64} : {:8} matches'.format(force_type.smirks, len(matches_for_this_type)))

        logger.info('{} matches identified'.format(len(matches)))
        return matches

    def assign_parameters(self, topology, system):
        """Assign parameters for the given Topology to the specified System object.

        Parameters
        ----------
        topology : openforcefield.topology.Topology
            The Topology for which parameters are to be assigned.
            Either a new Force will be created or parameters will be appended to an existing Force.
        system : simtk.openmm.System
            The OpenMM System object to add the Force (or append new parameters) to.
        """
        pass

    def postprocess_system(self, topology, system):
        """Allow the force to perform a a final post-processing pass on the System following parameter assignment, if needed.

        Parameters
        ----------
        topology : openforcefield.topology.Topology
            The Topology for which parameters are to be assigned.
            Either a new Force will be created or parameters will be appended to an existing Force.
        system : simtk.openmm.System
            The OpenMM System object to add the Force (or append new parameters) to.
        """
        pass


#=============================================================================================

class ConstraintHandler(ParameterHandler):
    """Handle SMIRNOFF ``<Constraints>`` tags

    ``ConstraintHandler`` must be applied before ``BondHandler`` and ``AngleHandler``,
    since those classes add constraints for which equilibrium geometries are needed from those tags.
    """
    class ConstraintType(ParameterType):
        """A SMIRNOFF constraint type"""
        def __init__(self, node, parent):
            super(ConstraintType, self).__init__(node, parent) # Base class handles ``smirks`` and ``id`` fields
            if 'distance' in node.attrib:
                self.distance = _extract_quantity_from_xml_element(node, parent, 'distance') # Constraint with specified distance will be added by ConstraintHandler
            else:
                self.distance = True # Constraint to equilibrium bond length will be added by HarmonicBondHandler

    _TAGNAME = 'Constraint'
    _VALENCE_TYPE = 'Bond' # ChemicalEnvironment valence type expected for SMARTS # TODO: Do we support more exotic types as well?
    _INFOTYPE = ConstraintType
    _OPENMMTYPE = None # don't create a corresponding OpenMM Force class
    _REQUIRED_UNITS = ['distance']

    def __init__(self, forcefield):
        super(ConstraintHandler, self).__init__(forcefield)

    def createForce(self, system, topology, **kwargs):
        constraints = self.getMatches(topology)
        for (atoms, constraint) in constraints.items():
            # Update constrained atom pairs in topology
            topology.add_constraint(*atoms, constraint.distance)
            # If a distance is specified (constraint.distance != True), add the constraint here.
            # Otherwise, the equilibrium bond length will be used to constrain the atoms in HarmonicBondHandler
            if constraint.distance is not True:
                system.addConstraint(*atoms, constraint.distance)

#=============================================================================================


class BondHandler(ParameterHandler):
    """Handle SMIRNOFF ``<BondForce>`` tags"""

    class BondType(ParameterType):
        """A SMIRNOFF Bond parameter type"""
        def __init__(self, node, parent):
            super(ConstraintType, self).__init__(node, parent) # Base class handles ``smirks`` and ``id`` fields

            # Determine if we are using fractional bond orders for this bond
            # First, check if this force uses fractional bond orders
            if 'fractional_bondorder_method' in parent.attrib:
                # If it does, see if this parameter line provides fractional bond order parameters
                if 'length_bondorder1' in node.attrib and 'k_bondorder1' in node.attrib:
                    # Store what interpolation scheme we're using
                    self.fractional_bondorder = parent.attrib['fractional_bondorder']
                    # Store bondorder1 and bondorder2 parameters
                    self.k = list()
                    self.length = list()
                    for ct in range(1,3):
                        self.length.append( _extract_quantity_from_xml_element(node, parent, 'length_bondorder%s' % ct, unit_name = 'length_unit') )
                        self.k.append( _extract_quantity_from_xml_element(node, parent, 'k_bondorder%s' % ct, unit_name = 'k_unit') )
                else:
                    self.fractional_bondorder = None
            else:
                self.fractional_bondorder = None

            # If no fractional bond orders, just get normal length and k
            if self.fractional_bondorder is None:
                self.length = _extract_quantity_from_xml_element(node, parent, 'length')
                self.k = _extract_quantity_from_xml_element(node, parent, 'k')

    _TAGNAME = 'BondForce' # SMIRNOFF tag name to process
    _VALENCE_TYPE = 'Bond' # ChemicalEnvironment valence type expected for SMARTS
    _INFOTYPE = BondType # class to hold force type info
    _OPENMMTYPE = openmm.HarmonicBondForce # OpenMM force class to create
    _DEPENDENCIES = [ConstraintHandler] # ConstraintHandler must be executed first

    def __init__(self, forcefield):
        super(HarmonicBondHandler, self).__init__(forcefield)

    def createForce(self, system, topology, **kwargs):
        # Create or retrieve existing OpenMM Force object
        force = super(BondHandler, self).createForce(system, topology, **kwargs)

        # Add all bonds to the system.
        bonds = self.getMatches(topology)
        skipped_constrained_bonds = 0 # keep track of how many bonds were constrained (and hence skipped)
        for (atoms, bond) in bonds.items():
            # Get corresponding particle indices in Topology
            particle_indices = tuple([ atom.particle_index for atom in atoms ])

            # Ensure atoms are actually bonded correct pattern in Topology
            topology.assert_bonded(atoms[0], atoms[1])

            # Compute equilibrium bond length and spring constant.
            if bond.fractional_bondorder is None:
                [k, length] = [bond.k, bond.length]
            else:
                # Interpolate using fractional bond orders
                # TODO: Do we really want to allow per-bond specification of interpolation schemes?
                order = topology.get_fractional_bond_order(*atoms)
                if bond.fractional_bondorder_interpolation == 'interpolate-linear':
                    k = bond.k[0] + (bond.k[1]-bond.k[0])*(order-1.)
                    length = bond.length[0] + (bond.length[1]-bond.length[0])*(order-1.)
                else:
                    raise Exception("Partial bondorder treatment {} is not implemented.".format(bond.fractional_bondorder))

            # Handle constraints.
            if topology.atom_pair_is_constrained(*atoms):
                # Atom pair is constrained; we don't need to add a bond term.
                skipped_constrained_bonds += 1
                # Check if we need to add the constraint here to the equilibrium bond length.
                if topology.atom_pair_is_constrained(*atoms) is True:
                    # Mark that we have now assigned a specific constraint distance to this constraint.
                    topology.add_constraint(*atoms, length)
                    # Add the constraint to the System.
                    system.addConstraint(*particle_indices, length)
                continue

            # Add harmonic bond to HarmonicBondForce
            force.addBond(*particle_indices, length, k)

        logger.info('{} bonds added ({} skipped due to constraints)'.format(len(bonds) - skipped_constrained_bonds, skipped_constrained_bonds))

        # Check that no topological bonds are missing force parameters
        _check_for_missing_valence_terms('BondForce', topology, bonds.keys(), topology.bonds)

#=============================================================================================


class AngleHandler(ParameterHandler):
    """Handle SMIRNOFF ``<AngleForce>`` tags"""

    class AngleType(ParameterType):
        """A SMIRNOFF angle type."""
        def __init__(self, node, parent):
            super(AngleType, self).__init__(node, parent) # base class handles ``smirks`` and ``id`` fields
            self.angle = _extract_quantity_from_xml_element(node, parent, 'angle')
            self.k = _extract_quantity_from_xml_element(node, parent, 'k')
            if 'fractional_bondorder' in parent.attrib:
                self.fractional_bondorder = parent.attrib['fractional_bondorder']
            else:
                self.fractional_bondorder = None

    _TAGNAME = 'AngleForce' # SMIRNOFF tag name to process
    _VALENCE_TYPE = 'Angle' # ChemicalEnvironment valence type expected for SMARTS
    _INFOTYPE = AngleType # class to hold force type info
    _OPENMMTYPE = openmm.HarmonicAngleForce # OpenMM force class to create

    def __init__(self, forcefield):
        super(AngleHandler, self).__init__(forcefield)

    def createForce(self, system, topology, **kwargs):
        force = super(AngleHandler, self).createForce(system, topology, **kwargs)

        # Add all angles to the system.
        angles = self.getMatches(topology)
        skipped_constrained_angles = 0 # keep track of how many angles were constrained (and hence skipped)
        for (atoms, angle) in angles.items():
            # Get corresponding particle indices in Topology
            particle_indices = tuple([ atom.particle_index for atom in atoms ])

            # Ensure atoms are actually bonded correct pattern in Topology
            for (i,j) in [ (0,1), (1,2) ]:
                topology.assert_bonded(atoms[i], atoms[j])

            if topology.is_constrained(atoms[0], atoms[1]) and topology.is_constrained(atoms[1], atoms[2]) and topology.is_constrained(atoms[0], atoms[2]):
                # Angle is constrained; we don't need to add an angle term.
                skipped_constrained_angles += 1
                continue

            force.addAngle(*particle_indices, angle.angle, angle.k)

        logger.info('{} angles added ({} skipped due to constraints)'.format(len(angles) - skipped_constrained_angles, skipped_constrained_angles))

        # Check that no topological angles are missing force parameters
        _check_for_missing_valence_terms('AngleForce', topology, angles.keys(), topology.angles())

#=============================================================================================


class ProperTorsionHandler(ParameterHandler):
    """Handle SMIRNOFF ``<ProperTorsionForce>`` tags"""

    class ProperTorsionType(ParameterType):
        """A SMIRNOFF torsion type for proper torsions."""
        def __init__(self, node, parent):
            super(ProperTorsionType, self).__init__(node, parent) # base class handles ``smirks`` and ``id`` fields
            self.periodicity = list()
            self.phase = list()
            self.k = list()

            # Check that the SMIRKS pattern matches the type it's supposed to
            try:
                chemenv = ChemicalEnvironment(self.smirks)
                thistype = chemenv.getType()
                if thistype != 'ProperTorsion':
                    raise Exception("Error: SMIRKS pattern %s (parameter %s) does not specify a %s torsion, but it is supposed to." % (self.smirks, self.pid, 'Proper'))
            except SMIRKSParsingError:
                print("Warning: Could not confirm whether smirks pattern %s is a valid %s torsion." % (self.smirks, self.torsiontype))


            # TODO: Fractional bond orders should be processed on the per-force basis instead of per-bond basis
            if 'fractional_bondorder_method' in parent.attrib:
                self.fractional_bondorder = parent.attrib['fractional_bondorder']
            else:
                self.fractional_bondorder = None

            # Store parameters.
            index = 1
            while 'phase%d'%index in node.attrib:
                self.periodicity.append( int(_extract_quantity_from_xml_element(node, parent, 'periodicity%d' % index)) )
                self.phase.append( _extract_quantity_from_xml_element(node, parent, 'phase%d' % index, unit_name='phase_unit') )
                self.k.append( _extract_quantity_from_xml_element(node, parent, 'k%d' % index, unit_name='k_unit') )
                # Optionally handle 'idivf', which divides the periodicity by the specified value
                if ('idivf%d' % index) in node.attrib:
                    idivf = _extract_quantity_from_xml_element(node, parent, 'idivf%d' % index)
                    self.k[-1] /= float(idivf)
                index += 1

            # Check for errors, i.e. 'phase' instead of 'phase1'
            # TODO: Can we raise a more useful error if there is no ``id``?
            if len(self.phase)==0:
               raise Exception("Error: Torsion with id %s has no parseable phase entries." % self.pid)

    _TAGNAME = 'ProperTorsionForce' # SMIRNOFF tag name to process
    _VALENCE_TYPE = 'ProperTorsion' # ChemicalEnvironment valence type expected for SMARTS
    _INFOTYPE = ProperTorsionType # info type to store
    _OPENMMTYPE = openmm.PeriodicTorsionForce # OpenMM force class to create

    def __init__(self, forcefield):
        super(ProperTorsionHandler, self).__init__(forcefield)

    def createForce(self, system, topology, **kwargs):
        force = super(ProperTorsionHandler, self).createForce(system, topology, **kwargs)

        # Add all proper torsions to the system.
        torsions = self.getMatches(topology)
        for (atoms, torsion) in torsions.items():
            # Ensure atoms are actually bonded correct pattern in Topology
            for (i,j) in [ (0,1), (1,2), (2,3) ]:
                topology.assert_bonded(atoms[i], atoms[j])

            for (periodicity, phase, k) in zip(torsion.periodicity, torsion.phase, torsion.k):
                force.addTorsion(atom_indices[0], atom_indices[1], atom_indices[2], atom_indices[3], periodicity, phase, k)

        logger.info('{} torsions added'.format(len(torsions)))

        # Check that no topological torsions are missing force parameters
        _check_for_missing_valence_terms('ProperTorsionForce', topology, torsions.keys(), topology.torsions())


class ImproperTorsionHandler(ParameterHandler):
    """Handle SMIRNOFF ``<ImproperTorsionForce>`` tags"""

    class ImproperTorsionType(ParameterType):
        """A SMIRNOFF torsion type for improper torsions."""
        def __init__(self, node, parent):
            super(ImproperTorsionType, self).__init__(node, parent) # base class handles ``smirks`` and ``id`` fields
            self.periodicity = list()
            self.phase = list()
            self.k = list()

            # Check that the SMIRKS pattern matches the type it's supposed to
            try:
                chemenv = ChemicalEnvironment(self.smirks)
                thistype = chemenv.getType()
                if thistype != 'ImproperTorsion':
                    raise Exception("Error: SMIRKS pattern %s (parameter %s) does not specify a %s torsion, but it is supposed to." % (self.smirks, self.pid, 'Improper'))
            except SMIRKSParsingError:
                print("Warning: Could not confirm whether smirks pattern %s is a valid %s torsion." % (self.smirks, self.torsiontype))

            if 'fractional_bondorder' in parent.attrib:
                self.fractional_bondorder = parent.attrib['fractional_bondorder']
            else:
                self.fractional_bondorder = None

            # Store parameters.
            index = 1
            while 'phase%d'%index in node.attrib:
                self.periodicity.append( int(_extract_quantity_from_xml_element(node, parent, 'periodicity%d' % index)) )
                self.phase.append( _extract_quantity_from_xml_element(node, parent, 'phase%d' % index, unit_name='phase_unit') )
                self.k.append( _extract_quantity_from_xml_element(node, parent, 'k%d' % index, unit_name='k_unit') )
                # Optionally handle 'idivf', which divides the periodicity by the specified value
                if ('idivf%d' % index) in node.attrib:
                    idivf = _extract_quantity_from_xml_element(node, parent, 'idivf%d' % index)
                    self.k[-1] /= float(idivf)
                index += 1
                # SMIRNOFF applies trefoil (three-fold, because of right-hand rule) impropers unlike AMBER
                # If it's an improper, divide by the factor of three internally
                if node.tag=='Improper':
                    self.k[-1] /= 3.
            # Check for errors, i.e. 'phase' instead of 'phase1'
            # TODO: What can we do if there is no ``id``?
            if len(self.phase)==0:
               raise Exception("Error: Torsion with id %s has no parseable phase entries." % self.pid)

    _TAGNAME = 'ImproperTorsionForce' # SMIRNOFF tag name to process
    _VALENCE_TYPE = 'ImproperTorsion' # ChemicalEnvironment valence type expected for SMARTS
    _INFOTYPE = ImproperTorsionType # info type to store
    _OPENMMTYPE = openmm.PeriodicTorsionForce # OpenMM force class to create

    def __init__(self, forcefield):
        super(ImproperTorsionHandler, self).__init__(forcefield)

    def createForce(self, system, topology, **kwargs):
        force = super(ImproperTorsionHandler, self).createForce(system, topology, **kwargs)

        # Add all improper torsions to the system
        torsions = self.getMatches(topology)
        for (atom_indices, improper) in impropers.items():
            # Ensure atoms are actually bonded correct pattern in Topology
            # For impropers, central atom is atom 1
            for (i,j) in [ (0,1), (1,2), (1,3) ]:
                topology.assert_bonded(atoms[i], atoms[j])

            # Impropers are applied in three paths around the trefoil having the same handedness
            for (periodicity, phase, k) in zip(improper.periodicity, improper.phase, improper.k):
                # Permute non-central atoms
                others = [ atom_indices[0], atom_indices[2], atom_indices[3] ]
                for p in [ (others[i], others[j], others[k]) for (i,j,k) in [(0,1,2), (1,2,0), (2,0,1)] ]
                    force.addTorsion(atom_indices[1], p[0], p[1], p[2], periodicity, phase, k)

        logger.info('{} impropers added, each applied in a six-fold trefoil' % (len(impropers)))

        # Check that no topological torsions are missing force parameters
        _check_for_missing_valence_terms('ImproperTorsionForce', topology, torsions.keys(), topology.impropers())

class vdWHandler(ParameterHandler):
    """Handle SMIRNOFF ``<vdW>`` tags"""

    _TAGNAME = 'vdW' # SMIRNOFF tag name to process
    _VALENCE_TYPE = 'Atom' # ChemicalEnvironment valence type expected for SMARTS
    _INFOTYPE = None # info type to store
    _OPENMMTYPE = openmm.NonbondedForce # OpenMM force class to create

    # TODO: Is this necessary
    SCALETOL = 1e-5

    # DEFAULTS
    _DEFAULTS = {
        'potential' : 'Lennard-Jones-12-6',
        'scale12' : 0.0,
        'scale13' : 0.0,
        'scale14' : 0.5,
        'scale15' : 1.0,
    }

    class vdWType(ParameterType):
        """A SMIRNOFF vdWForce type."""
        def __init__(self, node, parent):
            # NOTE: Currently we support radius definition via 'sigma' or 'rmin_half'.
            super(StericsType, self).__init__(node, parent) # base class handles ``smirks`` and ``id`` fields

            # Make sure we don't have BOTH rmin_half AND sigma
            try:
                a = _extract_quantity_from_xml_element(node, parent, 'sigma')
                a = _extract_quantity_from_xml_element(node, parent, 'rmin_half')
                raise Exception("Error: BOTH sigma and rmin_half cannot be specified simultaneously in the .offxml file.")
            except:
                pass

            # Handle Lennard-Jones sigma
            try:
                self.sigma = _extract_quantity_from_xml_element(node, parent, 'sigma')
            #Handle rmin_half, AMBER-style
            except:
                rmin_half = _extract_quantity_from_xml_element(node, parent, 'rmin_half', unit_name='sigma_unit')
                self.sigma = 2.*rmin_half/(2.**(1./6.))
            self.epsilon = _extract_quantity_from_xml_element(node, parent, 'epsilon')

    _TAGNAME = 'vdWForce' # SMIRNOFF tag name to process
    _INFOTYPE = vdWType # info type to store
    _OPENMMTYPE = openmm.NonbondedForce # OpenMM force class to create

    def __init__(self, forcefield):
        super(NonbondedParameterHandler, self).__init__(forcefield)


    # TODO: Handle the case where multiple <NonbondedForce> tags are found
    # if abs(Handler.coulomb14scale - float(element.attrib['coulomb14scale'])) > NonbondedHandler.SCALETOL or \
    #                 abs(Handler.lj14scale - float(element.attrib['lj14scale'])) > NonbondedHandler.SCALETOL:
    #            raise ValueError('Found multiple NonbondedForce tags with different 1-4 scales')
    #    for atom in element.findall('Atom'):
    #        Handler.registerAtom(atom, element)

    # TODO: nonbondedMethod and nonbondedCutoff should now be specified by StericsForce attributes
    def createForce(self, system, topology, nonbondedMethod=NoCutoff, nonbondedCutoff=0.9, **args):
        force = super(NonbondedParameterHandler, self).createForce(system, topology)

        methodMap = {NoCutoff:openmm.NonbondedForce.NoCutoff,
                     CutoffNonPeriodic:openmm.NonbondedForce.CutoffNonPeriodic,
                     CutoffPeriodic:openmm.NonbondedForce.CutoffPeriodic,
                     Ewald:openmm.NonbondedForce.Ewald,
                     PME:openmm.NonbondedForce.PME}
        if nonbondedMethod not in methodMap:
            raise ValueError('Illegal nonbonded method for NonbondedForce; method given was %s' % nonbondedMethod)
        force = openmm.NonbondedForce()
        force.setNonbondedMethod(methodMap[nonbondedMethod])
        force.setCutoffDistance(nonbondedCutoff)
        if 'ewaldErrorTolerance' in args:
            force.setEwaldErrorTolerance(args['ewaldErrorTolerance'])
        if 'useDispersionCorrection' in args:
            force.setUseDispersionCorrection(bool(args['useDispersionCorrection']))
        system.addForce(force)

        # Iterate over all defined Lennard-Jones types, allowing later matches to override earlier ones.
        atoms = self.getMatches(topology)

        # Create all particles.
        for particle in topology.particles:
            force.addParticle(0.0, 1.0, 0.0)

        # Set the particle Lennard-Jones terms.
        for (atoms, ljtype) in atoms.items():
            force.setParticleParameters(atoms[0].particle_index, 0.0, ljtype.sigma, ljtype.epsilon)

        # Check that no atoms are missing force parameters
        # QUESTION: Don't we want to allow atoms without force parameters? Or perhaps just *particles* without force parameters, but not atoms?
        _check_for_missing_valence_terms('NonbondedForce Lennard-Jones parameters', topology, atoms.keys(), topology.atoms)

        # Set the partial charges
        # TODO: We need to make sure we have already assigned partial charges to the Topology reference molecules
        for atom in topology.atoms:
            # Retrieve nonbonded parameters for reference atom (charge not set yet)
            _, sigma, epsilon = force.getParticleParameters(atom.particle_index)
            # Set the nonbonded force with the partial charge
            force.setParticleParameters(atom.particle_index, atom.charge, sigma, epsilon)

    def postprocessSystem(self, system, topology, **args):
        # Create exceptions based on bonds.
        # QUESTION: Will we want to do this for *all* cases, or would we ever want flexibility here?
        bond_particle_indices = [ (atom1.particle_index, atom2.particle_index) for (atom1, atom2) in topology.bonds ]
        for force in system.getForces():
            # TODO: Should we just store which `Force` object we are adding to and use that instead,
            # to prevent interference with other kinds of forces in the future?
            # TODO: Can we generalize this to allow for `CustomNonbondedForce` implementations too?
            if isinstance(force, openmm.NonbondedForce):
                nonbonded.createExceptionsFromBonds(bond_particle_indices, self.coulomb14scale, self.lj14scale)


class BondChargeCorrectionHandler(ParameterHandler):
    """Handle SMIRNOFF ``<BondChargeCorrection>`` tags"""

    class BondChargeCorrectionType(ParameterType):
        """A SMIRNOFF bond charge correction type."""
        def __init__(self, node, parent):
            super(BondChargeCorrectionHandler, self).__init__(node, parent) # base class handles ``smirks`` and ``id`` fields
            self.increment = _extract_quantity_from_xml_element(node, parent, 'increment')
            # If no units are specified, assume elementary charge
            if type(self.increment) == float:
                self.increment *= unit.elementary_charge

    _TAGNAME = 'BondChargeCorrection' # SMIRNOFF tag name to process
    _INFOTYPE = BondChargeCorrectionType # info type to store
    _OPENMMTYPE = openmm.NonbondedForce # OpenMM force class to create or utilize

    def __init__(self, forcefield):
        super(BondChargeCorrectionHandler, self).__init__(forcefield)

    #if element.attrib['method'] != existing[0]._initialChargeMethod:
    #raise Exception("Existing BondChargeCorrectionHandler uses initial charging method '%s' while new BondChargeCorrectionHandler requests '%s'" % (existing[0]._initialChargeMethod, element.attrib['method']))

    def createForce(self, system, topology, **args):
        # No forces are created by this Handler.
        pass

    # TODO: Move chargeModel and library residue charges to SMIRNOFF spec
    def postprocessSystem(self, system, topology, **kwargs):
        bonds = self.getMatches(topology)

        # Apply bond charge increments to all appropriate force groups
        # QUESTION: Should we instead apply this to the Topology in a preprocessing step, prior to spreading out charge onto virtual sites?
        for force in system.getForces():
            if force.__class__.__name__ in ['NonbondedForce']: # TODO: We need to apply this to all Force types that involve charges
                for (atoms, bond) in bonds.items():
                    # Get corresponding particle indices in Topology
                    particle_indices = tuple([ atom.particle_index for atom in atoms ])
                    # Retrieve parameters
                    [charge0, sigma0, epsilon0] = force.getParticleParameters(particle_indices[0])
                    [charge1, sigma1, epsilon1] = force.getParticleParameters(particle_indices[1])
                    # Apply bond charge increment
                    charge0 -= bond.increment
                    charge1 += bond.increment
                    # Update charges
                    force.setParticleParameters(particle_indices[0], charge0, sigma0, epsilon0)
                    force.setParticleParameters(particle_indices[1], charge1, sigma1, epsilon1)


class GBSAParameterHandler(ParameterHandler):
    """Handle SMIRNOFF ``<GBSAParameterHandler>`` tags"""
    # TODO: Differentiate between global and per-particle parameters for each model.

    # Global parameters for surface area (SA) component of model
    SA_expected_parameters = {
        'ACE' : ['surface_area_penalty', 'solvent_radius'],
        None : [],
    }

    # Per-particle parameters for generalized Born (GB) model
    GB_expected_parameters = {
        'HCT' : ['radius', 'scale'],
        'OBC1' : ['radius', 'scale'],
        'OBC2' : ['radius', 'scale'],
    }

    class GBSAType(ParameterType):
        """A SMIRNOFF GBSA type."""
        def __init__(self, node, parent):
            super(GBSAType, self).__init__(node, parent)

            # Store model parameters.
            gb_model = parent.attrib['gb_model']
            expected_parameters = GBSAParameterHandler.GB_expected_parameters[gb_model]
            provided_parameters = list()
            missing_parameters = list()
            for name in expected_parameters:
                if name in node.attrib:
                    provided_parameters.append(name)
                    value = _extract_quantity_from_xml_element(node, parent, name)
                    setattr(self, name, value)
                else:
                    missing_parameters.append(name)
            if len(missing_parameters) > 0:
                msg  = 'GBSAForce: missing per-atom parameters for tag %s' % str(node)
                msg += 'model "%s" requires specification of per-atom parameters %s\n' % (gb_model, str(expected_parameters))
                msg += 'provided parameters : %s\n' % str(provided_parameters)
                msg += 'missing parameters: %s' % str(missing_parameters)
                raise Exception(msg)

    def __init__(self, forcefield):
        super(GBSAParameterHandler, self).__init__(forcefield)

    # TODO: Fix this
    def parseElement(self):
        # Initialize GB model
        gb_model = element.attrib['gb_model']
        valid_GB_models = GBSAParameterHandler.GB_expected_parameters.keys()
        if not gb_model in valid_GB_models:
            raise Exception('Specified GBSAForce model "%s" not one of valid models: %s' % (gb_model, valid_GB_models))
        self.gb_model = gb_model

        # Initialize SA model
        sa_model = element.attrib['sa_model']
        valid_SA_models = GBSAParameterHandler.SA_expected_parameters.keys()
        if not sa_model in valid_SA_models:
            raise Exception('Specified GBSAForce SA_model "%s" not one of valid models: %s' % (sa_model, valid_SA_models))
        self.sa_model = sa_model

        # Store parameters for GB and SA models
        # TODO: Deep copy?
        self.parameters = element.attrib

    # TODO: Generalize this to allow forces to know when their OpenMM Force objects can be combined
    def checkCompatibility(self, Handler):
        """
        Check compatibility of this Handler with another Handlers.
        """
        Handler = existing[0]
        if (Handler.gb_model != self.gb_model):
            raise ValueError('Found multiple GBSAForce tags with different GB model specifications')
        if (Handler.sa_model != self.sa_model):
            raise ValueError('Found multiple GBSAForce tags with different SA model specifications')
        # TODO: Check other attributes (parameters of GB and SA models) automatically?

    def createForce(self, system, topology, **args):
        # TODO: Rework this
        from openforcefield.typing.engines.smirnoff import gbsaforces
        force_class = getattr(gbsaforces, self.gb_model)
        force = force_class(**self.parameters)
        system.addForce(force)

        # Add all GBSA terms to the system.
        expected_parameters = GBSAParameterHandler.GB_expected_parameters[self.gb_model]

        # Create all particles with parameters set to zero
        atoms = self.getMatches(topology)
        nparams = 1 + len(expected_parameters) # charge + GBSA parameters
        params = [ 0.0 for i in range(nparams) ]
        for particle in topology.particles():
            force.addParticle(params)
        # Set the GBSA parameters (keeping charges at zero for now)
        for (atoms, gbsa_type) in atoms.items():
            atom = atoms[0]
            # Set per-particle parameters for assigned parameters
            params = [atom.charge] + [ getattr(gbsa_type, name) for name in expected_parameters ]
            force.setParticleParameters(atom.particle_index, params)
