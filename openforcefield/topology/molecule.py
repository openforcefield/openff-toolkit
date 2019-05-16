#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================
"""
Molecular chemical entity representation and routines to interface with cheminformatics toolkits

.. todo::

   * Our main philosophy here is to keep the object contents of topology objects easily serializable/deserializable

   * Have ``Molecule`` raise an exception if loading/creating molecules with unspecified stereochemistry?
   * Create ``FrozenMolecule`` to represent immutable molecule
   * Make ``Atom`` and ``Bond`` an inner class of Molecule?
   * Add ``Molecule.from_smarts()`` or ``.from_tagged_smiles()`` to allow a tagged SMARTS string
     (where tags are zero-indexed atom indices) to be used to create a molecule with the given atom numbering.
   * How can we make the ``Molecule`` API more useful to codes like perses that modify molecules on the fly?
   * Use `attrs <http://www.attrs.org/>`_ for convenient class initialization?
   * JSON/BSON representations of objects?
   * Generalize Molecule infrastructure to provide "plug-in" support for cheminformatics toolkits
   * Do we need a way to write a bunch of molecules to a file, or serialize a set of molecules to a file?
     We currently don't have a way to do that through the ``Molecule`` API, even though there is a way to
     read multiple molecules via ``Molecules.from_file()``.
   * Should we allow the removal of atoms too?
   * Should invalidation of cached properties be handled via something like a tracked list?
   * Refactor toolkit encapsulation to generalize and provide only a few major toolkit methods and toolkit objects that can be queried for features
   * Speed up overall import time by putting non-global imports only where they are needed

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import numpy as np
from collections import OrderedDict
from copy import deepcopy

from simtk import unit
from simtk.openmm.app import element

import openforcefield
from openforcefield.utils import serialize_numpy, deserialize_numpy, quantity_to_string, string_to_quantity
from openforcefield.utils.toolkits import ToolkitRegistry, ToolkitWrapper, RDKitToolkitWrapper, OpenEyeToolkitWrapper,\
    InvalidToolkitError, GLOBAL_TOOLKIT_REGISTRY
from openforcefield.utils.toolkits import DEFAULT_AROMATICITY_MODEL
from openforcefield.utils.serialization import Serializable



#=============================================================================================
# GLOBAL PARAMETERS
#=============================================================================================

# TODO: Can we have the `ALLOWED_*_MODELS` list automatically appear in the docstrings below?
# TODO: Should `ALLOWED_*_MODELS` be objects instead of strings?
# TODO: Should these be imported from `openforcefield.cheminformatics.aromaticity_models` and `.bondorder_models`?

# TODO: Allow all OpenEye aromaticity models to be used with OpenEye names?
#       Only support OEAroModel_MDL in RDKit version?

#=============================================================================================
# PRIVATE SUBROUTINES
#=============================================================================================

#=============================================================================================
# Particle
#=============================================================================================


class Particle(Serializable):
    """
    Base class for all particles in a molecule.

    A particle object could be an ``Atom`` or a ``VirtualSite``.

    .. warning :: This API is experimental and subject to change.
    """

    @property
    def molecule(self):
        """
        The ``Molecule`` this atom is part of.

        .. todo::

           * Should we have a single unique ``Molecule`` for each molecule type in the system,
           or if we have multiple copies of the same molecule, should we have multiple ``Molecule``s?
        """
        return self._molecule

    @molecule.setter
    def molecule(self, molecule):
        """
        Set the particle's molecule pointer. Note that this will only work if the particle currently
        doesn't have a molecule
        """
        #TODO: Add informative exception here
        assert self._molecule is None
        self._molecule = molecule

    @property
    def molecule_particle_index(self):
        """
        Returns the index of this particle in its molecule
        """
        return self._molecule.particles.index(self)


    @property
    def name(self):
        """
        The name of the particle
        """
        return self._name

    def to_dict(self):
        """Convert to dictionary representation."""
        # Implement abstract method Serializable.to_dict()
        raise NotImplementedError()  # TODO

    @classmethod
    def from_dict(cls, d):
        """Static constructor from dictionary representation."""
        # Implement abstract method Serializable.to_dict()
        raise NotImplementedError()  # TODO


#=============================================================================================
# Atom
#=============================================================================================


class Atom(Particle):
    """
    A particle representing a chemical atom.

    Note that non-chemical virtual sites are represented by the ``VirtualSite`` object.

    .. todo::

       * Should ``Atom`` objects be immutable or mutable?
       * Do we want to support the addition of arbitrary additional properties,
         such as floating point quantities (e.g. ``charge``), integral quantities (such as ``id`` or ``serial`` index in a PDB file),
         or string labels (such as Lennard-Jones types)?

    .. todo :: Allow atoms to have associated properties.

    .. warning :: This API is experimental and subject to change.
    """

    def __init__(self,
                 atomic_number,
                 formal_charge,
                 is_aromatic,
                 name=None,
                 molecule=None,
                 stereochemistry=None):
        """
        Create an immutable Atom object.

        Object is serializable and immutable.

        .. todo :: Use attrs to validate?

        .. todo :: We can add setters if we need to.

        Parameters
        ----------
        atomic_number : int
            Atomic number of the atom
        formal_charge : int
            Formal charge of the atom
        is_aromatic : bool
            If True, atom is aromatic; if False, not aromatic
        stereochemistry : str, optional, default=None
            Either 'R' or 'S' for specified stereochemistry, or None for ambiguous stereochemistry
        name : str, optional, default=None
            An optional name to be associated with the atom

        Examples
        --------

        Create a non-aromatic carbon atom

        >>> atom = Atom(6, 0, False)

        Create a chiral carbon atom

        >>> atom = Atom(6, 0, False, stereochemistry='R', name='CT')

        """
        self._atomic_number = atomic_number
        self._formal_charge = formal_charge
        self._is_aromatic = is_aromatic
        self._stereochemistry = stereochemistry
        if name is None:
            name = ''
        self._name = name
        self._molecule = molecule
        ## From Jeff: I'm going to assume that this is implicit in the parent Molecule's ordering of atoms
        #self._molecule_atom_index = molecule_atom_index
        self._bonds = list()
        self._virtual_sites = list()

    # TODO: We can probably avoid an explicit call and determine this dynamically
    #   from self._molecule (maybe caching the result) to get rid of some bookkeeping.
    def add_bond(self, bond):
        """Adds a bond that this atom is involved in
        .. todo :: Is this how we want to keep records?

        Parameters
        ----------
        bond: an openforcefield.topology.molecule.Bond
            A bond involving this atom

"""

        self._bonds.append(bond)
        #self._stereochemistry = None

    def add_virtual_site(self, vsite):
        """Adds a bond that this atom is involved in
        .. todo :: Is this how we want to keep records?

        Parameters
        ----------
        bond: an openforcefield.topology.molecule.Bond
            A bond involving this atom

"""

        self._virtual_sites.append(vsite)

    def to_dict(self):
        """Return a dict representation of the atom.

"""
        # TODO
        atom_dict = OrderedDict()
        atom_dict['atomic_number'] = self._atomic_number
        atom_dict['formal_charge'] = self._formal_charge
        atom_dict['is_aromatic'] = self._is_aromatic
        atom_dict['stereochemistry'] = self._stereochemistry
        # TODO: Should we let atoms have names?
        atom_dict['name'] = self._name
        # TODO: Should this be implicit in the atom ordering when saved?
        #atom_dict['molecule_atom_index'] = self._molecule_atom_index
        return atom_dict

    @classmethod
    def from_dict(cls, atom_dict):
        """Create an Atom from a dict representation.

"""
        ## TODO: classmethod or static method? Classmethod is needed for Bond, so it have
        ## its _molecule set and then look up the Atom on each side of it by ID
        return cls.__init__(*atom_dict)

    @property
    def formal_charge(self):
        """
        The atom's formal charge
        """
        return self._formal_charge

    @property
    def partial_charge(self):
        """
        The partial charge of the atom, if any.

        Returns
        -------
        float or None
        """
        if self._molecule._partial_charges is None:
            return None
        else:
            index = self.molecule_atom_index
            return self._molecule._partial_charges[index]

    @property
    def is_aromatic(self):
        """
        The atom's is_aromatic flag
        """
        return self._is_aromatic

    @property
    def stereochemistry(self):
        """
        The atom's stereochemistry (if defined, otherwise None)
        """
        return self._stereochemistry

    @stereochemistry.setter
    def stereochemistry(self, value):
        """Set the atoms stereochemistry
        Parameters
        ----------
        value : str
            The stereochemistry around this atom, allowed values are "CW", "CCW", or None,
        """

        #if (value != 'CW') and (value != 'CCW') and not(value is None):
        #    raise Exception("Atom stereochemistry setter expected 'CW', 'CCW', or None. Received {} (type {})".format(value, type(value)))
        self._stereochemistry = value

    @property
    def element(self):
        """
        The element name

        """
        return element.Element.getByAtomicNumber(self._atomic_number)

    @property
    def atomic_number(self):
        """
        The integer atomic number of the atom.

        """
        return self._atomic_number

    @property
    def mass(self):
        """
        The standard atomic weight (abundance-weighted isotopic mass) of the atomic site.

        .. todo :: Should we discriminate between standard atomic weight and most abundant isotopic mass?
        TODO (from jeff): Are there atoms that have different chemical properties based on their isotopes?

        """
        return self.element.mass

    @property
    def name(self):
        """
        The name of this atom, if any
        """
        return self._name

    @name.setter
    def name(self, other):
        """

        Parameters
        ----------
        other : string
            The new name for this atom
        """
        if not (type(other) is str):
            raise Exception(
                "In setting atom name. Expected str, received {} (type {})".
                format(other, type(other)))
        self._name = other

    # TODO: How are we keeping track of bonds, angles, etc?

    @property
    def bonds(self):
        """
        The list of ``Bond`` objects this atom is involved in.

        """
        return self._bonds
        #for bond in self._bonds:
        #    yield bond

    @property
    #def bonded_to(self):
    def bonded_atoms(self):
        """
        The list of ``Atom`` objects this atom is involved in bonds with

        """
        for bond in self._bonds:
            for atom in bond.atoms:
                if not (atom == self):
                    # TODO: This seems dangerous. Ask John for a better way
                    yield atom

    def is_bonded_to(self, atom2):
        """
        Determine whether this atom is bound to another atom

        Parameters
        ----------
        atom2: openforcefield.topology.molecule.Atom
            a different atom in the same molecule

        Returns
        -------
        bool
            Whether this atom is bound to atom2
        """
        #TODO: Sanity check (check for same molecule?)
        assert self != atom2
        for bond in self._bonds:
            for bonded_atom in bond.atoms:
                if atom2 == bonded_atom:
                    return True
        return False

    @property
    def virtual_sites(self):
        """
        The list of ``VirtualSite`` objects this atom is involved in.

        """
        return self._virtual_sites
        #for vsite in self._vsites:
        #    yield vsite

    @property
    def molecule_atom_index(self):
        """
        The index of this Atom within the the list of atoms in ``Molecules``.
        Note that this can be different from ``molecule_particle_index``.

        """
        if self._molecule is None:
            raise ValueError('This Atom does not belong to a Molecule object')
        return self._molecule.atoms.index(self)

    @property
    def molecule_particle_index(self):
        """
        The index of this Atom within the the list of particles in the parent ``Molecule``.
        Note that this can be different from ``molecule_atom_index``.

        """
        if self._molecule is None:
            raise ValueError('This Atom does not belong to a Molecule object')
        return self._molecule.particles.index(self)

    # ## From Jeff: Not sure if we actually need this
    # @property
    # def topology_atom_index(self):
    #     """
    #     The index of this Atom within the the list of atoms in ``Topology``.
    #     Note that this can be different from ``particle_index``.
    #
    #     """
    #     if self._topology is None:
    #         raise ValueError('This Atom does not belong to a Topology object')
    #     # TODO: This will be slow; can we cache this and update it only when needed?
    #     #       Deleting atoms/molecules in the Topology would have to invalidate the cached index.
    #     return self._topology.atoms.index(self)

    def __repr__(self):
        # TODO: Also include particle_index and which molecule this atom belongs to?
        return "Atom(name={}, atomic number={})".format(
            self._name, self._atomic_number)

    def __str__(self):
        # TODO: Also include particle_index and which molecule this atom belongs to?
        return "<Atom name='{}' atomic number='{}'>".format(
            self._name, self._atomic_number)


#=============================================================================================
# VirtualSite
#=============================================================================================


class VirtualSite(Particle):
    """
    A particle representing a virtual site whose position is defined in terms of ``Atom`` positions.

    Note that chemical atoms are represented by the ``Atom``.


    .. warning :: This API is experimental and subject to change.

    .. todo::

       * Should a virtual site be able to belong to more than one Topology?
       * Should virtual sites be immutable or mutable?

    """

    def __init__(self,
                 atoms,
                 charge_increments=None,
                 epsilon=None,
                 sigma=None,
                 rmin_half=None,
                 name=None):
        """
        Base class for VirtualSites

        .. todo ::

           * This will need to be generalized for virtual sites to allow out-of-plane sites, which are not simply a linear combination of atomic positions
           * Add checking for allowed virtual site types
           * How do we want users to specify virtual site types?

        Parameters
        ----------
        atoms : list of Atom of shape [N]
            atoms[index] is the corresponding Atom for weights[index]
        charge_increments : list of floats of shape [N], optional, default=None
            The amount of charge to remove from the VirtualSite's atoms and put in the VirtualSite. Indexing in this list should match the ordering in the atoms list. Default is None.
        sigma : float, default=None
            Sigma term for VdW properties of virtual site. Default is None.
        epsilon : float
            Epsilon term for VdW properties of virtual site. Default is None.
        rmin_half : float
            Rmin_half term for VdW properties of virtual site. Default is None.
        name : string or None, default=None
            The name of this virtual site. Default is None.


        virtual_site_type : str
            Virtual site type.
        name : str or None, default=None
            The name of this virtual site. Default is None

"""

        # Ensure we have as many charge_increments as we do atoms
        if not (charge_increments is None):
            if not (len(charge_increments) == len(atoms)):
                raise Exception(
                    "VirtualSite definition must have same number of charge_increments ({}) and atoms({})"
                    .format(len(charge_increments), len(atoms)))

        # VdW parameters can either be epsilon+rmin_half or epsilon+sigma, but not both
        if not (epsilon is None):
            if ((rmin_half != None) and (sigma != None)):
                raise Exception(
                    "VirtualSite constructor given epsilon (value : {}), rmin_half (value : {}), and sigma (value : {}). If epsilon is nonzero, it should receive either rmin_half OR sigma"
                    .format(epsilon, rmin_half, sigma))
            if ((rmin_half is None) and (sigma is None)):
                raise Exception(
                    "VirtualSite constructor given epsilon (value : {}) but not given rmin_half (value : {}) or sigma (value : {})"
                    .format(epsilon, rmin_half, sigma))
            if (sigma is None):
                # TODO: Save the 6th root of 2 if this starts being slow.
                sigma = (2. * rmin_half) / (2.**(1. / 6))

        elif epsilon is None:
            if ((rmin_half != None) or (sigma != None)):
                raise Exception(
                    "VirtualSite constructor given rmin_half (value : {}) or sigma (value : {}), but not epsilon (value : {})"
                    .format(rmin_half, sigma, epsilon))

        # Perform type-checking
        for atom in atoms:
            assert isinstance(atom, Atom)
        for atom_index in range(len(atoms) - 1):
            assert atoms[atom_index].molecule is atoms[atom_index + 1].molecule
        assert isinstance(atoms[1].molecule, FrozenMolecule)

        if sigma is None:
            self._sigma = None
        else:
            assert hasattr(sigma, 'unit')
            assert unit.angstrom.is_compatible(sigma.unit)
            self._sigma = sigma.in_units_of(unit.angstrom)

        if epsilon is None:
            self._epsilon = None
        else:
            assert hasattr(epsilon, 'unit')
            assert (unit.kilojoule_per_mole).is_compatible(epsilon.unit)
            self._epsilon = epsilon.in_units_of(unit.kilojoule_per_mole)

        if charge_increments is None:
            self._charge_increments = None
        else:
            assert hasattr(charge_increments, 'unit')
            assert unit.elementary_charges.is_compatible(
                charge_increments.unit)
            self._charge_increments = charge_increments.in_units_of(
                unit.elementary_charges)

        self._atoms = list()
        for atom in atoms:
            atom.add_virtual_site(self)
            self._atoms.append(atom)
        self._molecule = atoms[0].molecule

        self._name = name

        # Subclassing makes _type unnecessary
        #self._type = None
        # TODO: Validate site types against allowed values

        #self._weights = np.array(weights) # make a copy and convert to array internally

    def to_dict(self):
        """Return a dict representation of the virtual site.

"""
        # Each subclass should have its own to_dict
        vsite_dict = OrderedDict()
        vsite_dict['name'] = self._name
        vsite_dict['atoms'] = tuple(
            [i.molecule_atom_index for i in self.atoms])
        vsite_dict['charge_increments'] = quantity_to_string(self._charge_increments)


        vsite_dict['epsilon'] = quantity_to_string(self._epsilon)


        vsite_dict['sigma'] = quantity_to_string(self._sigma)

        return vsite_dict


    @classmethod
    def from_dict(cls, vsite_dict):
        """Create a virtual site from a dict representation.

"""
        # Each subclass needs to have its own from_dict

        # Make a copy of the vsite_dict, where we'll unit-wrap the appropriate values
        vsite_dict_units = deepcopy(vsite_dict)

        # Attach units to epsilon term
        vsite_dict_units['epsilon'] = string_to_quantity(
            vsite_dict['epsilon'])
        vsite_dict_units['sigma'] = string_to_quantity(vsite_dict['sigma'])
        vsite_dict_units['charge_increments'] = string_to_quantity(
            vsite_dict['charge_increments'])

        return VirtualSite(**vsite_dict_units)

    @property
    def molecule_virtual_site_index(self):
        """
        The index of this VirtualSite within the list of virtual sites within ``Molecule``
        Note that this can be different from ``particle_index``.
        """
        #if self._topology is None:
        #    raise ValueError('This VirtualSite does not belong to a Topology object')
        # TODO: This will be slow; can we cache this and update it only when needed?
        #       Deleting atoms/molecules in the Topology would have to invalidate the cached index.
        return self._molecule.virtual_sites.index(self)

    @property
    def molecule_particle_index(self):
        """
        The index of this VirtualSite within the the list of particles in the parent ``Molecule``.
        Note that this can be different from ``molecule_virtual_site_index``.

        """
        if self._molecule is None:
            raise ValueError(
                'This VirtualSite does not belong to a Molecule object')
        return self._molecule.particles.index(self)

    @property
    def atoms(self):
        """
        Atoms on whose position this VirtualSite depends.
        """
        return self._atoms
        #for atom in self._atoms:
        #    yield atom

    @property
    def charge_increments(self):
        """
        Charges taken from this VirtualSite's atoms and given to the VirtualSite
        """
        return self._charge_increments

    @property
    def epsilon(self):
        """
        The VdW epsilon term of this VirtualSite
        """
        return self._epsilon

    @property
    def sigma(self):
        """
        The VdW sigma term of this VirtualSite
        """
        return self._sigma

    @property
    def rmin_half(self):
        """
        The VdW rmin_half term of this VirtualSite
        """
        rmin = 2.**(1. / 6) * self._sigma
        rmin_half = rmin / 2
        return rmin_half

    @property
    def name(self):
        """
        The name of this VirtualSite
        """
        return self._name

    @property
    def type(self):
        """The type of this VirtualSite (returns the class name as string)"""
        return self.__class__.__name__

    def __repr__(self):
        # TODO: Also include particle_index, which molecule this atom belongs to?
        return "VirtualSite(name={}, type={}, atoms={})".format(
            self.name, self.type, self.atoms)

    def __str__(self):
        # TODO: Also include particle_index, which molecule this atom belongs to?
        return "<VirtualSite name={} type={} atoms={}>".format(
            self.name, self.type, self.atoms)


class BondChargeVirtualSite(VirtualSite):
    """
    A particle representing a "Bond Charge"-type virtual site, in which the location of the charge is specified by the positions of two atoms. This supports placement of a virtual site S along a vector between two specified atoms, e.g. to allow for a sigma hole for halogens or similar contexts. With positive values of the distance, the virtual site lies outside the first indexed atom.

    .. warning :: This API is experimental and subject to change.
    """

    def __init__(self,
                 atoms,
                 distance,
                 charge_increments=None,
                 weights=None,
                 epsilon=None,
                 sigma=None,
                 rmin_half=None,
                 name=None):
        """
        Create a bond charge-type virtual site, in which the location of the charge is specified by the position of two atoms. This supports placement of a virtual site S along a vector between two specified atoms, e.g. to allow for a sigma hole for halogens or similar contexts. With positive values of the distance, the virtual site lies outside the first indexed atom.

        TODO: One of the examples on https://open-forcefield-toolkit.readthedocs.io/en/topology/smirnoff.html#virtualsites-virtual-sites-for-off-atom-charges has a BondCharge defined with three atoms -- How does that work?

        Parameters
        ----------
        atoms : list of openforcefield.topology.molecule.Atom objects of shape [N]
            The atoms defining the virtual site's position
        distance : float

        weights : list of floats of shape [N] or None, optional, default=None
            weights[index] is the weight of particles[index] contributing to the position of the virtual site. Default is None
        charge_increments : list of floats of shape [N], optional, default=None
            The amount of charge to remove from the VirtualSite's atoms and put in the VirtualSite. Indexing in this list should match the ordering in the atoms list. Default is None.
        epsilon : float
            Epsilon term for VdW properties of virtual site. Default is None.
        sigma : float, default=None
            Sigma term for VdW properties of virtual site. Default is None.
        rmin_half : float
            Rmin_half term for VdW properties of virtual site. Default is None.
        name : string or None, default=None
            The name of this virtual site. Default is None.
        """
        assert hasattr(distance, 'unit')
        assert unit.angstrom.is_compatible(distance.unit)

        super().__init__(
            atoms,
            charge_increments=charge_increments,
            epsilon=epsilon,
            sigma=sigma,
            rmin_half=rmin_half,
            name=name)
        self._distance = distance.in_units_of(unit.angstrom)

    def to_dict(self):
        vsite_dict = super().to_dict()
        vsite_dict['distance'] = quantity_to_string(self._distance)

        #type = self.type
        vsite_dict['vsite_type'] = self.type
        #vsite_dict['vsite_type'] = 'BondChargeVirtualSite'
        return vsite_dict

    @classmethod
    def from_dict(cls, vsite_dict):
        base_dict = deepcopy(vsite_dict)
        # Make sure it's the right type of virtual site
        assert vsite_dict['vsite_type'] == "BondChargeVirtualSite"
        base_dict.pop('vsite_type')
        base_dict.pop('distance')
        vsite = super().from_dict(**base_dict)
        vsite._distance = string_to_quantity(vsite_dict['distance'])
        return vsite

    @property
    def distance(self):
        """The distance parameter of the virtual site"""
        return self._distance


class _LonePairVirtualSite(VirtualSite):
    """Private base class for mono/di/trivalent lone pair virtual sites."""

    @classmethod
    def from_dict(cls, vsite_dict):
        base_dict = deepcopy(vsite_dict)

        # Make sure it's the right type of virtual site
        assert vsite_dict['vsite_type'] == cls.__name__
        base_dict.pop('vsite_type')
        base_dict.pop('distance')
        base_dict.pop('out_of_plane_angle')
        base_dict.pop('in_plane_angle')
        vsite = super().from_dict(**base_dict)
        vsite._distance = string_to_quantity(vsite_dict['distance'])
        vsite._in_plane_angle = string_to_quantity(
            vsite_dict['in_plane_angle'])
        vsite._out_of_plane_angle = string_to_quantity(
            vsite_dict['out_of_plane_angle'])
        return vsite


class MonovalentLonePairVirtualSite(VirtualSite):
    """
    A particle representing a "Monovalent Lone Pair"-type virtual site, in which the location of the charge is specified by the positions of three atoms. This is originally intended for situations like a carbonyl, and allows placement of a virtual site S at a specified distance d, in_plane_angle, and out_of_plane_angle relative to a central atom and two connected atoms.

    .. warning :: This API is experimental and subject to change.
   """

    def __init__(self,
                 atoms,
                 distance,
                 out_of_plane_angle,
                 in_plane_angle,
                 charge_increments=None,
                 weights=None,
                 epsilon=None,
                 sigma=None,
                 rmin_half=None,
                 name=None):
        """
        Create a bond charge-type virtual site, in which the location of the charge is specified by the position of three atoms.

        TODO : Do "weights" have any meaning here?

        Parameters
        ----------
        atoms : list of three openforcefield.topology.molecule.Atom objects
            The three atoms defining the virtual site's position
        distance : float

        out_of_plane_angle : float

        in_plane_angle : float

        epsilon : float
            Epsilon term for VdW properties of virtual site. Default is None.
        sigma : float, default=None
            Sigma term for VdW properties of virtual site. Default is None.
        rmin_half : float
            Rmin_half term for VdW properties of virtual site. Default is None.
        name : string or None, default=None
            The name of this virtual site. Default is None.


        """
        #assert isinstance(distance, unit.Quantity)
        # TODO: Check for proper number of atoms
        assert hasattr(distance, 'unit')
        assert unit.angstrom.is_compatible(distance.unit)
        assert hasattr(in_plane_angle, 'unit')
        assert unit.degree.is_compatible(in_plane_angle.unit)
        assert hasattr(out_of_plane_angle, 'unit')
        assert unit.degree.is_compatible(out_of_plane_angle.unit)

        assert len(atoms) == 3
        super().__init__(
            atoms,
            charge_increments=charge_increments,
            epsilon=epsilon,
            sigma=sigma,
            rmin_half=rmin_half,
            name=name)
        self._distance = distance.in_units_of(unit.angstrom)
        self._out_of_plane_angle = out_of_plane_angle.in_units_of(unit.degree)
        self._in_plane_angle = in_plane_angle.in_units_of(unit.degree)

    def to_dict(self):
        vsite_dict = super().to_dict()
        vsite_dict['distance'] = quantity_to_string(self._distance)
        vsite_dict['out_of_plane_angle'] = quantity_to_string(
            self._out_of_plane_angle)
        vsite_dict['in_plane_angle'] = quantity_to_string(self._in_plane_angle)
        #vsite_dict['vsite_type'] = self.type.fget()
        vsite_dict['vsite_type'] = self.type
        #vsite_dict['vsite_type'] = 'MonovalentLonePairVirtualSite'
        return vsite_dict

    @classmethod
    def from_dict(cls, vsite_dict):
        """
        Construct a new MonovalentLonePairVirtualSite from an serialized dictionary representation.

        Parameters
        ----------
        vsite_dict : dict
            The VirtualSite to deserialize.

        Returns
        -------
        The newly created MonovalentLonePairVirtualSite

        """
        # The function is overridden only to have a custom docstring.
        return super().from_dict(vsite_dict)

    @property
    def distance(self):
        """The distance parameter of the virtual site"""
        return self._distance

    @property
    def in_plane_angle(self):
        """The in_plane_angle parameter of the virtual site"""
        return self._in_plane_angle

    @property
    def out_of_plane_angle(self):
        """The out_of_plane_angle parameter of the virtual site"""
        return self._out_of_plane_angle


class DivalentLonePairVirtualSite(VirtualSite):
    """
    A particle representing a "Divalent Lone Pair"-type virtual site, in which the location of the charge is specified by the positions of three atoms. This is suitable for cases like four-point and five-point water models as well as pyrimidine; a charge site S lies a specified distance d from the central atom among three atoms along the bisector of the angle between the atoms (if out_of_plane_angle is zero) or out of the plane by the specified angle (if out_of_plane_angle is nonzero) with its projection along the bisector. For positive values of the distance d the virtual site lies outside the 2-1-3 angle and for negative values it lies inside.

    .. warning :: This API is experimental and subject to change.
    """

    def __init__(self,
                 atoms,
                 distance,
                 out_of_plane_angle,
                 in_plane_angle,
                 charge_increments=None,
                 weights=None,
                 epsilon=None,
                 sigma=None,
                 rmin_half=None,
                 name=None):
        """
        Create a divalent lone pair-type virtual site, in which the location of the charge is specified by the position of three atoms.

        TODO : Do "weights" have any meaning here?

        Parameters
        ----------
        atoms : list of 3 openforcefield.topology.molecule.Atom objects
            The three atoms defining the virtual site's position
        distance : float

        out_of_plane_angle : float

        in_plane_angle : float

        epsilon : float
            Epsilon term for VdW properties of virtual site. Default is None.
        sigma : float, default=None
            Sigma term for VdW properties of virtual site. Default is None.
        rmin_half : float
            Rmin_half term for VdW properties of virtual site. Default is None.
        name : string or None, default=None
            The name of this virtual site. Default is None.
        """
        #assert isinstance(distance, unit.Quantity)
        assert hasattr(distance, 'unit')
        assert unit.angstrom.is_compatible(distance.unit)
        assert hasattr(in_plane_angle, 'unit')
        assert unit.degree.is_compatible(in_plane_angle.unit)
        assert hasattr(out_of_plane_angle, 'unit')
        assert unit.degree.is_compatible(out_of_plane_angle.unit)
        assert len(atoms) == 3
        super().__init__(
            atoms,
            charge_increments=charge_increments,
            epsilon=epsilon,
            sigma=sigma,
            rmin_half=rmin_half,
            name=name)
        self._distance = distance.in_units_of(unit.angstrom)
        self._out_of_plane_angle = out_of_plane_angle.in_units_of(unit.degree)
        self._in_plane_angle = in_plane_angle.in_units_of(unit.degree)

    def to_dict(self):
        vsite_dict = super().to_dict()
        vsite_dict['distance'] = quantity_to_string(self._distance)
        vsite_dict['out_of_plane_angle'] = quantity_to_string(
            self._out_of_plane_angle)
        vsite_dict['in_plane_angle'] = quantity_to_string(self._in_plane_angle)
        vsite_dict['vsite_type'] = self.type
        return vsite_dict

    @classmethod
    def from_dict(cls, vsite_dict):
        """
        Construct a new DivalentLonePairVirtualSite from an serialized dictionary representation.

        Parameters
        ----------
        vsite_dict : dict
            The VirtualSite to deserialize.

        Returns
        -------
        The newly created DivalentLonePairVirtualSite

        """
        # The function is overridden only to have a custom docstring.
        return super().from_dict(vsite_dict)

    @property
    def distance(self):
        """The distance parameter of the virtual site"""
        return self._distance

    @property
    def in_plane_angle(self):
        """The in_plane_angle parameter of the virtual site"""
        return self._in_plane_angle

    @property
    def out_of_plane_angle(self):
        """The out_of_plane_angle parameter of the virtual site"""
        return self._out_of_plane_angle


class TrivalentLonePairVirtualSite(VirtualSite):
    """
    A particle representing a "Trivalent Lone Pair"-type virtual site, in which the location of the charge is specified by the positions of four atoms. This is suitable for planar or tetrahedral nitrogen lone pairs; a charge site S lies above the central atom (e.g. nitrogen a distance d along the vector perpendicular to the plane of the three connected atoms (2,3,4). With positive values of d the site lies above the nitrogen and with negative values it lies below the nitrogen.

    .. warning :: This API is experimental and subject to change.
    """

    def __init__(self,
                 atoms,
                 distance,
                 out_of_plane_angle,
                 in_plane_angle,
                 charge_increments=None,
                 weights=None,
                 epsilon=None,
                 sigma=None,
                 rmin_half=None,
                 name=None):
        """
        Create a trivalent lone pair-type virtual site, in which the location of the charge is specified by the position of four atoms.

        TODO : Do "weights" have any meaning here?

        Parameters
        ----------
        atoms : list of 4 openforcefield.topology.molecule.Atom objects
            The three atoms defining the virtual site's position
        distance : float

        out_of_plane_angle : float

        in_plane_angle : float

        epsilon : float
            Epsilon term for VdW properties of virtual site. Default is None.
        sigma : float, default=None
            Sigma term for VdW properties of virtual site. Default is None.
        rmin_half : float
            Rmin_half term for VdW properties of virtual site. Default is None.
        name : string or None, default=None
            The name of this virtual site. Default is None.

        """
        assert len(atoms) == 4
        assert hasattr(distance, 'unit')
        assert unit.angstrom.is_compatible(distance.unit)
        assert hasattr(in_plane_angle, 'unit')
        assert unit.degree.is_compatible(in_plane_angle.unit)
        assert hasattr(out_of_plane_angle, 'unit')
        assert unit.degree.is_compatible(out_of_plane_angle.unit)

        super().__init__(
            atoms,
            charge_increments=charge_increments,
            epsilon=epsilon,
            sigma=sigma,
            rmin_half=rmin_half,
            name=name)
        self._distance = distance.in_units_of(unit.angstrom)
        self._out_of_plane_angle = out_of_plane_angle.in_units_of(unit.degree)
        self._in_plane_angle = in_plane_angle.in_units_of(unit.degree)

    def to_dict(self):
        vsite_dict = super().to_dict()
        vsite_dict['distance'] = quantity_to_string(self._distance)
        vsite_dict['out_of_plane_angle'] = quantity_to_string(
            self._out_of_plane_angle)
        vsite_dict['in_plane_angle'] = quantity_to_string(self._in_plane_angle)
        vsite_dict['vsite_type'] = self.type
        return vsite_dict

    @classmethod
    def from_dict(cls, vsite_dict):
        """
        Construct a new TrivalentPairVirtualSite from an serialized dictionary representation.

        Parameters
        ----------
        vsite_dict : dict
            The VirtualSite to deserialize.


        Returns
        -------
        The newly created TrivalentLonePairVirtualSite

        """
        # The function is overridden only to have a custom docstring.
        return super().from_dict(vsite_dict)

    @property
    def distance(self):
        """The distance parameter of the virtual site"""
        return self._distance

    @property
    def in_plane_angle(self):
        """The in_plane_angle parameter of the virtual site"""
        return self._in_plane_angle

    @property
    def out_of_plane_angle(self):
        """The out_of_plane_angle parameter of the virtual site"""
        return self._out_of_plane_angle


# =============================================================================================
# Bond Stereochemistry
# =============================================================================================

#class BondStereochemistry(Serializable):
#"""
#Bond stereochemistry representation
#"""
#def __init__(self, stereo_type, neighbor1, neighbor2):
#    """
#
#    Parameters
#    ----------
#    stereo_type
#    neighbor1
#    neighbor2
#    """
#    assert isinstance(neighbor1, Atom)
#    assert isinstance(neighbor2, Atom)
#    # Use stereo_type @setter to check stereo type is a permitted value
#    self.stereo_type = stereo_type
#    self._neighbor1 = neighbor1
#    self._neighbor2 = neighbor2

#def to_dict(self):
#    bs_dict = OrderedDict()
#    bs_dict['stereo_type'] = self._stereo_type
#    bs_dict['neighbor1_index'] = self._neighbor1.molecule_atom_index
#    bs_dict['neighbor2_index'] = self._neighbor2.molecule_atom_index
#    return bs_dict

#classmethod
#def from_dict(cls, molecule, bs_dict):
#    neighbor1 = molecule.atoms[bs_dict['neighbor1_index']]
#    neighbor2 = molecule.atoms[bs_dict['neighbor2_index']]
#    return cls.__init__(bs_dict['stereo_type'], neighbor1, neighbor2)

#@property
#def stereo_type(self):
#    return self._stereo_type

#@stereo_type.setter
#def stereo_type(self, value):
#    assert (value == 'CIS') or (value == 'TRANS') or (value is None)
#    self._stereo_type = value

#@property
#def neighbor1(self):
#    return self._neighbor1

#@property
#def neighbor2(self):
#    return self._neighbor2

#@property
#def neighbors(self):
#    return (self._neighbor1, self._neighbor2)

# =============================================================================================
# Bond
# =============================================================================================


class Bond(Serializable):
    """
    Chemical bond representation.

    .. warning :: This API is experimental and subject to change.

    .. todo :: Allow bonds to have associated properties.

    Attributes
    ----------
    atom1, atom2 : openforcefield.topology.Atom
        Atoms involved in the bond
    bondtype : int
        Discrete bond type representation for the Open Forcefield aromaticity model
        TODO: Do we want to pin ourselves to a single standard aromaticity model?
    type : str
        String based bond type
    order : int
        Integral bond order
    fractional_bond_order : float, optional
        Fractional bond order, or None.


    .. warning :: This API is experimental and subject to change.
    """

    def __init__(self,
                 atom1,
                 atom2,
                 bond_order,
                 is_aromatic,
                 fractional_bond_order=None,
                 stereochemistry=None):
        """
        Create a new chemical bond.

        """
        assert type(atom1) == Atom
        assert type(atom2) == Atom
        assert atom1.molecule is atom2.molecule
        assert isinstance(atom1.molecule, FrozenMolecule)
        self._molecule = atom1.molecule

        self._atom1 = atom1
        self._atom2 = atom2

        atom1.add_bond(self)
        atom2.add_bond(self)
        # TODO: Check bondtype and fractional_bond_order are valid?
        # TODO: Dative bonds
        #self._type = bondtype
        self._fractional_bond_order = fractional_bond_order
        self._bond_order = bond_order
        self._is_aromatic = is_aromatic
        self._stereochemistry = stereochemistry

    def to_dict(self):
        """
        Return a dict representation of the bond.

        """
        bond_dict = OrderedDict()
        bond_dict['atom1'] = self.atom1.molecule_atom_index
        bond_dict['atom2'] = self.atom2.molecule_atom_index
        bond_dict['bond_order'] = self._bond_order
        bond_dict['is_aromatic'] = self._is_aromatic
        bond_dict['stereochemistry'] = self._stereochemistry
        bond_dict['fractional_bond_order'] = self._fractional_bond_order
        return bond_dict

    @classmethod
    def from_dict(cls, molecule, d):
        """Create a Bond from a dict representation."""
        # TODO
        d['molecule'] = molecule
        d['atom1'] = molecule.atoms[d['atom1']]
        d['atom2'] = molecule.atoms[d['atom2']]
        return cls(*d)

    @property
    def atom1(self):
        return self._atom1

    @property
    def atom2(self):
        return self._atom2

    @property
    def atom1_index(self):
        return self.molecule.atoms.index(self._atom1)

    @property
    def atom2_index(self):
        return self.molecule.atoms.index(self._atom2)

    @property
    def atoms(self):
        return (self._atom1, self._atom2)

    @property
    def bond_order(self):
        return self._bond_order

    @bond_order.setter
    def bond_order(self, value):
        self._bond_order = value

    @property
    def fractional_bond_order(self):
        return self._fractional_bond_order

    @fractional_bond_order.setter
    def fractional_bond_order(self, value):
        self._fractional_bond_order = value

    @property
    def stereochemistry(self):
        return self._stereochemistry

    @property
    def is_aromatic(self):
        return self._is_aromatic

    @property
    def molecule(self):
        return self._molecule

    @molecule.setter
    def molecule(self, value):
        """
        Sets the Bond's parent molecule. Can not be changed after assignment
        """
        assert self._molecule is None
        self._molecule = value

    @property
    def molecule_bond_index(self):
        """
        The index of this Bond within the the list of bonds in ``Molecules``.

        """
        if self._molecule is None:
            raise ValueError('This Atom does not belong to a Molecule object')
        return self._molecule.bonds.index(self)


#=============================================================================================
# Molecule
#=============================================================================================

# TODO: How do we automatically trigger invalidation of cached properties if an ``Atom``, ``Bond``, or ``VirtualSite`` is modified,
#       rather than added/deleted via the API? The simplest resolution is simply to make them immutable.


class FrozenMolecule(Serializable):
    """
    Immutable chemical representation of a molecule, such as a small molecule or biopolymer.

    .. todo :: What other API calls would be useful for supporting biopolymers as small molecules? Perhaps iterating over chains and residues?

    Examples
    --------

    Create a molecule from a sdf file

    >>> from openforcefield.utils import get_data_file_path
    >>> sdf_filepath = get_data_file_path('molecules/ethanol.sdf')
    >>> molecule = FrozenMolecule.from_file(sdf_filepath)

    Convert to OpenEye OEMol object

    >>> oemol = molecule.to_openeye()

    Create a molecule from an OpenEye molecule

    >>> molecule = FrozenMolecule.from_openeye(oemol)

    Convert to RDKit Mol object

    >>> rdmol = molecule.to_rdkit()

    Create a molecule from an RDKit molecule

    >>> molecule = FrozenMolecule.from_rdkit(rdmol)

    Create a molecule from IUPAC name (requires the OpenEye toolkit)

    >>> molecule = FrozenMolecule.from_iupac('imatinib')

    Create a molecule from SMILES

    >>> molecule = FrozenMolecule.from_smiles('Cc1ccccc1')

    .. warning :: This API is experimental and subject to change.


    """

    def __init__(self,
                 other=None,
                 file_format=None,
                 toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
                 allow_undefined_stereo=False):
        """
        Create a new FrozenMolecule object

        .. todo ::

           * If a filename or file-like object is specified but the file contains more than one molecule, what is the proper behavior?
           Read just the first molecule, or raise an exception if more than one molecule is found?

           * Should we also support SMILES strings or IUPAC names for ``other``?

        Parameters
        ----------
        other : optional, default=None
            If specified, attempt to construct a copy of the Molecule from the specified object.
            This can be any one of the following:

            * a :class:`Molecule` object
            * a file that can be used to construct a :class:`Molecule` object
            * an ``openeye.oechem.OEMol``
            * an ``rdkit.Chem.rdchem.Mol``
            * a serialized :class:`Molecule` object
        file_format : str, optional, default=None
            If providing a file-like object, you must specify the format
            of the data. If providing a file, the file format will attempt
            to be guessed from the suffix.
        toolkit_registry : a :class:`ToolkitRegistry` or :class:`ToolkitWrapper` object, optional, default=GLOBAL_TOOLKIT_REGISTRY
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for I/O operations
        allow_undefined_stereo : bool, default=False
            If loaded from a file and ``False``, raises an exception if
            undefined stereochemistry is detected during the molecule's
            construction.

        Examples
        --------

        Create an empty molecule:

        >>> empty_molecule = FrozenMolecule()

        Create a molecule from a file that can be used to construct a molecule,
        using either a filename or file-like object:

        >>> from openforcefield.utils import get_data_file_path
        >>> sdf_filepath = get_data_file_path('molecules/ethanol.sdf')
        >>> molecule = FrozenMolecule(sdf_filepath)
        >>> molecule = FrozenMolecule(open(sdf_filepath, 'r'), file_format='sdf')

        >>> import gzip
        >>> mol2_gz_filepath = get_data_file_path('molecules/toluene.mol2.gz')
        >>> molecule = FrozenMolecule(gzip.GzipFile(mol2_gz_filepath, 'r'), file_format='mol2')

        Create a molecule from another molecule:

        >>> molecule_copy = FrozenMolecule(molecule)

        Convert to OpenEye OEMol object

        >>> oemol = molecule.to_openeye()

        Create a molecule from an OpenEye molecule:

        >>> molecule = FrozenMolecule(oemol)

        Convert to RDKit Mol object

        >>> rdmol = molecule.to_rdkit()

        Create a molecule from an RDKit molecule:

        >>> molecule = FrozenMolecule(rdmol)

        Create a molecule from a serialized molecule object:

        >>> serialized_molecule = molecule.__getstate__()
        >>> molecule_copy = Molecule(serialized_molecule)

        """

        self._cached_smiles = None

        # Figure out if toolkit_registry is a whole registry, or just a single wrapper
        if isinstance(toolkit_registry, ToolkitRegistry):
            pass
        elif isinstance(toolkit_registry, ToolkitWrapper):
            toolkit = toolkit_registry
            toolkit_registry = ToolkitRegistry(toolkit_precedence=[])
            toolkit_registry.add_toolkit(toolkit)
        else:
            raise ValueError(
                "'toolkit_registry' must be either a ToolkitRegistry or a ToolkitWrapper"
            )

        if other is None:
            self._initialize()
        else:
            loaded = False
            if isinstance(
                    other,
                    openforcefield.topology.FrozenMolecule) and not (loaded):
                self._copy_initializer(other)
                loaded = True
            if isinstance(other,
                          openforcefield.topology.Molecule) and not (loaded):
                # TODO: This will need to be updated once FrozenMolecules and Molecules are significantly different
                self._copy_initializer(other)
                loaded = True
            if isinstance(other, OrderedDict) and not (loaded):
                self.__setstate__(other)
                loaded = True
            # Check through the toolkit registry to find a compatible wrapper for loading
            if not loaded:
                try:
                    result = toolkit_registry.call('from_object', other)
                except NotImplementedError:
                    pass
                else:
                    self._copy_initializer(result)
                    loaded = True
            # TODO: Make this compatible with file-like objects (I couldn't figure out how to make an oemolistream
            # from a fileIO object)
            if (isinstance(other, str)
                    or hasattr(other, 'read')) and not (loaded):
                mol = Molecule.from_file(
                    other,
                    file_format=file_format,
                    toolkit_registry=toolkit_registry,
                    allow_undefined_stereo=allow_undefined_stereo
                )  # returns a list only if multiple molecules are found
                if type(mol) == list:
                    raise ValueError(
                        'Specified file or file-like object must contain exactly one molecule'
                    )

                self._copy_initializer(mol)
                loaded = True
            if not (loaded):
                msg = 'Cannot construct openforcefield.topology.Molecule from {}\n'.format(
                    other)
                raise Exception(msg)

    ####################################################################################################
    # Safe serialization
    ####################################################################################################

    def to_dict(self):
        """
        Return a dictionary representation of the molecule.

        .. todo ::

           * Document the representation standard.
           * How do we do version control with this standard?

        Returns
        -------
        molecule_dict : OrderedDict
            A dictionary representation of the molecule.

        """
        molecule_dict = OrderedDict()
        molecule_dict['name'] = self._name
        ## From Jeff: If we go the properties-as-dict route, then _properties should, at
        ## the top level, be a dict. Should we go through recursively and ensure all values are dicts too?
        molecule_dict['atoms'] = [atom.to_dict() for atom in self._atoms]
        molecule_dict['virtual_sites'] = [
            vsite.to_dict() for vsite in self._virtual_sites
        ]
        molecule_dict['bonds'] = [bond.to_dict() for bond in self._bonds]
        # TODO: Charges
        # TODO: Properties
        # From Jeff: We could have the onerous requirement that all "properties" have to_dict() functions.
        # Or we could restrict properties to simple stuff (ints, strings, floats, and the like)
        # Or pickle anything unusual
        # Or not allow user-defined properties at all (just use our internal _cached_properties)
        #molecule_dict['properties'] = dict([(key, value._to_dict()) for key.value in self._properties])
        # TODO: Assuming "simple stuff" properties right now, figure out a better standard
        molecule_dict['properties'] = self._properties
        if hasattr(self, '_cached_properties'):
            molecule_dict['cached_properties'] = self._cached_properties
        # TODO: Conformers
        if self._conformers is None:
            molecule_dict['conformers'] = None
        else:
            molecule_dict['conformers'] = []
            molecule_dict[
                'conformers_unit'] = 'angstrom'  # Have this defined as a class variable?
            for conf in self._conformers:
                conf_unitless = (conf / unit.angstrom)
                conf_serialized, conf_shape = serialize_numpy((conf_unitless))
                molecule_dict['conformers'].append(conf_serialized)
        if self._partial_charges is None:
            molecule_dict['partial_charges'] = None
            molecule_dict['partial_charges_unit'] = None

        else:
            charges_unitless = self._partial_charges / unit.elementary_charge
            charges_serialized, charges_shape = serialize_numpy(
                charges_unitless)
            molecule_dict['partial_charges'] = charges_serialized
            molecule_dict['partial_charges_unit'] = 'elementary_charge'

        return molecule_dict

    def __hash__(self):
        """
        Returns a hash of this molecule. Used when checking molecule uniqueness in Topology creation.

        Returns
        -------
        string
        """
        return hash(self.to_smiles())

    @classmethod
    def from_dict(cls, molecule_dict):
        """
        Create a new Molecule from a dictionary representation

        Parameters
        ----------
        molecule_dict : OrderedDict
            A dictionary representation of the molecule.

        Returns
        -------
        molecule : Molecule
            A Molecule created from the dictionary representation

        """
        # This implementation is a compromise to let this remain as a classmethod
        mol = cls()
        mol._initialize_from_dict(molecule_dict)
        return mol

    def _initialize_from_dict(self, molecule_dict):
        """
        Initialize this Molecule from a dictionary representation

        Parameters
        ----------
        molecule_dict : OrderedDict
            A dictionary representation of the molecule.
        """
        # TODO: Provide useful exception messages if there are any failures

        self._initialize()
        self.name = molecule_dict['name']
        for atom_dict in molecule_dict['atoms']:
            self._add_atom(**atom_dict)

        # Handle virtual site unit reattachment and molecule tagging
        for vsite_dict in molecule_dict['virtual_sites']:
            vsite_dict_units = deepcopy(vsite_dict)

            # Attach units to epsilon term
            vsite_dict_units['epsilon'] = string_to_quantity(
                vsite_dict['epsilon'])
            vsite_dict_units['sigma'] = string_to_quantity(
                vsite_dict['sigma'])
            vsite_dict_units['charge_increments'] = string_to_quantity(
                vsite_dict['charge_increments'])

            # Call the correct molecule._add_X_virtual_site function, based on the stated type
            if vsite_dict_units['vsite_type'] == "BondChargeVirtualSite":
                del vsite_dict_units['vsite_type']
                vsite_dict_units['distance'] = string_to_quantity(
                    vsite_dict['distance'])
                self._add_bond_charge_virtual_site(**vsite_dict_units)

            elif vsite_dict_units[
                    'vsite_type'] == "MonovalentLonePairVirtualSite":
                del vsite_dict_units['vsite_type']
                vsite_dict_units['distance'] = string_to_quantity(
                    vsite_dict['distance'])
                vsite_dict_units['in_plane_angle'] = string_to_quantity(
                    vsite_dict['in_plane_angle'])
                vsite_dict_units['out_of_plane_angle'] = string_to_quantity(
                    vsite_dict['out_of_plane_angle'])
                self._add_monovalent_lone_pair_virtual_site(**vsite_dict_units)

            elif vsite_dict_units[
                    'vsite_type'] == "DivalentLonePairVirtualSite":
                del vsite_dict_units['vsite_type']
                vsite_dict_units['distance'] = string_to_quantity(
                    vsite_dict['distance'])
                vsite_dict_units['in_plane_angle'] = string_to_quantity(
                    vsite_dict['in_plane_angle'])
                vsite_dict_units['out_of_plane_angle'] = string_to_quantity(
                    vsite_dict['out_of_plane_angle'])
                self._add_divalent_lone_pair_virtual_site(**vsite_dict_units)

            elif vsite_dict_units[
                    'vsite_type'] == "TrivalentLonePairVirtualSite":
                del vsite_dict_units['vsite_type']
                vsite_dict_units['distance'] = string_to_quantity(
                    vsite_dict['distance'])
                vsite_dict_units['in_plane_angle'] = string_to_quantity(
                    vsite_dict['in_plane_angle'])
                vsite_dict_units['out_of_plane_angle'] = string_to_quantity(
                    vsite_dict['out_of_plane_angle'])
                self._add_trivalent_lone_pair_virtual_site(**vsite_dict_units)

            else:
                raise Exception("Vsite type {} not recognized".format(
                    vsite_dict['vsite_type']))

        for bond_dict in molecule_dict['bonds']:
            bond_dict['atom1'] = int(bond_dict['atom1'])
            bond_dict['atom2'] = int(bond_dict['atom2'])
            self._add_bond(**bond_dict)

        if molecule_dict['partial_charges'] is None:
            self._partial_charges = None
        else:
            charges_shape = (self.n_atoms, )
            partial_charges_unitless = deserialize_numpy(
                molecule_dict['partial_charges'], charges_shape)
            pc_unit = getattr(unit, molecule_dict['partial_charges_unit'])
            partial_charges = unit.Quantity(partial_charges_unitless, pc_unit)
            self._partial_charges = partial_charges

        if molecule_dict['conformers'] is None:
            self._conformers = None
        else:
            self._conformers = list()
            for ser_conf in molecule_dict['conformers']:
                # TODO: Update to use string_to_quantity
                conformers_shape = ((self.n_atoms, 3))
                conformer_unitless = deserialize_numpy(ser_conf,
                                                       conformers_shape)
                c_unit = getattr(unit, molecule_dict['conformers_unit'])
                conformer = unit.Quantity(conformer_unitless, c_unit)
                self._conformers.append(conformer)

        self._properties = molecule_dict['properties']

    def __repr__(self):
        """Return the SMILES of this molecule"""
        return "Molecule with name '{}' and SMILES '{}'".format(
            self.name, self.to_smiles())

    def __getstate__(self):
        return self.to_dict()

    def __setstate__(self, state):
        return self._initialize_from_dict(state)

    def _initialize(self):
        """
        Clear the contents of the current molecule.
        """
        self._name = ''
        self._atoms = list()
        self._virtual_sites = list()
        self._bonds = list()  # List of bonds between Atom objects
        self._properties = {}  # Attached properties to be preserved
        #self._cached_properties = None # Cached properties (such as partial charges) can be recomputed as needed
        self._partial_charges = None
        self._conformers = None  # Optional conformers

    def _copy_initializer(self, other):
        """
        Copy contents of the specified molecule

        .. todo :: Should this be a ``@staticmethod`` where we have an explicit copy constructor?

        Parameters
        ----------
        other : optional
            Overwrite the state of this FrozenMolecule with the specified FrozenMolecule object.
            A deep copy is made.

        """
        #assert isinstance(other, type(self)), "can only copy instances of {}".format(type(self))
        other_dict = other.to_dict()
        self._initialize_from_dict(other_dict)

    def __eq__(self, other):
        """Test two molecules for equality to see if they are the chemical species, but do not check other annotated properties.

        .. note ::

           Note that this method simply tests whether two molecules are identical chemical species using equivalence of their canonical isomeric SMILES.
           No effort is made to ensure that the atoms are in the same order or that any annotated properties are preserved.

        """
        return self.is_isomorphic(other)

    def to_smiles(self, toolkit_registry=GLOBAL_TOOLKIT_REGISTRY):
        """
        Return a canonical isomeric SMILES representation of the current molecule

        .. note :: RDKit and OpenEye versions will not necessarily return the same representation.

        Parameters
        ----------
        toolkit_registry : openforcefield.utils.toolkits.ToolRegistry or openforcefield.utils.toolkits.ToolkitWrapper, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for SMILES conversion

        Returns
        -------
        smiles : str
            Canonical isomeric explicit-hydrogen SMILES

        Examples
        --------

        >>> from openforcefield.utils import get_data_file_path
        >>> sdf_filepath = get_data_file_path('molecules/ethanol.sdf')
        >>> molecule = Molecule(sdf_filepath)
        >>> smiles = molecule.to_smiles()

        """
        smiles = self._cached_smiles

        if smiles is not None:
            return smiles

        if isinstance(toolkit_registry, ToolkitRegistry):
            smiles = toolkit_registry.call('to_smiles', self)
        elif isinstance(toolkit_registry, ToolkitWrapper):
            toolkit = toolkit_registry
            smiles = toolkit.to_smiles(self)
        else:
            raise Exception(
                'Invalid toolkit_registry passed to to_smiles. Expected ToolkitRegistry or ToolkitWrapper. Got  {}'
                .format(type(toolkit_registry)))

        self._cached_smiles = smiles
        return smiles

    @staticmethod
    def from_smiles(smiles, hydrogens_are_explicit=False, toolkit_registry=GLOBAL_TOOLKIT_REGISTRY):
        """
        Construct a Molecule from a SMILES representation

        Parameters
        ----------
        smiles : str
            The SMILES representation of the molecule.
        hydrogens_are_explicit : bool, default = False
            If False, the cheminformatics toolkit will perform hydrogen addition
        toolkit_registry : openforcefield.utils.toolkits.ToolRegistry or openforcefield.utils.toolkits.ToolkitWrapper, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for SMILES-to-molecule conversion

        Returns
        -------
        molecule : openforcefield.topology.Molecule

        Examples
        --------

        >>> molecule = Molecule.from_smiles('Cc1ccccc1')

        """
        if isinstance(toolkit_registry, ToolkitRegistry):
            return toolkit_registry.call('from_smiles', smiles)
        elif isinstance(toolkit_registry, ToolkitWrapper):
            toolkit = toolkit_registry
            return toolkit.from_smiles(smiles, hydrogens_are_explicit=hydrogens_are_explicit)
        else:
            raise Exception(
                'Invalid toolkit_registry passed to from_smiles. Expected ToolkitRegistry or ToolkitWrapper. Got  {}'
                .format(type(toolkit_registry)))

    def is_isomorphic(
            self, other,
            compare_atom_stereochemistry=True,
            compare_bond_stereochemistry=True,
    ):
        """
        Determines whether the molecules are isomorphic by comparing their graphs.

        Parameters
        ----------
        other : an openforcefield.topology.molecule.FrozenMolecule
            The molecule to test for isomorphism.
        compare_atom_stereochemistry : bool, optional
            If ``False``, atoms' stereochemistry is ignored for the
            purpose of determining equality. Default is ``True``.
        compare_bond_stereochemistry : bool, optional
            If ``False``, bonds' stereochemistry is ignored for the
            purpose of determining equality. Default is ``True``.

        Returns
        -------
        molecules_are_isomorphic : bool
        """
        import networkx as nx

        def node_match_func(x, y):
            is_equal = (
                (x['atomic_number'] == y['atomic_number']) and
                (x['is_aromatic'] == y['is_aromatic']) and
                (x['formal_charge'] == y['formal_charge'])
            )
            if compare_atom_stereochemistry:
                is_equal &= x['stereochemistry'] == y['stereochemistry']
            return is_equal

        def edge_match_func(x, y):
            # We don't need to check the exact bond order (which is 1 or 2)
            # if the bond is aromatic. This way we avoid missing a match only
            # if the alternate bond orders 1 and 2 are assigned differently.
            is_equal = x['is_aromatic'] == y['is_aromatic'] or x['bond_order'] == y['bond_order']
            if compare_bond_stereochemistry:
                is_equal &= x['stereochemistry'] == y['stereochemistry']
            return is_equal

        return nx.is_isomorphic(self.to_networkx(),
                                other.to_networkx(),
                                node_match=node_match_func,
                                edge_match=edge_match_func
                                )
        #if not (isinstance(other, FrozenMolecule)):
        #    other_fm = FrozenMolecule(other)
        #else:
        #    other_fm = other
        #self_smiles = self.to_smiles(toolkit_registry=toolkit_registry)
        #other_smiles = other_fm.to_smiles(toolkit_registry=toolkit_registry)
        #return self_smiles == other_smiles

    def generate_conformers(self,
                            toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
                            n_conformers=10,
                            clear_existing=True):
        """
        Generate conformers for this molecule using an underlying toolkit

        Parameters
        ----------
        toolkit_registry : openforcefield.utils.toolkits.ToolRegistry or openforcefield.utils.toolkits.ToolkitWrapper, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for SMILES-to-molecule conversion
        n_conformers : int, default=1
            The maximum number of conformers to produce
        clear_existing : bool, default=True
            Whether to overwrite existing conformers for the molecule

        Examples
        --------

        >>> molecule = Molecule.from_smiles('CCCCCC')
        >>> molecule.generate_conformers()

        Raises
        ------
        InvalidToolkitError
            If an invalid object is passed as the toolkit_registry parameter

        """
        if isinstance(toolkit_registry, ToolkitRegistry):
            return toolkit_registry.call('generate_conformers', self, n_conformers=n_conformers,
                                         clear_existing=clear_existing)
        elif isinstance(toolkit_registry, ToolkitWrapper):
            toolkit = toolkit_registry
            return toolkit.generate_conformers(self, n_conformers=n_conformers, clear_existing=clear_existing)
        else:
            raise InvalidToolkitError(
                'Invalid toolkit_registry passed to generate_conformers. Expected ToolkitRegistry or ToolkitWrapper. Got  {}'
                .format(type(toolkit_registry)))

    def compute_partial_charges_am1bcc(self, toolkit_registry=GLOBAL_TOOLKIT_REGISTRY):
        """
        Calculate partial atomic charges for this molecule using AM1-BCC run by an underlying toolkit

        Parameters
        ----------
        toolkit_registry : openforcefield.utils.toolkits.ToolRegistry or openforcefield.utils.toolkits.ToolkitWrapper, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for the calculation

        Examples
        --------

        >>> molecule = Molecule.from_smiles('CCCCCC')
        >>> molecule.generate_conformers()
        >>> molecule.compute_partial_charges_am1bcc()

        Raises
        ------
        InvalidToolkitError
            If an invalid object is passed as the toolkit_registry parameter

        """
        if isinstance(toolkit_registry, ToolkitRegistry):
            charges = toolkit_registry.call(
                      'compute_partial_charges_am1bcc',
                      self
            )
        elif isinstance(toolkit_registry, ToolkitWrapper):
            toolkit = toolkit_registry
            charges = toolkit.compute_partial_charges_am1bcc(self)
        else:
            raise InvalidToolkitError(
                'Invalid toolkit_registry passed to compute_partial_charges_am1bcc. Expected ToolkitRegistry or ToolkitWrapper. Got  {}'
                .format(type(toolkit_registry)))
        self.partial_charges = charges


    def compute_partial_charges(self,
                                #quantum_chemical_method='AM1-BCC',
                                #partial_charge_method='None',
                                toolkit_registry=GLOBAL_TOOLKIT_REGISTRY):
        """
        **Warning! Not Implemented!**
        Calculate partial atomic charges for this molecule using an underlying toolkit

        Parameters
        ----------
        quantum_chemical_method : string, default='AM1-BCC'
            The quantum chemical method to use for partial charge calculation.
        partial_charge_method : string, default='None'
            The partial charge calculation method to use for partial charge calculation.
        toolkit_registry : openforcefield.utils.toolkits.ToolRegistry or openforcefield.utils.toolkits.ToolkitWrapper, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for SMILES-to-molecule conversion

        Examples
        --------

        >>> molecule = Molecule.from_smiles('CCCCCC')
        >>> molecule.generate_conformers()

        Raises
        ------
        InvalidToolkitError
            If an invalid object is passed as the toolkit_registry parameter

        """
        raise NotImplementedError
        # TODO: Implement this in a way that's compliant with SMIRNOFF's <ChargeIncrementModel> tag when the spec gets finalized
        # if isinstance(toolkit_registry, ToolkitRegistry):
        #     charges = toolkit_registry.call(
        #               'compute_partial_charges_am1bcc',
        #               self,
        #     )
        # elif isinstance(toolkit_registry, ToolkitWrapper):
        #     toolkit = toolkit_registry
        #     charges = toolkit.compute_partial_charges_am1bcc(
        #         self,
        #         #quantum_chemical_method=quantum_chemical_method,
        #         #partial_charge_method=partial_charge_method
        #     )
        # else:
        #     raise InvalidToolkitError(
        #         'Invalid toolkit_registry passed to compute_partial_charges_am1bcc. Expected ToolkitRegistry or ToolkitWrapper. Got  {}'
        #         .format(type(toolkit_registry)))

    def compute_wiberg_bond_orders(self,
                                   charge_model=None,
                                   toolkit_registry=GLOBAL_TOOLKIT_REGISTRY):
        """
        Calculate wiberg bond orders for this molecule using an underlying toolkit

        Parameters
        ----------
        toolkit_registry : openforcefield.utils.toolkits.ToolRegistry or openforcefield.utils.toolkits.ToolkitWrapper, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for SMILES-to-molecule conversion
        charge_model : string, optional
            The charge model to use for partial charge calculation

        Examples
        --------

        >>> molecule = Molecule.from_smiles('CCCCCC')
        >>> molecule.generate_conformers()
        >>> molecule.compute_wiberg_bond_orders()

        Raises
        ------
        InvalidToolkitError
            If an invalid object is passed as the toolkit_registry parameter

        """
        if isinstance(toolkit_registry, ToolkitRegistry):
            return toolkit_registry.call(
                'compute_wiberg_bond_orders', self, charge_model=charge_model)
        elif isinstance(toolkit_registry, ToolkitWrapper):
            toolkit = toolkit_registry
            return toolkit.compute_wiberg_bond_orders(
                self, charge_model=charge_model)
        else:
            raise Exception(
                'Invalid toolkit_registry passed to compute_wiberg_bond_orders. Expected ToolkitRegistry or ToolkitWrapper. Got  {}'
                .format(type(toolkit_registry)))

    def _invalidate_cached_properties(self):
        """
        Indicate that the chemical entity has been altered.
        """
        #if hasattr(self, '_cached_properties'):
        #    delattr(self, '_cached_properties')
        self._conformers = None
        self._partial_charges = None
        self._propers = None
        self._impropers = None

        self._cached_smiles = None
        # TODO: Clear fractional bond orders

    def to_networkx(self):
        """Generate a NetworkX undirected graph from the Molecule.

        Nodes are Atoms labeled with particle indices and atomic elements (via the ``element`` node atrribute).
        Edges denote chemical bonds between Atoms.
        Virtual sites are not included, since they lack a concept of chemical connectivity.

        .. todo ::

           * Do we need a ``from_networkx()`` method? If so, what would the Graph be required to provide?
           * Should edges be labeled with discrete bond types in some aromaticity model?
           * Should edges be labeled with fractional bond order if a method is specified?
           * Should we add other per-atom and per-bond properties (e.g. partial charges) if present?
           * Can this encode bond/atom chirality?


        Returns
        -------
        graph : networkx.Graph
            The resulting graph, with nodes (atoms) labeled with atom indices, elements, stereochemistry and aromaticity
            flags and bonds with two atom indices, bond order, stereochemistry, and aromaticity flags

        Examples
        --------
        Retrieve the bond graph for imatinib (OpenEye toolkit required)

        >>> molecule = Molecule.from_iupac('imatinib')
        >>> nxgraph = molecule.to_networkx()

        """
        import networkx as nx
        G = nx.Graph()
        for atom in self.atoms:
            G.add_node(
                atom.molecule_atom_index, atomic_number=atom.atomic_number, is_aromatic=atom.is_aromatic,
                stereochemistry=atom.stereochemistry, formal_charge=atom.formal_charge)
            #G.add_node(atom.molecule_atom_index, attr_dict={'atomic_number': atom.atomic_number})
        for bond in self.bonds:
            G.add_edge(
                bond.atom1_index, bond.atom2_index, bond_order=bond.bond_order, is_aromatic=bond.is_aromatic,
                stereochemistry = bond.stereochemistry)
            #G.add_edge(bond.atom1_index, bond.atom2_index, attr_dict={'order':bond.bond_order})

        return G

    def _add_atom(self,
                  atomic_number,
                  formal_charge,
                  is_aromatic,
                  stereochemistry=None,
                  name=None):
        """
        Add an atom

        Parameters
        ----------
        atomic_number : int
            Atomic number of the atom
        formal_charge : int
            Formal charge of the atom
        is_aromatic : bool
            If True, atom is aromatic; if False, not aromatic
        stereochemistry : str, optional, default=None
            Either 'R' or 'S' for specified stereochemistry, or None if stereochemistry is irrelevant
        name : str, optional, default=None
            An optional name for the atom

        Returns
        -------
        index : int
            The index of the atom in the molecule

        Examples
        --------

        Define a methane molecule

        >>> molecule = Molecule()
        >>> molecule.name = 'methane'
        >>> C = molecule.add_atom(6, 0, False)
        >>> H1 = molecule.add_atom(1, 0, False)
        >>> H2 = molecule.add_atom(1, 0, False)
        >>> H3 = molecule.add_atom(1, 0, False)
        >>> H4 = molecule.add_atom(1, 0, False)
        >>> bond_idx = molecule.add_bond(C, H1, False, 1)
        >>> bond_idx = molecule.add_bond(C, H2, False, 1)
        >>> bond_idx = molecule.add_bond(C, H3, False, 1)
        >>> bond_idx = molecule.add_bond(C, H4, False, 1)

        """
        # Create an atom
        atom = Atom(
            atomic_number,
            formal_charge,
            is_aromatic,
            stereochemistry=stereochemistry,
            name=name,
            molecule=self)
        self._atoms.append(atom)
        #self._particles.append(atom)
        self._invalidate_cached_properties()
        return self._atoms.index(atom)

    def _add_bond_charge_virtual_site(self, atoms, distance, **kwargs):
        """
        Create a bond charge-type virtual site, in which the location of the charge is specified by the position of two
        atoms. This supports placement of a virtual site S along a vector between two specified atoms, e.g. to allow
        for a sigma hole for halogens or similar contexts. With positive values of the distance, the virtual site lies
        outside the first indexed atom.

        Parameters
        ----------
        atoms : list of openforcefield.topology.molecule.Atom objects of shape [N]
            The atoms defining the virtual site's position
        distance : float

        weights : list of floats of shape [N] or None, optional, default=None
            weights[index] is the weight of particles[index] contributing to the position of the virtual site. Default
            is None
        charge_increments : list of floats of shape [N], optional, default=None
            The amount of charge to remove from the VirtualSite's atoms and put in the VirtualSite. Indexing in this
            list should match the ordering in the atoms list. Default is None.
        epsilon : float
            Epsilon term for VdW properties of virtual site. Default is None.
        sigma : float, default=None
            Sigma term for VdW properties of virtual site. Default is None.
        rmin_half : float
            Rmin_half term for VdW properties of virtual site. Default is None.
        name : string or None, default=None
            The name of this virtual site. Default is None.

        Returns
        -------
        index : int
            The index of the newly-added virtual site in the molecule
        """
        # Check if function was passed list of atoms or atom indices
        if all([isinstance(atom, int) for atom in atoms]):
            atom_list = [self.atoms[atom_index] for atom_index in atoms]
        elif all([isinstance(atom, Atom) for atom in atoms]):
            atom_list = atoms
        else:
            raise Exception(
                'Invalid inputs to molecule._add_bond_charge_virtual_site.'
                ' Expected ints or Atoms. Received types {} '.format(
                    [type(i) for i in atoms]))
        # TODO: Check to make sure bond does not already exist
        vsite = BondChargeVirtualSite(atom_list, distance, **kwargs)
        self._virtual_sites.append(vsite)
        self._invalidate_cached_properties()
        return self._virtual_sites.index(vsite)

    def _add_monovalent_lone_pair_virtual_site(self, atoms, distance,
                                               out_of_plane_angle,
                                               in_plane_angle, **kwargs):
        """
        Create a bond charge-type virtual site, in which the location of the charge is specified by the position of
        three atoms.

        TODO : Do "weights" have any meaning here?

        Parameters
        ----------
        atoms : list of three openforcefield.topology.molecule.Atom objects
            The three atoms defining the virtual site's position
        distance : float

        out_of_plane_angle : float

        in_plane_angle : float

        epsilon : float
            Epsilon term for VdW properties of virtual site. Default is None.
        sigma : float, default=None
            Sigma term for VdW properties of virtual site. Default is None.
        rmin_half : float
            Rmin_half term for VdW properties of virtual site. Default is None.
        name : string or None, default=None
            The name of this virtual site. Default is None.

        Returns
        -------
        index : int
            The index of the newly-added virtual site in the molecule
        """
        # Check if function was passed list of atoms or atom indices
        if all([isinstance(atom, int) for atom in atoms]):
            atom_list = [self.atoms[atom_index] for atom_index in atoms]
        elif all([isinstance(atom, Atom) for atom in atoms]):
            atom_list = atoms
        else:
            raise Exception(
                'Invalid inputs to molecule._add_monovalent_lone_pair_virtual_site. Expected ints or Atoms.'
                ' Received types {} '.format([type(i) for i in atoms]))
        # TODO: Check to make sure bond does not already exist
        vsite = MonovalentLonePairVirtualSite(
            atom_list, distance, out_of_plane_angle, in_plane_angle, **kwargs)
        self._virtual_sites.append(vsite)
        self._invalidate_cached_properties()
        return self._virtual_sites.index(vsite)

    def _add_divalent_lone_pair_virtual_site(self, atoms, distance,
                                             out_of_plane_angle,
                                             in_plane_angle, **kwargs):
        """
        Create a divalent lone pair-type virtual site, in which the location of the charge is specified by the position
        of three atoms.

        TODO : Do "weights" have any meaning here?

        Parameters
        ----------
        atoms : list of 3 openforcefield.topology.molecule.Atom objects
            The three atoms defining the virtual site's position
        distance : float

        out_of_plane_angle : float

        in_plane_angle : float

        epsilon : float
            Epsilon term for VdW properties of virtual site. Default is None.
        sigma : float, default=None
            Sigma term for VdW properties of virtual site. Default is None.
        rmin_half : float
            Rmin_half term for VdW properties of virtual site. Default is None.
        name : string or None, default=None
            The name of this virtual site. Default is None.

        Returns
        -------
        index : int
            The index of the newly-added virtual site in the molecule
        """
        # Check if function was passed list of atoms or atom indices
        if all([isinstance(atom, int) for atom in atoms]):
            atom_list = [self.atoms[atom_index] for atom_index in atoms]
        elif all([isinstance(atom, Atom) for atom in atoms]):
            atom_list = atoms
        else:
            raise Exception(
                'Invalid inputs to molecule._add_divalent_lone_pair_virtual_site. Expected ints or Atoms. '
                'Received types {} '.format([type(i) for i in atoms]))
        # TODO: Check to make sure bond does not already exist
        vsite = DivalentLonePairVirtualSite(
            atom_list, distance, out_of_plane_angle, in_plane_angle, **kwargs)
        self._virtual_sites.append(vsite)
        self._invalidate_cached_properties()
        return self._virtual_sites.index(vsite)

    def _add_trivalent_lone_pair_virtual_site(self, atoms, distance,
                                              out_of_plane_angle,
                                              in_plane_angle, **kwargs):
        """
        Create a trivalent lone pair-type virtual site, in which the location of the charge is specified by the position
         of four atoms.

        TODO : Do "weights" have any meaning here?

        Parameters
        ----------
        atoms : list of 4 openforcefield.topology.molecule.Atom objects or atom indices
            The three atoms defining the virtual site's position
        distance : float

        out_of_plane_angle : float

        in_plane_angle : float

        epsilon : float
            Epsilon term for VdW properties of virtual site. Default is None.
        sigma : float, default=None
            Sigma term for VdW properties of virtual site. Default is None.
        rmin_half : float
            Rmin_half term for VdW properties of virtual site. Default is None.
        name : string or None, default=None
            The name of this virtual site. Default is None.
        """
        # Check if function was passed list of atoms or atom indices
        if all([isinstance(atom, int) for atom in atoms]):
            atom_list = [self.atoms[atom_index] for atom_index in atoms]
        elif all([isinstance(atom, Atom) for atom in atoms]):
            atom_list = atoms
        else:
            raise Exception(
                'Invalid inputs to molecule._add_trivalent_lone_pair_virtual_site. Expected ints or Atoms. Received types {} '
                .format([type(i) for i in atoms]))
        vsite = TrivalentLonePairVirtualSite(
            atom_list, distance, out_of_plane_angle, in_plane_angle, **kwargs)
        self._virtual_sites.append(vsite)
        self._invalidate_cached_properties()
        return self._virtual_sites.index(vsite)

    def _add_bond(self,
                  atom1,
                  atom2,
                  bond_order,
                  is_aromatic,
                  stereochemistry=None,
                  fractional_bond_order=None):
        """
        Add a bond between two specified atom indices

        Parameters
        ----------
        atom1 : int or openforcefield.topology.molecule.Atom
            Index of first atom or first atom
        atom2_index : int or openforcefield.topology.molecule.Atom
            Index of second atom or second atom
        bond_order : int
            Integral bond order of Kekulized form
        is_aromatic : bool
            True if this bond is aromatic, False otherwise
        stereochemistry : str, optional, default=None
            Either 'E' or 'Z' for specified stereochemistry, or None if stereochemistry is irrelevant
        fractional_bond_order : float, optional, default=None
            The fractional (eg. Wiberg) bond order
        Returns
        -------
        index : int
            The index of the bond in the molecule

        """
        if isinstance(atom1, int) and isinstance(atom2, int):
            atom1_atom = self.atoms[atom1]
            atom2_atom = self.atoms[atom2]
        elif isinstance(atom1, Atom) and isinstance(atom2, Atom):
            atom1_atom = atom1
            atom2_atom = atom2
        else:
            raise Exception(
                'Invalid inputs to molecule._add_bond. Expected ints or Atoms. '
                'Received {} (type {}) and {} (type {}) '.format(
                    atom1, type(atom1), atom2, type(atom2)))
        # TODO: Check to make sure bond does not already exist
        if atom1_atom.is_bonded_to(atom2_atom):
            raise Exception('Bond already exists between {} and {}'.format(
                atom1_atom, atom2_atom))
        bond = Bond(
            atom1_atom,
            atom2_atom,
            bond_order,
            is_aromatic,
            stereochemistry=stereochemistry,
            fractional_bond_order=fractional_bond_order)
        self._bonds.append(bond)
        self._invalidate_cached_properties()
        # TODO: This is a bad way to get bond index
        return self._bonds.index(bond)

    def _add_conformer(self, coordinates):
        """
        Add a conformation of the molecule

        Parameters
        ----------
        coordinates: A simtk vector wrapped unit quantity
            The coordinates of the conformer to add.

        Returns
        -------
        index: int
            The index of this conformer


        """
        new_conf = unit.Quantity(
            np.zeros((self.n_atoms, 3), np.float), unit.angstrom)
        if not (new_conf.shape == coordinates.shape):
            raise Exception(
                "molecule.add_conformer given input of the wrong shape: "
                "Given {}, expected {}".format(coordinates.shape,
                                               new_conf.shape))

        try:
            new_conf[:] = coordinates
        except AttributeError as e:
            print(e)
            raise Exception(
                'Coordinates passed to Molecule._add_conformer without units. Ensure that coordinates are '
                'of type simtk.units.Quantity')

        if self._conformers is None:
            self._conformers = []
        self._conformers.append(new_conf)
        return len(self._conformers)

    @property
    def partial_charges(self):
        """
        Returns the partial charges (if present) on the molecule

        Returns
        -------
        partial_charges : a simtk.unit.Quantity - wrapped numpy array [1 x n_atoms]
            The partial charges on this Molecule's atoms.
        """
        return self._partial_charges

    @partial_charges.setter
    def partial_charges(self, charges):
        """
        Set the atomic partial charges for this molecule

        Parameters
        ----------
        charges : a simtk.unit.Quantity - wrapped numpy array [1 x n_atoms]
            The partial charges to assign to the molecule. Must be in units compatible with simtk.unit.elementary_charge

        """
        assert hasattr(charges, 'unit')
        assert unit.elementary_charge.is_compatible(charges.unit)
        assert charges.shape == (self.n_atoms, )

        charges_ec = charges.in_units_of(unit.elementary_charge)
        self._partial_charges = charges_ec

    @property
    def n_particles(self):
        """
        The number of Particle objects, which corresponds to how many positions must be used.
        """
        return len(self._atoms) + len(self._virtual_sites)

    @property
    def n_atoms(self):
        """
        The number of Atom objects.
        """
        return len(self._atoms)

    @property
    def n_virtual_sites(self):
        """
        The number of VirtualSite objects.
        """
        return sum([1 for virtual_site in self.virtual_sites])

    @property
    def n_bonds(self):
        """
        The number of Bond objects.
        """
        return sum([1 for bond in self.bonds])

    @property
    def n_angles(self):
        """int: number of angles in the Molecule."""
        self._construct_angles()
        return len(self._angles)

    @property
    def n_propers(self):
        """int: number of proper torsions in the Molecule."""
        self._construct_torsions()
        return len(self._propers)

    @property
    def n_impropers(self):
        """int: number of improper torsions in the Molecule."""
        self._construct_torsions()
        return len(self._impropers)

    @property
    def particles(self):
        """
        Iterate over all Particle objects.
        """
        # TODO: Re-implement this when we see how it interfaces with Topology
        return self._atoms + self._virtual_sites

    @property
    def atoms(self):
        """
        Iterate over all Atom objects.
        """
        return self._atoms

    @property
    def conformers(self):
        """
        Iterate over all conformers in this molecule.
        """
        return self._conformers

    @property
    def n_conformers(self):
        """
        Iterate over all Atom objects.
        """
        if self._conformers is None:
            return 0
        return len(self._conformers)

    @property
    def virtual_sites(self):
        """
        Iterate over all VirtualSite objects.
        """
        return self._virtual_sites

    @property
    def bonds(self):
        """
        Iterate over all Bond objects.
        """
        return self._bonds

    @property
    def angles(self):
        """
        Get an iterator over all i-j-k angles.
        """
        self._construct_angles()
        return self._angles

    @property
    def torsions(self):
        """
        Get an iterator over all i-j-k-l torsions.
        Note that i-j-k-i torsions (cycles) are excluded.

        Returns
        -------
        torsions : iterable of 4-Atom tuples
        """
        self._construct_torsions()
        return self._torsions

    @property
    def propers(self):
        """
        Iterate over all proper torsions in the molecule

        .. todo::

           * Do we need to return a ``Torsion`` object that collects information about fractional bond orders?
        """
        self._construct_torsions()
        return self._propers

    @property
    def impropers(self):
        """
        Iterate over all proper torsions in the molecule

        .. todo::

           * Do we need to return a ``Torsion`` object that collects information about fractional bond orders?
        """
        self._construct_torsions()
        return self._impropers

    @property
    def total_charge(self):
        """
        Return the total charge on the molecule
        """
        return sum([atom.formal_charge for atom in self.atoms])

    @property
    def name(self):
        """
        The name (or title) of the molecule
        """
        return self._name

    @name.setter
    def name(self, other):
        """
        Set the name of this molecule
        """
        if other is None:
            self._name = ''
        elif type(other) is str:
            self._name = other
        else:
            raise Exception("Molecule name must be a string")

    @property
    def properties(self):
        """
        The properties dictionary of the molecule
        """
        return self._properties

    def chemical_environment_matches(self,
                                     query,
                                     toolkit_registry=GLOBAL_TOOLKIT_REGISTRY):
        """Retrieve all matches for a given chemical environment query.

        Parameters
        ----------
        query : str or ChemicalEnvironment
            SMARTS string (with one or more tagged atoms) or ``ChemicalEnvironment`` query
            Query will internally be resolved to SMIRKS using ``query.asSMIRKS()`` if it has an ``.asSMIRKS`` method.
        toolkit_registry : openforcefield.utils.toolkits.ToolRegistry or openforcefield.utils.toolkits.ToolkitWrapper, optional, default=GLOBAL_TOOLKIT_REGISTRY
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for chemical environment matches


        Returns
        -------
        matches : list of Atom tuples
            A list of all matching Atom tuples

        Examples
        --------
        Retrieve all the carbon-carbon bond matches in a molecule

        >>> molecule = Molecule.from_iupac('imatinib')
        >>> matches = molecule.chemical_environment_matches('[#6X3:1]~[#6X3:2]')

        .. todo ::

           * Do we want to generalize ``query`` to allow other kinds of queries, such as mdtraj DSL, pymol selections, atom index slices, etc?
             We could call it ``topology.matches(query)`` instead of ``chemical_environment_matches``

        """
        # Resolve to SMIRKS if needed
        # TODO: Update this to use updated ChemicalEnvironment API
        if hasattr(query, 'asSMIRKS'):
            smirks = query.asSMIRKS()
        elif type(query) == str:
            smirks = query
        else:
            raise ValueError(
                "'query' must be either a string or a ChemicalEnvironment")

        # Use specified cheminformatics toolkit to determine matches with specified aromaticity model
        # TODO: Simplify this by requiring a toolkit registry for the molecule?
        # TODO: Do we have to pass along an aromaticity model?
        if isinstance(toolkit_registry, ToolkitRegistry):
            matches = toolkit_registry.call('find_smarts_matches', self,
                                            smirks)
        elif isinstance(toolkit_registry, ToolkitWrapper):
            matches = toolkit_registry.find_smarts_matches(self, smirks)
        else:
            raise ValueError(
                "'toolkit_registry' must be either a ToolkitRegistry or a ToolkitWrapper"
            )

        return matches

    # TODO: Move OE-dependent parts of this to toolkits.py
    @classmethod
    @OpenEyeToolkitWrapper.requires_toolkit()
    def from_iupac(cls, iupac_name, **kwargs):
        """Generate a molecule from IUPAC or common name

        Parameters
        ----------
        iupac_name : str
            IUPAC name of molecule to be generated
        allow_undefined_stereo : bool, default=False
            If false, raises an exception if molecule contains undefined stereochemistry.

        Returns
        -------
        molecule : Molecule
            The resulting molecule with position

        .. note :: This method requires the OpenEye toolkit to be installed.

        Examples
        --------

        Create a molecule from a common name

        >>> molecule = Molecule.from_iupac('4-[(4-methylpiperazin-1-yl)methyl]-N-(4-methyl-3-{[4-(pyridin-3-yl)pyrimidin-2-yl]amino}phenyl)benzamide')

        Create a molecule from a common name

        >>> molecule = Molecule.from_iupac('imatinib')

        """
        from openeye import oechem, oeiupac
        oemol = oechem.OEMol()
        oeiupac.OEParseIUPACName(oemol, iupac_name)
        oechem.OETriposAtomNames(oemol)
        result = oechem.OEAddExplicitHydrogens(oemol)
        if result == False:
            raise Exception(
                "Addition of explicit hydrogens failed in from_iupac")
        return cls.from_openeye(oemol, **kwargs)

    # TODO: Move OE-dependent parts of this to toolkits.py
    @OpenEyeToolkitWrapper.requires_toolkit()
    def to_iupac(self):
        """Generate IUPAC name from Molecule

        Returns
        -------
        iupac_name : str
            IUPAC name of the molecule

        .. note :: This method requires the OpenEye toolkit to be installed.

        Examples
        --------

        >>> from openforcefield.utils import get_data_file_path
        >>> sdf_filepath = get_data_file_path('molecules/ethanol.sdf')
        >>> molecule = Molecule(sdf_filepath)
        >>> iupac_name = molecule.to_iupac()

        """
        from openeye import oeiupac
        return oeiupac.OECreateIUPACName(self.to_openeye())

    @staticmethod
    def from_topology(topology):
        """Return a Molecule representation of an openforcefield Topology containing a single Molecule object.

        Parameters
        ----------
        topology : openforcefield.topology.Topology
            The :class:`Topology` object containing a single :class:`Molecule` object.
            Note that OpenMM and MDTraj ``Topology`` objects are not supported.

        Returns
        -------
        molecule : openforcefield.topology.Molecule
            The Molecule object in the topology

        Raises
        ------
        ValueError
            If the topology does not contain exactly one molecule.

        Examples
        --------

        Create a molecule from a Topology object that contains exactly one molecule

        >>> molecule = Molecule.from_topology(topology)  # doctest: +SKIP

        """
        # TODO: Ensure we are dealing with an openforcefield Topology object
        if topology.n_topology_molecules != 1:
            raise ValueError('Topology must contain exactly one molecule')
        molecule = [i for i in topology.reference_molecules][0]
        return Molecule(molecule)

    def to_topology(self):
        """
        Return an openforcefield Topology representation containing one copy of this molecule

        Returns
        -------
        topology : openforcefield.topology.Topology
            A Topology representation of this molecule

        Examples
        --------

        >>> molecule = Molecule.from_iupac('imatinib')
        >>> topology = molecule.to_topology()

        """
        from openforcefield.topology import Topology
        return Topology.from_molecules(self)

    @staticmethod
    def from_file(file_path,
                  file_format=None,
                  toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
                  allow_undefined_stereo=False):
        """
        Create one or more molecules from a file

        .. todo::

           * Extend this to also include some form of .offmol Open Force Field Molecule format?
           * Generalize this to also include file-like objects?

        Parameters
        ----------
        file_path : str or file-like object
            The path to the file or file-like object to stream one or more molecules from.
        file_format : str, optional, default=None
            Format specifier, usually file suffix (eg. 'MOL2', 'SMI')
            Note that not all toolkits support all formats. Check ToolkitWrapper.toolkit_file_read_formats for your
            loaded toolkits for details.
        toolkit_registry : openforcefield.utils.toolkits.ToolRegistry or openforcefield.utils.toolkits.ToolkitWrapper,
        optional, default=GLOBAL_TOOLKIT_REGISTRY
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for file loading. If a Toolkit is passed, only
            the highest-precedence toolkit is used
        allow_undefined_stereo : bool, default=False
            If false, raises an exception if oemol contains undefined stereochemistry.

        Returns
        -------
        molecules : Molecule or list of Molecules
            If there is a single molecule in the file, a Molecule is returned;
            otherwise, a list of Molecule objects is returned.

        Examples
        --------
        >>> from openforcefield.tests.utils import get_monomer_mol2_file_path
        >>> mol2_file_path = get_monomer_mol2_file_path('cyclohexane')
        >>> molecule = Molecule.from_file(mol2_file_path)

        """

        if file_format is None:
            if not (isinstance(file_path, str)):
                raise Exception(
                    "If providing a file-like object for reading molecules, the format must be specified"
                )
            # Assume that files ending in ".gz" should use their second-to-last suffix for compatibility check
            # TODO: Will all cheminformatics packages be OK with gzipped files?
            if file_path[-3:] == '.gz':
                file_format = file_path.split('.')[-2]
            else:
                file_format = file_path.split('.')[-1]
        file_format = file_format.upper()

        # Determine which toolkit to use (highest priority that's compatible with input type)
        if isinstance(toolkit_registry, ToolkitRegistry):
            # TODO: Encapsulate this logic into ToolkitRegistry.call()?
            toolkit = None
            supported_read_formats = {}
            for query_toolkit in toolkit_registry.registered_toolkits:
                if file_format in query_toolkit.toolkit_file_read_formats:
                    toolkit = query_toolkit
                    break
                supported_read_formats[
                    query_toolkit.
                    toolkit_name] = query_toolkit.toolkit_file_read_formats
            if toolkit is None:
                raise NotImplementedError(
                    "No toolkits in registry can read file {} (format {}). Supported formats in the "
                    "provided ToolkitRegistry are {}".format(
                        file_path, file_format, supported_read_formats))

        elif isinstance(toolkit_registry, ToolkitWrapper):
            # TODO: Encapsulate this logic in ToolkitWrapper?
            toolkit = toolkit_registry
            if file_format not in toolkit.toolkit_file_read_formats:
                raise NotImplementedError(
                    "Toolkit {} can not read file {} (format {}). Supported formats for this toolkit "
                    "are {}".format(toolkit.toolkit_name, file_path,
                                    file_format,
                                    toolkit.toolkit_file_read_formats))
        else:
            raise ValueError(
                "'toolkit_registry' must be either a ToolkitRegistry or a ToolkitWrapper"
            )

        mols = list()

        if isinstance(file_path, str):
            mols = toolkit.from_file(
                file_path,
                file_format=file_format,
                allow_undefined_stereo=allow_undefined_stereo)
        elif hasattr(file_path, 'read'):
            file_obj = file_path
            mols = toolkit.from_file_obj(
                file_obj,
                file_format=file_format,
                allow_undefined_stereo=allow_undefined_stereo)

        if len(mols) == 0:
            raise Exception(
                'Unable to read molecule from file: {}'.format(file_path))
        elif len(mols) == 1:
            return mols[0]
        return mols

    def to_file(self,
                file_path,
                file_format,
                toolkit_registry=GLOBAL_TOOLKIT_REGISTRY):
        """Write the current molecule to a file or file-like object

        Parameters
        ----------
        file_path : str or file-like object
            A file-like object or the path to the file to be written.
        file_format : str
            Format specifier, one of ['MOL2', 'MOL2H', 'SDF', 'PDB', 'SMI', 'CAN', 'TDT']
            Note that not all toolkits support all formats
        toolkit_registry : openforcefield.utils.toolkits.ToolRegistry or openforcefield.utils.toolkits.ToolkitWrapper,
        optional, default=GLOBAL_TOOLKIT_REGISTRY
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for file writing. If a Toolkit is passed, only
            the highest-precedence toolkit is used

        Raises
        ------
        ValueError
            If the requested file_format is not supported by one of the installed cheminformatics toolkits

        Examples
        --------

        >>> molecule = Molecule.from_iupac('imatinib')
        >>> molecule.to_file('imatinib.mol2', file_format='mol2')  # doctest: +SKIP
        >>> molecule.to_file('imatinib.sdf', file_format='sdf')  # doctest: +SKIP
        >>> molecule.to_file('imatinib.pdb', file_format='pdb')  # doctest: +SKIP

        """

        if isinstance(toolkit_registry, ToolkitRegistry):
            pass
        elif isinstance(toolkit_registry, ToolkitWrapper):
            toolkit = toolkit_registry
            toolkit_registry = ToolkitRegistry(toolkit_precedence=[])
            toolkit_registry.add_toolkit(toolkit)
        else:
            raise ValueError(
                "'toolkit_registry' must be either a ToolkitRegistry or a ToolkitWrapper"
            )

        file_format = file_format.upper()

        # Take the first toolkit that can write the desired output format
        toolkit = None
        for query_toolkit in toolkit_registry.registered_toolkits:
            if file_format in query_toolkit.toolkit_file_write_formats:
                toolkit = query_toolkit
                break

        # Raise an exception if no toolkit was found to provide the requested file_format
        if toolkit is None:
            supported_formats = {}
            for toolkit in toolkit_registry.registered_toolkits:
                supported_formats[
                    toolkit.toolkit_name] = toolkit.toolkit_file_write_formats
            raise ValueError(
                'The requested file format ({}) is not available from any of the installed toolkits '
                '(supported formats: {})'.format(file_format,
                                                 supported_formats))

        # Write file
        if type(file_path) == str:
            # Open file for writing
            toolkit.to_file(self, file_path, file_format)
        else:
            toolkit.to_file_obj(self, file_path, file_format)

    @staticmethod
    @RDKitToolkitWrapper.requires_toolkit()
    def from_rdkit(rdmol, allow_undefined_stereo=False):
        """
        Create a Molecule from an RDKit molecule.

        Requires the RDKit to be installed.

        Parameters
        ----------
        rdmol : rkit.RDMol
            An RDKit molecule
        allow_undefined_stereo : bool, default=False
            If false, raises an exception if oemol contains undefined stereochemistry.

        Returns
        -------
        molecule : openforcefield.Molecule
            An openforcefield molecule

        Examples
        --------

        Create a molecule from an RDKit molecule

        >>> from rdkit import Chem
        >>> from openforcefield.tests.utils import get_data_file_path
        >>> rdmol = Chem.MolFromMolFile(get_data_file_path('systems/monomers/ethanol.sdf'))
        >>> molecule = Molecule.from_rdkit(rdmol)

        """
        toolkit = RDKitToolkitWrapper()
        return toolkit.from_rdkit(
            rdmol, allow_undefined_stereo=allow_undefined_stereo)

    @RDKitToolkitWrapper.requires_toolkit()
    def to_rdkit(self, aromaticity_model=DEFAULT_AROMATICITY_MODEL):
        """
        Create an RDKit molecule

        Requires the RDKit to be installed.

        Parameters
        ----------
        aromaticity_model : str, optional, default=DEFAULT_AROMATICITY_MODEL
            The aromaticity model to use

        Returns
        -------
        rdmol : rkit.RDMol
            An RDKit molecule

        Examples
        --------

        Convert a molecule to RDKit

        >>> from openforcefield.utils import get_data_file_path
        >>> sdf_filepath = get_data_file_path('molecules/ethanol.sdf')
        >>> molecule = Molecule(sdf_filepath)
        >>> rdmol = molecule.to_rdkit()

        """
        toolkit = RDKitToolkitWrapper()
        return toolkit.to_rdkit(self, aromaticity_model=aromaticity_model)

    @staticmethod
    @OpenEyeToolkitWrapper.requires_toolkit()
    def from_openeye(oemol, allow_undefined_stereo=False):
        """
        Create a Molecule from an OpenEye molecule.

        Requires the OpenEye toolkit to be installed.

        Parameters
        ----------
        oemol : openeye.oechem.OEMol
            An OpenEye molecule
        allow_undefined_stereo : bool, default=False
            If false, raises an exception if oemol contains undefined stereochemistry.

        Returns
        -------
        molecule : openforcefield.topology.Molecule
            An openforcefield molecule

        Examples
        --------

        Create a Molecule from an OpenEye OEMol

        >>> from openeye import oechem
        >>> from openforcefield.tests.utils import get_data_file_path
        >>> ifs = oechem.oemolistream(get_data_file_path('systems/monomers/ethanol.mol2'))
        >>> oemols = list(ifs.GetOEGraphMols())
        >>> molecule = Molecule.from_openeye(oemols[0])

        """
        toolkit = OpenEyeToolkitWrapper()
        return toolkit.from_openeye(
            oemol, allow_undefined_stereo=allow_undefined_stereo)

    @OpenEyeToolkitWrapper.requires_toolkit()
    def to_openeye(self, aromaticity_model=DEFAULT_AROMATICITY_MODEL):
        """
        Create an OpenEye molecule

        Requires the OpenEye toolkit to be installed.

        .. todo ::

           * Use stored conformer positions instead of an argument.
           * Should the aromaticity model be specified in some other way?

        Parameters
        ----------
        aromaticity_model : str, optional, default=DEFAULT_AROMATICITY_MODEL
            The aromaticity model to use

        Returns
        -------
        oemol : openeye.oechem.OEMol
            An OpenEye molecule

        Examples
        --------

        Create an OpenEye molecule from a Molecule

        >>> molecule = Molecule.from_smiles('CC')
        >>> oemol = molecule.to_openeye()

        """
        toolkit = OpenEyeToolkitWrapper()
        return toolkit.to_openeye(self, aromaticity_model=aromaticity_model)


    def get_fractional_bond_orders(self,
                                   method='Wiberg',
                                   toolkit_registry=GLOBAL_TOOLKIT_REGISTRY):
        """Get fractional bond orders.

        method : str, optional, default='Wiberg'
            The name of the charge method to use.
            Options are:
            * 'Wiberg' : Wiberg bond order
        toolkit_registry : openforcefield.utils.toolkits ToolkitRegistry
            The toolkit registry to use for molecule operations

        Examples
        --------

        Get fractional Wiberg bond orders

        >>> molecule = Molecule.from_iupac('imatinib')
        >>> molecule.generate_conformers()
        >>> fractional_bond_orders = molecule.get_fractional_bond_orders(method='Wiberg')


        .. todo::
            * Is it OK that the ``Molecule`` object does not store geometry, but will create it using ``openeye.omega`` or ``rdkit``?
            * Should this method assign fractional bond orders to the ``Bond``s in the molecule, a separate ``bond_orders`` molecule property,
              or just return the array of bond orders?
            * How do we add enough flexibility to specify the toolkit and optional parameters, such as:
              ``oequacpac.OEAssignPartialCharges(charged_copy, getattr(oequacpac, 'OECharges_AM1BCCSym'), False, False)``
            * Generalize to allow user to specify both QM method and bond order computation approach (e.g. ``AM1`` and ``Wiberg``)


        """
        # TODO: Let ToolkitRegistry handle this once compute_fractional_bond_orders will be added to the Wrappers API.
        if method != 'Wiberg':
            raise NotImplementedError('Only Wiberg bond order is currently implemented')
        # TODO: Use memoization to speed up subsequent calls; use decorator?
        fractional_bond_orders = toolkit_registry.call(
            'compute_wiberg_bond_orders', molecule=self)
        return fractional_bond_orders

    def _construct_angles(self):
        """
        Get an iterator over all i-j-k angles.
        """
        # TODO: Build Angle objects instead of tuple of atoms.
        if not hasattr(self, '_angles'):
            self._construct_bonded_atoms_list()
            self._angles = set()
            for atom1 in self._atoms:
                for atom2 in self._bondedAtoms[atom1]:
                    for atom3 in self._bondedAtoms[atom2]:
                        if atom1 == atom3:
                            continue
                        # TODO: Encapsulate this logic into an Angle class.
                        if atom1.molecule_atom_index < atom3.molecule_atom_index:
                            self._angles.add((atom1, atom2, atom3))
                        else:
                            self._angles.add((atom3, atom2, atom1))

    def _construct_torsions(self):
        """
        Construct sets containing the atoms improper and proper torsions
        """
        # TODO: Build Proper/ImproperTorsion objects instead of tuple of atoms.
        if not hasattr(self, '_torsions'):
            self._construct_bonded_atoms_list()

            #self._torsions = set()
            self._propers = set()
            self._impropers = set()
            for atom1 in self._atoms:
                for atom2 in self._bondedAtoms[atom1]:
                    for atom3 in self._bondedAtoms[atom2]:
                        if atom1 == atom3:
                            continue
                        for atom4 in self._bondedAtoms[atom3]:
                            if atom4 == atom2:
                                continue
                            # Exclude i-j-k-i
                            if atom1 == atom4:
                                continue

                            if atom1.molecule_atom_index < atom4.molecule_atom_index:
                                torsion = (atom1, atom2, atom3, atom4)
                            else:
                                torsion = (atom4, atom3, atom2, atom1)

                            self._propers.add(torsion)

                        for atom3i in self._bondedAtoms[atom2]:
                            if atom3i == atom3:
                                continue
                            if atom3i == atom1:
                                continue

                            improper = (atom1, atom2, atom3, atom3i)
                            self._impropers.add(improper)

            self._torsions = self._propers | self._impropers
        #return iter(self._torsions)

    def _construct_bonded_atoms_list(self):
        """
        Construct list of all atoms each atom is bonded to.

        """
        # TODO: Add this to cached_properties
        if not hasattr(self, '_bondedAtoms'):
            #self._atoms = [ atom for atom in self.atoms() ]
            self._bondedAtoms = dict()
            for atom in self._atoms:
                self._bondedAtoms[atom] = set()
            for bond in self._bonds:
                atom1 = self.atoms[bond.atom1_index]
                atom2 = self.atoms[bond.atom2_index]
                self._bondedAtoms[atom1].add(atom2)
                self._bondedAtoms[atom2].add(atom1)

    def _is_bonded(self, atom_index_1, atom_index_2):
        """Return True if atoms are bonded, False if not.

        Parameters
        ----------
        atom_index_1 : int
        atom_index_2 : int
            Atom indices

        Returns
        -------
        is_bonded : bool
            True if atoms are bonded, False otherwise


        """
        self._construct_bonded_atoms_list()
        atom1 = self._atoms[atom_index_1]
        atom2 = self._atoms[atom_index_2]
        return atom2 in self._bondedAtoms[atom1]


class Molecule(FrozenMolecule):
    """
    Mutable chemical representation of a molecule, such as a small molecule or biopolymer.

    .. todo :: What other API calls would be useful for supporting biopolymers as small molecules? Perhaps iterating over chains and residues?

    Examples
    --------

    Create a molecule from an sdf file

    >>> from openforcefield.utils import get_data_file_path
    >>> sdf_filepath = get_data_file_path('molecules/ethanol.sdf')
    >>> molecule = Molecule(sdf_filepath)

    Convert to OpenEye OEMol object

    >>> oemol = molecule.to_openeye()

    Create a molecule from an OpenEye molecule

    >>> molecule = Molecule.from_openeye(oemol)

    Convert to RDKit Mol object

    >>> rdmol = molecule.to_rdkit()

    Create a molecule from an RDKit molecule

    >>> molecule = Molecule.from_rdkit(rdmol)

    Create a molecule from IUPAC name (requires the OpenEye toolkit)

    >>> molecule = Molecule.from_iupac('imatinib')

    Create a molecule from SMILES

    >>> molecule = Molecule.from_smiles('Cc1ccccc1')

    .. warning :: This API is experimental and subject to change.

    """

    def __init__(self, *args, **kwargs):
        """
        Create a new Molecule object

        Parameters
        ----------
        other : optional, default=None
            If specified, attempt to construct a copy of the Molecule from the specified object.
            This can be any one of the following:

            * a :class:`Molecule` object
            * a file that can be used to construct a :class:`Molecule` object
            * an ``openeye.oechem.OEMol``
            * an ``rdkit.Chem.rdchem.Mol``
            * a serialized :class:`Molecule` object

        Examples
        --------

        Create an empty molecule:

        >>> empty_molecule = Molecule()

        Create a molecule from a file that can be used to construct a molecule,
        using either a filename or file-like object:

        >>> from openforcefield.utils import get_data_file_path
        >>> sdf_filepath = get_data_file_path('molecules/ethanol.sdf')
        >>> molecule = Molecule(sdf_filepath)
        >>> molecule = Molecule(open(sdf_filepath, 'r'), file_format='sdf')

        >>> import gzip
        >>> mol2_gz_filepath = get_data_file_path('molecules/toluene.mol2.gz')
        >>> molecule = Molecule(gzip.GzipFile(mol2_gz_filepath, 'r'), file_format='mol2')

        Create a molecule from another molecule:

        >>> molecule_copy = Molecule(molecule)

        Convert to OpenEye OEMol object

        >>> oemol = molecule.to_openeye()

        Create a molecule from an OpenEye molecule:

        >>> molecule = Molecule(oemol)

        Convert to RDKit Mol object

        >>> rdmol = molecule.to_rdkit()

        Create a molecule from an RDKit molecule:

        >>> molecule = Molecule(rdmol)

        Create a molecule from a serialized molecule object:

        >>> serialized_molecule = molecule.__getstate__()
        >>> molecule_copy = Molecule(serialized_molecule)

        .. todo ::

           * If a filename or file-like object is specified but the file contains more than one molecule, what is the
           proper behavior? Read just the first molecule, or raise an exception if more than one molecule is found?

           * Should we also support SMILES strings or IUPAC names for ``other``?

        """
        #super(self, Molecule).__init__(*args, **kwargs)
        super(Molecule, self).__init__(*args, **kwargs)

    # TODO: Change this to add_atom(Atom) to improve encapsulation and extensibility?
    def add_atom(self,
                 atomic_number,
                 formal_charge,
                 is_aromatic,
                 stereochemistry=None,
                 name=None):
        """
        Add an atom

        Parameters
        ----------
        atomic_number : int
            Atomic number of the atom
        formal_charge : int
            Formal charge of the atom
        is_aromatic : bool
            If True, atom is aromatic; if False, not aromatic
        stereochemistry : str, optional, default=None
            Either 'R' or 'S' for specified stereochemistry, or None if stereochemistry is irrelevant
        name : str, optional, default=None
            An optional name for the atom

        Returns
        -------
        index : int
            The index of the atom in the molecule

        Examples
        --------

        Define a methane molecule

        >>> molecule = Molecule()
        >>> molecule.name = 'methane'
        >>> C = molecule.add_atom(6, 0, False)
        >>> H1 = molecule.add_atom(1, 0, False)
        >>> H2 = molecule.add_atom(1, 0, False)
        >>> H3 = molecule.add_atom(1, 0, False)
        >>> H4 = molecule.add_atom(1, 0, False)
        >>> bond_idx = molecule.add_bond(C, H1, False, 1)
        >>> bond_idx = molecule.add_bond(C, H2, False, 1)
        >>> bond_idx = molecule.add_bond(C, H3, False, 1)
        >>> bond_idx = molecule.add_bond(C, H4, False, 1)

        """
        atom_index = self._add_atom(
            atomic_number,
            formal_charge,
            is_aromatic,
            stereochemistry=stereochemistry,
            name=name)
        return atom_index

    def add_bond_charge_virtual_site(self,
                                     atoms,
                                     distance,
                                     charge_increments=None,
                                     weights=None,
                                     epsilon=None,
                                     sigma=None,
                                     rmin_half=None,
                                     name=''):
        """
        Create a bond charge-type virtual site, in which the location of the charge is specified by the position of two atoms. This supports placement of a virtual site S along a vector between two specified atoms, e.g. to allow for a sigma hole for halogens or similar contexts. With positive values of the distance, the virtual site lies outside the first indexed atom.
        Parameters
        ----------
        atoms : list of openforcefield.topology.molecule.Atom objects or ints of shape [N
            The atoms defining the virtual site's position or their indices
        distance : float

        weights : list of floats of shape [N] or None, optional, default=None
            weights[index] is the weight of particles[index] contributing to the position of the virtual site. Default is None
        charge_increments : list of floats of shape [N], optional, default=None
            The amount of charge to remove from the VirtualSite's atoms and put in the VirtualSite. Indexing in this list should match the ordering in the atoms list. Default is None.
        epsilon : float
            Epsilon term for VdW properties of virtual site. Default is None.
        sigma : float, default=None
            Sigma term for VdW properties of virtual site. Default is None.
        rmin_half : float
            Rmin_half term for VdW properties of virtual site. Default is None.
        name : string or None, default=''
            The name of this virtual site. Default is ''.

        Returns
        -------
        index : int
            The index of the newly-added virtual site in the molecule

        """

        vsite_index = self._add_bond_charge_virtual_site(
            atoms,
            distance,
            weights=weights,
            charge_increments=charge_increments,
            epsilon=epsilon,
            sigma=sigma,
            rmin_half=rmin_half,
            name=name)
        return vsite_index

    def add_monovalent_lone_pair_virtual_site(self, atoms, distance,
                                              out_of_plane_angle,
                                              in_plane_angle, **kwargs):
        """
        Create a bond charge-type virtual site, in which the location of the charge is specified by the position of three atoms.

        TODO : Do "weights" have any meaning here?

        Parameters
        ----------
        atoms : list of three openforcefield.topology.molecule.Atom objects or ints
            The three atoms defining the virtual site's position or their molecule atom indices
        distance : float

        out_of_plane_angle : float

        in_plane_angle : float

        epsilon : float
            Epsilon term for VdW properties of virtual site. Default is None.
        sigma : float, default=None
            Sigma term for VdW properties of virtual site. Default is None.
        rmin_half : float
            Rmin_half term for VdW properties of virtual site. Default is None.
        name : string or None, default=''
            The name of this virtual site. Default is ''.

        Returns
        -------
        index : int
            The index of the newly-added virtual site in the molecule


        """

        #vsite_index = self._add_monovalent_lone_pair_virtual_site(self, atoms, distance, out_of_plane_angle, in_plane_angle, charge_increments=charge_increments, weights=weights, epsilon=epsilon, sigma=sigma, rmin_half=rmin_half, name=name)
        vsite_index = self._add_monovalent_lone_pair_virtual_site(
            atoms, distance, out_of_plane_angle, in_plane_angle, **kwargs)
        return vsite_index

    #def add_divalent_lone_pair_virtual_site(self, atoms, distance, out_of_plane_angle, in_plane_angle, charge_increments=None, weights=None, epsilon=None, sigma=None, rmin_half=None, name=None):
    def add_divalent_lone_pair_virtual_site(self, atoms, distance,
                                            out_of_plane_angle, in_plane_angle,
                                            **kwargs):
        """
        Create a divalent lone pair-type virtual site, in which the location of the charge is specified by the position of three atoms.

        TODO : Do "weights" have any meaning here?

        Parameters
        ----------
        atoms : list of 3 openforcefield.topology.molecule.Atom objects or ints
            The three atoms defining the virtual site's position or their molecule atom indices
        distance : float

        out_of_plane_angle : float

        in_plane_angle : float

        epsilon : float
            Epsilon term for VdW properties of virtual site. Default is None.
        sigma : float, default=None
            Sigma term for VdW properties of virtual site. Default is None.
        rmin_half : float
            Rmin_half term for VdW properties of virtual site. Default is None.
        name : string or None, default=''
            The name of this virtual site. Default is ''.

        Returns
        -------
        index : int
            The index of the newly-added virtual site in the molecule

        """
        #vsite_index = self._add_divalent_lone_pair_virtual_site(self, atoms, distance, out_of_plane_angle, in_plane_angle, charge_increments=charge_increments, weights=weights, epsilon=epsilon, sigma=sigma, rmin_half=rmin_half, name=name)
        vsite_index = self._add_divalent_lone_pair_virtual_site(
            atoms, distance, out_of_plane_angle, in_plane_angle, **kwargs)
        return vsite_index

    def add_trivalent_lone_pair_virtual_site(self, atoms, distance,
                                             out_of_plane_angle,
                                             in_plane_angle, **kwargs):
        """
        Create a trivalent lone pair-type virtual site, in which the location of the charge is specified by the position of four atoms.

        TODO : Do "weights" have any meaning here?

        Parameters
        ----------
        atoms : list of 4 openforcefield.topology.molecule.Atom objects or ints
            The three atoms defining the virtual site's position or their molecule atom indices
        distance : float

        out_of_plane_angle : float

        in_plane_angle : float

        epsilon : float
            Epsilon term for VdW properties of virtual site. Default is None.
        sigma : float, default=None
            Sigma term for VdW properties of virtual site. Default is None.
        rmin_half : float
            Rmin_half term for VdW properties of virtual site. Default is None.
        name : string or None, default=''
            The name of this virtual site. Default is ''.

        Returns
        -------
        index : int
            The index of the newly-added virtual site in the molecule

        """
        vsite_index = self._add_trivalent_lone_pair_virtual_site(
            atoms, distance, out_of_plane_angle, in_plane_angle, **kwargs)
        return vsite_index

    def add_bond(self,
                 atom1,
                 atom2,
                 bond_order,
                 is_aromatic,
                 stereochemistry=None,
                 fractional_bond_order=None):
        """
        Add a bond between two specified atom indices


        Parameters
        ----------
        atom1 : int or openforcefield.topology.molecule.Atom
            Index of first atom
        atom2 : int or openforcefield.topology.molecule.Atom
            Index of second atom
        bond_order : int
            Integral bond order of Kekulized form
        is_aromatic : bool
            True if this bond is aromatic, False otherwise
        stereochemistry : str, optional, default=None
            Either 'E' or 'Z' for specified stereochemistry, or None if stereochemistry is irrelevant
        fractional_bond_order : float, optional, default=None
            The fractional (eg. Wiberg) bond order

        Returns
        -------
        index: int
            Index of the bond in this molecule

"""
        bond_index = self._add_bond(
            atom1,
            atom2,
            bond_order,
            is_aromatic,
            stereochemistry=stereochemistry,
            fractional_bond_order=fractional_bond_order)
        return bond_index

    def add_conformer(self, coordinates):
        """
        # TODO: Should this not be public?
        Adds a conformer of the molecule

        Parameters
        ----------
        coordinates: simtk.unit.Quantity(np.array) with shape (n_atoms, 3)
            Coordinates of the new conformer, with the first dimension of the array corresponding to the atom index in
            the Molecule's indexing system.
        Returns
        -------
        index: int
            Index of the conformer in the Molecule


"""
        return self._add_conformer(coordinates)
