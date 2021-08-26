"""
Components of molecular topologies other than Topology and Molecule
"""
import operator
import warnings
from abc import abstractmethod
from collections import OrderedDict
from copy import deepcopy
from typing import Optional, Union

import numpy as np
from simtk import unit
from simtk.openmm.app import Element, element

import openff.toolkit
from openff.toolkit.utils import quantity_to_string, string_to_quantity
from openff.toolkit.utils.exceptions import (
    InvalidConformerError,
    NotAttachedToMoleculeError,
    SmilesParsingError,
)
from openff.toolkit.utils.serialization import Serializable
from openff.toolkit.utils.toolkits import (
    DEFAULT_AROMATICITY_MODEL,
    GLOBAL_TOOLKIT_REGISTRY,
    InvalidToolkitRegistryError,
    OpenEyeToolkitWrapper,
    RDKitToolkitWrapper,
    ToolkitRegistry,
    ToolkitWrapper,
    UndefinedStereochemistryError,
)
from openff.toolkit.utils.utils import MissingDependencyError, requires_package


class Particle(Serializable):
    """
    Base class for all particles in a molecule.

    A particle object could be an ``Atom`` or a ``VirtualSite``.

    .. warning :: This API is experimental and subject to change.
    """

    @property
    def molecule(self):
        """
        The ``Molecule`` this particle is part of.

        .. todo::

            * Should we have a single unique ``Molecule`` for each molecule
              type in the system, or if we have multiple copies of the same
              molecule, should we have multiple ``Molecule``\ s?

        """
        return self._molecule

    @molecule.setter
    def molecule(self, molecule):
        """
        Set the particle's molecule pointer. Note that this will only work if the particle currently
        doesn't have a molecule
        """
        err = f"{type(self).__name__} already has an associated molecule"
        assert self._molecule is None, err
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


# =============================================================================================
# Atom
# =============================================================================================


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

    def __init__(
        self,
        atomic_number,
        formal_charge,
        is_aromatic,
        name=None,
        molecule=None,
        stereochemistry=None,
    ):
        """
        Create an immutable Atom object.

        Object is serializable and immutable.

        .. todo :: Use attrs to validate?

        .. todo :: We can add setters if we need to.

        Parameters
        ----------
        atomic_number : int
            Atomic number of the atom
        formal_charge : int or simtk.unit.Quantity-wrapped int with dimension "charge"
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
        # Use the setter here, since it will handle either ints or Quantities
        self.formal_charge = formal_charge
        self._is_aromatic = is_aromatic
        self._stereochemistry = stereochemistry
        if name is None:
            name = ""
        self._name = name
        self._molecule = molecule
        ## From Jeff: I'm going to assume that this is implicit in the parent Molecule's ordering of atoms
        # self._molecule_atom_index = molecule_atom_index
        self._bonds = list()
        self._virtual_sites = list()

    # TODO: We can probably avoid an explicit call and determine this dynamically
    #   from self._molecule (maybe caching the result) to get rid of some bookkeeping.
    def add_bond(self, bond):
        """Adds a bond that this atom is involved in
        .. todo :: Is this how we want to keep records?

        Parameters
        ----------
        bond: an openff.toolkit.topology.molecule.Bond
            A bond involving this atom
        """

        self._bonds.append(bond)
        # self._stereochemistry = None

    def add_virtual_site(self, vsite):
        """Adds a bond that this atom is involved in
        .. todo :: Is this how we want to keep records?

        Parameters
        ----------
        bond: an openff.toolkit.topology.molecule.Bond
            A bond involving this atom
        """

        self._virtual_sites.append(vsite)

    def to_dict(self):
        """Return a dict representation of the atom."""
        # TODO
        atom_dict = OrderedDict()
        atom_dict["atomic_number"] = self._atomic_number
        atom_dict["formal_charge"] = self._formal_charge / unit.elementary_charge
        atom_dict["is_aromatic"] = self._is_aromatic
        atom_dict["stereochemistry"] = self._stereochemistry
        # TODO: Should we let atoms have names?
        atom_dict["name"] = self._name
        # TODO: Should this be implicit in the atom ordering when saved?
        # atom_dict['molecule_atom_index'] = self._molecule_atom_index
        return atom_dict

    @classmethod
    def from_dict(cls, atom_dict):
        """Create an Atom from a dict representation."""
        ## TODO: classmethod or static method? Classmethod is needed for Bond, so it have
        ## its _molecule set and then look up the Atom on each side of it by ID
        return cls.__init__(*atom_dict)

    @property
    def formal_charge(self):
        """
        The atom's formal charge
        """
        return self._formal_charge

    @formal_charge.setter
    def formal_charge(self, other):
        from openff.toolkit.utils.utils import check_units_are_compatible

        """
        Set the atom's formal charge. Accepts either ints or simtk.unit.Quantity-wrapped ints with units of charge.
        """
        if isinstance(other, int):
            self._formal_charge = other * unit.elementary_charge
        else:
            check_units_are_compatible("formal charge", other, unit.elementary_charge)
            self._formal_charge = other

    @property
    def partial_charge(self):
        """
        The partial charge of the atom, if any.

        Returns
        -------
        simtk.unit.Quantity with dimension of atomic charge, or None if no charge has been specified
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

        # if (value != 'CW') and (value != 'CCW') and not(value is None):
        #    raise Exception("Atom stereochemistry setter expected 'CW', 'CCW', or None. Received {} (type {})".format(value, type(value)))
        self._stereochemistry = value

    @property
    def element(self):
        """
        The element of this atom.

        Returns
        -------
        simtk.openmm.app.element.Element
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
                "In setting atom name. Expected str, received {} (type {})".format(
                    other, type(other)
                )
            )
        self._name = other

    # TODO: How are we keeping track of bonds, angles, etc?

    @property
    def bonds(self):
        """
        The list of ``Bond`` objects this atom is involved in.

        """
        return self._bonds
        # for bond in self._bonds:
        #    yield bond

    @property
    # def bonded_to(self):
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
        atom2: openff.toolkit.topology.molecule.Atom
            a different atom in the same molecule

        Returns
        -------
        bool
            Whether this atom is bound to atom2
        """
        # TODO: Sanity check (check for same molecule?)
        assert self != atom2
        for bond in self._bonds:
            for bonded_atom in bond.atoms:
                if atom2 == bonded_atom:
                    return True
        return False

    @property
    def is_in_ring(self):
        """
        Return whether or not this atom is in a ring(s) (of any size)

        """
        if self._molecule is None:
            raise NotAttachedToMoleculeError(
                "This Atom does not belong to a Molecule object"
            )

        return any([self.molecule_atom_index in ring for ring in self._molecule.rings])

    @property
    def virtual_sites(self):
        """
        The list of ``VirtualSite`` objects this atom is involved in.

        """
        return self._virtual_sites
        # for vsite in self._vsites:
        #    yield vsite

    @property
    def molecule_atom_index(self):
        """
        The index of this Atom within the the list of atoms in ``Molecules``.
        Note that this can be different from ``molecule_particle_index``.

        """
        if self._molecule is None:
            raise ValueError("This Atom does not belong to a Molecule object")
        return self._molecule.atoms.index(self)

    @property
    def molecule_particle_index(self):
        """
        The index of this Atom within the the list of particles in the parent ``Molecule``.
        Note that this can be different from ``molecule_atom_index``.

        """
        if self._molecule is None:
            raise ValueError("This Atom does not belong to a Molecule object")
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
        return "Atom(name={}, atomic number={})".format(self._name, self._atomic_number)

    def __str__(self):
        # TODO: Also include particle_index and which molecule this atom belongs to?
        return "<Atom name='{}' atomic number='{}'>".format(
            self._name, self._atomic_number
        )


# =============================================================================================
# VirtualParticle
# =============================================================================================


class VirtualParticle(Particle):
    """
    A single particle owned by a VirtualSite

    .. warning :: This API is experimental and subject to change.
    """

    def __init__(self, vsite, orientation, name=None):
        """
        A single particle owned by a VirtualSite

        Parameters
        ----------
        vsite : openff.toolkit.topology.VirtualSite
            The parent VirtualSite of this VirtualParticle
        orientation : tuple of int
            Molecule atom indices of parent atoms
        name : str, optional
            The name of the particle

        """
        self._virtual_site = vsite
        self._molecule = vsite.molecule
        self._orientation = orientation
        self._name = name

    @property
    def virtual_site(self):
        return self._virtual_site

    @property
    def orientation(self):
        return self._orientation

    @property
    def virtual_site_particle_index(self):
        """
        The index of the particle relative to its owning virtual site. Normally
        this should either be 0 or 1.
        """
        return self.virtual_site.orientations.index(self.orientation)

    def _position(self, atom_positions):
        """
        Calculations the position of a virtual particle, as defined by the OpenMM
        :class:`simtk.openmm.openmm.LocalCoordinatesSite` definition.

        The frame is first constructed using the input atoms, where the weights defined
        by each virtual site are used. The virtual particle positions are then
        determined by setting the displacements, also determined uniquely by each
        virtual site definition.

        Note that, following the definition of the OpenMM LocalCoordinatesSite, the
        frame is forced to be orthogonal. This is first enforced by only allowing the
        x- and y-axis to be defined, since the z-axis must be normal to this plane.
        Then, y is then reset to be normal to the zx plane. This should ensure that the
        frame is orthonormal (after normalization).

        Note that this returns a 1D flat list as it is meant to be appended into a
        (M, 3) array via the public interface.

        Parameters
        ----------
        atom_positions: iterable of int
            The indices of the atoms, relative to the indices defined by the owning
            molecule. This is necessary since this particle has a certain orientation,
            so the input atoms must be in the original input ordering which was used to
            define the orientation.

        Returns
        -------
        :class:`simtk.unit.Quantity` of dimension [Length] in unit Angstrom wrapping an
        numpy.ndarray
        """

        # Although the above docstring claims that we fully implement
        # the OpenMM behavior, it has not been compared to OpenMM
        # at the source code level. If positions seem to be inconsistent,
        # please submit a bug report! We have tests to make sure our
        # implemented types are correct, so we are interested in cases
        # where custom virtual sites cause breakage.

        originwt, xdir, ydir = self.virtual_site.local_frame_weights
        disp = self.virtual_site.local_frame_position
        _unit = disp.unit
        x, y, z = disp / _unit

        # this pulls the correct ordering of the atoms
        pos = []
        for atom in self._orientation:
            pos.append(atom_positions[atom])

        atom_positions = pos

        originwt = np.atleast_2d(originwt)
        atom_positions = np.atleast_2d(atom_positions)

        origin = np.dot(originwt, atom_positions).sum(axis=0)

        xaxis, yaxis = np.dot(np.vstack((xdir, ydir)), atom_positions)

        zaxis = np.cross(xaxis, yaxis)
        yaxis = np.cross(zaxis, xaxis)

        def _normalize(axis):
            L = np.linalg.norm(axis)
            if L > 0.0:
                axis /= L
            return axis

        xaxis, yaxis, zaxis = map(_normalize, (xaxis, yaxis, zaxis))

        position = origin + x * xaxis + y * yaxis + z * zaxis

        return unit.Quantity(position, unit=_unit)

    def _extract_position_from_conformer(self, conformation):

        indices = [atom.molecule_atom_index for atom in self.virtual_site.atoms]

        atom_positions = [conformation[i] for i in indices]

        return atom_positions

    def _get_conformer(self, conformer_idx):

        assert self.molecule
        assert len(self.molecule.conformers) > 0

        conformer = self.molecule.conformers[conformer_idx]

        return conformer

    def compute_position_from_conformer(self, conformer_idx):
        """
        Compute the position of this virtual particle given an existing
        conformer owned by the parent molecule/virtual site.

        Parameters
        ----------
        conformer_idx : int
            The index of the conformer in the owning molecule.

        Returns
        -------
        :class:`simtk.unit.Quantity` of dimension [Length] in unit Angstroms wrapping a
        numpy.ndarray
            The positions of the virtual particles belonging to this virtual site.
            The array is the size (M, 3) where M is the number of virtual particles
            belonging to this virtual site.

        """

        atom_positions = self._get_conformer(conformer_idx)

        return self.compute_position_from_atom_positions(atom_positions)

    def compute_position_from_atom_positions(self, atom_positions):
        """
        Compute the position of this virtual site particle given a set of coordinates.

        Parameters
        ----------
        atom_positions : :class:`simtk.unit.Quantity` of dimension [Length] wrapping a
        numpy.ndarray
            The positions of all atoms in the molecule. The array is the size (N, 3)
            where N is the number of atoms in the molecule.

        Returns
        -------
        :class:`simtk.unit.Quantity` of dimension [Length] in unit Angstroms wrapping a
        numpy.ndarray
            The positions of the virtual particles belonging to this virtual site.
            The array is the size (M, 3) where M is the number of virtual particles
            belonging to this virtual site.

        """

        return self._position(atom_positions)


# =============================================================================================
# VirtualSite
# =============================================================================================


class VirtualSite(Particle):
    """
    A container representing one or more virtual particles whose positions are
    defined in terms of ``Atom`` positions. This container enables the coupling
    of particles that are symmetric about some axis/plane of the underlying
    atoms. For example, a single virtual site can represent two lone pairs of a
    water molecule, where the angle and distance parameters are expected to stay
    coupled, and are reflections across the plane of symmetry.

    Note that chemical atoms are represented by the ``Atom``.


    .. warning :: This API is experimental and subject to change.

    .. todo::

       * Should a virtual site be able to belong to more than one Topology?
       * Should virtual sites be immutable or mutable?

    """

    def __init__(
        self,
        atoms,
        charge_increments=None,
        epsilon=None,
        sigma=None,
        rmin_half=None,
        name=None,
        orientations=None,
    ):
        """
        Base class for VirtualSites

        .. todo ::

           * change sigma/epsilon/rmin_half to have units

        Parameters
        ----------
        atoms : list of Atom of shape [N]
            atoms[index] is the corresponding Atom
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
        orientation : list of int tuples or None, default=None
            The ordering of the atoms used to define the frame of the virtual site.
        """

        # Ensure we have as many charge_increments as we do atoms
        if not (charge_increments is None):
            if not (len(charge_increments) == len(atoms)):
                raise Exception(
                    "VirtualSite definition must have same number of charge_increments ({}) and atoms({})".format(
                        len(charge_increments), len(atoms)
                    )
                )
        else:
            charge_increments = ([0.0] * len(atoms)) * unit.elementary_charge

        # set sane defaults for OpenMM
        if epsilon is None and rmin_half is None:
            epsilon = 0.0 * unit.kilocalorie_per_mole
        if sigma is None and rmin_half is None:
            sigma = 0.0 * unit.angstrom

        # VdW parameters can either be epsilon+rmin_half or epsilon+sigma, but not both
        if not (epsilon is None):
            if (rmin_half is not None) and (sigma is not None):
                raise Exception(
                    "VirtualSite constructor given epsilon (value : {}), rmin_half (value : {}), and sigma (value : {}). If epsilon is nonzero, it should receive either rmin_half OR sigma".format(
                        epsilon, rmin_half, sigma
                    )
                )
            if (rmin_half is None) and (sigma is None):
                raise Exception(
                    "VirtualSite constructor given epsilon (value : {}) but not given rmin_half (value : {}) or sigma (value : {})".format(
                        epsilon, rmin_half, sigma
                    )
                )
            if sigma is None:
                # TODO: Save the 6th root of 2 if this starts being slow.
                sigma = (2.0 * rmin_half) / (2.0 ** (1.0 / 6))

        elif epsilon is None:
            if (rmin_half is not None) or (sigma is not None):
                raise Exception(
                    "VirtualSite constructor given rmin_half (value : {}) or sigma (value : {}), but not epsilon (value : {})".format(
                        rmin_half, sigma, epsilon
                    )
                )

        # Perform type-checking
        # for atom in atoms:
        #     assert isinstance(atom, Atom)
        # for atom_index in range(len(atoms) - 1):
        #     assert atoms[atom_index].molecule is atoms[atom_index + 1].molecule
        # assert isinstance(atoms[1].molecule, FrozenMolecule)

        if sigma is None:
            self._sigma = None
        else:
            assert hasattr(sigma, "unit")
            assert unit.angstrom.is_compatible(sigma.unit)
            self._sigma = sigma.in_units_of(unit.angstrom)

        if epsilon is None:
            self._epsilon = None
        else:
            assert hasattr(epsilon, "unit")
            assert (unit.kilojoule_per_mole).is_compatible(epsilon.unit)
            self._epsilon = epsilon.in_units_of(unit.kilojoule_per_mole)

        if charge_increments is None:
            self._charge_increments = None
        else:
            for ci in charge_increments:
                assert hasattr(ci, "unit")
                assert unit.elementary_charges.is_compatible(ci.unit)
            self._charge_increments = [
                ci.value_in_unit(unit.elementary_charges) for ci in charge_increments
            ] * unit.elementary_charges

        self._atoms = list()

        for atom in atoms:
            atom.add_virtual_site(self)
            self._atoms.append(atom)
        self._molecule = atoms[0].molecule

        self._name = name

        if orientations is None:
            ornt = [tuple(atom.molecule_atom_index for atom in atoms)]
            self._orientations = ornt
            self._particles = {ornt[0]: VirtualParticle(self, ornt[0])}
        else:
            ornt = None
            if type(orientations[0]) is int:
                ornt = [tuple(orientations)]
            else:
                ornt = [tuple(x) for x in orientations]
            self._orientations = ornt
            self._particles = dict(
                {order: VirtualParticle(self, order) for order in ornt}
            )

        # Subclassing makes _type unnecessary
        # self._type = None
        # TODO: Validate site types against allowed values

        # self._weights = np.array(weights) # make a copy and convert to array internally

    def __eq__(self, other):
        if not issubclass(type(other), VirtualSite):
            return False
        if self.type != other.type:
            return False
        same_name = self.name == other.name
        same_indices = self.atoms == other.atoms
        same_mol = self.molecule is other.molecule
        same_vsite = same_name and same_indices and same_mol
        return same_vsite

    def to_dict(self):
        """
        Return a dict representation of the virtual site.

        """
        # Each subclass should have its own to_dict
        vsite_dict = OrderedDict()
        vsite_dict["name"] = self._name
        vsite_dict["atoms"] = tuple([i.molecule_atom_index for i in self.atoms])
        vsite_dict["charge_increments"] = quantity_to_string(self._charge_increments)

        vsite_dict["epsilon"] = quantity_to_string(self._epsilon)

        vsite_dict["sigma"] = quantity_to_string(self._sigma)
        vsite_dict["orientations"] = self._orientations

        # skip packing the particles; they are created dynamically

        return vsite_dict

    @classmethod
    def from_dict(cls, vsite_dict):
        """Create a virtual site from a dict representation."""
        # Each subclass needs to have its own from_dict

        # Make a copy of the vsite_dict, where we'll unit-wrap the appropriate values
        vsite_dict_units = deepcopy(vsite_dict)

        # Attach units to epsilon term
        vsite_dict_units["epsilon"] = string_to_quantity(vsite_dict["epsilon"])
        vsite_dict_units["sigma"] = string_to_quantity(vsite_dict["sigma"])
        vsite_dict_units["charge_increments"] = string_to_quantity(
            vsite_dict["charge_increments"]
        )

        vsite_dict_units["orientation"] = cls._orientation

        return VirtualSite(**vsite_dict_units)

    def index_of_orientation(self, virtual_particle):
        """
        Return the orientation used by the given virtual particle.

        Parameters
        ----------
        virtual_particle : VirtualParticle
            The virtual particle contained in this virual site

        Returns
        -------
        A tuple of atom indices
        """

        for i, vp in enumerate(self.particles):
            if vp.orientation == virtual_particle.orientation:
                return i
        assert ValueError(
            "The given virtual particle was not found in this Virtual Site"
        )

    @property
    def orientations(self):
        """
        The orientations used by the virtual site particles.

        Orientations are an implementation to allow generation and coupling of multiple
        particles using the same physical definition. We can do this by allowing each
        particle to use a specific ordering of bases when calculating the positions.
        This is similar to improper torsion angles: the angle you find depends on the
        atom ordering used in the calculation.

        Before the positions are constructed, the parent atoms are reordered according
        to the particle's orientation. Each virtual particle has exactly one
        orientation. Since the frame of the virtual site is defined by a static list of
        weights and masks, we are able to influence how the local frame is constructed
        by crafting specific ordering the parent atoms.

        As a concrete example, we could define a TIP5 water by using one virtual site,
        and the particles have orientations (0, 1, 2) and (2, 1, 0). This means that,
        given that we are using a right-handed coordinate system, the z-axis will
        point in opposite directions for each particle. Using the same
        ``out_of_plane_angle`` and ``distance`` will therefore result in two unique
        particle positions.

        Using the toolkit API allows arbitrary selection of orientations. The SMIRNOFF
        specification, via the offxml file format, the orientations are controlled
        bondtype the "match" attribute. In this case, only the keywords "once" and
        "all_permuations" are allowed, meaning only the first orientation or all
        possible orientations are generated.

        The virtual site adders via :class:`Molecule` simplify this by optionally using
        a ``symmetric`` kwarg, which is the equivalent to the XML ``match`` keyword
        described above. However, the symmetric kwarg is not available for sites
        which symmetry is not possible, e.g. :class:`TrivalentLonePairVirtualSite`,
        provided a layer of sanity checking.  For the TIP5 example above, setting
        ``symmetric=True`` (the default) should automatically produce both particles.

        Parameters
        ----------

        Returns
        -------
        List of tuples of ints specifying the ordering of the parent atoms.
        """
        return self._orientations

    @property
    def particles(self):
        """
        Particles owned by this VirtualSite
        """
        for vp in self._particles.values():
            yield vp

    @property
    def n_particles(self):
        """
        The number of particles that the virtual site represents
        """
        # Virtual sites can represent multiple particles in a system
        # Assume a 1 to 1 mapping of orientations to particles for now
        # This means a virtualsite can only represent a single physical set
        # of parameters (distance, angle, etc)
        return len(self._particles)

    @property
    def molecule_virtual_site_index(self):
        """
        The index of this VirtualSite within the list of virtual sites within ``Molecule``
        Note that this can be different from ``particle_index``.
        """
        # if self._topology is None:
        #    raise ValueError('This VirtualSite does not belong to a Topology object')
        # TODO: This will be slow; can we cache this and update it only when needed?
        #       Deleting atoms/molecules in the Topology would have to invalidate the cached index.
        return self._molecule.virtual_sites.index(self)

    # @property
    # def molecule_particle_index(self):
    #     """
    #     The index of this VirtualSite within the the list of particles in the parent ``Molecule``.
    #     Note that this can be different from ``molecule_virtual_site_index``.

    #     """
    #     if self._molecule is None:
    #         raise ValueError(
    #             'This VirtualSite does not belong to a Molecule object')
    #     return self._molecule.particles.index(self)

    @property
    def atoms(self):
        """
        Atoms on whose position this VirtualSite depends.
        """
        return self._atoms
        # for atom in self._atoms:
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
        rmin = 2.0 ** (1.0 / 6) * self._sigma
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

    @property
    @abstractmethod
    def local_frame_weights(self):
        """
        The per-atom weights used to define the virtual site frame.

        The SMIRNOFF virtual sites use the definition of
        :class:`simtk.openmm.openmm.LocalCoordinatesSite` implemented by OpenMM.
        As such, the weights are used to determine the origin and the x and y axes of
        the local frame. Since the frame is an orthogonal bases, the z axis is not
        specified as it is assumed to be the cross of the x and y axes (using a
        right-handed coordinates).

        The weights defined refer to the weights of each atom's positions. For the
        origin, the weights must sum to 1. For the x and y axes, the weights much each
        sum to 0. For example, for a custom bond charge virtual site with two atoms:

        - Origin: [.5, .5] The origin of the frame is always in between atom 1 and
          atom 2. The calculation is 0.5 * atom1.xyz + 0.5 * atom2.xyz
        - X-Axis: [-1, 1] The x-axis points from atom 1 to atom 2. Positive
          displacements of this axis are closer to atom 2.
        - Y-Axis: [0, 0] This axis must be defined, so here we set it to the null
          space. Any displacements along y are sent to 0. Because of this, the z-axis
          will also be 0.

        The displacements along the axes defined here are defined/returned by
        :attr:`VirtualSite.local_frame_position`.

        To implement a new virtual site type (using a LocalCoordinatesSite
        definition), override this function.

        Parameters
        ----------

        Returns
        -------
        Tuple of list of weights used to define the origin, x-axis, and y-axis
        """

    @property
    @abstractmethod
    def local_frame_position(self):
        """
        The displacements of the virtual site relative to the local frame.

        The SMIRNOFF virtual sites use the definition of
        :class:`simtk.openmm.openmm.LocalCoordinatesSite` as implemented by OpenMM.
        As such, the frame positions refer to positions as defined by the frame, or the
        local axes defined by the owning atoms (see
        :attr:`VirtualSite.local_frame_weights`).

        To implement a new virtual site type (using a LocalCoordinatesSite
        definition), override this function.

        Parameters
        ----------

        Returns
        -------
        :class:`simtk.unit.Quantity` of dimension [Length] wrapping a list of
        displacements in the local frame for the x, y, and z directions.
        """

    def __repr__(self):
        # TODO: Also include particle_index, which molecule this atom belongs to?
        return "VirtualSite(name={}, type={}, atoms={})".format(
            self.name, self.type, self.atoms
        )

    def __str__(self):
        # TODO: Also include particle_index, which molecule this atom belongs to?
        return "<VirtualSite name={} type={} atoms={} particles={}>".format(
            self.name, self.type, self.atoms, self.n_particles
        )

    def _openmm_virtual_site(self, atoms):

        from simtk.openmm import LocalCoordinatesSite

        originwt, xdir, ydir = self.local_frame_weights
        pos = self.local_frame_position

        return LocalCoordinatesSite(atoms, originwt, xdir, ydir, pos)

    def compute_positions_from_conformer(self, conformer_idx):
        """
        Compute the position of the virtual site particles given an existing conformer
        owned by the parent molecule.

        Parameters
        ----------
        conformer_idx : int
            The index of the conformer in the owning molecule.

        Returns
        -------
        :class:`simtk.unit.Quantity` of dimension [Length] in unit Angstroms wrapping a
        numpy.ndarray
            The positions of the virtual particles belonging to this virtual site.
            The array is the size (M, 3) where M is the number of virtual particles
            belonging to this virtual site.
        """

        positions = []
        for vp in self.particles:
            vp_pos = vp.compute_position_from_conformer(conformer_idx)
            positions.append(vp_pos.value_in_unit(unit.angstrom))

        return unit.Quantity(np.array(positions).reshape(-1, 3), unit=unit.angstrom)

    def compute_positions_from_atom_positions(self, atom_positions):
        """
        Compute the positions of the virtual site particles given a set of coordinates.

        Parameters
        ----------
        atom_positions : :class:`simtk.unit.Quantity` of dimension [Length] wrapping a
        numpy.ndarray
            The positions of all atoms in the molecule. The array is the size (N, 3)
            where N is the number of atoms in the molecule.

        Returns
        -------
        :class:`simtk.unit.Quantity` of dimension [Length] in unit Angstroms wrapping a
        numpy.ndarray
            The positions of the virtual particles belonging to this virtual site.
            The array is the size (M, 3) where M is the number of virtual particles
            belonging to this virtual site.
        """

        positions = []
        for vp in self.particles:
            vp_pos = vp.compute_position_from_atom_positions(atom_positions)
            positions.extend(vp_pos.value_in_unit(unit.angstrom))

        return unit.Quantity(np.array(positions).reshape(-1, 3), unit=unit.angstrom)


class BondChargeVirtualSite(VirtualSite):
    """
    A particle representing a "Bond Charge"-type virtual site, in which the location of the charge is specified by the positions of two atoms. This supports placement of a virtual site S along a vector between two specified atoms, e.g. to allow for a sigma hole for halogens or similar contexts. With positive values of the distance, the virtual site lies outside the first indexed atom.

    .. warning :: This API is experimental and subject to change.
    """

    def __init__(
        self,
        atoms,
        distance,
        charge_increments=None,
        epsilon=None,
        sigma=None,
        rmin_half=None,
        name=None,
        orientations=None,
    ):
        """
        Create a bond charge-type virtual site, in which the location of the charge is specified by the position of two atoms. This supports placement of a virtual site S along a vector between two specified atoms, e.g. to allow for a sigma hole for halogens or similar contexts. With positive values of the distance, the virtual site lies outside the first indexed atom.

        TODO: One of the examples on https://open-forcefield-toolkit.readthedocs.io/en/topology/smirnoff.html#virtualsites-virtual-sites-for-off-atom-charges has a BondCharge defined with three atoms -- How does that work?

        Parameters
        ----------
        atoms : list of openff.toolkit.topology.molecule.Atom objects of shape [N]
            The atoms defining the virtual site's position

        distance : :class:`simtk.unit.Quantity` of dimension [Length] wrapping a scalar

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
        orientations : list of tuples of 3 Atoms or ints
            The permutations of the matched atoms that should be used to define
            the orientation of each virtual site particle
        """
        assert hasattr(distance, "unit")
        assert unit.angstrom.is_compatible(distance.unit)

        super().__init__(
            atoms,
            charge_increments=charge_increments,
            epsilon=epsilon,
            sigma=sigma,
            rmin_half=rmin_half,
            name=name,
            orientations=orientations,
        )
        self._distance = distance.in_units_of(unit.angstrom)

    def __eq__(self, other):
        return super().__eq__(other)

    def to_dict(self):
        vsite_dict = super().to_dict()
        vsite_dict["distance"] = quantity_to_string(self._distance)

        vsite_dict["vsite_type"] = self.type
        vsite_dict["orientations"] = self._orientations

        return vsite_dict

    @classmethod
    def from_dict(cls, vsite_dict):
        base_dict = deepcopy(vsite_dict)
        # Make sure it's the right type of virtual site
        assert vsite_dict["vsite_type"] == "BondChargeVirtualSite"
        base_dict.pop("vsite_type")
        base_dict.pop("distance")
        vsite = super().from_dict(**base_dict)
        vsite._distance = string_to_quantity(vsite_dict["distance"])
        return vsite

    @property
    def distance(self):
        """The distance parameter of the virtual site"""
        return self._distance

    @property
    def local_frame_weights(self):
        """
        Returns the local frame weights used to calculate the particle positions.
        See :attr:`VirtualSite.local_frame_weights` for a general description.

        Bond charge virtual sites are defined by the axis defined by the two
        atoms that define the bond. Since the virtual site position is defined
        solely by this axis, the other y-axis is defined but not used.

        Parameters
        ----------

        Returns
        -------
        Tuple of list of weights used to define the origin, x-axis, and y-axis.
        """

        originwt = [1.0, 0.0]  # first atom is origin

        xdir = [-1.0, 1.0]

        ydir = [-1.0, 1.0]

        return originwt, xdir, ydir

    @property
    def local_frame_position(self):
        """
        The displacements of the virtual site relative to the local frame.
        See :attr:`VirtualSite.local_frame_position` for a general description.

        Parameters
        ----------

        Returns
        -------
        :class:`simtk.unit.Quantity` of dimension [Length] wrapping a list of
        displacements in the local frame for the x, y, and z directions.
        """

        # since the origin is atom 1, and xdir is a unit vector pointing
        # towards the center of the other atoms, we want the
        # vsite to point away from the unit vector to achieve the desired
        # distance
        _unit = self._distance.unit
        pos = _unit * [-self._distance / _unit, 0.0, 0.0]

        return pos

    def get_openmm_virtual_site(self, atoms):
        """
        Returns the OpenMM virtual site corresponding to this BondChargeVirtualSite.

        Parameters
        ----------
        atoms : iterable of int
            The indices of the atoms involved in this virtual site.

        Returns
        -------
        :class:`simtk.openmm.openmm.LocalCoordinatesSite`
        """
        assert len(atoms) >= 2
        return self._openmm_virtual_site(atoms)


class MonovalentLonePairVirtualSite(VirtualSite):
    """
    A particle representing a "Monovalent Lone Pair"-type virtual site, in which the location of the charge is specified by the positions of three atoms. This is originally intended for situations like a carbonyl, and allows placement of a virtual site S at a specified distance d, in_plane_angle, and out_of_plane_angle relative to a central atom and two connected atoms.

    .. warning :: This API is experimental and subject to change.
    """

    def __init__(
        self,
        atoms,
        distance,
        out_of_plane_angle,
        in_plane_angle,
        charge_increments=None,
        epsilon=None,
        sigma=None,
        rmin_half=None,
        name=None,
        orientations=None,
    ):
        """
        Create a bond charge-type virtual site, in which the location of the charge is specified by the position of three atoms.

        Parameters
        ----------
        atoms : list of three openff.toolkit.topology.molecule.Atom objects
            The three atoms defining the virtual site's position

        distance : :class:`simtk.unit.Quantity` of dimension [Length] wrapping a scalar

        out_of_plane_angle : :class:`simtk.unit.Quantity` of dimension [Angle] wrapping
        a scalar

        in_plane_angle : :class:`simtk.unit.Quantity` of dimension [Angle] wrapping a
        scalar

        epsilon : float
            Epsilon term for VdW properties of virtual site. Default is None.
        sigma : float, default=None
            Sigma term for VdW properties of virtual site. Default is None.
        rmin_half : float
            Rmin_half term for VdW properties of virtual site. Default is None.
        name : string or None, default=None
            The name of this virtual site. Default is None.
        orientations : list of tuples of 3 Atoms or ints
            The permutations of the matched atoms that should be used to define
            the orientation of each virtual site particle
        """
        # assert isinstance(distance, unit.Quantity)
        # TODO: Check for proper number of atoms
        assert hasattr(distance, "unit")
        assert unit.angstrom.is_compatible(distance.unit)
        assert hasattr(in_plane_angle, "unit")
        assert unit.degree.is_compatible(in_plane_angle.unit)
        assert hasattr(out_of_plane_angle, "unit")
        assert unit.degree.is_compatible(out_of_plane_angle.unit)

        assert len(atoms) == 3
        super().__init__(
            atoms,
            charge_increments=charge_increments,
            epsilon=epsilon,
            sigma=sigma,
            rmin_half=rmin_half,
            name=name,
            orientations=orientations,
        )
        self._distance = distance.in_units_of(unit.angstrom)
        self._out_of_plane_angle = out_of_plane_angle.in_units_of(unit.degree)
        self._in_plane_angle = in_plane_angle.in_units_of(unit.degree)

    def to_dict(self):
        vsite_dict = super().to_dict()
        vsite_dict["distance"] = quantity_to_string(self._distance)
        vsite_dict["out_of_plane_angle"] = quantity_to_string(self._out_of_plane_angle)
        vsite_dict["in_plane_angle"] = quantity_to_string(self._in_plane_angle)
        vsite_dict["vsite_type"] = self.type
        return vsite_dict

    def __eq__(self, other):
        return super().__eq__(other)

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
        vsite = super().from_dict(vsite_dict)
        vsite._out_of_plane_angle = string_to_quantity(vsite_dict["out_of_plane_angle"])
        vsite._in_plane_angle = string_to_quantity(vsite_dict["in_plane_angle"])
        return vsite

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

    @property
    def local_frame_weights(self):
        """
        Returns the local frame weights used to calculate the particle positions.
        See :attr:`VirtualSite.local_frame_weights` for a general description.

        Parameters
        ----------

        Returns
        -------
        Tuple of list of weights used to define the origin, x-axis, and y-axis.
        """

        originwt = [1.0, 0.0, 0.0]

        xdir = [-1.0, 1.0, 0.0]
        ydir = [-1.0, 0.0, 1.0]

        return originwt, xdir, ydir

    @property
    def local_frame_position(self):
        """
        The displacements of the virtual site relative to the local frame.
        See :attr:`VirtualSite.local_frame_position` for a general description.

        Parameters
        ----------

        Returns
        -------
        :class:`simtk.unit.Quantity` of dimension [Length] wrapping a list of displacements
        in the local frame for the x, y, and z directions.
        """

        theta = self._in_plane_angle.value_in_unit(unit.radians)
        psi = self._out_of_plane_angle.value_in_unit(unit.radians)

        _unit = self._distance.unit
        pos = unit.Quantity(
            [
                self._distance / _unit * np.cos(theta) * np.cos(psi),
                self._distance / _unit * np.sin(theta) * np.cos(psi),
                self._distance / _unit * np.sin(psi),
            ],
            unit=_unit,
        )

        return pos

    def get_openmm_virtual_site(self, atoms):
        """
        Returns the OpenMM virtual site corresponding to this
        MonovalentLonePairVirtualSite.

        Parameters
        ----------
        atoms : iterable of int
            The indices of the atoms involved in this virtual site.

        Returns
        -------
        :class:`simtk.openmm.openmm.LocalCoordinatesSite`
        """

        assert len(atoms) >= 3
        return self._openmm_virtual_site(atoms)


class DivalentLonePairVirtualSite(VirtualSite):
    """
    A particle representing a "Divalent Lone Pair"-type virtual site, in which the location of the charge is specified by the positions of three atoms. This is suitable for cases like four-point and five-point water models as well as pyrimidine; a charge site S lies a specified distance d from the central atom among three atoms along the bisector of the angle between the atoms (if out_of_plane_angle is zero) or out of the plane by the specified angle (if out_of_plane_angle is nonzero) with its projection along the bisector. For positive values of the distance d the virtual site lies outside the 2-1-3 angle and for negative values it lies inside.
    """

    def __init__(
        self,
        atoms,
        distance,
        out_of_plane_angle,
        charge_increments=None,
        epsilon=None,
        sigma=None,
        rmin_half=None,
        name=None,
        orientations=None,
    ):
        """
        Create a divalent lone pair-type virtual site, in which the location of the charge is specified by the position of three atoms.

        Parameters
        ----------
        atoms : list of 3 openff.toolkit.topology.molecule.Atom objects
            The three atoms defining the virtual site's position

        distance : :class:`simtk.unit.Quantity` of dimension [Length] wrapping a scalar

        out_of_plane_angle : :class:`simtk.unit.Quantity` of dimension [Angle] wrapping
        a scalar

        epsilon : float
            Epsilon term for VdW properties of virtual site. Default is None.
        sigma : float, default=None
            Sigma term for VdW properties of virtual site. Default is None.
        rmin_half : float
            Rmin_half term for VdW properties of virtual site. Default is None.
        name : string or None, default=None
            The name of this virtual site. Default is None.
        orientations : list of tuples of 3 Atoms or ints
            The permutations of the matched atoms that should be used to define
            the orientation of each virtual site particle
        """
        assert hasattr(distance, "unit")
        assert unit.angstrom.is_compatible(distance.unit)

        assert hasattr(out_of_plane_angle, "unit")
        assert unit.degree.is_compatible(out_of_plane_angle.unit)

        assert len(atoms) == 3
        super().__init__(
            atoms,
            charge_increments=charge_increments,
            epsilon=epsilon,
            sigma=sigma,
            rmin_half=rmin_half,
            name=name,
            orientations=orientations,
        )
        self._distance = distance.in_units_of(unit.angstrom)
        self._out_of_plane_angle = out_of_plane_angle.in_units_of(unit.degree)

    def __eq__(self, other):
        return super().__eq__(other)

    def to_dict(self):
        vsite_dict = super().to_dict()
        vsite_dict["distance"] = quantity_to_string(self._distance)
        vsite_dict["out_of_plane_angle"] = quantity_to_string(self._out_of_plane_angle)
        vsite_dict["vsite_type"] = self.type
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
        vsite = super().from_dict(vsite_dict)
        vsite._out_of_plane_angle = string_to_quantity(vsite_dict["out_of_plane_angle"])
        return vsite

    @property
    def distance(self):
        """The distance parameter of the virtual site"""
        return self._distance

    @property
    def out_of_plane_angle(self):
        """The out_of_plane_angle parameter of the virtual site"""
        return self._out_of_plane_angle

    @property
    def local_frame_weights(self):
        """
        Returns the local frame weights used to calculate the particle positions.
        See :attr:`VirtualSite.local_frame_weights` for a general description.

        Parameters
        ----------

        Returns
        -------
        Tuple of list of weights used to define the origin, x-axis, and y-axis.
        """

        originwt = [0.0, 1.0, 0.0]

        xdir = [0.5, -1.0, 0.5]
        ydir = [1.0, -1.0, 0.0]

        return originwt, xdir, ydir

    @property
    def local_frame_position(self):
        """
        The displacements of the virtual site relative to the local frame.
        See :attr:`VirtualSite.local_frame_position` for a general description.

        Parameters
        ----------

        Returns
        -------
        :class:`simtk.unit.Quantity` of dimension [Length] wrapping a list of
        displacements in the local frame for the x, y, and z directions.
        """

        theta = self._out_of_plane_angle.value_in_unit(unit.radians)

        _unit = self._distance.unit

        pos = _unit * [
            -self._distance / _unit * np.cos(theta),
            0.0,
            self._distance / _unit * np.sin(theta),
        ]  # pos of the vsite in local crds
        return pos

    def get_openmm_virtual_site(self, atoms):
        """
        Returns the OpenMM virtual site corresponding to this
        DivalentLonePairVirtualSite.

        Parameters
        ----------
        atoms : iterable of int
            The indices of the atoms involved in this virtual site.

        Returns
        -------
        :class:`simtk.openmm.openmm.LocalCoordinatesSite`
        """

        assert len(atoms) >= 3
        return self._openmm_virtual_site(atoms)


class TrivalentLonePairVirtualSite(VirtualSite):
    """
    A particle representing a "Trivalent Lone Pair"-type virtual site, in which the location of the charge is specified by the positions of four atoms. This is suitable for planar or tetrahedral nitrogen lone pairs; a charge site S lies above the central atom (e.g. nitrogen a distance d along the vector perpendicular to the plane of the three connected atoms (2,3,4). With positive values of d the site lies above the nitrogen and with negative values it lies below the nitrogen.

    .. warning :: This API is experimental and subject to change.
    """

    def __init__(
        self,
        atoms,
        distance,
        charge_increments=None,
        epsilon=None,
        sigma=None,
        rmin_half=None,
        name=None,
        orientations=None,
    ):
        """
        Create a trivalent lone pair-type virtual site, in which the location of the charge is specified by the position of four atoms.

        Parameters
        ----------
        atoms : list of 4 openff.toolkit.topology.molecule.Atom objects
            The three atoms defining the virtual site's position

        distance : :class:`simtk.unit.Quantity` of dimension [Length] wrapping a scalar

        epsilon : float
            Epsilon term for VdW properties of virtual site. Default is None.
        sigma : float, default=None
            Sigma term for VdW properties of virtual site. Default is None.
        rmin_half : float
            Rmin_half term for VdW properties of virtual site. Default is None.
        name : string or None, default=None
            The name of this virtual site. Default is None.
        orientations : list of tuples of 3 Atoms or ints
            The permutations of the matched atoms that should be used to define
            the orientation of each virtual site particle
        """
        assert len(atoms) == 4

        assert hasattr(distance, "unit")
        assert unit.angstrom.is_compatible(distance.unit)

        super().__init__(
            atoms,
            charge_increments=charge_increments,
            epsilon=epsilon,
            sigma=sigma,
            rmin_half=rmin_half,
            name=name,
            orientations=orientations,
        )
        self._distance = distance.in_units_of(unit.angstrom)

    def __eq__(self, other):
        return super().__eq__(other)

    def to_dict(self):
        vsite_dict = super().to_dict()
        vsite_dict["distance"] = quantity_to_string(self._distance)
        vsite_dict["vsite_type"] = self.type
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
    def local_frame_weights(self):
        """
        Returns the local frame weights used to calculate the particle positions.
        See :attr:`VirtualSite.local_frame_weights` for a general description.

        Parameters
        ----------

        Returns
        -------
        Tuple of list of weights used to define the origin, x-axis, and y-axis.
        """

        originwt = [0.0, 1.0, 0.0, 0.0]

        xdir = [1 / 3, -1.0, 1 / 3, 1 / 3]

        # ydir does not matter
        ydir = [1.0, -1.0, 0.0, 0.0]

        return originwt, xdir, ydir

    @property
    def local_frame_position(self):
        """
        The displacements of the virtual site relative to the local frame.
        See :attr:`VirtualSite.local_frame_position` for a general description.

        Parameters
        ----------

        Returns
        -------
        :class:`simtk.unit.Quantity` of dimension [Length] wrapping a list of
        displacements in the local frame for the x, y, and z directions.
        """

        _unit = self._distance.unit
        pos = unit.Quantity([-self._distance / _unit, 0.0, 0.0], unit=_unit)

        return pos

    def get_openmm_virtual_site(self, atoms):
        """
        Returns the OpenMM virtual site corresponding to this
        TrivalentLonePairVirtualSite.

        Parameters
        ----------
        atoms : iterable of int
            The indices of the atoms involved in this virtual site.

        Returns
        -------
        :class:`simtk.openmm.openmm.LocalCoordinatesSite`
        """

        assert len(atoms) >= 4
        return self._openmm_virtual_site(atoms)


# =============================================================================================
# Bond Stereochemistry
# =============================================================================================

# class BondStereochemistry(Serializable):
# """
# Bond stereochemistry representation
# """
# def __init__(self, stereo_type, neighbor1, neighbor2):
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

# def to_dict(self):
#    bs_dict = OrderedDict()
#    bs_dict['stereo_type'] = self._stereo_type
#    bs_dict['neighbor1_index'] = self._neighbor1.molecule_atom_index
#    bs_dict['neighbor2_index'] = self._neighbor2.molecule_atom_index
#    return bs_dict

# classmethod
# def from_dict(cls, molecule, bs_dict):
#    neighbor1 = molecule.atoms[bs_dict['neighbor1_index']]
#    neighbor2 = molecule.atoms[bs_dict['neighbor2_index']]
#    return cls.__init__(bs_dict['stereo_type'], neighbor1, neighbor2)

# @property
# def stereo_type(self):
#    return self._stereo_type

# @stereo_type.setter
# def stereo_type(self, value):
#    assert (value == 'CIS') or (value == 'TRANS') or (value is None)
#    self._stereo_type = value

# @property
# def neighbor1(self):
#    return self._neighbor1

# @property
# def neighbor2(self):
#    return self._neighbor2

# @property
# def neighbors(self):
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
    atom1, atom2 : openff.toolkit.topology.Atom
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

    def __init__(
        self,
        atom1,
        atom2,
        bond_order,
        is_aromatic,
        fractional_bond_order=None,
        stereochemistry=None,
    ):
        """
        Create a new chemical bond.

        """
        from openff.toolkit.topology.molecule import FrozenMolecule

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
        # self._type = bondtype
        self._fractional_bond_order = fractional_bond_order
        self._bond_order = bond_order
        self._is_aromatic = is_aromatic
        self._stereochemistry = stereochemistry

    def to_dict(self):
        """
        Return a dict representation of the bond.

        """
        bond_dict = OrderedDict()
        bond_dict["atom1"] = self.atom1.molecule_atom_index
        bond_dict["atom2"] = self.atom2.molecule_atom_index
        bond_dict["bond_order"] = self._bond_order
        bond_dict["is_aromatic"] = self._is_aromatic
        bond_dict["stereochemistry"] = self._stereochemistry
        bond_dict["fractional_bond_order"] = self._fractional_bond_order
        return bond_dict

    @classmethod
    def from_dict(cls, molecule, d):
        """Create a Bond from a dict representation."""
        # TODO
        d["molecule"] = molecule
        d["atom1"] = molecule.atoms[d["atom1"]]
        d["atom2"] = molecule.atoms[d["atom2"]]
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
            raise ValueError("This Atom does not belong to a Molecule object")
        return self._molecule.bonds.index(self)

    @property
    def is_in_ring(self):
        """
        Return whether or not this bond is in a ring(s) (of any size)

        """
        if self._molecule is None:
            raise NotAttachedToMoleculeError(
                "This Bond does not belong to a Molecule object"
            )

        for ring in self._molecule.rings:
            if self.atom1.molecule_atom_index in ring:
                if self.atom2.molecule_atom_index in ring:
                    return True
        return False

    def __repr__(self):
        return f"Bond(atom1 index={self.atom1_index}, atom2 index={self.atom2_index})"

    def __str__(self):
        return (
            f"<Bond atom1 index='{self.atom1_index}', atom2 index='{self.atom2_index}'>"
        )
