#!/usr/bin/env python

# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================
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
from openff.toolkit.topology.components import (
    Atom,
    Bond,
    BondChargeVirtualSite,
    DivalentLonePairVirtualSite,
    MonovalentLonePairVirtualSite,
    TrivalentLonePairVirtualSite,
)
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

# =============================================================================================
# Moclecule
# =============================================================================================

# TODO: How do we automatically trigger invalidation of cached properties if an ``Atom``, ``Bond``, or ``VirtualSite`` is modified,
#       rather than added/deleted via the API? The simplest resolution is simply to make them immutable.


class FrozenMolecule(Serializable):
    """
    Immutable chemical representation of a molecule, such as a small molecule or biopolymer.

    .. todo :: What other API calls would be useful for supporting biopolymers
               as small molecules? Perhaps iterating over chains and residues?

    Examples
    --------

    Create a molecule from a sdf file

    >>> from openff.toolkit.utils import get_data_file_path
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

    def __init__(
        self,
        other=None,
        file_format=None,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
        allow_undefined_stereo=False,
    ):
        """
        Create a new FrozenMolecule object

        .. todo ::

           * If a filename or file-like object is specified but the file
             contains more than one molecule, what is the proper behavior?
             Read just the first molecule, or raise an exception if more
             than one molecule is found?

           * Should we also support SMILES strings or IUPAC names for
             ``other``\ ?

        Parameters
        ----------
        other : optional, default=None
            If specified, attempt to construct a copy of the Molecule from
            the specified object. This can be any one of the following:

            * a :class:`Molecule` object
            * a file that can be used to construct a :class:`Molecule` object
            * an ``openeye.oechem.OEMol``
            * an ``rdkit.Chem.rdchem.Mol``
            * a serialized :class:`Molecule` object

        file_format : str, optional, default=None
            If providing a file-like object, you must specify the format
            of the data. If providing a file, the file format will attempt
            to be guessed from the suffix.
        toolkit_registry : a :class:`ToolkitRegistry` or
            :class:`ToolkitWrapper` object, optional,
            default=GLOBAL_TOOLKIT_REGISTRY :class:`ToolkitRegistry`
            or :class:`ToolkitWrapper` to use for I/O operations
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

        >>> from openff.toolkit.utils import get_data_file_path
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
            raise InvalidToolkitRegistryError(
                "'toolkit_registry' must be either a ToolkitRegistry or a ToolkitWrapper"
            )

        if other is None:
            self._initialize()
        else:
            loaded = False
            # Start a list of the ValueErrors the following logic encounters, so we can print it out
            # if there turned out to be no way to load this input
            value_errors = list()

            if isinstance(other, openff.toolkit.topology.FrozenMolecule) and not (
                loaded
            ):
                self._copy_initializer(other)
                loaded = True
            if isinstance(other, openff.toolkit.topology.Molecule) and not (loaded):
                # TODO: This will need to be updated once FrozenMolecules and Molecules are significantly different
                self._copy_initializer(other)
                loaded = True
            if isinstance(other, OrderedDict) and not (loaded):
                self.__setstate__(other)
                loaded = True

            # Check through the toolkit registry to find a compatible wrapper for loading
            if not loaded:
                try:
                    # Each ToolkitWrapper may provide a from_object method, which turns some particular type(s)
                    # of object into OFFMols. For example, RDKitToolkitWrapper's from_object method will
                    # return an OFFMol if provided with an RDMol, or raise a ValueError if it is provided
                    # an OEMol (or anything else). This makes the assumption that any non-ValueError errors raised
                    # by the toolkit _really are_ bad and should be raised immediately, which may be a bad assumption.
                    result = toolkit_registry.call(
                        "from_object",
                        other,
                        allow_undefined_stereo=allow_undefined_stereo,
                        raise_exception_types=[UndefinedStereochemistryError],
                        _cls=self.__class__,
                    )
                # NotImplementedError should never be raised... Only from_file and from_file_obj are provided
                # in the base ToolkitWrapper class and require overwriting, so from_object should be excluded
                # except NotImplementedError as e:
                #    raise e
                # The toolkit registry will aggregate all errors except UndefinedStereochemistryErrors into a single
                # ValueError, which we should catch and and store that here.
                except ValueError as e:
                    value_errors.append(e)
                else:
                    self._copy_initializer(result)
                    loaded = True
            # TODO: Make this compatible with file-like objects (I couldn't figure out how to make an oemolistream
            # from a fileIO object)
            if (isinstance(other, str) or hasattr(other, "read")) and not (loaded):
                try:
                    mol = Molecule.from_file(
                        other,
                        file_format=file_format,
                        toolkit_registry=toolkit_registry,
                        allow_undefined_stereo=allow_undefined_stereo,
                    )  # returns a list only if multiple molecules are found
                    if type(mol) == list:
                        raise ValueError(
                            "Specified file or file-like object must contain exactly one molecule"
                        )
                except ValueError as e:
                    value_errors.append(e)
                else:
                    self._copy_initializer(mol)
                    loaded = True

            # If none of the above methods worked, raise a ValueError summarizing the
            # errors from the different loading attempts

            if not (loaded):
                msg = "Cannot construct openff.toolkit.topology.Molecule from {}\n".format(
                    other
                )
                for value_error in value_errors:
                    msg += str(value_error)
                raise ValueError(msg)

    @property
    def has_unique_atom_names(self):
        """True if the molecule has unique atom names, False otherwise."""
        unique_atom_names = set([atom.name for atom in self.atoms])
        if len(unique_atom_names) < self.n_atoms:
            return False
        return True

    def generate_unique_atom_names(self):
        """
        Generate unique atom names using element name and number of times that element has occurred
        e.g. 'C1', 'H1', 'O1', 'C2', ...

        """
        from collections import defaultdict

        element_counts = defaultdict(int)
        for atom in self.atoms:
            symbol = atom.element.symbol
            element_counts[symbol] += 1
            atom.name = symbol + str(element_counts[symbol])

    def _validate(self):
        """
        Validate the molecule, ensuring it has unique atom names

        """
        if not self.has_unique_atom_names:
            self.generate_unique_atom_names()

    def strip_atom_stereochemistry(
        self, smarts, toolkit_registry=GLOBAL_TOOLKIT_REGISTRY
    ):
        """Delete stereochemistry information for certain atoms, if it is present.
        This method can be used to "normalize" molecules imported from different cheminformatics
        toolkits, which differ in which atom centers are considered stereogenic.

        Parameters
        ----------
        smarts: str or ChemicalEnvironment
            Tagged SMARTS with a single atom with index 1. Any matches for this atom will have any assigned
            stereocheistry information removed.
        toolkit_registry : a :class:`ToolkitRegistry` or :class:`ToolkitWrapper` object, optional, default=GLOBAL_TOOLKIT_REGISTRY
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for I/O operations

        """
        from openff.toolkit.typing.chemistry.environment import AtomChemicalEnvironment

        chem_env = AtomChemicalEnvironment(smarts)
        matches = self.chemical_environment_matches(
            chem_env, toolkit_registry=toolkit_registry
        )

        for match in set(matches):
            atom_idx = match[0]
            self.atoms[atom_idx].stereochemistry = None

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
        from openff.toolkit.utils.utils import serialize_numpy

        molecule_dict = OrderedDict()
        molecule_dict["name"] = self._name
        ## From Jeff: If we go the properties-as-dict route, then _properties should, at
        ## the top level, be a dict. Should we go through recursively and ensure all values are dicts too?
        molecule_dict["atoms"] = [atom.to_dict() for atom in self._atoms]
        molecule_dict["virtual_sites"] = [
            vsite.to_dict() for vsite in self._virtual_sites
        ]
        molecule_dict["bonds"] = [bond.to_dict() for bond in self._bonds]
        # TODO: Charges
        # TODO: Properties
        # From Jeff: We could have the onerous requirement that all "properties" have to_dict() functions.
        # Or we could restrict properties to simple stuff (ints, strings, floats, and the like)
        # Or pickle anything unusual
        # Or not allow user-defined properties at all (just use our internal _cached_properties)
        # molecule_dict['properties'] = dict([(key, value._to_dict()) for key.value in self._properties])
        # TODO: Assuming "simple stuff" properties right now, figure out a better standard
        molecule_dict["properties"] = self._properties
        if hasattr(self, "_cached_properties"):
            molecule_dict["cached_properties"] = self._cached_properties
        # TODO: Conformers
        if self._conformers is None:
            molecule_dict["conformers"] = None
        else:
            molecule_dict["conformers"] = []
            molecule_dict[
                "conformers_unit"
            ] = "angstrom"  # Have this defined as a class variable?
            for conf in self._conformers:
                conf_unitless = conf / unit.angstrom
                conf_serialized, conf_shape = serialize_numpy((conf_unitless))
                molecule_dict["conformers"].append(conf_serialized)
        if self._partial_charges is None:
            molecule_dict["partial_charges"] = None
            molecule_dict["partial_charges_unit"] = None

        else:
            charges_unitless = self._partial_charges / unit.elementary_charge
            charges_serialized, charges_shape = serialize_numpy(charges_unitless)
            molecule_dict["partial_charges"] = charges_serialized
            molecule_dict["partial_charges_unit"] = "elementary_charge"

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
        from openff.toolkit.utils.utils import deserialize_numpy

        self._initialize()
        self.name = molecule_dict["name"]
        for atom_dict in molecule_dict["atoms"]:
            self._add_atom(**atom_dict)

        # Handle virtual site unit reattachment and molecule tagging
        for vsite_dict in molecule_dict["virtual_sites"]:
            vsite_dict_units = deepcopy(vsite_dict)

            # Attach units to epsilon term
            vsite_dict_units["epsilon"] = string_to_quantity(vsite_dict["epsilon"])
            vsite_dict_units["sigma"] = string_to_quantity(vsite_dict["sigma"])
            vsite_dict_units["charge_increments"] = string_to_quantity(
                vsite_dict["charge_increments"]
            )
            vsite_dict_units["orientations"] = vsite_dict["orientations"]

            # Call the correct molecule._add_X_virtual_site function, based on the stated type
            # Also generate the Atom objects from the atom indices
            atoms = [self._atoms[i] for i in vsite_dict_units["atoms"]]
            vsite_dict_units["atoms"] = atoms
            if vsite_dict_units["vsite_type"] == "BondChargeVirtualSite":
                del vsite_dict_units["vsite_type"]
                vsite_dict_units["distance"] = string_to_quantity(
                    vsite_dict["distance"]
                )
                self._add_bond_charge_virtual_site(**vsite_dict_units)

            elif vsite_dict_units["vsite_type"] == "MonovalentLonePairVirtualSite":
                del vsite_dict_units["vsite_type"]
                vsite_dict_units["distance"] = string_to_quantity(
                    vsite_dict["distance"]
                )
                vsite_dict_units["in_plane_angle"] = string_to_quantity(
                    vsite_dict["in_plane_angle"]
                )
                vsite_dict_units["out_of_plane_angle"] = string_to_quantity(
                    vsite_dict["out_of_plane_angle"]
                )
                self._add_monovalent_lone_pair_virtual_site(**vsite_dict_units)

            elif vsite_dict_units["vsite_type"] == "DivalentLonePairVirtualSite":
                del vsite_dict_units["vsite_type"]
                vsite_dict_units["distance"] = string_to_quantity(
                    vsite_dict["distance"]
                )
                vsite_dict_units["out_of_plane_angle"] = string_to_quantity(
                    vsite_dict["out_of_plane_angle"]
                )
                self._add_divalent_lone_pair_virtual_site(**vsite_dict_units)

            elif vsite_dict_units["vsite_type"] == "TrivalentLonePairVirtualSite":
                del vsite_dict_units["vsite_type"]
                vsite_dict_units["distance"] = string_to_quantity(
                    vsite_dict["distance"]
                )
                self._add_trivalent_lone_pair_virtual_site(**vsite_dict_units)

            else:
                raise Exception(
                    "Vsite type {} not recognized".format(vsite_dict["vsite_type"])
                )

        for bond_dict in molecule_dict["bonds"]:
            bond_dict["atom1"] = int(bond_dict["atom1"])
            bond_dict["atom2"] = int(bond_dict["atom2"])
            self._add_bond(**bond_dict)

        if molecule_dict["partial_charges"] is None:
            self._partial_charges = None
        else:
            charges_shape = (self.n_atoms,)
            partial_charges_unitless = deserialize_numpy(
                molecule_dict["partial_charges"], charges_shape
            )
            pc_unit = getattr(unit, molecule_dict["partial_charges_unit"])
            partial_charges = unit.Quantity(partial_charges_unitless, pc_unit)
            self._partial_charges = partial_charges

        if molecule_dict["conformers"] is None:
            self._conformers = None
        else:
            self._conformers = list()
            for ser_conf in molecule_dict["conformers"]:
                # TODO: Update to use string_to_quantity
                conformers_shape = (self.n_atoms, 3)
                conformer_unitless = deserialize_numpy(ser_conf, conformers_shape)
                c_unit = getattr(unit, molecule_dict["conformers_unit"])
                conformer = unit.Quantity(conformer_unitless, c_unit)
                self._conformers.append(conformer)

        self._properties = molecule_dict["properties"]

    def __repr__(self):
        """Return the SMILES of this molecule"""
        return "Molecule with name '{}' and SMILES '{}'".format(
            self.name, self.to_smiles()
        )

    def __getstate__(self):
        return self.to_dict()

    def __setstate__(self, state):
        return self._initialize_from_dict(state)

    def _initialize(self):
        """
        Clear the contents of the current molecule.
        """
        self._name = ""
        self._atoms = list()
        self._virtual_sites = list()
        self._bonds = list()  # List of bonds between Atom objects
        self._properties = {}  # Attached properties to be preserved
        # self._cached_properties = None # Cached properties (such as partial charges) can be recomputed as needed
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
        # assert isinstance(other, type(self)), "can only copy instances of {}".format(type(self))

        # Run a deepcopy here so that items that were _always_ dict (like other.properties) will
        # not have any references to the old molecule
        other_dict = deepcopy(other.to_dict())
        self._initialize_from_dict(other_dict)

    def __eq__(self, other):
        """Test two molecules for equality to see if they are the chemical species, but do not check other annotated properties.

        .. note ::

           Note that this method simply tests whether two molecules are identical chemical species using equivalence of their canonical isomeric SMILES.
           No effort is made to ensure that the atoms are in the same order or that any annotated properties are preserved.

        """
        # updated to use the new isomorphic checking method, with full matching
        # TODO the doc string did not match the previous function what matching should this method do?
        return Molecule.are_isomorphic(self, other, return_atom_map=False)[0]

    def to_smiles(
        self,
        isomeric=True,
        explicit_hydrogens=True,
        mapped=False,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
    ):
        """
        Return a canonical isomeric SMILES representation of the current molecule.
        A partially mapped smiles can also be generated for atoms of interest by supplying an `atom_map` to the
        properties dictionary.

        .. note :: RDKit and OpenEye versions will not necessarily return the same representation.

        Parameters
        ----------
        isomeric: bool optional, default= True
            return an isomeric smiles
        explicit_hydrogens: bool optional, default=True
            return a smiles string containing all hydrogens explicitly
        mapped: bool optional, default=False
            return a explicit hydrogen mapped smiles, the atoms to be mapped can be controlled by supplying an
            atom map into the properties dictionary. If no mapping is passed all atoms will be mapped in order, else
            an atom map dictionary from the current atom index to the map id should be supplied with no duplicates.
            The map ids (values) should start from 0 or 1.
        toolkit_registry : openff.toolkit.utils.toolkits.ToolkitRegistry or openff.toolkit.utils.toolkits.ToolkitWrapper, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for SMILES conversion

        Returns
        -------
        smiles : str
            Canonical isomeric explicit-hydrogen SMILES

        Examples
        --------

        >>> from openff.toolkit.utils import get_data_file_path
        >>> sdf_filepath = get_data_file_path('molecules/ethanol.sdf')
        >>> molecule = Molecule(sdf_filepath)
        >>> smiles = molecule.to_smiles()

        """
        # Initialize cached_smiles dict for this molecule if none exists
        if self._cached_smiles is None:
            self._cached_smiles = {}

        # Figure out which toolkit should be used to create the SMILES
        if isinstance(toolkit_registry, ToolkitRegistry):
            to_smiles_method = toolkit_registry.resolve("to_smiles")
        elif isinstance(toolkit_registry, ToolkitWrapper):
            to_smiles_method = toolkit_registry.to_smiles
        else:
            raise InvalidToolkitRegistryError(
                "Invalid toolkit_registry passed to to_smiles. Expected ToolkitRegistry or ToolkitWrapper. Got  {}".format(
                    type(toolkit_registry)
                )
            )

        # Get a string representation of the function containing the toolkit name so we can check
        # if a SMILES was already cached for this molecule. This will return, for example
        # "RDKitToolkitWrapper.to_smiles"
        smiles_hash = (
            to_smiles_method.__qualname__
            + str(isomeric)
            + str(explicit_hydrogens)
            + str(mapped)
        )
        smiles_hash += str(self._properties.get("atom_map", None))
        # Check to see if a SMILES for this molecule was already cached using this method
        if smiles_hash in self._cached_smiles:
            return self._cached_smiles[smiles_hash]
        else:
            smiles = to_smiles_method(self, isomeric, explicit_hydrogens, mapped)
            self._cached_smiles[smiles_hash] = smiles
            return smiles

    @classmethod
    def from_inchi(
        cls,
        inchi,
        allow_undefined_stereo=False,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
    ):
        """
        Construct a Molecule from a InChI representation

        Parameters
        ----------
        inchi : str
            The InChI representation of the molecule.

        allow_undefined_stereo : bool, default=False
            Whether to accept InChI with undefined stereochemistry. If False,
            an exception will be raised if a InChI with undefined stereochemistry
            is passed into this function.

        toolkit_registry : openff.toolkit.utils.toolkits.ToolRegistry or openff.toolkit.utils.toolkits.ToolkitWrapper, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for InChI-to-molecule conversion


        Returns
        -------
        molecule : openff.toolkit.topology.Molecule

        Examples
        --------
        Make cis-1,2-Dichloroethene:

        >>> molecule = Molecule.from_inchi('InChI=1S/C2H2Cl2/c3-1-2-4/h1-2H/b2-1-')
        """

        if isinstance(toolkit_registry, ToolkitRegistry):
            molecule = toolkit_registry.call(
                "from_inchi",
                inchi,
                _cls=cls,
                allow_undefined_stereo=allow_undefined_stereo,
            )
        elif isinstance(toolkit_registry, ToolkitWrapper):
            toolkit = toolkit_registry
            molecule = toolkit.from_inchi(
                inchi, _cls=cls, allow_undefined_stereo=allow_undefined_stereo
            )
        else:
            raise InvalidToolkitRegistryError(
                "Invalid toolkit_registry passed to from_inchi. Expected ToolkitRegistry or ToolkitWrapper. Got  {}".format(
                    type(toolkit_registry)
                )
            )

        return molecule

    def to_inchi(self, fixed_hydrogens=False, toolkit_registry=GLOBAL_TOOLKIT_REGISTRY):
        """
        Create an InChI string for the molecule using the requested toolkit backend.
        InChI is a standardised representation that does not capture tautomers unless specified using the fixed hydrogen
        layer.

        For information on InChi see here https://iupac.org/who-we-are/divisions/division-details/inchi/

        Parameters
        ----------
        fixed_hydrogens: bool, default=False
            If a fixed hydrogen layer should be added to the InChI, if `True` this will produce a non standard specific
            InChI string of the molecule.

        toolkit_registry : openff.toolkit.utils.toolkits.ToolRegistry or openff.toolkit.utils.toolkits.ToolkitWrapper, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for molecule-to-InChI conversion

        Returns
        --------
        inchi: str
            The InChI string of the molecule.

        Raises
        -------
        InvalidToolkitRegistryError
             If an invalid object is passed as the toolkit_registry parameter
        """

        if isinstance(toolkit_registry, ToolkitRegistry):
            inchi = toolkit_registry.call(
                "to_inchi", self, fixed_hydrogens=fixed_hydrogens
            )
        elif isinstance(toolkit_registry, ToolkitWrapper):
            toolkit = toolkit_registry
            inchi = toolkit.to_inchi(self, fixed_hydrogens=fixed_hydrogens)
        else:
            raise InvalidToolkitRegistryError(
                "Invalid toolkit_registry passed to to_inchi. Expected ToolkitRegistry or ToolkitWrapper. Got  {}".format(
                    type(toolkit_registry)
                )
            )

        return inchi

    def to_inchikey(
        self, fixed_hydrogens=False, toolkit_registry=GLOBAL_TOOLKIT_REGISTRY
    ):
        """
        Create an InChIKey for the molecule using the requested toolkit backend.
        InChIKey is a standardised representation that does not capture tautomers unless specified using the fixed hydrogen
        layer.

        For information on InChi see here https://iupac.org/who-we-are/divisions/division-details/inchi/

        Parameters
        ----------
        fixed_hydrogens: bool, default=False
            If a fixed hydrogen layer should be added to the InChI, if `True` this will produce a non standard specific
            InChI string of the molecule.

        toolkit_registry : openff.toolkit.utils.toolkits.ToolRegistry or openff.toolkit.utils.toolkits.ToolkitWrapper, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for molecule-to-InChIKey conversion

        Returns
        --------
        inchi_key: str
            The InChIKey representation of the molecule.

        Raises
        -------
        InvalidToolkitRegistryError
             If an invalid object is passed as the toolkit_registry parameter
        """

        if isinstance(toolkit_registry, ToolkitRegistry):
            inchi_key = toolkit_registry.call(
                "to_inchikey", self, fixed_hydrogens=fixed_hydrogens
            )
        elif isinstance(toolkit_registry, ToolkitWrapper):
            toolkit = toolkit_registry
            inchi_key = toolkit.to_inchikey(self, fixed_hydrogens=fixed_hydrogens)
        else:
            raise InvalidToolkitRegistryError(
                "Invalid toolkit_registry passed to to_inchikey. Expected ToolkitRegistry or ToolkitWrapper. Got  {}".format(
                    type(toolkit_registry)
                )
            )

        return inchi_key

    @classmethod
    def from_smiles(
        cls,
        smiles,
        hydrogens_are_explicit=False,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
        allow_undefined_stereo=False,
    ):
        """
        Construct a Molecule from a SMILES representation

        Parameters
        ----------
        smiles : str
            The SMILES representation of the molecule.
        hydrogens_are_explicit : bool, default = False
            If False, the cheminformatics toolkit will perform hydrogen addition
        toolkit_registry : openff.toolkit.utils.toolkits.ToolkitRegistry or openff.toolkit.utils.toolkits.ToolkitWrapper, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for SMILES-to-molecule conversion
        allow_undefined_stereo : bool, default=False
            Whether to accept SMILES with undefined stereochemistry. If False,
            an exception will be raised if a SMILES with undefined stereochemistry
            is passed into this function.

        Returns
        -------
        molecule : openff.toolkit.topology.Molecule

        Examples
        --------

        >>> molecule = Molecule.from_smiles('Cc1ccccc1')

        """
        if isinstance(toolkit_registry, ToolkitRegistry):
            molecule = toolkit_registry.call(
                "from_smiles",
                smiles,
                hydrogens_are_explicit=hydrogens_are_explicit,
                allow_undefined_stereo=allow_undefined_stereo,
                _cls=cls,
            )
        elif isinstance(toolkit_registry, ToolkitWrapper):
            toolkit = toolkit_registry
            molecule = toolkit.from_smiles(
                smiles,
                hydrogens_are_explicit=hydrogens_are_explicit,
                allow_undefined_stereo=allow_undefined_stereo,
                _cls=cls,
            )
        else:
            raise InvalidToolkitRegistryError(
                "Invalid toolkit_registry passed to from_smiles. Expected ToolkitRegistry or ToolkitWrapper. Got  {}".format(
                    type(toolkit_registry)
                )
            )

        return molecule

    @staticmethod
    def are_isomorphic(
        mol1,
        mol2,
        return_atom_map=False,
        aromatic_matching=True,
        formal_charge_matching=True,
        bond_order_matching=True,
        atom_stereochemistry_matching=True,
        bond_stereochemistry_matching=True,
        strip_pyrimidal_n_atom_stereo=True,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
    ):
        """
        Determines whether the two molecules are isomorphic by comparing their graph representations and the chosen
        node/edge attributes. Minimally connections and atomic_number are checked.

        If nx.Graphs() are given they must at least have atomic_number attributes on nodes.
        other optional attributes for nodes are: is_aromatic, formal_charge and stereochemistry.
        optional attributes for edges are: is_aromatic, bond_order and stereochemistry.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        mol1 : an openff.toolkit.topology.molecule.FrozenMolecule or TopologyMolecule or nx.Graph()
        mol2 : an openff.toolkit.topology.molecule.FrozenMolecule or TopologyMolecule or nx.Graph()
            The molecule to test for isomorphism.

        return_atom_map: bool, default=False, optional
            will return an optional dict containing the atomic mapping.

        aromatic_matching: bool, default=True, optional
            compare the aromatic attributes of bonds and atoms.

        formal_charge_matching: bool, default=True, optional
            compare the formal charges attributes of the atoms.

        bond_order_matching: bool, deafult=True, optional
            compare the bond order on attributes of the bonds.

        atom_stereochemistry_matching : bool, default=True, optional
            If ``False``, atoms' stereochemistry is ignored for the
            purpose of determining equality.

        bond_stereochemistry_matching : bool, default=True, optional
            If ``False``, bonds' stereochemistry is ignored for the
            purpose of determining equality.

        strip_pyrimidal_n_atom_stereo: bool, default=True, optional
            If ``True``, any stereochemistry defined around pyrimidal
            nitrogen stereocenters will be disregarded in the isomorphism
            check.

        toolkit_registry : openff.toolkit.utils.toolkits.ToolkitRegistry or openff.toolkit.utils.toolkits.ToolkitWrapper, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for
            removing stereochemistry from pyrimidal nitrogens.

        Returns
        -------
        molecules_are_isomorphic : bool

        atom_map : default=None, Optional,
            [Dict[int,int]] ordered by mol1 indexing {mol1_index: mol2_index}
            If molecules are not isomorphic given input arguments, will return None instead of dict.
        """

        # Do a quick hill formula check first
        if Molecule.to_hill_formula(mol1) != Molecule.to_hill_formula(mol2):
            return False, None

        # Build the user defined matching functions
        def node_match_func(x, y):
            # always match by atleast atomic number
            is_equal = x["atomic_number"] == y["atomic_number"]
            if aromatic_matching:
                is_equal &= x["is_aromatic"] == y["is_aromatic"]
            if formal_charge_matching:
                is_equal &= x["formal_charge"] == y["formal_charge"]
            if atom_stereochemistry_matching:
                is_equal &= x["stereochemistry"] == y["stereochemistry"]
            return is_equal

        # check if we want to do any bond matching if not the function is None
        if aromatic_matching or bond_order_matching or bond_stereochemistry_matching:

            def edge_match_func(x, y):
                # We don't need to check the exact bond order (which is 1 or 2)
                # if the bond is aromatic. This way we avoid missing a match only
                # if the alternate bond orders 1 and 2 are assigned differently.
                if aromatic_matching and bond_order_matching:
                    is_equal = (x["is_aromatic"] == y["is_aromatic"]) or (
                        x["bond_order"] == y["bond_order"]
                    )
                elif aromatic_matching:
                    is_equal = x["is_aromatic"] == y["is_aromatic"]
                elif bond_order_matching:
                    is_equal = x["bond_order"] == y["bond_order"]
                else:
                    is_equal = None
                if bond_stereochemistry_matching:
                    if is_equal is None:
                        is_equal = x["stereochemistry"] == y["stereochemistry"]
                    else:
                        is_equal &= x["stereochemistry"] == y["stereochemistry"]

                return is_equal

        else:
            edge_match_func = None

        # Here we should work out what data type we have, also deal with lists?
        def to_networkx(data):
            """For the given data type, return the networkx graph"""
            import networkx as nx

            from openff.toolkit.topology import TopologyMolecule

            if strip_pyrimidal_n_atom_stereo:
                SMARTS = "[N+0X3:1](-[*])(-[*])(-[*])"

            if isinstance(data, FrozenMolecule):
                # Molecule class instance
                if strip_pyrimidal_n_atom_stereo:
                    # Make a copy of the molecule so we don't modify the original
                    data = deepcopy(data)
                    data.strip_atom_stereochemistry(
                        SMARTS, toolkit_registry=toolkit_registry
                    )
                return data.to_networkx()
            elif isinstance(data, TopologyMolecule):
                # TopologyMolecule class instance
                if strip_pyrimidal_n_atom_stereo:
                    # Make a copy of the molecule so we don't modify the original
                    ref_mol = deepcopy(data.reference_molecule)
                    ref_mol.strip_atom_stereochemistry(
                        SMARTS, toolkit_registry=toolkit_registry
                    )
                return ref_mol.to_networkx()

            elif isinstance(data, nx.Graph):
                return data

            else:
                raise NotImplementedError(
                    f"The input type {type(data)} is not supported,"
                    f"please supply an openff.toolkit.topology.molecule.Molecule,"
                    f"openff.toolkit.topology.topology.TopologyMolecule or networkx.Graph "
                    f"representation of the molecule."
                )

        mol1_netx = to_networkx(mol1)
        mol2_netx = to_networkx(mol2)
        from networkx.algorithms.isomorphism import GraphMatcher

        GM = GraphMatcher(
            mol1_netx, mol2_netx, node_match=node_match_func, edge_match=edge_match_func
        )
        isomorphic = GM.is_isomorphic()

        if isomorphic and return_atom_map:
            topology_atom_map = GM.mapping

            # reorder the mapping by keys
            sorted_mapping = {}
            for key in sorted(topology_atom_map.keys()):
                sorted_mapping[key] = topology_atom_map[key]

            return isomorphic, sorted_mapping

        else:
            return isomorphic, None

    def is_isomorphic_with(self, other, **kwargs):
        """
        Check if the molecule is isomorphic with the other molecule which can be an openff.toolkit.topology.Molecule,
        or TopologyMolecule or nx.Graph(). Full matching is done using the options described bellow.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        other: openff.toolkit.topology.Molecule or TopologyMolecule or nx.Graph()

        return_atom_map: bool, default=False, optional
            will return an optional dict containing the atomic mapping.

        aromatic_matching: bool, default=True, optional
        compare the aromatic attributes of bonds and atoms.

        formal_charge_matching: bool, default=True, optional
        compare the formal charges attributes of the atoms.

        bond_order_matching: bool, deafult=True, optional
        compare the bond order on attributes of the bonds.

        atom_stereochemistry_matching : bool, default=True, optional
            If ``False``, atoms' stereochemistry is ignored for the
            purpose of determining equality.

        bond_stereochemistry_matching : bool, default=True, optional
            If ``False``, bonds' stereochemistry is ignored for the
            purpose of determining equality.

        strip_pyrimidal_n_atom_stereo: bool, default=True, optional
            If ``True``, any stereochemistry defined around pyrimidal
            nitrogen stereocenters will be disregarded in the isomorphism
            check.

        toolkit_registry : openff.toolkit.utils.toolkits.ToolkitRegistry or openff.toolkit.utils.toolkits.ToolkitWrapper, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for
            removing stereochemistry from pyrimidal nitrogens.

        Returns
        -------
        isomorphic : bool
        """

        return Molecule.are_isomorphic(
            self,
            other,
            return_atom_map=False,
            aromatic_matching=kwargs.get("aromatic_matching", True),
            formal_charge_matching=kwargs.get("formal_charge_matching", True),
            bond_order_matching=kwargs.get("bond_order_matching", True),
            atom_stereochemistry_matching=kwargs.get(
                "atom_stereochemistry_matching", True
            ),
            bond_stereochemistry_matching=kwargs.get(
                "bond_stereochemistry_matching", True
            ),
            strip_pyrimidal_n_atom_stereo=kwargs.get(
                "strip_pyrimidal_n_atom_stereo", True
            ),
            toolkit_registry=kwargs.get("toolkit_registry", GLOBAL_TOOLKIT_REGISTRY),
        )[0]

    def generate_conformers(
        self,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
        n_conformers=10,
        rms_cutoff=None,
        clear_existing=True,
    ):
        """
        Generate conformers for this molecule using an underlying toolkit.

        If ``n_conformers=0``, no toolkit wrapper will be called. If ``n_conformers=0``
        and ``clear_existing=True``, ``molecule.conformers`` will be set to ``None``.

        Parameters
        ----------
        toolkit_registry : openff.toolkit.utils.toolkits.ToolkitRegistry or openff.toolkit.utils.toolkits.ToolkitWrapper, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for SMILES-to-molecule conversion
        n_conformers : int, default=1
            The maximum number of conformers to produce
        rms_cutoff : simtk.Quantity-wrapped float, in units of distance, optional, default=None
            The minimum RMS value at which two conformers are considered redundant and one is deleted. Precise
            implementation of this cutoff may be toolkit-dependent. If ``None``, the cutoff is set to be the default value
            for each ``ToolkitWrapper`` (generally 1 Angstrom).
        clear_existing : bool, default=True
            Whether to overwrite existing conformers for the molecule

        Examples
        --------

        >>> molecule = Molecule.from_smiles('CCCCCC')
        >>> molecule.generate_conformers()

        Raises
        ------
        InvalidToolkitRegistryError
            If an invalid object is passed as the toolkit_registry parameter

        """
        # If no conformers are requested, do not call to a ToolkitWrapper at all
        if n_conformers == 0:
            if clear_existing:
                self._conformers = None
            return

        if isinstance(toolkit_registry, ToolkitRegistry):
            return toolkit_registry.call(
                "generate_conformers",
                self,
                n_conformers=n_conformers,
                rms_cutoff=rms_cutoff,
                clear_existing=clear_existing,
            )
        elif isinstance(toolkit_registry, ToolkitWrapper):
            toolkit = toolkit_registry
            return toolkit.generate_conformers(
                self,
                n_conformers=n_conformers,
                rms_cutoff=rms_cutoff,
                clear_existing=clear_existing,
            )
        else:
            raise InvalidToolkitRegistryError(
                "Invalid toolkit_registry passed to generate_conformers. Expected ToolkitRegistry or ToolkitWrapper. Got  {}".format(
                    type(toolkit_registry)
                )
            )

    def compute_virtual_site_positions_from_conformer(self, conformer_idx):
        """
        Compute the position of all virtual sites given an existing
        conformer specified by its index.

        Parameters
        ----------
        conformer_idx : int
            The index of the conformer.

        Returns
        -------
        :class:`simtk.unit.Quantity` of dimension [Length] in unit Angstroms wrapping a
        numpy.ndarray
            The positions of the virtual particles belonging to this virtual site.
            The array is the size (M, 3) where M is the number of virtual particles
            belonging to this virtual site.
        """

        atom_positions = self._conformers[conformer_idx]
        return self.compute_virtual_site_positions_from_atom_positions(atom_positions)

    def compute_virtual_site_positions_from_atom_positions(self, atom_positions):
        """
        Compute the positions of the virtual sites in this molecule given a set of
        external coordinates. The coordinates do not need come from an internal
        conformer, but are assumed to have the same shape and be in the same order.

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

        for vsite in self.virtual_sites:
            vsite_pos = vsite.compute_positions_from_atom_positions(atom_positions)
            positions.append(vsite_pos.value_in_unit(unit.angstrom))

        return unit.Quantity(np.array(positions).reshape(-1, 3), unit=unit.angstrom)

    def apply_elf_conformer_selection(
        self,
        percentage: float = 2.0,
        limit: int = 10,
        toolkit_registry: Optional[
            Union[ToolkitRegistry, ToolkitWrapper]
        ] = GLOBAL_TOOLKIT_REGISTRY,
        **kwargs,
    ):
        """Applies the `ELF method
        <https://docs.eyesopen.com/toolkits/python/quacpactk/molchargetheory.html#elf-conformer-selection>`_
        to select a set of diverse conformers which have minimal
        electrostatically strongly interacting functional groups from a
        molecules conformers.

        Notes
        -----
        * The input molecule should have a large set of conformers already
          generated to select the ELF conformers from.
        * The selected conformers will be retained in the `conformers` list
          while unselected conformers will be discarded.

        See Also
        --------
        OpenEyeToolkitWrapper.apply_elf_conformer_selection
        RDKitToolkitWrapper.apply_elf_conformer_selection

        Parameters
        ----------
        toolkit_registry
            The underlying toolkit to use to select the ELF conformers.
        percentage
            The percentage of conformers with the lowest electrostatic
            interaction energies to greedily select from.
        limit
            The maximum number of conformers to select.
        """
        if isinstance(toolkit_registry, ToolkitRegistry):
            toolkit_registry.call(
                "apply_elf_conformer_selection",
                molecule=self,
                percentage=percentage,
                limit=limit,
                **kwargs,
            )
        elif isinstance(toolkit_registry, ToolkitWrapper):
            toolkit = toolkit_registry
            toolkit.apply_elf_conformer_selection(
                self, molecule=self, percentage=percentage, limit=limit, **kwargs
            )
        else:
            raise InvalidToolkitRegistryError(
                f"Invalid toolkit_registry passed to apply_elf_conformer_selection."
                f"Expected ToolkitRegistry or ToolkitWrapper. Got "
                f"{type(toolkit_registry)}"
            )

    def compute_partial_charges_am1bcc(
        self,
        use_conformers=None,
        strict_n_conformers=False,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
    ):
        """
        Calculate partial atomic charges for this molecule using AM1-BCC run by an underlying toolkit
        and assign them to this molecule's ``partial_charges`` attribute.

        Parameters
        ----------
        strict_n_conformers : bool, default=False
            Whether to raise an exception if an invalid number of conformers is provided for the given charge method.
            If this is False and an invalid number of conformers is found, a warning will be raised.
        use_conformers : iterable of simtk.unit.Quantity-wrapped numpy arrays, each with shape (n_atoms, 3) and dimension of distance. Optional, default=None
            Coordinates to use for partial charge calculation.
            If None, an appropriate number of conformers for the given charge method will be generated.
        toolkit_registry : openff.toolkit.utils.toolkits.ToolkitRegistry or openff.toolkit.utils.toolkits.ToolkitWrapper, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for the calculation

        Examples
        --------

        >>> molecule = Molecule.from_smiles('CCCCCC')
        >>> molecule.generate_conformers()
        >>> molecule.compute_partial_charges_am1bcc()

        Raises
        ------
        InvalidToolkitRegistryError
            If an invalid object is passed as the toolkit_registry parameter

        """
        self.assign_partial_charges(
            partial_charge_method="am1bcc",
            use_conformers=use_conformers,
            strict_n_conformers=strict_n_conformers,
            toolkit_registry=toolkit_registry,
        )

    def assign_partial_charges(
        self,
        partial_charge_method,
        strict_n_conformers=False,
        use_conformers=None,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
        normalize_partial_charges=True,
    ):
        """
        Calculate partial atomic charges for this molecule using an underlying toolkit, and assign
        the new values to the partial_charges attribute.

        Parameters
        ----------
        partial_charge_method : string
            The partial charge calculation method to use for partial charge calculation.
        strict_n_conformers : bool, default=False
            Whether to raise an exception if an invalid number of conformers is provided for the given charge method.
            If this is False and an invalid number of conformers is found, a warning will be raised.
        use_conformers : iterable of simtk.unit.Quantity-wrapped numpy arrays, each with shape (n_atoms, 3) and dimension of distance. Optional, default=None
            Coordinates to use for partial charge calculation. If None, an appropriate number of conformers will be generated.
        toolkit_registry : openff.toolkit.utils.toolkits.ToolkitRegistry or openff.toolkit.utils.toolkits.ToolkitWrapper, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for the calculation.
        normalize_partial_charges : bool, default=True
            Whether to offset partial charges so that they sum to the total formal charge of the molecule.
            This is used to prevent accumulation of rounding errors when the partial charge assignment method
            returns values at limited precision.

        Examples
        --------

        >>> molecule = Molecule.from_smiles('CCCCCC')
        >>> molecule.assign_partial_charges('am1-mulliken')

        Raises
        ------
        InvalidToolkitRegistryError
            If an invalid object is passed as the toolkit_registry parameter

        """
        if isinstance(toolkit_registry, ToolkitRegistry):
            # We may need to try several toolkitwrappers to find one
            # that supports the desired partial charge method, so we
            # tell the ToolkitRegistry to continue trying ToolkitWrappers
            # if one raises an error (raise_exception_types=[])
            toolkit_registry.call(
                "assign_partial_charges",
                molecule=self,
                partial_charge_method=partial_charge_method,
                use_conformers=use_conformers,
                strict_n_conformers=strict_n_conformers,
                normalize_partial_charges=normalize_partial_charges,
                raise_exception_types=[],
                _cls=self.__class__,
            )
        elif isinstance(toolkit_registry, ToolkitWrapper):
            toolkit = toolkit_registry
            toolkit.assign_partial_charges(
                self,
                partial_charge_method=partial_charge_method,
                use_conformers=use_conformers,
                strict_n_conformers=strict_n_conformers,
                normalize_partial_charges=normalize_partial_charges,
                _cls=self.__class__,
            )
        else:
            raise InvalidToolkitRegistryError(
                f"Invalid toolkit_registry passed to assign_partial_charges."
                f"Expected ToolkitRegistry or ToolkitWrapper. Got  {type(toolkit_registry)}"
            )

    def _normalize_partial_charges(self):
        """
        Add offsets to each partial charge to ensure that they sum to the formal charge of the molecule,
        to the limit of a python float's precision. Modifies the partial charges in-place.
        """
        expected_charge = self.total_charge

        current_charge = 0.0 * unit.elementary_charge
        for pc in self.partial_charges:
            current_charge += pc

        charge_offset = (expected_charge - current_charge) / self.n_atoms

        self.partial_charges += charge_offset

    def assign_fractional_bond_orders(
        self,
        bond_order_model=None,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
        use_conformers=None,
    ):
        """
        Update and store list of bond orders this molecule. Bond orders are stored on each
        bond, in the ``bond.fractional_bond_order`` attribute.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        toolkit_registry : openff.toolkit.utils.toolkits.ToolkitRegistry or openff.toolkit.utils.toolkits.ToolkitWrapper, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for SMILES-to-molecule conversion
        bond_order_model : string, optional. Default=None
            The bond order model to use for fractional bond order calculation. If ``None``, "am1-wiberg" will be used.
        use_conformers : iterable of simtk.unit.Quantity(np.array) with shape (n_atoms, 3) and dimension of distance, optional, default=None
            The conformers to use for fractional bond order calculation. If ``None``, an appropriate number
            of conformers will be generated by an available ToolkitWrapper.

        Examples
        --------

        >>> molecule = Molecule.from_smiles('CCCCCC')
        >>> molecule.assign_fractional_bond_orders()

        Raises
        ------
        InvalidToolkitRegistryError
            If an invalid object is passed as the toolkit_registry parameter

        """
        bond_order_model = bond_order_model.lower()

        if isinstance(toolkit_registry, ToolkitRegistry):
            return toolkit_registry.call(
                "assign_fractional_bond_orders",
                self,
                bond_order_model=bond_order_model,
                use_conformers=use_conformers,
            )
        elif isinstance(toolkit_registry, ToolkitWrapper):
            toolkit = toolkit_registry
            return toolkit.assign_fractional_bond_orders(
                self, bond_order_model=bond_order_model, use_conformers=use_conformers
            )
        else:
            raise InvalidToolkitRegistryError(
                f"Invalid toolkit_registry passed to assign_fractional_bond_orders. "
                f"Expected ToolkitRegistry or ToolkitWrapper. Got {type(toolkit_registry)}."
            )

    def _invalidate_cached_properties(self):
        """
        Indicate that the chemical entity has been altered.
        """
        # if hasattr(self, '_cached_properties'):
        #    delattr(self, '_cached_properties')
        self._conformers = None
        self._partial_charges = None
        self._propers = None
        self._impropers = None

        self._cached_smiles = None
        # TODO: Clear fractional bond orders
        self._rings = None

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
                atom.molecule_atom_index,
                atomic_number=atom.atomic_number,
                is_aromatic=atom.is_aromatic,
                stereochemistry=atom.stereochemistry,
                formal_charge=atom.formal_charge,
            )
            # G.add_node(atom.molecule_atom_index, attr_dict={'atomic_number': atom.atomic_number})
        for bond in self.bonds:
            G.add_edge(
                bond.atom1_index,
                bond.atom2_index,
                bond_order=bond.bond_order,
                is_aromatic=bond.is_aromatic,
                stereochemistry=bond.stereochemistry,
            )
            # G.add_edge(bond.atom1_index, bond.atom2_index, attr_dict={'order':bond.bond_order})

        return G

    def find_rotatable_bonds(
        self, ignore_functional_groups=None, toolkit_registry=GLOBAL_TOOLKIT_REGISTRY
    ):
        """
        Find all bonds classed as rotatable ignoring any matched to the ``ignore_functional_groups`` list.

        Parameters
        ----------
        ignore_functional_groups: optional, List[str], default=None,
            A list of bond SMARTS patterns to be ignored when finding rotatable bonds.

        toolkit_registry: openff.toolkit.utils.toolkits.ToolkitRegistry or openff.toolkit.utils.toolkits.ToolkitWrapper, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for SMARTS matching

        Returns
        -------
        bonds: List[openff.toolkit.topology.molecule.Bond]
            The list of openff.toolkit.topology.molecule.Bond instances which are rotatable.
        """

        # general rotatable bond smarts taken from RDKit
        # https://github.com/rdkit/rdkit/blob/1bf6ef3d65f5c7b06b56862b3fb9116a3839b229/rdkit/Chem/Lipinski.py#L47%3E
        rotatable_bond_smarts = "[!$(*#*)&!D1:1]-&!@[!$(*#*)&!D1:2]"

        # get all of the general matches
        general_matches = self.chemical_environment_matches(
            query=rotatable_bond_smarts, toolkit_registry=toolkit_registry
        )

        # this will give all forwards and backwards matches, so condense them down with this function
        def condense_matches(matches):
            condensed_matches = set()
            for m in matches:
                condensed_matches.add(tuple(sorted(m)))
            return condensed_matches

        general_bonds = condense_matches(general_matches)

        # now refine the list using the ignore groups
        if ignore_functional_groups is not None:
            matches_to_ignore = set()

            # make ignore_functional_groups an iterable object
            if isinstance(ignore_functional_groups, str):
                ignore_functional_groups = [ignore_functional_groups]
            else:
                try:
                    iter(ignore_functional_groups)
                except TypeError:
                    ignore_functional_groups = [ignore_functional_groups]

            # find the functional groups to remove
            for functional_group in ignore_functional_groups:
                # note I run the searches through this function so they have to be SMIRKS?
                ignore_matches = self.chemical_environment_matches(
                    query=functional_group, toolkit_registry=toolkit_registry
                )
                ignore_matches = condense_matches(ignore_matches)
                # add the new matches to the matches to ignore
                matches_to_ignore.update(ignore_matches)

            # now remove all the matches
            for match in matches_to_ignore:
                try:
                    general_bonds.remove(match)
                # if the key is not in the list, the ignore pattern was not valid
                except KeyError:
                    continue

        # gather a list of bond instances to return
        rotatable_bonds = [self.get_bond_between(*bond) for bond in general_bonds]
        return rotatable_bonds

    def _add_atom(
        self, atomic_number, formal_charge, is_aromatic, stereochemistry=None, name=None
    ):
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
            molecule=self,
        )
        self._atoms.append(atom)
        # self._particles.append(atom)
        self._invalidate_cached_properties()
        return self._atoms.index(atom)

    def _add_virtual_site(self, vsite, replace=False):
        replaced = False
        for i, existing_vsite in enumerate(self._virtual_sites):
            same_vsite = existing_vsite == vsite
            if same_vsite:
                if replace:
                    self._virtual_sites[i] = vsite
                    replaced = True
                    break
                else:
                    error_msg = (
                        "Attempted to add the new virtual site:\n{}\n"
                        + "to molecule: \n{}\nAnother vsite with the same type "
                        + "already exists and replace=False. Existing vsite "
                        + "is:\n{}\n"
                    ).format(vsite, self, existing_vsite)
                    raise Exception(error_msg)
        if not replaced:
            self._virtual_sites.append(vsite)
        return self._virtual_sites.index(vsite)

    def _add_bond_charge_virtual_site(self, atoms, distance, **kwargs):
        """
        Create a bond charge-type virtual site, in which the location of the charge is specified by the position of two
        atoms. This supports placement of a virtual site S along a vector between two specified atoms, e.g. to allow
        for a sigma hole for halogens or similar contexts. With positive values of the distance, the virtual site lies
        outside the first indexed atom.

        Parameters
        ----------
        atoms : list of openff.toolkit.topology.molecule.Atom objects of shape [N]
            The atoms defining the virtual site's position

        distance : :class:`simtk.unit.Quantity` of dimension [Length] wrapping a scalar

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
        orientations : list of int 2-tuples
            The orientations that should be used to create the virtual site.
            Each orientation corresponds to an individual virtual particle.

        Returns
        -------
        index : int
            The index of the newly-added virtual site in the molecule
        """

        replace = kwargs.pop("replace", False)

        vsite = BondChargeVirtualSite(atoms, distance, **kwargs)

        self._add_virtual_site(vsite, replace=replace)
        return self._virtual_sites.index(vsite)

    def _add_monovalent_lone_pair_virtual_site(
        self, atoms, distance, out_of_plane_angle, in_plane_angle, **kwargs
    ):
        """
        Create a bond charge-type virtual site, in which the location of the charge is specified by the position of
        three atoms.

        Parameters
        ----------
        atoms : list of three :class:`openff.toolkit.topology.molecule.Atom` objects
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
        orientations : list of int 3-tuples
            The orientations that should be used to create the virtual site.
            Each orientation corresponds to an individual virtual particle.

        Returns
        -------
        index : int
            The index of the newly-added virtual site in the molecule
        """

        replace = kwargs.pop("replace", False)

        vsite = MonovalentLonePairVirtualSite(
            atoms, distance, out_of_plane_angle, in_plane_angle, **kwargs
        )

        self._add_virtual_site(vsite, replace=replace)
        return self._virtual_sites.index(vsite)

    def _add_divalent_lone_pair_virtual_site(
        self, atoms, distance, out_of_plane_angle, **kwargs
    ):
        """
        Create a divalent lone pair-type virtual site, in which the location of the charge is specified by the position
        of three atoms.

        Parameters
        ----------
        atoms : list of three :class:`openff.toolkit.topology.molecule.Atom` objects
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
        orientations : list of int 3-tuples
            The orientations that should be used to create the virtual site.
            Each orientation corresponds to an individual virtual particle.

        Returns
        -------
        index : int
            The index of the newly-added virtual site in the molecule
        """

        replace = kwargs.pop("replace", False)

        vsite = DivalentLonePairVirtualSite(
            atoms, distance, out_of_plane_angle, **kwargs
        )

        self._add_virtual_site(vsite, replace=replace)
        return self._virtual_sites.index(vsite)

    def _add_trivalent_lone_pair_virtual_site(self, atoms, distance, **kwargs):
        """
        Create a trivalent lone pair-type virtual site, in which the location of the charge is specified by the position
         of four atoms.

        Parameters
        ----------
        atoms : list of 4 :class:`openff.toolkit.topology.molecule.Atom` objects
            The four atoms defining the virtual site's position

        distance : :class:`simtk.unit.Quantity` of dimension [Length] wrapping a scalar

        epsilon : float
            Epsilon term for VdW properties of virtual site. Default is None.
        sigma : float, default=None
            Sigma term for VdW properties of virtual site. Default is None.
        rmin_half : float
            Rmin_half term for VdW properties of virtual site. Default is None.
        name : string or None, default=None
            The name of this virtual site. Default is None.
        """

        replace = kwargs.pop("replace", False)

        vsite = TrivalentLonePairVirtualSite(atoms, distance, **kwargs)

        self._add_virtual_site(vsite, replace=replace)
        return self._virtual_sites.index(vsite)

    def _add_bond(
        self,
        atom1,
        atom2,
        bond_order,
        is_aromatic,
        stereochemistry=None,
        fractional_bond_order=None,
    ):
        """
        Add a bond between two specified atom indices

        Parameters
        ----------
        atom1 : int or openff.toolkit.topology.molecule.Atom
            Index of first atom or first atom
        atom2_index : int or openff.toolkit.topology.molecule.Atom
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
                "Invalid inputs to molecule._add_bond. Expected ints or Atoms. "
                "Received {} (type {}) and {} (type {}) ".format(
                    atom1, type(atom1), atom2, type(atom2)
                )
            )
        # TODO: Check to make sure bond does not already exist
        if atom1_atom.is_bonded_to(atom2_atom):
            raise Exception(
                "Bond already exists between {} and {}".format(atom1_atom, atom2_atom)
            )
        bond = Bond(
            atom1_atom,
            atom2_atom,
            bond_order,
            is_aromatic,
            stereochemistry=stereochemistry,
            fractional_bond_order=fractional_bond_order,
        )
        self._bonds.append(bond)
        self._invalidate_cached_properties()
        # TODO: This is a bad way to get bond index
        return self._bonds.index(bond)

    def _add_conformer(self, coordinates):
        """
        Add a conformation of the molecule

        Parameters
        ----------
        coordinates: simtk.unit.Quantity(np.array) with shape (n_atoms, 3) and dimension of distance
            Coordinates of the new conformer, with the first dimension of the array corresponding to the atom index in
            the Molecule's indexing system.

        Returns
        -------
        index: int
            The index of this conformer
        """
        new_conf = unit.Quantity(
            np.zeros(shape=(self.n_atoms, 3), dtype=float), unit.angstrom
        )
        if not (new_conf.shape == coordinates.shape):
            raise Exception(
                "molecule.add_conformer given input of the wrong shape: "
                "Given {}, expected {}".format(coordinates.shape, new_conf.shape)
            )

        try:
            new_conf[:] = coordinates
        except AttributeError as e:
            print(e)
            raise Exception(
                "Coordinates passed to Molecule._add_conformer without units. Ensure that coordinates are "
                "of type simtk.units.Quantity"
            )

        if self._conformers is None:
            # TODO should we checking that the exact same conformer is not in the list already?
            self._conformers = []
        self._conformers.append(new_conf)
        return len(self._conformers)

    @property
    def partial_charges(self):
        """
        Returns the partial charges (if present) on the molecule.

        Returns
        -------
        partial_charges : a simtk.unit.Quantity - wrapped numpy array [1 x n_atoms] or None
            The partial charges on this Molecule's atoms. Returns None if no charges have been specified.
        """
        return self._partial_charges

    @partial_charges.setter
    def partial_charges(self, charges):
        """
        Set the atomic partial charges for this molecule.

        Parameters
        ----------
        charges : None or a simtk.unit.Quantity - wrapped numpy array [1 x n_atoms]
            The partial charges to assign to the molecule. If not None, must be in units compatible with simtk.unit.elementary_charge

        """
        if charges is None:
            self._partial_charges = None
        else:
            assert hasattr(charges, "unit")
            assert unit.elementary_charge.is_compatible(charges.unit)
            assert charges.shape == (self.n_atoms,)

            charges_ec = charges.in_units_of(unit.elementary_charge)
            self._partial_charges = charges_ec

    @property
    def n_particles(self):
        """
        The number of Particle objects, which corresponds to how many positions must be used.
        """
        return len(self._atoms) + sum(
            vsite.n_particles for vsite in self._virtual_sites
        )

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
        return len(self._virtual_sites)

    @property
    def n_virtual_particles(self):
        """
        The number of VirtualParticle objects.
        """
        return sum(vsite.n_particles for vsite in self._virtual_sites)

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
        """int: number of possible improper torsions in the Molecule."""
        self._construct_torsions()
        return len(self._impropers)

    @property
    def n_rings(self):
        """Return the number of rings found in the Molecule

        Requires the RDKit to be installed.

        .. note ::

            For systems containing some special cases of connected rings, this
            function may not be well-behaved and may report a different number
            rings than expected. Some problematic cases include networks of many
            (5+) rings or bicyclic moieties (i.e. norbornane).

        """
        return len(self.rings)

    @property
    def particles(self):
        """
        Iterate over all Particle objects.
        """

        return self._atoms + [
            ptl for vsite in self._virtual_sites for ptl in vsite.particles
        ]

    @property
    def atoms(self):
        """
        Iterate over all Atom objects.
        """
        return self._atoms

    @property
    def conformers(self):
        """
        Returns the list of conformers for this molecule. This returns a list of simtk.unit.Quantity-wrapped numpy
        arrays, of shape (3 x n_atoms) and with dimensions of distance. The return value is the actual list of
        conformers, and changes to the contents affect the original FrozenMolecule.

        """
        return self._conformers

    @property
    def n_conformers(self):
        """
        Returns the number of conformers for this molecule.
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
        Iterate over all improper torsions in the molecule.

        .. todo ::
           * Do we need to return a ``Torsion`` object that collects information about fractional bond orders?

        Returns
        -------
        impropers : set of tuple
            An iterator of tuples, each containing the indices of atoms making
            up a possible improper torsion.

        See Also
        --------
        smirnoff_impropers, amber_impropers
        """
        self._construct_torsions()
        return self._impropers

    @property
    def smirnoff_impropers(self):
        """
        Iterate over improper torsions in the molecule, but only those with
        trivalent centers, reporting the central atom second in each improper.

        Note that it's possible that a trivalent center will not have an improper assigned.
        This will depend on the force field that is used.

        Also note that this will return 6 possible atom orderings around each improper
        center. In current SMIRNOFF parameterization, three of these six
        orderings will be used for the actual assignment of the improper term
        and measurement of the angles. These three orderings capture the three unique
        angles that could be calculated around the improper center, therefore the sum
        of these three terms will always return a consistent energy.

        The exact three orderings that will be applied during parameterization can not be
        determined in this method, since it requires sorting the particle indices, and
        those indices may change when this molecule is added to a Topology.

        For more details on the use of three-fold ('trefoil') impropers, see
        https://open-forcefield-toolkit.readthedocs.io/en/latest/smirnoff.html#impropertorsions

        Returns
        -------
        impropers : set of tuple
            An iterator of tuples, each containing the indices of atoms making
            up a possible improper torsion. The central atom is listed second
            in each tuple.

        See Also
        --------
        impropers, amber_impropers

        """
        # TODO: Replace with non-cheminformatics-toolkit method
        #       (ie. just looping over all atoms and finding ones that have 3 bonds?)

        smirnoff_improper_smarts = "[*:1]~[X3:2](~[*:3])~[*:4]"
        improper_idxs = self.chemical_environment_matches(smirnoff_improper_smarts)
        smirnoff_impropers = {
            tuple(self.atoms[idx] for idx in imp) for imp in improper_idxs
        }
        return smirnoff_impropers

    @property
    def amber_impropers(self):
        """
        Iterate over improper torsions in the molecule, but only those with
        trivalent centers, reporting the central atom first in each improper.

        Note that it's possible that a trivalent center will not have an improper assigned.
        This will depend on the force field that is used.

        Also note that this will return 6 possible atom orderings around each improper
        center. In current AMBER parameterization, one of these six
        orderings will be used for the actual assignment of the improper term
        and measurement of the angle. This method does not encode the logic to
        determine which of the six orderings AMBER would use.

        Returns
        -------
        impropers : set of tuple
            An iterator of tuples, each containing the indices of atoms making
            up a possible improper torsion. The central atom is listed first in
            each tuple.

        See Also
        --------
        impropers, smirnoff_impropers

        """
        # TODO: Replace with non-cheminformatics-toolkit method
        #       (ie. just looping over all atoms and finding ones that have 3 bonds?)
        amber_improper_smarts = "[X3:1](~[*:2])(~[*:3])~[*:4]"
        improper_idxs = self.chemical_environment_matches(amber_improper_smarts)
        amber_impropers = {
            tuple(self.atoms[idx] for idx in imp) for imp in improper_idxs
        }
        return amber_impropers

    def _nth_degree_neighbors(self, n_degrees):
        import networkx as nx

        mol_graph = self.to_networkx()

        for node_i in mol_graph.nodes:
            for node_j in mol_graph.nodes:
                if node_i == node_j:
                    continue

                path_length = nx.shortest_path_length(mol_graph, node_i, node_j)

                if path_length == n_degrees:
                    if node_i > node_j:
                        continue
                    yield (self.atoms[node_i], self.atoms[node_j])

    def nth_degree_neighbors(self, n_degrees):
        """
        Return canonicalized pairs of atoms whose shortest separation is `exactly` n bonds.
        Only pairs with increasing atom indices are returned.

        Parameters
        ----------
        n: int
            The number of bonds separating atoms in each pair

        Returns
        -------
        neighbors: iterator of tuple of Atom
            Tuples (len 2) of atom that are separated by ``n`` bonds.

        Notes
        -----

        The criteria used here relies on minimum distances; when there are multiple valid
        paths between atoms, such as atoms in rings, the shortest path is considered.
        For example, two atoms in "meta" positions with respect to each other in a benzene
        are separated by two paths, one length 2 bonds and the other length 4 bonds. This
        function would consider them to be 2 apart and would not include them if ``n=4`` was
        passed.

        """
        if n_degrees <= 0:
            raise ValueError(
                "Cannot consider neighbors separated by 0 or fewer atoms. Asked to consider "
                f"path lengths of {n_degrees}."
            )
        else:
            return self._nth_degree_neighbors(n_degrees=n_degrees)

    @property
    def total_charge(self):
        """
        Return the total charge on the molecule
        """
        charge_sum = 0.0 * unit.elementary_charge
        for atom in self.atoms:
            charge_sum += atom.formal_charge
        return charge_sum

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
            self._name = ""
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

    @property
    def hill_formula(self):
        """
        Get the Hill formula of the molecule
        """
        return Molecule.to_hill_formula(self)

    @staticmethod
    def to_hill_formula(molecule):
        """
        Generate the Hill formula from either a :class:`FrozenMolecule`, :class:`TopologyMolecule` or
        ``nx.Graph()`` of the molecule

        .. TODO: Set up Intersphinx for nx.Graph(), above

        Parameters
        -----------
        molecule : FrozenMolecule, TopologyMolecule or nx.Graph()

        Returns
        ----------
        formula : the Hill formula of the molecule

        Raises
        -----------
        NotImplementedError : if the molecule is not of one of the specified types.
        """

        import networkx as nx

        from openff.toolkit.topology import TopologyMolecule

        # check for networkx then assuming we have a Molecule or TopologyMolecule instance just try and
        # extract the info. Note we do not type check the TopologyMolecule due to cyclic dependencies
        if isinstance(molecule, nx.Graph):
            atom_nums = list(
                dict(molecule.nodes(data="atomic_number", default=1)).values()
            )

        elif isinstance(molecule, TopologyMolecule):
            atom_nums = [atom.atomic_number for atom in molecule.atoms]

        elif isinstance(molecule, FrozenMolecule):
            atom_nums = [atom.atomic_number for atom in molecule.atoms]

        else:
            raise NotImplementedError(
                f"The input type {type(molecule)} is not supported,"
                f"please supply an openff.toolkit.topology.molecule.Molecule,"
                f"openff.toolkit.topology.topology.TopologyMolecule or networkx representaion "
                f"of the molecule."
            )

        # make a correct hill formula representation following this guide
        # https://en.wikipedia.org/wiki/Chemical_formula#Hill_system

        # create the counter dictionary using chemical symbols
        from collections import Counter

        atom_symbol_counts = Counter(
            Element.getByAtomicNumber(atom_num).symbol for atom_num in atom_nums
        )

        formula = []
        # Check for C and H first, to make a correct hill formula
        for el in ["C", "H"]:
            if el in atom_symbol_counts:
                count = atom_symbol_counts.pop(el)
                formula.append(el)
                if count > 1:
                    formula.append(str(count))

        # now get the rest of the elements in alphabetical ordering
        for el in sorted(atom_symbol_counts.keys()):
            count = atom_symbol_counts.pop(el)
            formula.append(el)
            if count > 1:
                formula.append(str(count))

        return "".join(formula)

    def chemical_environment_matches(
        self,
        query,
        unique=False,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
    ):
        """Retrieve all matches for a given chemical environment query.

        Parameters
        ----------
        query : str or ChemicalEnvironment
            SMARTS string (with one or more tagged atoms) or ``ChemicalEnvironment`` query.
            Query will internally be resolved to SMIRKS using ``query.asSMIRKS()`` if it has an ``.asSMIRKS`` method.
        toolkit_registry : openff.toolkit.utils.toolkits.ToolkitRegistry or openff.toolkit.utils.toolkits.ToolkitWrapper, optional, default=GLOBAL_TOOLKIT_REGISTRY
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for chemical environment matches


        Returns
        -------
        matches : list of atom index tuples
            A list of tuples, containing the indices of the matching atoms.

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
        if hasattr(query, "smirks"):
            smirks = query.smirks
        elif type(query) == str:
            smirks = query
        else:
            raise ValueError("'query' must be either a string or a ChemicalEnvironment")

        # Use specified cheminformatics toolkit to determine matches with specified aromaticity model
        # TODO: Simplify this by requiring a toolkit registry for the molecule?
        # TODO: Do we have to pass along an aromaticity model?
        if isinstance(toolkit_registry, ToolkitRegistry):
            matches = toolkit_registry.call(
                "find_smarts_matches",
                self,
                smirks,
                unique=unique,
            )
        elif isinstance(toolkit_registry, ToolkitWrapper):
            matches = toolkit_registry.find_smarts_matches(
                self,
                smirks,
                unique=unique,
            )
        else:
            raise InvalidToolkitRegistryError(
                "'toolkit_registry' must be either a ToolkitRegistry or a ToolkitWrapper"
            )

        return matches

    @classmethod
    def from_iupac(
        cls,
        iupac_name,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
        allow_undefined_stereo=False,
        **kwargs,
    ):
        """Generate a molecule from IUPAC or common name

        Parameters
        ----------
        iupac_name : str
            IUPAC name of molecule to be generated
        toolkit_registry : openff.toolkit.utils.toolkits.ToolkitRegistry or openff.toolkit.utils.toolkits.ToolkitWrapper, optional, default=GLOBAL_TOOLKIT_REGISTRY
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for chemical environment matches
        allow_undefined_stereo : bool, default=False
            If false, raises an exception if molecule contains undefined stereochemistry.

        Returns
        -------
        molecule : Molecule
            The resulting molecule with position

        .. note :: This method requires the OpenEye toolkit to be installed.

        Examples
        --------

        Create a molecule from an IUPAC name

        >>> molecule = Molecule.from_iupac('4-[(4-methylpiperazin-1-yl)methyl]-N-(4-methyl-3-{[4-(pyridin-3-yl)pyrimidin-2-yl]amino}phenyl)benzamide')

        Create a molecule from a common name

        >>> molecule = Molecule.from_iupac('imatinib')

        """
        if isinstance(toolkit_registry, ToolkitRegistry):
            molecule = toolkit_registry.call(
                "from_iupac",
                iupac_name,
                allow_undefined_stereo=allow_undefined_stereo,
                _cls=cls,
                **kwargs,
            )
        elif isinstance(toolkit_registry, ToolkitWrapper):
            toolkit = toolkit_registry
            molecule = toolkit.from_iupac(
                iupac_name,
                allow_undefined_stereo=allow_undefined_stereo,
                _cls=cls ** kwargs,
            )
        else:
            raise Exception(
                "Invalid toolkit_registry passed to from_iupac. Expected ToolkitRegistry or ToolkitWrapper. Got  {}".format(
                    type(toolkit_registry)
                )
            )

        return molecule

    def to_iupac(self, toolkit_registry=GLOBAL_TOOLKIT_REGISTRY):
        """Generate IUPAC name from Molecule

        Returns
        -------
        iupac_name : str
            IUPAC name of the molecule

        .. note :: This method requires the OpenEye toolkit to be installed.

        Examples
        --------

        >>> from openff.toolkit.utils import get_data_file_path
        >>> sdf_filepath = get_data_file_path('molecules/ethanol.sdf')
        >>> molecule = Molecule(sdf_filepath)
        >>> iupac_name = molecule.to_iupac()

        """
        if isinstance(toolkit_registry, ToolkitRegistry):
            to_iupac_method = toolkit_registry.resolve("to_iupac")
        elif isinstance(toolkit_registry, ToolkitWrapper):
            to_iupac_method = toolkit_registry.to_iupac
        else:
            raise Exception(
                "Invalid toolkit_registry passed to to_iupac. Expected ToolkitRegistry or ToolkitWrapper. Got  {}".format(
                    type(toolkit_registry)
                )
            )

        # TODO: Can `to_iupac` fail if given a well-behaved OFFMol/OEMol?
        result = to_iupac_method(self)
        return result

    @classmethod
    def from_topology(cls, topology):
        """Return a Molecule representation of an OpenFF Topology containing a single Molecule object.

        Parameters
        ----------
        topology : openff.toolkit.topology.Topology
            The :class:`Topology` object containing a single :class:`Molecule` object.
            Note that OpenMM and MDTraj ``Topology`` objects are not supported.

        Returns
        -------
        molecule : openff.toolkit.topology.Molecule
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
        # TODO: Ensure we are dealing with an OpenFF Topology object
        if topology.n_topology_molecules != 1:
            raise ValueError("Topology must contain exactly one molecule")
        molecule = [i for i in topology.reference_molecules][0]
        return cls(molecule)

    def to_topology(self):
        """
        Return an OpenFF Topology representation containing one copy of this molecule

        Returns
        -------
        topology : openff.toolkit.topology.Topology
            A Topology representation of this molecule

        Examples
        --------

        >>> molecule = Molecule.from_iupac('imatinib')
        >>> topology = molecule.to_topology()

        """
        from openff.toolkit.topology import Topology

        return Topology.from_molecules(self)

    @classmethod
    def from_file(
        cls,
        file_path,
        file_format=None,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
        allow_undefined_stereo=False,
    ):
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
        toolkit_registry : openff.toolkit.utils.toolkits.ToolkitRegistry or openff.toolkit.utils.toolkits.ToolkitWrapper,
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
        >>> from openff.toolkit.tests.utils import get_monomer_mol2_file_path
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
            if file_path[-3:] == ".gz":
                file_format = file_path.split(".")[-2]
            else:
                file_format = file_path.split(".")[-1]
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
                    query_toolkit.toolkit_name
                ] = query_toolkit.toolkit_file_read_formats
            if toolkit is None:
                msg = (
                    f"No toolkits in registry can read file {file_path} (format {file_format}). Supported "
                    f"formats in the provided ToolkitRegistry are {supported_read_formats}. "
                )
                # Per issue #407, not allowing RDKit to read mol2 has confused a lot of people. Here we add text
                # to the error message that will hopefully reduce this confusion.
                if file_format == "MOL2" and RDKitToolkitWrapper.is_available():
                    msg += (
                        f"RDKit does not fully support input of molecules from mol2 format unless they "
                        f"have Corina atom types, and this is not common in the simulation community. For this "
                        f"reason, the Open Force Field Toolkit does not use "
                        f"RDKit to read .mol2. Consider reading from SDF instead. If you would like to attempt "
                        f"to use RDKit to read mol2 anyway, you can load the molecule of interest into an RDKit "
                        f"molecule and use openff.toolkit.topology.Molecule.from_rdkit, but we do not recommend this."
                    )
                elif file_format == "PDB" and RDKitToolkitWrapper.is_available():
                    msg += (
                        "RDKit can not safely read PDBs on their own. Information about bond order and aromaticity "
                        "is likely to be lost. PDBs can be used along with a valid smiles string with RDKit using "
                        "the constructor Molecule.from_pdb_and_smiles(file_path, smiles)"
                    )
                raise NotImplementedError(msg)

        elif isinstance(toolkit_registry, ToolkitWrapper):
            # TODO: Encapsulate this logic in ToolkitWrapper?
            toolkit = toolkit_registry
            if file_format not in toolkit.toolkit_file_read_formats:
                msg = (
                    f"Toolkit {toolkit.toolkit_name} can not read file {file_path} (format {file_format}). Supported "
                    f"formats for this toolkit are {toolkit.toolkit_file_read_formats}."
                )
                if toolkit.toolkit_name == "The RDKit" and file_format == "PDB":
                    msg += (
                        "RDKit can however read PDBs with a valid smiles string using the "
                        "Molecule.from_pdb_and_smiles(file_path, smiles) constructor"
                    )
                raise NotImplementedError(msg)
        else:
            raise InvalidToolkitRegistryError(
                "'toolkit_registry' must be either a ToolkitRegistry or a ToolkitWrapper"
            )

        mols = list()

        if isinstance(file_path, str):
            mols = toolkit.from_file(
                file_path,
                file_format=file_format,
                allow_undefined_stereo=allow_undefined_stereo,
                _cls=cls,
            )
        elif hasattr(file_path, "read"):
            file_obj = file_path
            mols = toolkit.from_file_obj(
                file_obj,
                file_format=file_format,
                allow_undefined_stereo=allow_undefined_stereo,
                _cls=cls,
            )

        if len(mols) == 0:
            raise Exception("Unable to read molecule from file: {}".format(file_path))
        elif len(mols) == 1:
            return mols[0]

        return mols

    def _to_xyz_file(self, file_path):
        """
        Write the current molecule and its conformers to a multiframe xyz file, if the molecule
        has no current coordinates all atoms will be set to 0,0,0 in keeping with the behaviour of the
        backend toolkits.

        Information on the type of XYZ file written can be found here <http://openbabel.org/wiki/XYZ_(format)>.

        Parameters
        ----------
        file_path : str or file-like object
            A file-like object or the path to the file to be written.
        """

        # If we do not have a conformer make one with all zeros
        if self.n_conformers == 0:
            conformers = [
                unit.Quantity(np.zeros((self.n_atoms, 3), dtype=float), unit.angstrom)
            ]

        else:
            conformers = self._conformers

        if len(conformers) == 1:
            end = ""
            title = (
                lambda frame: f'{self.name if self.name != "" else self.hill_formula}{frame}\n'
            )
        else:
            end = 1
            title = (
                lambda frame: f'{self.name if self.name != "" else self.hill_formula} Frame {frame}\n'
            )

        # check if we have a file path or an open file object
        if isinstance(file_path, str):
            xyz_data = open(file_path, "w")
        else:
            xyz_data = file_path

        # add the data to the xyz_data list
        for i, geometry in enumerate(conformers, 1):
            xyz_data.write(f"{self.n_atoms}\n" + title(end))
            for j, atom_coords in enumerate(geometry.in_units_of(unit.angstrom)):
                x, y, z = atom_coords._value
                xyz_data.write(
                    f"{self.atoms[j].element.symbol}       {x: .10f}   {y: .10f}   {z: .10f}\n"
                )

            # now we up the frame count
            end = i + 1

        # now close the file
        xyz_data.close()

    def to_file(self, file_path, file_format, toolkit_registry=GLOBAL_TOOLKIT_REGISTRY):
        """Write the current molecule to a file or file-like object

        Parameters
        ----------
        file_path : str or file-like object
            A file-like object or the path to the file to be written.
        file_format : str
            Format specifier, one of ['MOL2', 'MOL2H', 'SDF', 'PDB', 'SMI', 'CAN', 'TDT']
            Note that not all toolkits support all formats
        toolkit_registry : openff.toolkit.utils.toolkits.ToolkitRegistry or openff.toolkit.utils.toolkits.ToolkitWrapper,
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
            raise InvalidToolkitRegistryError(
                "'toolkit_registry' must be either a ToolkitRegistry or a ToolkitWrapper"
            )

        file_format = file_format.upper()
        # check if xyz, use the toolkit independent method.
        if file_format == "XYZ":
            return self._to_xyz_file(file_path=file_path)

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
                    toolkit.toolkit_name
                ] = toolkit.toolkit_file_write_formats
            raise ValueError(
                "The requested file format ({}) is not available from any of the installed toolkits "
                "(supported formats: {})".format(file_format, supported_formats)
            )

        # Write file
        if type(file_path) == str:
            # Open file for writing
            toolkit.to_file(self, file_path, file_format)
        else:
            toolkit.to_file_obj(self, file_path, file_format)

    def enumerate_tautomers(
        self, max_states=20, toolkit_registry=GLOBAL_TOOLKIT_REGISTRY
    ):
        """
        Enumerate the possible tautomers of the current molecule

        Parameters
        ----------
        max_states: int optional, default=20
            The maximum amount of molecules that should be returned

        toolkit_registry: openff.toolkit.utils.toolkits.ToolkitRegistry or openff.toolkit.utils.toolkits.ToolkitWrapper, default=GLOBAL_TOOLKIT_REGISTRY
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use to enumerate the tautomers.

        Returns
        -------
        molecules: List[openff.toolkit.topology.Molecule]
            A list of openff.toolkit.topology.Molecule instances not including the input molecule.
        """

        if isinstance(toolkit_registry, ToolkitRegistry):
            molecules = toolkit_registry.call(
                "enumerate_tautomers", molecule=self, max_states=max_states
            )

        elif isinstance(toolkit_registry, ToolkitWrapper):
            molecules = toolkit_registry.enumerate_tautomers(
                self, max_states=max_states
            )

        else:
            raise InvalidToolkitRegistryError(
                "'toolkit_registry' must be either a ToolkitRegistry or a ToolkitWrapper"
            )

        return molecules

    def enumerate_stereoisomers(
        self,
        undefined_only=False,
        max_isomers=20,
        rationalise=True,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
    ):
        """
        Enumerate the stereocenters and bonds of the current molecule.

        Parameters
        ----------
        undefined_only: bool optional, default=False
            If we should enumerate all stereocenters and bonds or only those with undefined stereochemistry

        max_isomers: int optional, default=20
            The maximum amount of molecules that should be returned

        rationalise: bool optional, default=True
            If we should try to build and rationalise the molecule to ensure it can exist

        toolkit_registry: openff.toolkit.utils.toolkits.ToolkitRegistry or openff.toolkit.utils.toolkits.ToolkitWrapper, default=GLOBAL_TOOLKIT_REGISTRY
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use to enumerate the stereoisomers.

        Returns
        --------
        molecules: List[openff.toolkit.topology.Molecule]
            A list of :class:`Molecule` instances not including the input molecule.

        """

        if isinstance(toolkit_registry, ToolkitRegistry):
            molecules = toolkit_registry.call(
                "enumerate_stereoisomers",
                molecule=self,
                undefined_only=undefined_only,
                max_isomers=max_isomers,
                rationalise=rationalise,
            )

        elif isinstance(toolkit_registry, ToolkitWrapper):
            molecules = toolkit_registry.enumerate_stereoisomers(
                self,
                undefined_only=undefined_only,
                max_isomers=max_isomers,
                rationalise=rationalise,
            )

        else:
            raise InvalidToolkitRegistryError(
                "'toolkit_registry' must be either a ToolkitRegistry or a ToolkitWrapper"
            )

        return molecules

    @OpenEyeToolkitWrapper.requires_toolkit()
    def enumerate_protomers(self, max_states=10):
        """
        Enumerate the formal charges of a molecule to generate different protomoers.

        Parameters
        ----------
        max_states: int optional, default=10,
            The maximum number of protomer states to be returned.

        Returns
        -------
        molecules: List[openff.toolkit.topology.Molecule],
            A list of the protomers of the input molecules not including the input.
        """

        toolkit = OpenEyeToolkitWrapper()
        molecules = toolkit.enumerate_protomers(molecule=self, max_states=max_states)

        return molecules

    @classmethod
    @RDKitToolkitWrapper.requires_toolkit()
    def from_rdkit(
        cls, rdmol, allow_undefined_stereo=False, hydrogens_are_explicit=False
    ):
        """
        Create a Molecule from an RDKit molecule.

        Requires the RDKit to be installed.

        Parameters
        ----------
        rdmol : rkit.RDMol
            An RDKit molecule
        allow_undefined_stereo : bool, default=False
            If ``False``, raises an exception if ``rdmol`` contains undefined stereochemistry.
        hydrogens_are_explicit : bool, default=False
            If ``False``, RDKit will perform hydrogen addition using ``Chem.AddHs``

        Returns
        -------
        molecule : openff.toolkit.topology.Molecule
            An OpenFF molecule

        Examples
        --------

        Create a molecule from an RDKit molecule

        >>> from rdkit import Chem
        >>> from openff.toolkit.tests.utils import get_data_file_path
        >>> rdmol = Chem.MolFromMolFile(get_data_file_path('systems/monomers/ethanol.sdf'))
        >>> molecule = Molecule.from_rdkit(rdmol)

        """
        toolkit = RDKitToolkitWrapper()
        molecule = toolkit.from_rdkit(
            rdmol,
            allow_undefined_stereo=allow_undefined_stereo,
            hydrogens_are_explicit=hydrogens_are_explicit,
            _cls=cls,
        )
        return molecule

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
        rdmol : rdkit.RDMol
            An RDKit molecule

        Examples
        --------

        Convert a molecule to RDKit

        >>> from openff.toolkit.utils import get_data_file_path
        >>> sdf_filepath = get_data_file_path('molecules/ethanol.sdf')
        >>> molecule = Molecule(sdf_filepath)
        >>> rdmol = molecule.to_rdkit()

        """
        toolkit = RDKitToolkitWrapper()
        return toolkit.to_rdkit(self, aromaticity_model=aromaticity_model)

    @classmethod
    @OpenEyeToolkitWrapper.requires_toolkit()
    def from_openeye(cls, oemol, allow_undefined_stereo=False):
        """
        Create a ``Molecule`` from an OpenEye molecule.

        Requires the OpenEye toolkit to be installed.

        Parameters
        ----------
        oemol : openeye.oechem.OEMol
            An OpenEye molecule
        allow_undefined_stereo : bool, default=False
            If ``False``, raises an exception if oemol contains undefined stereochemistry.

        Returns
        -------
        molecule : openff.toolkit.topology.Molecule
            An OpenFF molecule

        Examples
        --------

        Create a ``Molecule`` from an OpenEye OEMol

        >>> from openeye import oechem
        >>> from openff.toolkit.tests.utils import get_data_file_path
        >>> ifs = oechem.oemolistream(get_data_file_path('systems/monomers/ethanol.mol2'))
        >>> oemols = list(ifs.GetOEGraphMols())
        >>> molecule = Molecule.from_openeye(oemols[0])

        """
        toolkit = OpenEyeToolkitWrapper()
        molecule = toolkit.from_openeye(
            oemol, allow_undefined_stereo=allow_undefined_stereo, _cls=cls
        )
        return molecule

    @requires_package("qcelemental")
    def to_qcschema(self, multiplicity=1, conformer=0, extras=None):
        """
        Create a QCElemental Molecule.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        multiplicity : int, default=1,
            The multiplicity of the molecule;
            sets ``molecular_multiplicity`` field for QCElemental Molecule.

        conformer : int, default=0,
            The index of the conformer to use for the QCElemental Molecule geometry.

        extras : dict, default=None
            A dictionary that should be included in the ``extras`` field on the QCElemental Molecule.
            This can be used to include extra information, such as a smiles representation.

        Returns
        ---------
        qcelemental.models.Molecule
            A validated QCElemental Molecule.

        Examples
        --------

        Create a QCElemental Molecule:

        >>> import qcelemental as qcel
        >>> mol = Molecule.from_smiles('CC')
        >>> mol.generate_conformers(n_conformers=1)
        >>> qcemol = mol.to_qcschema()

        Raises
        --------
        MissingDependencyError
            If qcelemental is not installed, the qcschema can not be validated.
        InvalidConformerError
            No conformer found at the given index.

        """

        import qcelemental as qcel

        # get/ check the geometry
        try:
            geometry = self.conformers[conformer].in_units_of(unit.bohr)
        except (IndexError, TypeError):
            raise InvalidConformerError(
                "The molecule must have a conformation to produce a valid qcschema; "
                f"no conformer was found at index {conformer}."
            )

        # Gather the required qcschema data
        charge = self.total_charge / unit.elementary_charge
        connectivity = [
            (bond.atom1_index, bond.atom2_index, bond.bond_order) for bond in self.bonds
        ]
        symbols = [
            Element.getByAtomicNumber(atom.atomic_number).symbol for atom in self.atoms
        ]
        if extras is not None:
            extras[
                "canonical_isomeric_explicit_hydrogen_mapped_smiles"
            ] = self.to_smiles(mapped=True)
        else:
            extras = {
                "canonical_isomeric_explicit_hydrogen_mapped_smiles": self.to_smiles(
                    mapped=True
                )
            }

        schema_dict = {
            "symbols": symbols,
            "geometry": geometry,
            "connectivity": connectivity,
            "molecular_charge": charge,
            "molecular_multiplicity": multiplicity,
            "extras": extras,
        }

        return qcel.models.Molecule.from_data(schema_dict, validate=True)

    @classmethod
    def from_mapped_smiles(
        cls,
        mapped_smiles,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
        allow_undefined_stereo=False,
    ):
        """
        Create an :class:`Molecule` from a mapped SMILES made with cmiles.
        The molecule will be in the order of the indexing in the mapped smiles string.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        mapped_smiles: str,
            A CMILES-style mapped smiles string with explicit hydrogens.

        toolkit_registry : openff.toolkit.utils.toolkits.ToolkitRegistry or openff.toolkit.utils.toolkits.ToolkitWrapper, optional
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for SMILES-to-molecule conversion

        allow_undefined_stereo : bool, default=False
            If false, raises an exception if oemol contains undefined stereochemistry.

        Returns
        ----------
        offmol : openff.toolkit.topology.molecule.Molecule
            An OpenFF molecule instance.

        Raises
        --------
        SmilesParsingError
            If the given SMILES had no indexing picked up by the toolkits.
        """

        # create the molecule from the smiles and check we have the right number of indexes
        # in the mapped SMILES
        offmol = cls.from_smiles(
            mapped_smiles,
            hydrogens_are_explicit=True,
            toolkit_registry=toolkit_registry,
            allow_undefined_stereo=allow_undefined_stereo,
        )

        # check we found some mapping and remove it as we do not want to expose atom maps
        try:
            mapping = offmol._properties.pop("atom_map")
        except KeyError:
            raise SmilesParsingError(
                "The given SMILES has no indexing, please generate a valid explicit hydrogen "
                "mapped SMILES using cmiles."
            )

        if len(mapping) != offmol.n_atoms:
            raise SmilesParsingError(
                "The mapped smiles does not contain enough indexes to remap the molecule."
            )

        # remap the molecule using the atom map found in the smiles
        # the order is mapping = Dict[current_index: new_index]
        # first renumber the mapping dict indexed from 0, currently from 1 as 0 indicates no mapping in toolkits
        adjusted_mapping = dict((current, new - 1) for current, new in mapping.items())

        return offmol.remap(adjusted_mapping, current_to_new=True)

    @classmethod
    @requires_package("qcelemental")
    def from_qcschema(
        cls,
        qca_record,
        client=None,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
        allow_undefined_stereo=False,
    ):
        """
        Create a Molecule from a QCArchive molecule record or dataset entry
        based on attached cmiles information.

        For a molecule record, a conformer will be set from its geometry.

        For a dataset entry, if a corresponding client instance is provided,
        the starting geometry for that entry will be used as a conformer.

        A QCElemental Molecule produced from ``Molecule.to_qcschema`` can be round-tripped
        through this method to produce a new, valid Molecule.

        Parameters
        ----------
        qca_record : dict
            A QCArchive molecule record or dataset entry.

        client : optional, default=None,
            A qcportal.FractalClient instance to use for fetching an initial geometry.
            Only used if ``qca_record`` is a dataset entry.

        toolkit_registry : openff.toolkit.utils.toolkits.ToolkitRegistry or openff.toolkit.utils.toolkits.ToolkitWrapper, optional
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for SMILES-to-molecule conversion

        allow_undefined_stereo : bool, default=False
            If false, raises an exception if qca_record contains undefined stereochemistry.

        Returns
        -------
        molecule : openff.toolkit.topology.Molecule
            An OpenFF molecule instance.

        Examples
        --------
        Get Molecule from a QCArchive molecule record:

        >>> from qcportal import FractalClient
        >>> client = FractalClient()
        >>> offmol = Molecule.from_qcschema(client.query_molecules(molecular_formula="C16H20N3O5")[0])

        Get Molecule from a QCArchive optimization entry:

        >>> from qcportal import FractalClient
        >>> client = FractalClient()
        >>> optds = client.get_collection("OptimizationDataset",
                                          "SMIRNOFF Coverage Set 1")
        >>> offmol = Molecule.from_qcschema(optds.get_entry('coc(o)oc-0'))

        Same as above, but with conformer(s) from initial molecule(s) by providing client to database:

        >>> offmol = Molecule.from_qcschema(optds.get_entry('coc(o)oc-0'), client=client)

        Raises
        -------
        AttributeError
            - If the record dict can not be made from ``qca_record``.
            - If a ``client`` is passed and it could not retrieve the initial molecule.

        KeyError
            If the dict does not contain the ``canonical_isomeric_explicit_hydrogen_mapped_smiles``.

        InvalidConformerError
            Silent error, if the conformer could not be attached.
        """

        # We can accept the Dataset entry record or the dict with JSON encoding
        # lets get it all in the dict rep
        if not isinstance(qca_record, dict):
            try:
                qca_record = qca_record.dict(encoding="json")
            except AttributeError:
                raise AttributeError(
                    "The object passed could not be converted to a dict with json encoding"
                )

        # identify if this is a dataset entry
        if "attributes" in qca_record:
            mapped_smiles = qca_record["attributes"][
                "canonical_isomeric_explicit_hydrogen_mapped_smiles"
            ]
            if client is not None:
                # try and find the initial molecule conformations and attach them
                # collect the input molecules
                try:
                    input_mols = client.query_molecules(
                        id=qca_record["initial_molecules"]
                    )
                except KeyError:
                    # this must be an optimisation record
                    input_mols = client.query_molecules(
                        id=qca_record["initial_molecule"]
                    )
                except AttributeError:
                    raise AttributeError(
                        "The provided client can not query molecules, make sure it is an instance of"
                        "qcportal.client.FractalClient() with the correct address."
                    )
            else:
                input_mols = []

        # identify if this is a molecule record
        elif "extras" in qca_record:
            mapped_smiles = qca_record["extras"][
                "canonical_isomeric_explicit_hydrogen_mapped_smiles"
            ]
            input_mols = [qca_record]
        else:
            raise KeyError(
                "The record must contain the hydrogen mapped smiles to be safely made from the archive. "
                "It is not present in either 'attributes' or 'extras' on the provided `qca_record`"
            )

        # make a new molecule that has been reordered to match the cmiles mapping
        offmol = cls.from_mapped_smiles(
            mapped_smiles,
            toolkit_registry=toolkit_registry,
            allow_undefined_stereo=allow_undefined_stereo,
        )

        # now for each molecule convert and attach the input geometry
        initial_ids = {}
        for molecule in input_mols:
            if not isinstance(molecule, dict):
                mol = molecule.dict(encoding="json")
            else:
                mol = molecule

            geometry = unit.Quantity(
                np.array(mol["geometry"], float).reshape(-1, 3), unit.bohr
            )
            try:
                offmol._add_conformer(geometry.in_units_of(unit.angstrom))
                # in case this molecule didn't come from a server at all
                if "id" in mol:
                    initial_ids[mol["id"]] = offmol.n_conformers - 1
            except InvalidConformerError:
                print(
                    "Invalid conformer for this molecule, the geometry could not be attached."
                )

        # attach a dict that has the initial molecule ids and the number of the conformer it is stored in
        # if it's empty, don't bother
        if initial_ids:
            offmol._properties["initial_molecules"] = initial_ids

        return offmol

    @classmethod
    @RDKitToolkitWrapper.requires_toolkit()
    def from_pdb_and_smiles(cls, file_path, smiles, allow_undefined_stereo=False):
        """
        Create a Molecule from a pdb file and a SMILES string using RDKit.

        Requires RDKit to be installed.

        .. warning :: This API is experimental and subject to change.

        The molecule is created and sanitised based on the SMILES string, we then find a mapping
        between this molecule and one from the PDB based only on atomic number and connections.
        The SMILES molecule is then reindexed to match the PDB, the conformer is attached, and the
        molecule returned.

        Note that any stereochemistry in the molecule is set by the SMILES, and not the coordinates
        of the PDB.

        Parameters
        ----------
        file_path: str
            PDB file path
        smiles : str
            a valid smiles string for the pdb, used for stereochemistry, formal charges, and bond order
        allow_undefined_stereo : bool, default=False
            If false, raises an exception if SMILES contains undefined stereochemistry.

        Returns
        --------
        molecule : openff.toolkit.Molecule
            An OFFMol instance with ordering the same as used in the PDB file.

        Raises
        ------
        InvalidConformerError
            If the SMILES and PDB molecules are not isomorphic.
        """

        toolkit = RDKitToolkitWrapper()
        return toolkit.from_pdb_and_smiles(
            file_path, smiles, allow_undefined_stereo, _cls=cls
        )

    def canonical_order_atoms(self, toolkit_registry=GLOBAL_TOOLKIT_REGISTRY):
        """
        Canonical order the atoms in a copy of the molecule using a toolkit, returns a new copy.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        toolkit_registry : openff.toolkit.utils.toolkits.ToolkitRegistry or openff.toolkit.utils.toolkits.ToolkitWrapper, optional
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for SMILES-to-molecule conversion

         Returns
        -------
        molecule : openff.toolkit.topology.Molecule
            An new OpenFF style molecule with atoms in the canonical order.
        """

        if isinstance(toolkit_registry, ToolkitRegistry):
            return toolkit_registry.call("canonical_order_atoms", self)
        elif isinstance(toolkit_registry, ToolkitWrapper):
            toolkit = toolkit_registry
            return toolkit.canonical_order_atoms(self)
        else:
            raise InvalidToolkitRegistryError(
                "Invalid toolkit_registry passed to from_smiles. Expected ToolkitRegistry or ToolkitWrapper. Got  {}".format(
                    type(toolkit_registry)
                )
            )

    def remap(self, mapping_dict, current_to_new=True):
        """
        Remap all of the indexes in the molecule to match the given mapping dict

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        mapping_dict : dict,
            A dictionary of the mapping between indexes, this should start from 0.
        current_to_new : bool, default=True
            If this is ``True``, then ``mapping_dict`` is of the form ``{current_index: new_index}``;
            otherwise, it is of the form ``{new_index: current_index}``

        Returns
        -------
        new_molecule :  openff.toolkit.topology.molecule.Molecule
            An openff.toolkit.Molecule instance with all attributes transferred, in the PDB order.
        """

        if self.n_virtual_sites != 0:
            raise NotImplementedError("We can not remap virtual sites yet!")

        # make sure the size of the mapping matches the current molecule
        if len(mapping_dict) != self.n_atoms:
            raise ValueError(
                f"The number of mapping indices({len(mapping_dict)}) does not match the number of"
                f"atoms in this molecule({self.n_atoms})"
            )

        # make two mapping dicts we need new to old for atoms
        # and old to new for bonds
        if current_to_new:
            cur_to_new = mapping_dict
            new_to_cur = dict(zip(mapping_dict.values(), mapping_dict.keys()))
        else:
            new_to_cur = mapping_dict
            cur_to_new = dict(zip(mapping_dict.values(), mapping_dict.keys()))

        new_molecule = self.__class__()
        new_molecule.name = self.name

        try:
            # add the atoms list
            for i in range(self.n_atoms):
                # get the old atom info
                old_atom = self._atoms[new_to_cur[i]]
                new_molecule._add_atom(**old_atom.to_dict())
        # this is the first time we access the mapping; catch an index error here corresponding to mapping that starts
        # from 0 or higher
        except (KeyError, IndexError):
            raise IndexError(
                f"The mapping supplied is missing a relation corresponding to atom({i})"
            )

        # add the bonds but with atom indexes in a sorted ascending order
        for bond in self._bonds:
            atoms = sorted([cur_to_new[bond.atom1_index], cur_to_new[bond.atom2_index]])
            bond_dict = bond.to_dict()
            bond_dict["atom1"] = atoms[0]
            bond_dict["atom2"] = atoms[1]
            new_molecule._add_bond(**bond_dict)

        # we can now resort the bonds
        sorted_bonds = sorted(
            new_molecule.bonds, key=operator.attrgetter("atom1_index", "atom2_index")
        )
        new_molecule._bonds = sorted_bonds

        # remap the charges
        if self.partial_charges is not None:
            new_charges = np.zeros(self.n_atoms)
            for i in range(self.n_atoms):
                new_charges[i] = self.partial_charges[new_to_cur[i]].value_in_unit(
                    unit.elementary_charge
                )
            new_molecule.partial_charges = new_charges * unit.elementary_charge

        # remap the conformers there can be more than one
        if self.conformers is not None:
            for conformer in self.conformers:
                new_conformer = np.zeros((self.n_atoms, 3))
                for i in range(self.n_atoms):
                    new_conformer[i] = conformer[new_to_cur[i]].value_in_unit(
                        unit.angstrom
                    )
                new_molecule._add_conformer(new_conformer * unit.angstrom)

        # move any properties across
        new_molecule._properties = self._properties

        return new_molecule

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

    def _construct_angles(self):
        """
        Get an iterator over all i-j-k angles.
        """
        # TODO: Build Angle objects instead of tuple of atoms.
        if not hasattr(self, "_angles"):
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
        if not hasattr(self, "_torsions"):
            self._construct_bonded_atoms_list()

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

    def _construct_bonded_atoms_list(self):
        """
        Construct list of all atoms each atom is bonded to.

        """
        # TODO: Add this to cached_properties
        if not hasattr(self, "_bondedAtoms"):
            # self._atoms = [ atom for atom in self.atoms() ]
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

    def get_bond_between(self, i, j):
        """Returns the bond between two atoms

        Parameters
        ----------
        i, j : int or Atom
            Atoms or atom indices to check

        Returns
        -------
        bond : Bond
            The bond between i and j.

        """
        if isinstance(i, int) and isinstance(j, int):
            atom_i = self._atoms[i]
            atom_j = self._atoms[j]
        elif isinstance(i, Atom) and isinstance(j, Atom):
            atom_i = i
            atom_j = j
        else:
            raise TypeError(
                "Invalid input passed to get_bond_between(). Expected ints or Atoms, "
                "got {} and {}".format(i, j)
            )

        for bond in atom_i.bonds:

            for atom in bond.atoms:

                if atom == atom_i:
                    continue

                if atom == atom_j:
                    return bond

        from openff.toolkit.topology import NotBondedError

        raise NotBondedError("No bond between atom {} and {}".format(i, j))

    @property
    def rings(self):
        """Return the number of rings in this molecule.

        Requires the RDKit to be installed.

        .. note ::

            For systems containing some special cases of connected rings, this
            function may not be well-behaved and may report a different number
            rings than expected. Some problematic cases include networks of many
            (5+) rings or bicyclic moieties (i.e. norbornane).

        """
        if self._rings is None:
            self._get_rings()
        return self._rings

    @RDKitToolkitWrapper.requires_toolkit()
    def _get_rings(self):
        """
        Call out to RDKitToolkitWrapper methods to find the rings in this molecule.

        Requires the RDKit to be installed.

        .. note ::

            For systems containing some special cases of connected rings, this
            function may not be well-behaved and may report a different number
            rings than expected. Some problematic cases include networks of many
            (5+) rings or bicyclic moieties (i.e. norbornane).

        .. todo :: This could be refactored to use ToolkitWrapper.call() to flexibly
            access other toolkits, if find_rings is implemented.

        Returns
        -------
        rings : tuple of tuple of int
            A nested tuple with one subtuple per ring and each subtuple containing
            a tuple of the indices of atoms containing with it. If no rings are
            found, a single empty tuple is returned.

        """
        toolkit = RDKitToolkitWrapper()
        rings = toolkit.find_rings(self)
        self._rings = rings


class Molecule(FrozenMolecule):
    """
    Mutable chemical representation of a molecule, such as a small molecule or biopolymer.

    .. todo :: What other API calls would be useful for supporting biopolymers as small molecules? Perhaps iterating over chains and residues?

    Examples
    --------

    Create a molecule from an sdf file

    >>> from openff.toolkit.utils import get_data_file_path
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
            If specified, attempt to construct a copy of the Molecule from the
            specified object. This can be any one of the following:

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

        >>> from openff.toolkit.utils import get_data_file_path
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

           * If a filename or file-like object is specified but the file
             contains more than one molecule, what is the proper behavior?
             Read just the first molecule, or raise an exception if more
             than one molecule is found?

           * Should we also support SMILES strings or IUPAC names for
             ``other``?

        """
        # super(self, Molecule).__init__(*args, **kwargs)
        super(Molecule, self).__init__(*args, **kwargs)

    # TODO: Change this to add_atom(Atom) to improve encapsulation and extensibility?
    def add_atom(
        self, atomic_number, formal_charge, is_aromatic, stereochemistry=None, name=None
    ):
        """
        Add an atom

        Parameters
        ----------
        atomic_number : int
            Atomic number of the atom
        formal_charge : int
            Formal charge of the atom
        is_aromatic : bool
            If ``True``, atom is aromatic; if ``False``, not aromatic
        stereochemistry : str, optional, default=None
            Either ``'R'`` or ``'S'`` for specified stereochemistry, or ``None`` if stereochemistry is irrelevant
        name : str, optional
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
            name=name,
        )
        return atom_index

    def add_bond_charge_virtual_site(self, atoms, distance, **kwargs):
        """
        Add a virtual site representing the charge on a bond.

        Create a bond charge-type virtual site, in which the location of the
        charge is specified by the position of two atoms. This supports
        placement of a virtual site :math:`S` along a vector between two specified
        atoms, e.g. to allow for a sigma hole for halogens or similar contexts.
        With positive values of the distance, the virtual site lies outside the
        first indexed atom.

        Parameters
        ----------
        atoms : list of :class:`openff.toolkit.topology.molecule.Atom` objects
            The atoms defining the virtual site's position

        distance : :class:`simtk.unit.Quantity` of dimension [Length] wrapping a scalar

        charge_increments : list of floats of shape [N], optional, default=None
            The amount of charge to remove from the VirtualSite's atoms and put
            in the VirtualSite. Indexing in this list should match the ordering
            in the atoms list. Default is None.
        epsilon : float
            Epsilon term for VdW properties of virtual site. Default is ``None``.
        sigma : float, default=None
            Sigma term for VdW properties of virtual site. Default is ``None``.
        rmin_half : float
            Rmin_half term for VdW properties of virtual site. Default is ``None``.
        name : string or None, default=''
            The name of this virtual site. Default is ''.
        symmetric : bool, default=True
            Whether to make virtual site symmetric by creating two particles
            instead of just one. As an example, for N_2 this should be set to
            True to model both lone pairs with the same parameters.

        Returns
        -------
        index : int
            The index of the newly-added virtual site in the molecule

        """
        if kwargs.get("symmetric", True) is True:
            kwargs["orientations"] = [(0, 1), (1, 0)]
        else:
            kwargs["orientations"] = [(0, 1)]
        kwargs.pop("symmetric", None)

        vsite_index = self._add_bond_charge_virtual_site(atoms, distance, **kwargs)
        return vsite_index

    def add_monovalent_lone_pair_virtual_site(
        self, atoms, distance, out_of_plane_angle, in_plane_angle, **kwargs
    ):
        """
        Create a bond charge-type virtual site, in which the location of the charge is specified by the position of three atoms.

        Parameters
        ----------
        atoms : list of three :class:`openff.toolkit.topology.molecule.Atom` objects
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
        name : string or None, default=''
            The name of this virtual site. Default is ''.
        symmetric : bool, default=False
            Whether to make virtual site symmetric by creating two particles
            instead of just one. Note that because this site is defined is placed
            on the noncentral atom, setting this to True will place one particle
            on atom1, and the other on atom3.

        Returns
        -------
        index : int
            The index of the newly-added virtual site in the molecule


        """

        if kwargs.get("symmetric", False) is True:
            kwargs["orientations"] = [(0, 1, 2), (2, 1, 0)]
        else:
            kwargs["orientations"] = [(0, 1, 2)]
        kwargs.pop("symmetric", None)

        vsite_index = self._add_monovalent_lone_pair_virtual_site(
            atoms, distance, out_of_plane_angle, in_plane_angle, **kwargs
        )
        return vsite_index

    # def add_divalent_lone_pair_virtual_site(self, atoms, distance, out_of_plane_angle, in_plane_angle, charge_increments=None, weights=None, epsilon=None, sigma=None, rmin_half=None, name=None):
    def add_divalent_lone_pair_virtual_site(
        self, atoms, distance, out_of_plane_angle, **kwargs
    ):
        """
        Create a divalent lone pair-type virtual site, in which the location of the charge is specified by the position of three atoms.

        Parameters
        ----------
        atoms : list of three :class:`openff.toolkit.topology.molecule.Atom` objects
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
        name : string or None, default=''
            The name of this virtual site. Default is ''.
        symmetric : bool, default=True
            Whether to make virtual site symmetric by creating two particles
            instead of just one. As an example, for TIP5 should be set to True
            to model both lone pairs with the same parameters.


        Returns
        -------
        index : int
            The index of the newly-added virtual site in the molecule

        """
        if kwargs.get("symmetric", True) is True:
            kwargs["orientations"] = [(0, 1, 2), (2, 1, 0)]
        else:
            kwargs["orientations"] = [(0, 1, 2)]
        kwargs.pop("symmetric", None)

        vsite_index = self._add_divalent_lone_pair_virtual_site(
            atoms, distance, out_of_plane_angle, **kwargs
        )
        return vsite_index

    def add_trivalent_lone_pair_virtual_site(self, atoms, distance, **kwargs):
        """
        Create a trivalent lone pair-type virtual site, in which the location of the charge is specified by the position of four atoms.

        Parameters
        ----------
        atoms : list of four :class:`openff.toolkit.topology.molecule.Atom` objects
            The four atoms defining the virtual site's position

        distance : simtk.unit.Quantity of dimension [Length] wrapping a scalar

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

        # This virtual site only makes sense with a single orientation

        kwargs["orientations"] = [(0, 1, 2, 3)]
        kwargs.pop("symmetric", None)

        vsite_index = self._add_trivalent_lone_pair_virtual_site(
            atoms, distance, **kwargs
        )
        return vsite_index

    def add_bond(
        self,
        atom1,
        atom2,
        bond_order,
        is_aromatic,
        stereochemistry=None,
        fractional_bond_order=None,
    ):
        """
        Add a bond between two specified atom indices


        Parameters
        ----------
        atom1 : int or openff.toolkit.topology.molecule.Atom
            Index of first atom
        atom2 : int or openff.toolkit.topology.molecule.Atom
            Index of second atom
        bond_order : int
            Integral bond order of Kekulized form
        is_aromatic : bool
            True if this bond is aromatic, False otherwise
        stereochemistry : str, optional, default=None
            Either ``'E'`` or ``'Z'`` for specified stereochemistry, or ``None`` if stereochemistry is irrelevant
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
            fractional_bond_order=fractional_bond_order,
        )
        return bond_index

    def add_conformer(self, coordinates):
        """
        Add a conformation of the molecule

        Parameters
        ----------
        coordinates: simtk.unit.Quantity(np.array) with shape (n_atoms, 3) and dimension of distance
            Coordinates of the new conformer, with the first dimension of the array corresponding to the atom index in
            the Molecule's indexing system.

        Returns
        -------
        index: int
            The index of this conformer
        """

        # TODO how can be check that a set of coords and no connections
        #   is a conformation that does not change connectivity?

        return self._add_conformer(coordinates)

    def visualize(self, backend="rdkit", width=500, height=300):
        """
        Render a visualization of the molecule in Jupyter

        Parameters
        ----------
        backend : str, optional, default='rdkit'
            Which visualization engine to use. Choose from:

            - rdkit
            - openeye
            - nglview (conformers needed)

        width : int, optional, default=500
            Width of the generated representation (only applicable to
            ``backend=openeye``)
        height : int, optional, default=300
            Width of the generated representation (only applicable to
            ``backend=openeye``)

        Returns
        -------
        object
            Depending on the backend chosen:

            - rdkit, openeye → IPython.display.Image
            - nglview → nglview.NGLWidget

        """
        from openff.toolkit.utils.toolkits import OPENEYE_AVAILABLE, RDKIT_AVAILABLE

        backend = backend.lower()

        if backend == "nglview":
            try:
                import nglview as nv
            except ImportError:
                raise MissingDependencyError("nglview")
            if self.conformers:
                from openff.toolkit.utils.viz import _OFFTrajectoryNGLView

                trajectory_like = _OFFTrajectoryNGLView(self)
                widget = nv.NGLWidget(trajectory_like)
                return widget
            else:
                raise ValueError(
                    "Visualizing with NGLview requires that the molecule has "
                    "conformers."
                )
        if backend == "rdkit":
            if RDKIT_AVAILABLE:
                from rdkit.Chem.Draw import IPythonConsole

                return self.to_rdkit()
            else:
                warnings.warn(
                    "RDKit was requested as a visualization backend but "
                    "it was not found to be installed. Falling back to "
                    "trying to using OpenEye for visualization."
                )
                backend = "openeye"
        if backend == "openeye":
            if OPENEYE_AVAILABLE:
                from IPython.display import Image
                from openeye import oedepict

                oemol = self.to_openeye()

                opts = oedepict.OE2DMolDisplayOptions(
                    width, height, oedepict.OEScale_AutoScale
                )

                oedepict.OEPrepareDepiction(oemol)
                img = oedepict.OEImage(width, height)
                display = oedepict.OE2DMolDisplay(oemol, opts)
                oedepict.OERenderMolecule(img, display)
                png = oedepict.OEWriteImageToString("png", img)
                return Image(png)

        raise ValueError("Could not find an appropriate backend")

    def _ipython_display_(self):
        from IPython.display import display

        try:
            return display(self.visualize(backend="nglview"))
        except (ImportError, ValueError):
            pass

        try:
            return display(self.visualize(backend="rdkit"))
        except ValueError:
            pass

        try:
            return display(self.visualize(backend="openeye"))
        except ValueError:
            pass
