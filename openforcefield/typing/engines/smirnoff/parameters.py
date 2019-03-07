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
from enum import Enum
import logging
import itertools

from collections import OrderedDict

import numpy as np

import lxml.etree as etree

from simtk import openmm, unit
from simtk.openmm.app import element as elem

from openforcefield.utils import detach_units, attach_units, unit_to_string, string_to_unit, extract_serialized_units_from_dict
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
# and returns the corresponding Handler object; 2) a create_force() method that constructs the Force object and adds it
# to the System; and 3) a labelForce() method that provides access to which
# terms are applied to which atoms in specified oemols.
# The static method should be added to the parsers map.
#=============================================================================================


class SMIRNOFFSpecError(Exception):
    """
    Exception for when a non-spec keyword is read and cosmetic attributes are not allowed.
    """

    def __init__(self, msg):
        super().__init__(self, msg)
        self.msg = msg

class IncompatibleUnitError(Exception):
    """
    Exception for when a parameter is in the wrong units for a ParameterHandler's unit system
    """

    def __init__(self, msg):
        super().__init__(self, msg)
        self.msg = msg


class IncompatibleParameterError(Exception):
    """
    Exception for when a set of parameters is scientifically incompatible with another
    """

    def __init__(self, msg):
        super().__init__(self, msg)
        self.msg = msg


class NonbondedMethod(Enum):
    """
    An enumeration of the nonbonded methods
    """
    NoCutoff = 0
    CutoffPeriodic = 1
    CutoffNonPeriodic = 2
    Ewald = 3
    PME = 4


# We can't actually make this derive from dict, because it's possible for the user to change SMIRKS
# of parameters already in the list, which would cause the ParameterType object's SMIRKS and
# the dictionary key's SMIRKS to be out of sync.
class ParameterList(list):
    """
    Parameter list that also supports accessing items by SMARTS string. Remembers the
    most recent parameter that was added.
    """

    # TODO: Make this faster by caching SMARTS -> index lookup?

    # TODO: Override __del__ to make sure we don't remove root atom type

    # TODO: Allow retrieval by `id` as well

    def __init__(self, input_parameter_list=None):
        """
        Initialize a new ParameterList, optionally providing a list of ParameterType objects
        to initially populate it.

        Parameters
        ----------
        input_parameter_list: list[ParameterType], default=None
            A pre-existing list of ParameterType-based objects. If None, this ParameterList
            will be initialized empty.
        """
        super().__init__()

        # We keep track of the last parameter added as this will be
        # used to set output units during serialization
        self._last_added_param = None
        input_parameter_list = input_parameter_list or []
        for input_parameter in input_parameter_list:
            self.append(input_parameter)
            self._last_added_param = input_parameter_list

    @property
    def last_added_parameter(self):
        """
        Get a copy of the last parameter added to this ParameterList. Important
        for serializing as the last parameter added determines the units used for
        serializing all other parameters of this type.

        Returns
        -------
        parameter : a ParameterType-derived object
            The last parameter added to this ParameterList
        """
        return self._last_added_param

    # def __setitem__(self, key, value):
    #     """
    #     Add a new entry to the ParameterList.
    #     Parameters
    #     ----------
    #     key
    #     value
    #     """
    #     super().__setitem__(key, value)
    #     self._last_added_param = value

    def append(self, parameter):
        """
        Add a ParameterType object to the end of the ParameterList

        Parameters
        ----------
        parameter : a ParameterType-derived object

        """
        super().append(parameter)
        self._last_added_param = parameter

    def extend(self, other):
        """
        Add a ParameterType object to the end of the ParameterList

        Parameters
        ----------
        parameter : a ParameterType-derived object

        """
        if not isinstance(other, ParameterList):
            msg = 'ParameterList.extend(other) expected instance of ParameterList, ' \
                  'but received {} (type {}) instead'.format(other, type(other))
            raise TypeError(other)
        super().extend(other)
        if len(other) > 0:
            self._last_added_param = other[-1]


    def insert(self, index, parameter):
        """
        Add a ParameterType object as if this were a list

        Parameters
        ----------
        parameter : a ParameterType-derived object

        """
        super().insert(index, parameter)
        self._last_added_param = parameter

    def __delitem__(self, item):
        """
        Delete item by index or SMIRKS
        """
        if type(item) is str:
            # Try to find by SMIRKS
            for parameter in self:
                if parameter.smirks == item:
                    self.remove(parameter)
                    return

        # Try numerical index access. This will grab the item to remove by index, and
        # then call __delitem__ again on its own SMIRKS, finishing in the "if" statement above
        item_to_remove = self[item].smirks
        del self[item_to_remove]


    def __getitem__(self, item):
        """Retrieve item by index or SMIRKS
        """
        if type(item) == str:
            # Try to retrieve by SMIRKS
            for result in self:
                if result.smirks == item:
                    return result

        # Try traditional access
        result = super().__getitem__(item)
        return result

    # TODO: Override __setitem__ and __del__ to ensure we can slice by SMIRKS as well

    def __contains__(self, item):
        """Check to see if either Parameter or SMIRKS is contained in parameter list.
        """
        if type(item) == str:
            # Special case for SMIRKS strings
            if item in [result.smirks for result in self]:
                return True
        # Fall back to traditional access
        return list.__contains__(self, item)

    def to_list(self, return_cosmetic_attributes=False):
        """
        Render this ParameterList to a normal list, serializing each ParameterType-derived object in it.

        Parameters
        ----------

        return_cosmetic_attributes : bool, optional. default = False
            Whether to return non-spec attributes of each ParameterType-derived object.

        Returns
        -------
        parameter_list : List[dict]
            A serialized representation of a parameter type, with all unit-bearing quantities reduced to a unitless form.
        attached_units : dict['X_unit' : str]
            A dict for converting each serialized quantity back to its unit-bearing form.
        """


        parameter_list = list()

        # If output_units is None, initialize it as an empty dict.
        # It will be populated as ParameterTypes are serialized.
        #if output_units is None:
        #    output_units = dict()

        # Iterate backwards over the list so that we get units from the last-read
        # parameters first, and convert subsequent parameters to those units.
        for parameter in self:
            parameter_dict = parameter.to_dict(return_cosmetic_attributes=return_cosmetic_attributes)
            # parameter_dict_unitless, attached_units = detach_units(parameter_dict, output_units)
            parameter_list.append(parameter_dict)
            # output_units.update(attached_units)

        # Reverse the order of the list (since it was serialized backwards)
        #smirnoff_list.reverse()
        return parameter_list



# TODO: Rename to better reflect role as parameter base class?
class ParameterType(object):
    """
    Base class for SMIRNOFF parameter types.

    """
    # This list will be used for validating input and detecting cosmetic attributes.
    _VALENCE_TYPE = None  # ChemicalEnvironment valence type string expected by SMARTS string for this Handler
    _SMIRNOFF_ATTRIBS = ['smirks'] # Attributes expected per the SMIRNOFF spec.
    _REQUIRE_UNITS = {} # A dict of attribs which will be checked for unit compatibility
    # TODO: Make sure that this doesn't get defined at the class level
    _OPTIONAL_ATTRIBS = ['id', 'parent_id'] # Attributes in the SMIRNOFF spec that may
                                            # be present but have no impact on performance
    _INDEXED_ATTRIBS = []  # list of attribs that will have consecutive numerical suffixes starting at 1
    _ATTRIBS_TO_TYPE = {}  # dict of attributes that need to be cast to a type (like int or float) to be interpreted




    # TODO: Can we provide some shared tools for returning settable/gettable attributes, and checking unit-bearing attributes?

    def __init__(self, smirks=None, permit_cosmetic_attributes=False, **kwargs):
        """
        Create a ParameterType

        Parameters
        ----------
        smirks : str
            The SMIRKS match for the provided parameter type.
        permit_cosmetic_attributes : bool optional. Default = False
            Whether to store non-spec kwargs as "cosmetic attributes", which can be accessed and written out.

        """
        import openforcefield.utils.toolkits

        self._COSMETIC_ATTRIBS = []  # A list that may be populated to record the cosmetic
        # attributes read from a SMIRNOFF data source

        if smirks is None:
            raise ValueError("'smirks' must be specified")

        # TODO: Make better switch for toolkit registry
        if openforcefield.utils.toolkits.OPENEYE_AVAILABLE:
            toolkit = 'openeye'
        elif openforcefield.utils.toolkits.RDKIT_AVAILABLE:
            toolkit = 'rdkit'
        ChemicalEnvironment.validate(
            smirks, ensure_valence_type=self._VALENCE_TYPE, toolkit=toolkit)

        self._smirks = smirks


        # First look for indexed attribs, removing them from kwargs as they're found
        for unidx_key in self._INDEXED_ATTRIBS:
            index = 1
            idx_key = unidx_key+str(index)
            if idx_key in kwargs:
                setattr(self, unidx_key, list())
            while idx_key in kwargs:
                val = kwargs[idx_key]
                if unidx_key in self._REQUIRE_UNITS:
                    if not val.unit.is_compatible(self._REQUIRE_UNITS[unidx_key]):

                        msg = "{} constructor received kwarg {} with value {}, " \
                              "which is incompatible with expected unit {}".format(self.__class__,
                                                                                   idx_key,
                                                                                   val,
                                                                                   self._REQUIRE_UNITS[unidx_key])
                        raise SMIRNOFFSpecError(msg)
                if unidx_key in self._ATTRIBS_TO_TYPE:
                    type_to_cast = self._ATTRIBS_TO_TYPE[unidx_key]
                    val = type_to_cast(val)
                getattr(self, unidx_key).append(val)
                del kwargs[idx_key]
                index += 1
                idx_key = unidx_key+str(index)



        for key, val in kwargs.items():
            if key in self._REQUIRE_UNITS:
                # TODO: Add dynamic property-getter/setter for each thing in self._REQUIRE_UNITS
                if not val.unit.is_compatible(self._REQUIRE_UNITS[key]):
                    msg = "{} constructor received kwarg {} with value {}, " \
                          "which is incompatible with expected unit {}".format(self.__class__,
                                                                               key,
                                                                               val,
                                                                               self._REQUIRE_UNITS[key])
                    raise SMIRNOFFSpecError(msg)

                if key in self._ATTRIBS_TO_TYPE:
                    type_to_cast = self._ATTRIBS_TO_TYPE[key]
                    val = type_to_cast(val)
            # Iterate through
            #attr_name = '_' + key
            if key in self._SMIRNOFF_ATTRIBS:
                setattr(self, key, val)

            # If it's an optional attrib,
            elif key in self._OPTIONAL_ATTRIBS:
                setattr(self, key, val)

            # Handle all unknown kwargs as cosmetic so we can write them back out
            elif permit_cosmetic_attributes:
                self._COSMETIC_ATTRIBS.append(key)
                setattr(self, key, val)
            else:
                raise SMIRNOFFSpecError("Incompatible kwarg {} passed to {} constructor. If this is "
                                        "a desired cosmetic attribute, consider setting "
                                        "'permit_cosmetic_attributes=True'".format({key: val}, self.__class__))



    @property
    def smirks(self):
        return self._smirks

    @smirks.setter
    def smirks(self, smirks):
        # Validate the SMIRKS string to ensure it matches the expected parameter type,
        # raising an exception if it is invalid or doesn't tag a valid set of atoms
        # TODO: Add check to make sure we can't make tree non-hierarchical
        #       This would require parameter type knows which ParameterList it belongs to
        ChemicalEnvironment.validate(
            smirks, ensure_valence_type=self._VALENCE_TYPE)
        self._smirks = smirks

    # TODO: Can we automatically check unit compatibilities for other parameters we create?
    # For example, if we have a parameter with units energy/distance**2, can we check to make
    # sure the dimensionality is preserved when the parameter is modified?

    def to_dict(self, return_cosmetic_attributes=False):
        """
        Convert this ParameterType-derived object to dict. A unit-bearing attribute ('X') will be converted to two dict
        entries, one (['X'] containing the unitless value, and another (['X_unit']) containing a string representation
        of its unit.

        Parameters
        ----------
        return_cosmetic_attributes : bool, optional. default = False
            Whether to return non-spec attributes of this ParameterType


        Returns
        -------
        smirnoff_dict : dict
            The SMIRNOFF-compliant dict representation of this ParameterType-derived object.
        output_units : dict[str: simtk.unit.Unit]
            A mapping from each simtk.unit.Quanitity-valued ParameterType attribute
            to the unit it was converted to during serialization.

        """
        #from simtk import unit


        #output_units = {}

        # Make a list of all attribs that should be included in the
        # returned dict (call list() to make a copy)
        attribs_to_return = list(self._SMIRNOFF_ATTRIBS)
        attribs_to_return += [opt_attrib for opt_attrib in self._OPTIONAL_ATTRIBS if hasattr(self, opt_attrib)]
        if return_cosmetic_attributes:
            attribs_to_return += self._COSMETIC_ATTRIBS

        # Start populating a dict of the attribs
        smirnoff_dict = OrderedDict()
        # If attribs_to_return is ordered here, that will effectively be an informal output ordering
        for attrib_name in attribs_to_return:
            attrib_value = self.__getattribute__(attrib_name)
            #if isinstance(attrib_value, unit.Quantity):
            #    # If the user specified a preferred output unit for this attrib
            #    if attrib_name in output_units:
            #        # convert attrib_val to the desired unit
            #        output_unit = output_units[attrib_name]
            #        # ser_result is a dict of {'unitless_value': val, 'unit': simtk.unit.Unit}
            #        ser_result = serialize_quantity(attrib_value, output_unit=output_unit)
            #
            #    # If the user didn't specify a preferred output unit for this attrib
            #    else:
            #        # Returns a dict of {'unitless_value': val, 'unit': simtk.unit.Unit}
            #        ser_result = serialize_quantity(attrib_value)
            #        # Since a preferred output unit wasn't specified, make this it.
            #        output_units[attrib_name + '_unit'] = ser_result['unit']
            #    smirnoff_dict[attrib_name] = ser_result['unitless_value']
            #
            ## If it's not a Quantity, just add the raw value to the dict
            #else:
            if type(attrib_value) is list:
                for idx, val in enumerate(attrib_value):
                    smirnoff_dict[attrib_name + str(idx+1)] = val
            else:
                smirnoff_dict[attrib_name] = attrib_value

        return smirnoff_dict

# TODO: Should we have a parameter handler registry?


class ParameterHandler(object):
    """Virtual base class for parameter handlers.

    Parameter handlers are configured with some global parameters for a given section, and

    .. warning

       Parameter handler objects can only belong to a single :class:`ForceField` object.
       If you need to create a copy to attach to a different :class:`ForceField` object, use ``create_copy()``.

    """

    _TAGNAME = None  # str of section type handled by this ParameterHandler (XML element name for SMIRNOFF XML representation)
    _INFOTYPE = None  # container class with type information that will be stored in self._types
    _OPENMMTYPE = None  # OpenMM Force class (or None if no equivalent)
    _DEPENDENCIES = None  # list of ParameterHandler classes that must precede this, or None
    _REQUIRED_SPEC_ATTRIBS = []
    _DEFAULT_SPEC_ATTRIBS = {}  # dict of tag-level attributes and their default values
    _OPTIONAL_SPEC_ATTRIBS = []  # list of non-required attributes that can be defined on initialization
    _INDEXED_ATTRIBS = []  # list of parameter attribs that will have consecutive numerical suffixes starting at 1
    _REQUIRE_UNITS = {}  # dict of {header attrib : unit } for input checking
    _ATTRIBS_TO_TYPE = {} # dict of attribs that must be cast to a specific type to be interpreted correctly
    _KWARGS = [] # Kwargs to catch when create_force is called
    _SMIRNOFF_VERSION_INTRODUCED = 0.0  # the earliest version of SMIRNOFF spec that supports this ParameterHandler
    _SMIRNOFF_VERSION_DEPRECATED = None  # if deprecated, the first SMIRNOFF version number it is no longer used


    def __init__(self, permit_cosmetic_attributes=False, **kwargs):
        """
        Initialize a ParameterHandler, optionally with a list of parameters and other kwargs.

        Parameters
        ----------
        permit_cosmetic_attributes : bool
            Whether to accept non-spec kwargs
        **kwargs : dict
            The dict representation of the SMIRNOFF data source

        """

        self._COSMETIC_ATTRIBS = []  # list of cosmetic header attributes to
        # remember and optionally write out

        # Ensure that all required attribs are present
        for reqd_attrib in self._REQUIRED_SPEC_ATTRIBS:
            if not reqd_attrib in kwargs:
                msg = "{} requires {} as a parameter during initialization, however this is not " \
                      "provided. Defined kwargs are {}".format(self.__class__,
                                                               reqd_attrib,
                                                               list(kwargs.keys()))
                raise SMIRNOFFSpecError(msg)

        # list of ParameterType objects # TODO: Change to method accessor so we can access as list or dict
        self._parameters = ParameterList()

        # Handle all the unknown kwargs as cosmetic so we can write them back out
        allowed_header_attribs = self._REQUIRED_SPEC_ATTRIBS + \
                                 list(self._DEFAULT_SPEC_ATTRIBS.keys()) + \
                                 self._OPTIONAL_SPEC_ATTRIBS

        # Check for attribs that need to be casted to specific types
        for attrib, type_to_cast in self._ATTRIBS_TO_TYPE.items():
            if attrib in kwargs:
                kwargs[attrib] = type_to_cast(kwargs[attrib])

        # Check for indexed attribs
        for attrib_basename in self._INDEXED_ATTRIBS:
            attrib_unit_key = attrib_basename + '_unit'

            index = 1
            attrib_w_index = '{}{}'.format(attrib_basename, index)
            while attrib_w_index in kwargs:
                # As long as we keep finding higher-indexed entries for
                # this attrib, add them to the expected arguments
                allowed_header_attribs.append(attrib_w_index)

                # If there's a unit for this attrib, copy unit entries for each index instance
                if attrib_unit_key in kwargs:
                    kwargs[attrib_w_index+'_unit'] = kwargs[attrib_unit_key]

        # Attach units to the handler kwargs, if applicable
        unitless_kwargs, attached_units = extract_serialized_units_from_dict(kwargs)
        smirnoff_data = attach_units(unitless_kwargs, attached_units)

        # Add default values to smirnoff_data if they're not already there
        for default_key, default_val in self._DEFAULT_SPEC_ATTRIBS.items():
            if not (default_key in kwargs):
                smirnoff_data[default_key] = default_val

        # Perform unit compatibility checks
        for key, val in smirnoff_data.items():
            if key in self._REQUIRE_UNITS:
                if not val.unit.is_compatible(self._REQUIRE_UNITS[key]):
                    msg = "{} constructor received kwarg {} with value {}, " \
                          "which is incompatible with expected unit {}".format(self.__class__,
                                                                               key,
                                                                               val,
                                                                               self._REQUIRE_UNITS[key])
                    raise SMIRNOFFSpecError(msg)

        for key, val in smirnoff_data.items():
            # If we're reading the parameter list, iterate through and attach units to
            # each parameter_dict, then use it to initialize a ParameterType
            if key == self._TAGNAME:
                for unitless_param_dict in val:
                    param_dict = attach_units(unitless_param_dict, attached_units)
                    new_parameter = self._INFOTYPE(**param_dict,
                                                   permit_cosmetic_attributes=permit_cosmetic_attributes)
                    self._parameters.append(new_parameter)
            elif key in allowed_header_attribs:
                attr_name = '_' + key
                # TODO: create @property.setter here if attrib requires unit
                setattr(self, attr_name, val)
            elif permit_cosmetic_attributes:
                self._COSMETIC_ATTRIBS.append(key)
                attr_name = '_' + key
                setattr(self, attr_name, val)





    # TODO: Do we need to return these, or can we handle this internally
    @property
    def known_kwargs(self):
        """List of kwargs that can be parsed by the function.
        """
        # TODO: Should we use introspection to inspect the function signature instead?
        return set(self._KWARGS)

    #@classmethod
    def check_parameter_compatibility(self, parameter_kwargs):
        """
        Check to make sure that the fields requiring defined units are compatible with the required units for the
        Parameters handled by this ParameterHandler

        Parameters
        ----------
        parameter_kwargs: dict
            The dict that will be used to construct the ParameterType

        Raises
        ------
        Raises a ValueError if the parameters are incompatible.
        """
        for key in parameter_kwargs:
            if key in self._REQUIRE_UNITS:
                reqd_unit = self._REQUIRE_UNITS[key]
                #if arg in cls._REQUIRE_UNITS:
                #    raise Exception(cls)
                #    reqd_unit = cls._REQUIRE_UNITS[arg]
                val = parameter_kwargs[key]
                if not (reqd_unit.is_compatible(val.unit)):
                    raise IncompatibleUnitError(
                        "Input unit {} is not compatible with ParameterHandler unit {}"
                        .format(val.unit, reqd_unit))

    def check_handler_compatibility(self, handler_kwargs):
        """
        Checks if a set of kwargs used to create a ParameterHandler are compatible with this ParameterHandler. This is
        called if a second handler is attempted to be initialized for the same tag.

        Parameters
        ----------
        handler_kwargs : dict
            The kwargs that would be used to construct

        Raises
        ------
        IncompatibleParameterError if handler_kwargs are incompatible with existing parameters.
        """
        pass

    # TODO: Can we ensure SMIRKS and other parameters remain valid after manipulation?
    def add_parameter(self, parameter_kwargs):
        """Add a parameter to the forcefield, ensuring all parameters are valid.

        Parameters
        ----------
        parameter : dict
            The kwargs to pass to the ParameterHandler.INFOTYPE (a ParameterType) constructor
        """
        #if not(isinstance(parameter, ParameterType)):
        #    raise TypeError("Inappropriate object type passed to ParameterHandler.add_parameter(): {}".format(parameter))
        # TODO: Do we need to check for incompatibility with existing parameters?

        # Perform unit compatibility checks
        self.check_parameter_compatibility(parameter_kwargs)
        # Check for correct SMIRKS valence

        new_parameter = self._INFOTYPE(**parameter_kwargs)
        self._parameters.append(new_parameter)

    def get_parameter(self, parameter_attrs):
        """
        Return the parameters in this ParameterHandler that match the parameter_attrs argument

        Parameters
        ----------
        parameter_attrs : dict of {attr: value}
            The attrs mapped to desired values (for example {"smirks": "[*:1]~[#16:2]=,:[#6:3]~[*:4]", "id": "t105"} )


        Returns
        -------
        list of ParameterType-derived objects
            A list of matching ParameterType-derived objects
        """
        # TODO: This is a necessary API point for Lee-Ping's ForceBalance

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
        logger.info(self.__class__.__name__)  # TODO: Overhaul logging
        matches = ValenceDict()
        for force_type in self._parameters:
            matches_for_this_type = {}
            #atom_top_indexes = [()]
            for atoms in entity.chemical_environment_matches(
                    force_type.smirks):
                atom_top_indexes = tuple(
                    [atom.topology_particle_index for atom in atoms])
                matches_for_this_type[atom_top_indexes] = force_type
            #matches_for_this_type = { atoms : force_type for atoms in entity.chemical_environment_matches(force_type.smirks }
            matches.update(matches_for_this_type)
            logger.info('{:64} : {:8} matches'.format(
                force_type.smirks, len(matches_for_this_type)))

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

    def postprocess_system(self, topology, system, **kwargs):
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

    def to_dict(self, output_units=None, return_cosmetic_attributes=False):
        """
        Convert this ParameterHandler to an OrderedDict, compliant with the SMIRNOFF data spec.

        Parameters
        ----------
        output_units : dict[str : simtk.unit.Unit], optional. Default = None
            A mapping from the ParameterType attribute name to the output unit its value should be converted to.
        return_cosmetic_attributes : bool, optional. Default = False.
            Whether to return non-spec parameter and header attributes in this ParameterHandler.

        Returns
        -------
        smirnoff_data : OrderedDict
            SMIRNOFF-spec compliant representation of this ParameterHandler and its internal ParameterList.
        """
        smirnoff_data = OrderedDict()

        # Set default output units to those from the last parameter added to the ParameterList
        if (output_units is None):
            if (self._parameters.last_added_parameter is not None):
                _, output_units = detach_units(self._parameters.last_added_parameter.to_dict())
            else:
                output_units = dict()

        # Populate parameter list
        parameter_list = self._parameters.to_list(return_cosmetic_attributes=return_cosmetic_attributes)
        unitless_parameter_list = list()

        # Detach units into a separate dict.
        for parameter_dict in parameter_list:
            unitless_parameter_dict, attached_units = detach_units(parameter_dict, output_units=output_units)
            unitless_parameter_list.append(unitless_parameter_dict)
            output_units.update(attached_units)

        # Collapse down indexed attribute units
        # (eg. {'k1_unit': angstrom, 'k2_unit': angstrom} --> {'k_unit': angstrom})
        for attrib_key in self._INDEXED_ATTRIBS:
            index = 1
            # Store a variable that is 'k1_unit'
            idxed_attrib_unit_key = attrib_key + str(index) + '_unit'
            # See if 'k1_unit' is in output_units
            if idxed_attrib_unit_key in output_units:
                # If so, define 'k_unit' and add it to the output_units dict
                attrib_unit_key = attrib_key + '_unit'
                output_units[attrib_unit_key] = output_units[idxed_attrib_unit_key]
            # Increment the 'kN_unit' value, checking that each is the same as the
            # 'k1_unit' value, and deleting them from output_units
            while idxed_attrib_unit_key in output_units:
                # Ensure that no different units are defined for higher indexes of this attrib
                assert output_units[attrib_unit_key] == output_units[idxed_attrib_unit_key]
                del output_units[idxed_attrib_unit_key]
                index += 1
                idxed_attrib_unit_key = attrib_key + str(index) + '_unit'


        smirnoff_data[self._TAGNAME] = unitless_parameter_list


        # Collect the names of handler attributes to return
        header_attribs_to_return = self._REQUIRED_SPEC_ATTRIBS + list(self._DEFAULT_SPEC_ATTRIBS.keys())

        # Check whether the optional attribs are defined, and add them if so
        for key in self._OPTIONAL_SPEC_ATTRIBS:
            if hasattr(self, key):
                header_attribs_to_return.append(key)
        # Add the cosmetic attributes if requested
        if return_cosmetic_attributes:
            header_attribs_to_return += self._COSMETIC_ATTRIBS


        # Go through the attribs of this ParameterHandler and collect the appropriate values to return
        header_attribute_dict = {}
        for header_attribute in header_attribs_to_return:
            value = getattr(self, '_' + header_attribute)
            header_attribute_dict[header_attribute] = value

        # Detach all units from the header attribs
        unitless_header_attribute_dict, attached_header_units = detach_units(header_attribute_dict)

        # Convert all header attrib units (eg. {'length_unit': simtk.unit.angstrom}) to strings (eg.
        # {'length_unit': 'angstrom'}) and add them to the header attribute dict
        # TODO: Should I check for collisions between parameter "_unit" keys and header "_unit" keys?
        output_units.update(attached_header_units)
        for key, value in output_units.items():
            value_str = unit_to_string(value)
            # Made a note to add this to the smirnoff spec
            output_units[key] = value_str

        # Reattach the attached units here
        smirnoff_data.update(unitless_header_attribute_dict)
        smirnoff_data.update(output_units)
        return smirnoff_data
#=============================================================================================


class ConstraintHandler(ParameterHandler):
    """Handle SMIRNOFF ``<Constraints>`` tags

    ``ConstraintHandler`` must be applied before ``BondHandler`` and ``AngleHandler``,
    since those classes add constraints for which equilibrium geometries are needed from those tags.
    """

    class ConstraintType(ParameterType):
        """A SMIRNOFF constraint type"""
        _VALENCE_TYPE = 'Bond'
        _SMIRNOFF_ATTRIBS = ['smirks']  # Attributes expected per the SMIRNOFF spec.
        _OPTIONAL_ATTRIBS = ['distance', 'id', 'parent_id']
        _REQUIRE_UNITS = {'distance': unit.angstrom}
        def __init__(self, **kwargs):
            super().__init__(**kwargs)
            # if 'distance' in node.attrib:
            #     self.distance = _extract_quantity_from_xml_element(
            #         node, parent, 'distance'
            #     )  # Constraint with specified distance will be added by ConstraintHandler
            # else:
            #     self.distance = True  # Constraint to equilibrium bond length will be added by HarmonicBondHandler

    _TAGNAME = 'Constraint'
    _INFOTYPE = ConstraintType
    _OPENMMTYPE = None  # don't create a corresponding OpenMM Force class


    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def create_force(self, system, topology, **kwargs):
        constraints = self.get_matches(topology)
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
        _VALENCE_TYPE = 'Bond' # ChemicalEnvironment valence type string expected by SMARTS string for this Handler
        _SMIRNOFF_ATTRIBS = ['smirks', 'length', 'k']  # Attributes expected per the SMIRNOFF spec.
        _REQUIRE_UNITS = {'length' : unit.angstrom,
                          'k' : unit.kilocalorie_per_mole / unit.angstrom**2}
        _INDEXED_ATTRIBS = ['k']  # May be indexed (by integer bond order) if fractional bond orders are used

        def __init__(self, **kwargs):
            super().__init__(**kwargs)  # Base class handles ``smirks`` and ``id`` fields


    _TAGNAME = 'Bonds'  # SMIRNOFF tag name to process
    _INFOTYPE = BondType  # class to hold force type info
    _OPENMMTYPE = openmm.HarmonicBondForce  # OpenMM force class to create
    _DEPENDENCIES = [ConstraintHandler]  # ConstraintHandler must be executed first
    _DEFAULT_SPEC_ATTRIBS = {'potential': 'harmonic',
                             'fractional_bondorder_method': None,
                             'fractional_bondorder_interpolation': 'linear'}
    _INDEXED_ATTRIBS = ['k'] # May be indexed (by integer bond order) if fractional bond orders are used

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def create_force(self, system, topology, **kwargs):
        # Create or retrieve existing OpenMM Force object
        # TODO: The commented line below should replace the system.getForce search
        #force = super(BondHandler, self).create_force(system, topology, **kwargs)
        existing = [system.getForce(i) for i in range(system.getNumForces())]
        existing = [f for f in existing if type(f) == self._OPENMMTYPE]
        if len(existing) == 0:
            force = self._OPENMMTYPE()
            system.addForce(force)
        else:
            force = existing[0]

        # Add all bonds to the system.
        bonds = self.get_matches(topology)
        skipped_constrained_bonds = 0  # keep track of how many bonds were constrained (and hence skipped)
        for (atoms, bond_params) in bonds.items():
            # Get corresponding particle indices in Topology
            #particle_indices = tuple([ atom.particle_index for atom in atoms ])

            # Ensure atoms are actually bonded correct pattern in Topology
            topology.assert_bonded(atoms[0], atoms[1])

            # Compute equilibrium bond length and spring constant.
            topology_bond = topology.get_bond_between(atoms[0], atoms[1])
            if topology_bond.bond.fractional_bond_order is None:
                [k, length] = [bond_params.k, bond_params.length]
            else:
                # Interpolate using fractional bond orders
                # TODO: Do we really want to allow per-bond specification of interpolation schemes?
                order = topology_bond.bond.fractional_bond_order
                if self.fractional_bondorder_interpolation == 'interpolate-linear':
                    k = bond_params.k[0] + (bond_params.k[1] - bond_params.k[0]) * (order - 1.)
                    length = bond_params.length[0] + (
                        bond_params.length[1] - bond_params.length[0]) * (order - 1.)
                else:
                    raise Exception(
                        "Partial bondorder treatment {} is not implemented.".
                        format(self.fractional_bondorder_method))

            # Handle constraints.
            # TODO: I don't understand why there are two if statements checking the same thing here.
            if topology.is_constrained(*atoms):
                # Atom pair is constrained; we don't need to add a bond term.
                skipped_constrained_bonds += 1
                # Check if we need to add the constraint here to the equilibrium bond length.
                if topology.is_constrained(*atoms) is True: # Note: This could have a value of the constraint length
                    # Mark that we have now assigned a specific constraint distance to this constraint.
                    topology.add_constraint(*atoms, length)
                    # Add the constraint to the System.
                system.addConstraint(*atoms, length)
                #system.addConstraint(*particle_indices, length)
                continue

            # Add harmonic bond to HarmonicBondForce
            force.addBond(*atoms, length, k)
            #force.addBond(*particle_indices, length, k)

        logger.info('{} bonds added ({} skipped due to constraints)'.format(
            len(bonds) - skipped_constrained_bonds, skipped_constrained_bonds))

        # Check that no topological bonds are missing force parameters
        #_check_for_missing_valence_terms('BondForce', topology, bonds.keys(), topology.bonds)


#=============================================================================================


class AngleHandler(ParameterHandler):
    """Handle SMIRNOFF ``<AngleForce>`` tags"""

    class AngleType(ParameterType):
        """A SMIRNOFF angle type."""
        _VALENCE_TYPE = 'Angle'  # ChemicalEnvironment valence type string expected by SMARTS string for this Handler
        _SMIRNOFF_ATTRIBS = ['smirks', 'angle', 'k']  # Attributes expected per the SMIRNOFF spec.
        _REQUIRE_UNITS = {'angle': unit.degree,
                          'k': unit.kilocalorie_per_mole / unit.degree**2}


        def __init__(self, **kwargs):
            super().__init__(**kwargs)  # base class handles ``smirks`` and ``id`` fields


    _TAGNAME = 'Angles'  # SMIRNOFF tag name to process
    _INFOTYPE = AngleType  # class to hold force type info
    _OPENMMTYPE = openmm.HarmonicAngleForce  # OpenMM force class to create
    _DEFAULT_SPEC_ATTRIBS = {'potential': 'harmonic'}

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def create_force(self, system, topology, **kwargs):
        #force = super(AngleHandler, self).create_force(system, topology, **kwargs)
        existing = [system.getForce(i) for i in range(system.getNumForces())]
        existing = [f for f in existing if type(f) == self._OPENMMTYPE]
        if len(existing) == 0:
            force = self._OPENMMTYPE()
            system.addForce(force)
        else:
            force = existing[0]

        # Add all angles to the system.
        angles = self.get_matches(topology)
        skipped_constrained_angles = 0  # keep track of how many angles were constrained (and hence skipped)
        for (atoms, angle) in angles.items():
            # Get corresponding particle indices in Topology
            #particle_indices = tuple([ atom.particle_index for atom in atoms ])

            # Ensure atoms are actually bonded correct pattern in Topology
            for (i, j) in [(0, 1), (1, 2)]:
                topology.assert_bonded(atoms[i], atoms[j])

            if topology.is_constrained(
                    atoms[0], atoms[1]) and topology.is_constrained(
                        atoms[1], atoms[2]) and topology.is_constrained(
                            atoms[0], atoms[2]):
                # Angle is constrained; we don't need to add an angle term.
                skipped_constrained_angles += 1
                continue

            force.addAngle(*atoms, angle.angle, angle.k)

        logger.info('{} angles added ({} skipped due to constraints)'.format(
            len(angles) - skipped_constrained_angles,
            skipped_constrained_angles))

        # Check that no topological angles are missing force parameters
        #_check_for_missing_valence_terms('AngleForce', topology, angles.keys(), topology.angles())


#=============================================================================================


class ProperTorsionHandler(ParameterHandler):
    """Handle SMIRNOFF ``<ProperTorsionForce>`` tags"""

    class ProperTorsionType(ParameterType):
        """A SMIRNOFF torsion type for proper torsions."""

        _VALENCE_TYPE = 'ProperTorsion'
        _SMIRNOFF_ATTRIBS = ['smirks', 'periodicity', 'phase', 'k']  # Attributes expected per the SMIRNOFF spec.
        _REQUIRE_UNITS = {'k': unit.kilocalorie_per_mole,
                          'phase': unit.degree}
        _OPTIONAL_ATTRIBS = ['id', 'parent_id', 'idivf']
        _INDEXED_ATTRIBS = ['k', 'phase', 'periodicity', 'idivf']
        _ATTRIBS_TO_TYPE = {'periodicity': int,
                            'idivf': float}

        def __init__(self, **kwargs):
            super().__init__(**kwargs)  # base class handles ``smirks`` and ``id`` fields


    _TAGNAME = 'ProperTorsions'  # SMIRNOFF tag name to process
    _INFOTYPE = ProperTorsionType  # info type to store
    _OPENMMTYPE = openmm.PeriodicTorsionForce  # OpenMM force class to create
    _DEFAULT_SPEC_ATTRIBS = {'potential': 'charmm',
                             'default_idivf': 'auto'}
    _INDEXED_ATTRIBS = ['k', 'phase', 'periodicity', 'idivf']

    def __init__(self, potential=None, **kwargs):

        # NOTE: We do not want to overwrite idivf values here! If they're missing from the ParameterType
        # dictionary, that means they should be set to defualt _AT SYSTEM CREATION TIME_. The user may
        # change that default to a different value than it is now. The solution here will be to leave
        # those idivfX values uninitialized and deal with it during system creation

        super().__init__(**kwargs)


    def create_force(self, system, topology, **kwargs):
        #force = super(ProperTorsionHandler, self).create_force(system, topology, **kwargs)
        existing = [system.getForce(i) for i in range(system.getNumForces())]
        existing = [f for f in existing if type(f) == self._OPENMMTYPE]
        if len(existing) == 0:
            force = self._OPENMMTYPE()
            system.addForce(force)
        else:
            force = existing[0]
        # Add all proper torsions to the system.
        torsions = self.get_matches(topology)
        for (atom_indices, torsion) in torsions.items():
            # Ensure atoms are actually bonded correct pattern in Topology
            for (i, j) in [(0, 1), (1, 2), (2, 3)]:
                topology.assert_bonded(atom_indices[i], atom_indices[j])
            for (periodicity, phase, k, idivf) in zip(torsion.periodicity,
                                               torsion.phase, torsion.k, torsion.idivf):
                if idivf == 'auto':
                    # TODO: Implement correct "auto" behavior
                    raise NotImplementedError("The OpenForceField toolkit hasn't implemented "
                                              "support for the torsion `idivf` value of 'auto'")

                force.addTorsion(atom_indices[0], atom_indices[1],
                                 atom_indices[2], atom_indices[3], periodicity,
                                 phase, k/idivf)

        logger.info('{} torsions added'.format(len(torsions)))

        # Check that no topological torsions are missing force parameters
        #_check_for_missing_valence_terms('ProperTorsionForce', topology, torsions.keys(), topology.torsions())


class ImproperTorsionHandler(ParameterHandler):
    """Handle SMIRNOFF ``<ImproperTorsionForce>`` tags"""

    class ImproperTorsionType(ParameterType):
        """A SMIRNOFF torsion type for improper torsions."""
        _VALENCE_TYPE = 'ImproperTorsion'
        _SMIRNOFF_ATTRIBS = ['smirks', 'periodicity', 'phase', 'k']  # Attributes expected per the SMIRNOFF spec.
        _REQUIRE_UNITS = {'k': unit.kilocalorie_per_mole,
                          'phase': unit.degree}
        _OPTIONAL_ATTRIBS = ['id', 'parent_id', 'idivf']
        _INDEXED_ATTRIBS = ['k', 'phase', 'periodicity', 'idivf']
        _ATTRIBS_TO_TYPE = {'periodicity': int,
                            'idivf': float}

        def __init__(self, **kwargs):
            super().__init__( **kwargs)


    _TAGNAME = 'ImproperTorsions'  # SMIRNOFF tag name to process
    _INFOTYPE = ImproperTorsionType  # info type to store
    _OPENMMTYPE = openmm.PeriodicTorsionForce  # OpenMM force class to create
    _HANDLER_DEFAULTS = {'potential': 'charmm',
                         'default_idivf': 'auto'}
    _INDEXED_ATTRIBS = ['k', 'phase', 'periodicity', 'idivf']

    def __init__(self, potential=None, **kwargs):
        # Cast necessary kwargs to int
        super().__init__(**kwargs)



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
        logger.info(self.__class__.__name__)  # TODO: Overhaul logging
        matches = ImproperDict()
        for force_type in self._parameters:
            matches_for_this_type = {}
            #atom_top_indexes = [()]
            for atoms in entity.chemical_environment_matches(
                    force_type.smirks):
                atom_top_indexes = tuple(
                    [atom.topology_particle_index for atom in atoms])
                matches_for_this_type[atom_top_indexes] = force_type
            #matches_for_this_type = { atoms : force_type for atoms in entity.chemical_environment_matches(force_type.smirks }
            matches.update(matches_for_this_type)
            logger.info('{:64} : {:8} matches'.format(
                force_type.smirks, len(matches_for_this_type)))

        logger.info('{} matches identified'.format(len(matches)))
        return matches


    def create_force(self, system, topology, **kwargs):
        #force = super(ImproperTorsionHandler, self).create_force(system, topology, **kwargs)
        #force = super().create_force(system, topology, **kwargs)
        existing = [system.getForce(i) for i in range(system.getNumForces())]
        existing = [
            f for f in existing if type(f) == openmm.PeriodicTorsionForce
        ]
        if len(existing) == 0:
            force = openmm.PeriodicTorsionForce()
            system.addForce(force)
        else:
            force = existing[0]

        # Add all improper torsions to the system
        impropers = self.get_matches(topology)
        for (atom_indices, improper) in impropers.items():
            # Ensure atoms are actually bonded correct pattern in Topology
            # For impropers, central atom is atom 1
            for (i, j) in [(0, 1), (1, 2), (1, 3)]:
                topology.assert_bonded(atom_indices[i], atom_indices[j])
                #topology.assert_bonded(topology.atom(atom_indices[i]), topology.atom(atom_indices[j]))

            # Impropers are applied in three paths around the trefoil having the same handedness
            for (improper_periodicity, improper_phase, improper_k, improper_idivf) in zip(improper.periodicity,
                                               improper.phase, improper.k, improper.idivf):
                # TODO: Implement correct "auto" behavior
                if improper_idivf == 'auto':
                    improper_idivf = 3
                    logger.warning("The OpenForceField toolkit hasn't implemented "
                                   "support for the torsion `idivf` value of 'auto'."
                                   "Currently assuming a value of '3' for impropers.")
                # Permute non-central atoms
                others = [atom_indices[0], atom_indices[2], atom_indices[3]]
                # ((0, 1, 2), (1, 2, 0), and (2, 0, 1)) are the three paths around the trefoil
                for p in [(others[i], others[j], others[k]) for (i, j, k) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]]:
                    # The torsion force gets added three times, since the k is divided by three
                    force.addTorsion(atom_indices[1], p[0], p[1], p[2],
                                     improper_periodicity, improper_phase, improper_k/improper_idivf)
        logger.info(
            '{} impropers added, each applied in a six-fold trefoil'.format(
                len(impropers)))

        # Check that no topological torsions are missing force parameters
        #_check_for_missing_valence_terms('ImproperTorsionForce', topology, torsions.keys(), topology.impropers())


class vdWHandler(ParameterHandler):
    """Handle SMIRNOFF ``<vdW>`` tags"""

    class vdWType(ParameterType):
        """A SMIRNOFF vdWForce type."""
        _VALENCE_TYPE = 'Atom'  # ChemicalEnvironment valence type expected for SMARTS
        _SMIRNOFF_ATTRIBS = ['smirks', 'epsilon'] # Attributes expected per the SMIRNOFF spec.
        _OPTIONAL_ATTRIBS = ['id', 'parent_id', 'sigma', 'rmin_half']
        _REQUIRE_UNITS = {
            'epsilon': unit.kilocalorie_per_mole,
            'sigma': unit.angstrom,
            'rmin_half': unit.angstrom
        }

        def __init__(self, **kwargs):



            sigma = kwargs.get('sigma', None)
            rmin_half = kwargs.get('rmin_half', None)
            if (sigma is None) and (rmin_half is None):
                raise SMIRNOFFSpecError("Either sigma or rmin_half must be specified.")
            if (sigma is not None) and (rmin_half is not None):
                raise SMIRNOFFSpecError(
                    "BOTH sigma and rmin_half cannot be specified simultaneously."
                )

            # TODO: Is it necessary to force everything to be sigma? We could handle having either when create_force runs
            if (rmin_half is not None):
                kwargs['sigma'] = 2. * rmin_half / (2.**(1. / 6.))
                del kwargs['rmin_half']

            super().__init__(**kwargs)


        @property
        def attrib(self):
            """Return all storable attributes as a dict.
            """
            names = ['smirks', 'sigma', 'epsilon']
            return {
                name: getattr(self, name)
                for name in names if hasattr(self, name)
            }

    _TAGNAME = 'vdW'  # SMIRNOFF tag name to process
    _INFOTYPE = vdWType  # info type to store
    _OPENMMTYPE = openmm.NonbondedForce  # OpenMM force class to create
    _KWARGS = ['ewaldErrorTolerance', 'useDispersionCorrection'] # Kwargs to catch when create_force is called
    _REQUIRE_UNITS = {'switch': unit.angstrom,
                      'cutoff': unit.angstrom}
    _DEFAULT_SPEC_ATTRIBS = {
        'potential': 'Lennard-Jones-12-6',
        'combining_rules': 'Loentz-Berthelot',
        'scale12': 0.0,
        'scale13': 0.0,
        'scale14': 0.5,
        'scale15': 1.0,
        'switch': 8.0 * unit.angstroms,
        'cutoff': 9.0 * unit.angstroms,
        'long_range_dispersion': 'isotropic',
        'nonbonded_method': NonbondedMethod.NoCutoff
    }
    _ATTRIBS_TO_TYPE = {'scale12': float,
                        'scale13': float,
                        'scale14': float,
                        'scale15': float}

    # TODO: Is this necessary
    _SCALETOL = 1e-5

    _NONBOND_METHOD_MAP = {
        NonbondedMethod.NoCutoff:
        openmm.NonbondedForce.NoCutoff,
        NonbondedMethod.CutoffPeriodic:
        openmm.NonbondedForce.CutoffPeriodic,
        NonbondedMethod.CutoffNonPeriodic:
        openmm.NonbondedForce.CutoffNonPeriodic,
        NonbondedMethod.Ewald:
        openmm.NonbondedForce.Ewald,
        NonbondedMethod.PME:
        openmm.NonbondedForce.PME
    }

    def __init__(self,
                 # forcefield,
                 # scale12=None,
                 # scale13=None,
                 # scale14=None,
                 # scale15=None,
                 # potential=None,
                 # switch=None,
                 # cutoff=None,
                 # long_range_dispersion=None,
                 # combining_rules=None,
                 # nonbonded_method=None,
                 **kwargs):
        super().__init__(**kwargs)

        # TODO: Find a better way to set defaults
        # TODO: Validate these values against the supported output types (openMM force kwargs?)
        # TODO: Add conditional logic to assign NonbondedMethod and check compatibility

        # # Set the nonbonded method
        # if nonbonded_method is None:
        #     self._nonbonded_method = NonbondedMethod.NoCutoff
        # else:
        #     # If it's a string that's the name of a nonbonded method
        #     if type(nonbonded_method) is str:
        #         self._nonbonded_method = NonbondedMethod[nonbonded_method]
        #
        #     # If it's an enum'ed value of NonbondedMethod
        #     elif nonbonded_method in NonbondedMethod:
        #         self._nonbonded_method = nonbonded_method
        #     # If it's an openMM nonbonded method, reverse it back to a package-independent enum
        #     elif nonbonded_method in self._NONBOND_METHOD_MAP.values():
        #         for key, val in self._NONBOND_METHOD_MAP.items():
        #             if nonbonded_method == val:
        #                 self._nonbonded_method = key
        #                 break



    def check_handler_compatibility(self,
                                    handler_kwargs,
                                    assume_missing_is_default=True):
        """
               Checks if a set of kwargs used to create a ParameterHandler are compatible with this ParameterHandler. This is
               called if a second handler is attempted to be initialized for the same tag. If no value is given for a field, it
               will be assumed to expect the ParameterHandler class default.

               Parameters
               ----------
               handler_kwargs : dict
                   The kwargs that would be used to construct a ParameterHandler
               assume_missing_is_default : bool
                   If True, will assume that parameters not specified in handler_kwargs would have been set to the default.
                   Therefore, an exception is raised if the ParameterHandler is incompatible with the default value for a
                   unspecified field.

               Raises
               ------
               IncompatibleParameterError if handler_kwargs are incompatible with existing parameters.
               """
        compare_attr_to_kwargs = {
            self._scale12: 'scale12',
            self._scale13: 'scale13',
            self._scale14: 'scale14',
            self._scale15: 'scale15'
        }
        for attr, kwarg_key in compare_attr_to_kwargs.items():
            kwarg_val = handler_kwargs.get(kwarg_key,
                                           self._DEFAULTS[kwarg_key])
            if abs(kwarg_val - attr) > self._SCALETOL:
                raise IncompatibleParameterError(
                    "Difference between '{}' values is beyond allowed tolerance {}. "
                    "(handler value: {}, incompatible valie: {}".format(
                        kwarg_key, self._SCALETOL, attr, kwarg_val))

        # TODO: Test for other possible incompatibilities here -- Probably just check for string equality for now,
        # detailed check will require some openMM/MD expertise)
        #self._potential: 'potential',
        #self._combining_rules: 'combining_rules',
        #self._switch: 'switch',
        #self._cutoff: 'cutoff',
        #self._long_range_dispersion:'long_range_dispersion'
        #}

    # TODO: nonbondedMethod and nonbondedCutoff should now be specified by StericsForce attributes
    def create_force(self, system, topology, **kwargs):

        force = openmm.NonbondedForce()
        nonbonded_method = self._NONBOND_METHOD_MAP[self._nonbonded_method]
        force.setNonbondedMethod(nonbonded_method)
        force.setCutoffDistance(self._cutoff.in_units_of(unit.nanometer))
        if 'ewaldErrorTolerance' in kwargs:
            force.setEwaldErrorTolerance(kwargs['ewaldErrorTolerance'])
        if 'useDispersionCorrection' in kwargs:
            force.setUseDispersionCorrection(
                bool(kwargs['useDispersionCorrection']))
        system.addForce(force)

        # Iterate over all defined Lennard-Jones types, allowing later matches to override earlier ones.
        atoms = self.get_matches(topology)

        # Create all particles.
        for particle in topology.topology_particles:
            force.addParticle(0.0, 1.0, 0.0)

        # Set the particle Lennard-Jones terms.
        for (atoms, ljtype) in atoms.items():
            force.setParticleParameters(atoms[0], 0.0, ljtype.sigma,
                                        ljtype.epsilon)

        # Check that no atoms are missing force parameters
        # QUESTION: Don't we want to allow atoms without force parameters? Or perhaps just *particles* without force parameters, but not atoms?
        # TODO: Enable this check
        #_check_for_missing_valence_terms('NonbondedForce Lennard-Jones parameters', topology, atoms.keys(), topology.atoms)


    # TODO: Can we express separate constraints for postprocessing and normal processing?
    def postprocess_system(self, system, topology, **kwargs):
        # Create exceptions based on bonds.
        # TODO: This postprocessing must occur after the ChargeIncrementModelHandler
        # QUESTION: Will we want to do this for *all* cases, or would we ever want flexibility here?
        bond_particle_indices = []
        for bond in topology.topology_bonds:
            topology_atoms = [atom for atom in bond.atoms]
            bond_particle_indices.append(
                (topology_atoms[0].topology_particle_index,
                 topology_atoms[1].topology_particle_index))
        for force in system.getForces():
            # TODO: Should we just store which `Force` object we are adding to and use that instead,
            # to prevent interference with other kinds of forces in the future?
            # TODO: Can we generalize this to allow for `CustomNonbondedForce` implementations too?
            if isinstance(force, openmm.NonbondedForce):
                #nonbonded.createExceptionsFromBonds(bond_particle_indices, self.coulomb14scale, self.lj14scale)

                # TODO: Don't mess with electrostatic scaling here. Have a separate electrostatics handler.
                force.createExceptionsFromBonds(bond_particle_indices, 0.83333,
                                                self._scale14)
                #force.createExceptionsFromBonds(bond_particle_indices, self.coulomb14scale, self._scale14)

class ToolkitAM1BCCHandler(ParameterHandler):
    """Handle SMIRNOFF ``<ToolkitAM1BCC>`` tags"""

    _TAGNAME = 'ToolkitAM1BCC'  # SMIRNOFF tag name to process
    _OPENMMTYPE = openmm.NonbondedForce  # OpenMM force class to create or utilize
    _DEPENDENCIES = [vdWHandler] # vdWHandler must first run NonBondedForce.addParticle for each particle in the topology
    _KWARGS = ['charge_from_molecules', 'toolkit_registry'] # Kwargs to catch when create_force is called



    def __init__(self, **kwargs):
        super().__init__(**kwargs)



    def check_handler_compatibility(self, handler_kwargs, assume_missing_is_default=True):
        """
        Checks if a set of kwargs used to create a ParameterHandler are compatible with this ParameterHandler. This is
        called if a second handler is attempted to be initialized for the same tag. If no value is given for a field, it
        will be assumed to expect the ParameterHandler class default.

        Parameters
        ----------
        handler_kwargs : dict
            The kwargs that would be used to construct a ParameterHandler
        assume_missing_is_default : bool
            If True, will assume that parameters not specified in handler_kwargs would have been set to the default.
            Therefore, an exception is raised if the ParameterHandler is incompatible with the default value for a
            unspecified field.

        Raises
        ------
        IncompatibleParameterError if handler_kwargs are incompatible with existing parameters.
        """
        pass

    def assign_charge_from_molecules(self, molecule, charge_mols):
        """
        Given an input molecule, checks against a list of molecules for an isomorphic match. If found, assigns
        partial charges from the match to the input molecule.

        Parameters
        ----------
        molecule : an openforcefield.topology.FrozenMolecule
            The molecule to have partial charges assigned if a match is found.
        charge_mols : list of [openforcefield.topology.FrozenMolecule]
            A list of molecules with charges already assigned.

        Returns
        -------
        match_found : bool
            Whether a match was found. If True, the input molecule will have been modified in-place.
        """

        from networkx.algorithms.isomorphism import GraphMatcher
        import simtk.unit
        # Define the node/edge attributes that we will use to match the atoms/bonds during molecule comparison
        node_match_func = lambda x, y: ((x['atomic_number'] == y['atomic_number']) and
                                        (x['stereochemistry'] == y['stereochemistry']) and
                                        (x['is_aromatic'] == y['is_aromatic'])
                                        )
        edge_match_func = lambda x, y: ((x['bond_order'] == y['bond_order']) and
                                        (x['stereochemistry'] == y['stereochemistry']) and
                                        (x['is_aromatic'] == y['is_aromatic'])
                                        )
        # Check each charge_mol for whether it's isomorphic to the input molecule
        for charge_mol in charge_mols:
            if molecule.is_isomorphic(charge_mol):
                # Take the first valid atom indexing map
                ref_mol_G = molecule.to_networkx()
                charge_mol_G = charge_mol.to_networkx()
                GM = GraphMatcher(
                    charge_mol_G,
                    ref_mol_G,
                    node_match=node_match_func,
                    edge_match=edge_match_func)
                for mapping in GM.isomorphisms_iter():
                    topology_atom_map = mapping
                    break
                # Set the partial charges

                # Get the partial charges
                # Make a copy of the charge molecule's charges array (this way it's the right shape)
                temp_mol_charges = simtk.unit.Quantity(charge_mol.partial_charges)
                for charge_idx, ref_idx in topology_atom_map.items():
                    temp_mol_charges[ref_idx] = charge_mol.partial_charges[charge_idx]
                molecule.partial_charges = temp_mol_charges
                return True

        # If no match was found, return False
        return False

    def create_force(self, system, topology, **kwargs):

        from openforcefield.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY
        from openforcefield.topology import FrozenMolecule, TopologyAtom, TopologyVirtualSite

        existing = [system.getForce(i) for i in range(system.getNumForces())]
        existing = [f for f in existing if type(f) == self._OPENMMTYPE]
        if len(existing) == 0:
            force = self._OPENMMTYPE()
            system.addForce(force)
        else:
            force = existing[0]

        for ref_mol in topology.reference_molecules:

            # Make a temporary copy of ref_mol to assign charges from charge_mol
            temp_mol = FrozenMolecule(ref_mol)

            # First, check whether any of the reference molecules in the topology are in the charge_from_mol list
            charges_from_charge_mol = False
            if 'charge_from_molecules' in kwargs:
                charges_from_charge_mol = self.assign_charge_from_molecules(temp_mol, kwargs['charge_from_molecules'])

            # If the molecule wasn't assigned parameters from a manually-input charge_mol, calculate them here
            if not(charges_from_charge_mol):
                toolkit_registry = kwargs.get('toolkit_registry', GLOBAL_TOOLKIT_REGISTRY)
                temp_mol.generate_conformers(num_conformers=10, toolkit_registry=toolkit_registry)
                #temp_mol.compute_partial_charges(quantum_chemical_method=self._quantum_chemical_method,
                #                                 partial_charge_method=self._partial_charge_method)
                temp_mol.compute_partial_charges_am1bcc(toolkit_registry=toolkit_registry)
            # Assign charges to relevant atoms
            for topology_molecule in topology._reference_molecule_to_topology_molecules[ref_mol]:
                for topology_particle in topology_molecule.particles:
                    topology_particle_index = topology_particle.topology_particle_index
                    if type(topology_particle) is TopologyAtom:
                        ref_mol_particle_index = topology_particle.atom.molecule_particle_index
                    if type(topology_particle) is TopologyVirtualSite:
                        ref_mol_particle_index = topology_particle.virtual_site.molecule_particle_index
                    particle_charge = temp_mol._partial_charges[ref_mol_particle_index]

                    # Retrieve nonbonded parameters for reference atom (charge not set yet)
                    _, sigma, epsilon = force.getParticleParameters(topology_particle_index)
                    # Set the nonbonded force with the partial charge
                    force.setParticleParameters(topology_particle_index,
                                                particle_charge, sigma,
                                                epsilon)



    # TODO: Move chargeModel and library residue charges to SMIRNOFF spec
    def postprocess_system(self, system, topology, **kwargs):
        bonds = self.get_matches(topology)

        # Apply bond charge increments to all appropriate force groups
        # QUESTION: Should we instead apply this to the Topology in a preprocessing step, prior to spreading out charge onto virtual sites?
        for force in system.getForces():
            if force.__class__.__name__ in [
                    'NonbondedForce'
            ]:  # TODO: We need to apply this to all Force types that involve charges, such as (Custom)GBSA forces and CustomNonbondedForce
                for (atoms, bond) in bonds.items():
                    # Get corresponding particle indices in Topology
                    particle_indices = tuple(
                        [atom.particle_index for atom in atoms])
                    # Retrieve parameters
                    [charge0, sigma0, epsilon0] = force.getParticleParameters(
                        particle_indices[0])
                    [charge1, sigma1, epsilon1] = force.getParticleParameters(
                        particle_indices[1])
                    # Apply bond charge increment
                    charge0 -= bond.increment
                    charge1 += bond.increment
                    # Update charges
                    force.setParticleParameters(particle_indices[0], charge0,
                                                sigma0, epsilon0)
                    force.setParticleParameters(particle_indices[1], charge1,
                                                sigma1, epsilon1)
                    # TODO: Calculate exceptions


class ChargeIncrementModelHandler(ParameterHandler):
    """Handle SMIRNOFF ``<ChargeIncrementModel>`` tags"""

    class ChargeIncrementType(ParameterType):
        """A SMIRNOFF bond charge correction type."""
        _VALENCE_TYPE = 'Bond'  # ChemicalEnvironment valence type expected for SMARTS
        _SMIRNOFF_ATTRIBS = ['smirks', 'chargeIncrement']
        _REQUIRE_UNITS = {
            'chargeIncrement': unit.elementary_charge
        }
        _INDEXED_ATTRIBS = ['chargeIncrement']

        def __init__(self, node, parent):
            super().__init__(**kwargs)

    _TAGNAME = 'ChargeIncrementModel'  # SMIRNOFF tag name to process
    _INFOTYPE = ChargeIncrementType  # info type to store
    _OPENMMTYPE = openmm.NonbondedForce  # OpenMM force class to create or utilize
    # TODO: The structure of this is still undecided
    _KWARGS = ['charge_from_molecules']
    _DEFAULTS = {'number_of_conformers': 10,
                 'quantum_chemical_method': 'AM1',
                 'partial_charge_method': 'CM2'}
    _ALLOWED_VALUES = {'quantum_chemical_method': ['AM1'],
                       'partial_charge_method': ['CM2']}



    def __init__(self, **kwargs):
        raise NotImplementedError("ChangeIncrementHandler is not yet implemented, pending finalization of the "
                                  "SMIRNOFF spec")
        super().__init__(**kwargs)

        if number_of_conformers is None:
            self._number_of_conformers = self._DEFAULTS['number_of_conformers']
        elif type(number_of_conformers) is str:
            self._number_of_conformers = int(number_of_conformers)
        else:
            self._number_of_conformers = number_of_conformers

        if quantum_chemical_method is None:
            self._quantum_chemical_method = self._DEFAULTS['quantum_chemical_method']
        elif number_of_conformers in self._ALLOWED_VALUES['quantum_chemical_method']:
            self._number_of_conformers = number_of_conformers

        if partial_charge_method is None:
            self._partial_charge_method = self._DEFAULTS['partial_charge_method']
        elif partial_charge_method in self._ALLOWED_VALUES['partial_charge_method']:
            self._partial_charge_method = partial_charge_method



    def check_handler_compatibility(self, handler_kwargs, assume_missing_is_default=True):
        """
        Checks if a set of kwargs used to create a ParameterHandler are compatible with this ParameterHandler. This is
        called if a second handler is attempted to be initialized for the same tag. If no value is given for a field, it
        will be assumed to expect the ParameterHandler class default.

        Parameters
        ----------
        handler_kwargs : dict
            The kwargs that would be used to construct a ParameterHandler
        assume_missing_is_default : bool
            If True, will assume that parameters not specified in handler_kwargs would have been set to the default.
            Therefore, an exception is raised if the ParameterHandler is incompatible with the default value for a
            unspecified field.

        Raises
        ------
        IncompatibleParameterError if handler_kwargs are incompatible with existing parameters.
        """
        compare_kwarg_to_attr = {
            'number_of_conformers': self._number_of_conformers,
            'quantum_chemical_method': self._quantum_chemical_method,
            'partial_charge_method': self._partial_charge_method,
        }

        for kwarg_key, attr in compare_kwarg_to_attr.items():
            # Skip this comparison if the kwarg isn't in handler_kwargs and we're not comparing against defaults
            if not(assume_missing_is_default) and not(kwarg_key in handler_kwargs.keys()):
                continue

            kwarg_val = handler_kwargs.get(kwarg_key, self._DEFAULTS[kwarg_key])
            if kwarg_val != attr:
                raise IncompatibleParameterError(
                    "Incompatible '{}' values found during handler compatibility check."
                    "(existing handler value: {}, new existing value: {}".format(kwarg_key, attr, kwarg_val))

    def assign_charge_from_molecules(self, molecule, charge_mols):
        """
        Given an input molecule, checks against a list of molecules for an isomorphic match. If found, assigns
        partial charges from the match to the input molecule.

        Parameters
        ----------
        molecule : an openforcefield.topology.FrozenMolecule
            The molecule to have partial charges assigned if a match is found.
        charge_mols : list of [openforcefield.topology.FrozenMolecule]
            A list of molecules with charges already assigned.

        Returns
        -------
        match_found : bool
            Whether a match was found. If True, the input molecule will have been modified in-place.
        """

        from networkx.algorithms.isomorphism import GraphMatcher
        # Define the node/edge attributes that we will use to match the atoms/bonds during molecule comparison
        node_match_func = lambda x, y: ((x['atomic_number'] == y['atomic_number']) and
                                        (x['stereochemistry'] == y['stereochemistry']) and
                                        (x['is_aromatic'] == y['is_aromatic'])
                                        )
        edge_match_func = lambda x, y: ((x['order'] == y['order']) and
                                        (x['stereochemistry'] == y['stereochemistry']) and
                                        (x['is_aromatic'] == y['is_aromatic'])
                                        )
        # Check each charge_mol for whether it's isomorphic to the input molecule
        for charge_mol in charge_mols:
            if molecule.is_isomorphic(charge_mol):
                # Take the first valid atom indexing map
                ref_mol_G = molecule.to_networkx()
                charge_mol_G = charge_mol.to_networkX()
                GM = GraphMatcher(
                    charge_mol_G,
                    ref_mol_G,
                    node_match=node_match_func,
                    edge_match=edge_match_func)
                for mapping in GM.isomorphisms_iter():
                    topology_atom_map = mapping
                    break
                # Set the partial charges
                charge_mol_charges = charge_mol.get_partial_charges()
                temp_mol_charges = charge_mol_charges.copy()
                for charge_idx, ref_idx in topology_atom_map:
                    temp_mol_charges[ref_idx] = charge_mol_charges[charge_idx]
                molecule.set_partial_charges(temp_mol_charges)
                return True

        # If no match was found, return False
        return False

    def create_force(self, system, topology, **kwargs):


        from openforcefield.topology import FrozenMolecule, TopologyAtom, TopologyVirtualSite

        existing = [system.getForce(i) for i in range(system.getNumForces())]
        existing = [f for f in existing if type(f) == self._OPENMMTYPE]
        if len(existing) == 0:
            force = self._OPENMMTYPE()
            system.addForce(force)
        else:
            force = existing[0]

        for ref_mol in topology.reference_molecules:

            # Make a temporary copy of ref_mol to assign charges from charge_mol
            temp_mol = FrozenMolecule(ref_mol)

            # First, check whether any of the reference molecules in the topology are in the charge_from_mol list
            charges_from_charge_mol = False
            if 'charge_from_mol' in kwargs:
                charges_from_charge_mol = self.assign_charge_from_molecules(temp_mol, kwargs['charge_from_mol'])

            # If the molecule wasn't assigned parameters from a manually-input charge_mol, calculate them here
            if not(charges_from_charge_mol):
                temp_mol.generate_conformers(num_conformers=10)
                temp_mol.compute_partial_charges(quantum_chemical_method=self._quantum_chemical_method,
                                                 partial_charge_method=self._partial_charge_method)

            # Assign charges to relevant atoms
            for topology_molecule in topology._reference_molecule_to_topology_molecules[ref_mol]:
                for topology_particle in topology_molecule.particles:
                    topology_particle_index = topology_particle.topology_particle_index
                    if type(topology_particle) is TopologyAtom:
                        ref_mol_particle_index = topology_particle.atom.molecule_particle_index
                    if type(topology_particle) is TopologyVirtualSite:
                        ref_mol_particle_index = topology_particle.virtual_site.molecule_particle_index
                    particle_charge = temp_mol._partial_charges[ref_mol_particle_index]

                    # Retrieve nonbonded parameters for reference atom (charge not set yet)
                    _, sigma, epsilon = force.getParticleParameters(topology_particle_index)
                    # Set the nonbonded force with the partial charge
                    force.setParticleParameters(topology_particle_index,
                                                particle_charge, sigma,
                                                epsilon)



    # TODO: Move chargeModel and library residue charges to SMIRNOFF spec
    def postprocess_system(self, system, topology, **kwargs):
        bonds = self.get_matches(topology)

        # Apply bond charge increments to all appropriate force groups
        # QUESTION: Should we instead apply this to the Topology in a preprocessing step, prior to spreading out charge onto virtual sites?
        for force in system.getForces():
            if force.__class__.__name__ in [
                    'NonbondedForce'
            ]:  # TODO: We need to apply this to all Force types that involve charges, such as (Custom)GBSA forces and CustomNonbondedForce
                for (atoms, bond) in bonds.items():
                    # Get corresponding particle indices in Topology
                    particle_indices = tuple(
                        [atom.particle_index for atom in atoms])
                    # Retrieve parameters
                    [charge0, sigma0, epsilon0] = force.getParticleParameters(
                        particle_indices[0])
                    [charge1, sigma1, epsilon1] = force.getParticleParameters(
                        particle_indices[1])
                    # Apply bond charge increment
                    charge0 -= bond.increment
                    charge1 += bond.increment
                    # Update charges
                    force.setParticleParameters(particle_indices[0], charge0,
                                                sigma0, epsilon0)
                    force.setParticleParameters(particle_indices[1], charge1,
                                                sigma1, epsilon1)
                    # TODO: Calculate exceptions


class GBSAParameterHandler(ParameterHandler):
    """Handle SMIRNOFF ``<GBSAParameterHandler>`` tags"""
    # TODO: Differentiate between global and per-particle parameters for each model.

    # Global parameters for surface area (SA) component of model
    SA_expected_parameters = {
        'ACE': ['surface_area_penalty', 'solvent_radius'],
        None: [],
    }

    # Per-particle parameters for generalized Born (GB) model
    GB_expected_parameters = {
        'HCT': ['radius', 'scale'],
        'OBC1': ['radius', 'scale'],
        'OBC2': ['radius', 'scale'],
    }

    class GBSAType(ParameterType):
        """A SMIRNOFF GBSA type."""
        _VALENCE_TYPE = 'Atom'
        _SMIRNOFF_ATTRIBS = ['smirks', 'radius', 'scale']
        _REQUIRE_UNITS = {'radius': unit.angstrom}
        _ATTRIBS_TO_TYPE = {'scale': float}

        def __init__(self, **kwargs):
            super().__init__(**kwargs)

            # # Store model parameters.
            # gb_model = parent.attrib['gb_model']
            # expected_parameters = GBSAParameterHandler.GB_expected_parameters[
            #     gb_model]
            # provided_parameters = list()
            # missing_parameters = list()
            # for name in expected_parameters:
            #     if name in node.attrib:
            #         provided_parameters.append(name)
            #         value = _extract_quantity_from_xml_element(
            #             node, parent, name)
            #         setattr(self, name, value)
            #     else:
            #         missing_parameters.append(name)
            # if len(missing_parameters) > 0:
            #     msg = 'GBSAForce: missing per-atom parameters for tag %s' % str(
            #         node)
            #     msg += 'model "%s" requires specification of per-atom parameters %s\n' % (
            #         gb_model, str(expected_parameters))
            #     msg += 'provided parameters : %s\n' % str(provided_parameters)
            #     msg += 'missing parameters: %s' % str(missing_parameters)
            #     raise Exception(msg)

    # TODO: Finish this
    _TAGNAME = 'GBSA'
    _INFOTYPE = GBSAType
    #_OPENMMTYPE =

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

    # TODO: Fix this
    def parseElement(self):
        # Initialize GB model
        gb_model = element.attrib['gb_model']
        valid_GB_models = GBSAParameterHandler.GB_expected_parameters.keys()
        if not gb_model in valid_GB_models:
            raise Exception(
                'Specified GBSAForce model "%s" not one of valid models: %s' %
                (gb_model, valid_GB_models))
        self.gb_model = gb_model

        # Initialize SA model
        sa_model = element.attrib['sa_model']
        valid_SA_models = GBSAParameterHandler.SA_expected_parameters.keys()
        if not sa_model in valid_SA_models:
            raise Exception(
                'Specified GBSAForce SA_model "%s" not one of valid models: %s'
                % (sa_model, valid_SA_models))
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
            raise ValueError(
                'Found multiple GBSAForce tags with different GB model specifications'
            )
        if (Handler.sa_model != self.sa_model):
            raise ValueError(
                'Found multiple GBSAForce tags with different SA model specifications'
            )
        # TODO: Check other attributes (parameters of GB and SA models) automatically?

    def create_force(self, system, topology, **args):
        # TODO: Rework this
        from openforcefield.typing.engines.smirnoff import gbsaforces
        force_class = getattr(gbsaforces, self.gb_model)
        force = force_class(**self.parameters)
        system.addForce(force)

        # Add all GBSA terms to the system.
        expected_parameters = GBSAParameterHandler.GB_expected_parameters[
            self.gb_model]

        # Create all particles with parameters set to zero
        atoms = self.getMatches(topology)
        nparams = 1 + len(expected_parameters)  # charge + GBSA parameters
        params = [0.0 for i in range(nparams)]
        for particle in topology.topology_particles():
            force.addParticle(params)
        # Set the GBSA parameters (keeping charges at zero for now)
        for (atoms, gbsa_type) in atoms.items():
            atom = atoms[0]
            # Set per-particle parameters for assigned parameters
            params = [atom.charge] + [
                getattr(gbsa_type, name) for name in expected_parameters
            ]
            force.setParticleParameters(atom.particle_index, params)
