"""
Parameter handlers for the SMIRNOFF force field engine

This file contains standard parameter handlers for the SMIRNOFF force field engine.
These classes implement the object model for self-contained parameter assignment.
New pluggable handlers can be created by creating subclasses of :class:`ParameterHandler`.

"""

__all__ = [
    "DuplicateParameterError",
    "DuplicateVirtualSiteTypeException",
    "FractionalBondOrderInterpolationMethodUnsupportedError",
    "IncompatibleParameterError",
    "NonintegralMoleculeChargeException",
    "NotEnoughPointsForInterpolationError",
    "ParameterLookupError",
    "SMIRNOFFSpecError",
    "SMIRNOFFSpecUnimplementedError",
    "UnassignedAngleParameterException",
    "UnassignedBondParameterException",
    "UnassignedMoleculeChargeException",
    "UnassignedProperTorsionParameterException",
    "UnassignedValenceParameterException",
    "NonbondedMethod",
    "ParameterList",
    "ParameterType",
    "ParameterHandler",
    "ParameterAttribute",
    "MappedParameterAttribute",
    "IndexedParameterAttribute",
    "IndexedMappedParameterAttribute",
    "ConstraintHandler",
    "BondHandler",
    "AngleHandler",
    "ProperTorsionHandler",
    "ImproperTorsionHandler",
    "ElectrostaticsHandler",
    "LibraryChargeHandler",
    "vdWHandler",
    "GBSAHandler",
    "ToolkitAM1BCCHandler",
    "VirtualSiteHandler",
]

import abc
import copy
import functools
import inspect
import logging
import re
from collections import OrderedDict, defaultdict
from enum import Enum
from itertools import combinations

try:
    import openmm
    from openmm import unit
except ImportError:
    from simtk import openmm, unit

from openff.toolkit.topology import (
    ImproperDict,
    TagSortedDict,
    Topology,
    UnsortedDict,
    ValenceDict,
)
from openff.toolkit.topology.molecule import Molecule
from openff.toolkit.topology.topology import NotBondedError
from openff.toolkit.typing.chemistry import ChemicalEnvironment
from openff.toolkit.utils.collections import ValidatedDict, ValidatedList
from openff.toolkit.utils.exceptions import (
    DuplicateParameterError,
    DuplicateVirtualSiteTypeException,
    FractionalBondOrderInterpolationMethodUnsupportedError,
    IncompatibleParameterError,
    NonintegralMoleculeChargeException,
    NotEnoughPointsForInterpolationError,
    ParameterLookupError,
    SMIRNOFFSpecError,
    SMIRNOFFSpecUnimplementedError,
    UnassignedAngleParameterException,
    UnassignedBondParameterException,
    UnassignedMoleculeChargeException,
    UnassignedProperTorsionParameterException,
    UnassignedValenceParameterException,
)
from openff.toolkit.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY
from openff.toolkit.utils.utils import (
    IncompatibleUnitError,
    all_subclasses,
    attach_units,
    extract_serialized_units_from_dict,
    object_to_quantity,
)

logger = logging.getLogger(__name__)


class NonbondedMethod(Enum):
    """
    An enumeration of the nonbonded methods
    """

    NoCutoff = 0
    CutoffPeriodic = 1
    CutoffNonPeriodic = 2
    Ewald = 3
    PME = 4


def _linear_inter_or_extrapolate(points_dict, x_query):
    """
    Linearly interpolate or extrapolate based on a piecewise linear function defined by a set of points.
    This function is designed to work with key:value pairs where the value may be a openmm.unit.Quantity.

    Parameters
    ----------
    points_dict : dict{float: float or float-valued openmm.unit.Quantity}
        A dictionary with each item representing a point, where the key is the X value and the value is the Y value.
    x_query : float
        The X value of the point to interpolate or extrapolate.

    Returns
    -------
    y_value : float or float-valued openmm.unit.Quantity
        The result of interpolation/extrapolation.
    """

    # pre-empt case where no interpolation is necessary
    if x_query in points_dict:
        return points_dict[x_query]

    if len(points_dict) < 2:
        raise NotEnoughPointsForInterpolationError(
            f"Unable to perform interpolation with less than two points. "
            f"points_dict: {points_dict}   x_query: {x_query}"
        )
    # TODO: error out for nonsensical fractional bond orders

    # find the nearest point beneath our queried x value
    try:
        below = max(bo for bo in points_dict if bo < x_query)
    except ValueError:
        below = None
    # find the nearest point above our queried x value
    try:
        above = min(bo for bo in points_dict if bo > x_query)
    except ValueError:
        above = None

    # handle case where we can clearly interpolate
    if (above is not None) and (below is not None):
        return points_dict[below] + (points_dict[above] - points_dict[below]) * (
            (x_query - below) / (above - below)
        )

    # error if we can't hope to interpolate at all
    elif (above is None) and (below is None):
        raise NotImplementedError(
            f"Failed to find interpolation references for "
            f"`x_query` '{x_query}', "
            f"with `points_dict` '{points_dict}'"
        )

    # extrapolate for fractional bond orders below our lowest defined bond order
    elif below is None:
        bond_orders = sorted(points_dict)
        y0 = points_dict[bond_orders[0]]
        y_diff = points_dict[bond_orders[1]] - points_dict[bond_orders[0]]
        x_diff = bond_orders[1] - bond_orders[0]
        k = y0 - y_diff / x_diff * (bond_orders[0] - x_query)
        return k

    # extrapolate for fractional bond orders above our highest defined bond order
    elif above is None:
        bond_orders = sorted(points_dict)
        y0 = points_dict[bond_orders[-1]]
        y_diff = points_dict[bond_orders[-1]] - points_dict[bond_orders[-2]]
        x_diff = bond_orders[-1] - bond_orders[-2]
        k = y0 + y_diff / x_diff * (x_query - bond_orders[-1])
        return k


# TODO: This is technically a validator, not a converter, but ParameterAttribute doesn't support them yet (it'll be easy if we switch to use the attrs library).
def _allow_only(allowed_values):
    """A converter that checks the new value is only in a set."""
    allowed_values = frozenset(allowed_values)

    def _value_checker(instance, attr, new_value):
        # This statement means that, in the "SMIRNOFF Data Dict" format, the string "None"
        # and the Python None are the same thing
        if new_value == "None":
            new_value = None

        # Ensure that the new value is in the list of allowed values
        if new_value not in allowed_values:

            err_msg = (
                f"Attempted to set {instance.__class__.__name__}.{attr.name} "
                f"to {new_value}. Currently, only the following values "
                f"are supported: {sorted(allowed_values)}."
            )
            raise SMIRNOFFSpecError(err_msg)
        return new_value

    return _value_checker


# TODO: Think about adding attrs to the dependencies and inherit from attr.ib
class ParameterAttribute:
    """A descriptor for ``ParameterType`` attributes.

    The descriptors allows associating to the parameter a default value,
    which makes the attribute optional, a unit, and a custom converter.

    Because we may want to have ``None`` as a default value, required
    attributes have the ``default`` set to the special type ``UNDEFINED``.

    Converters can be both static or instance functions/methods with
    respective signatures::

        converter(value): -> converted_value
        converter(instance, parameter_attribute, value): -> converted_value

    A decorator syntax is available (see example below).

    Parameters
    ----------
    default : object, optional
        When specified, the descriptor makes this attribute optional by
        attaching a default value to it.
    unit : openmm.unit.Quantity, optional
        When specified, only quantities with compatible units are allowed
        to be set, and string expressions are automatically parsed into a
        ``Quantity``.
    converter : callable, optional
        An optional function that can be used to convert values before
        setting the attribute.

    See Also
    --------
    IndexedParameterAttribute
        A parameter attribute with multiple terms.

    Examples
    --------

    Create a parameter type with an optional and a required attribute.

    >>> class MyParameter:
    ...     attr_required = ParameterAttribute()
    ...     attr_optional = ParameterAttribute(default=2)
    ...
    >>> my_par = MyParameter()

    Even without explicit assignment, the default value is returned.

    >>> my_par.attr_optional
    2

    If you try to access an attribute without setting it first, an
    exception is raised.

    >>> my_par.attr_required
    Traceback (most recent call last):
    ...
    AttributeError: 'MyParameter' object has no attribute '_attr_required'

    The attribute allow automatic conversion and validation of units.

    >>> from openmm import unit
    >>> class MyParameter:
    ...     attr_quantity = ParameterAttribute(unit=unit.angstrom)
    ...
    >>> my_par = MyParameter()
    >>> my_par.attr_quantity = '1.0 * nanometer'
    >>> my_par.attr_quantity
    Quantity(value=1.0, unit=nanometer)
    >>> my_par.attr_quantity = 3.0
    Traceback (most recent call last):
    ...
    openff.toolkit.utils.utils.IncompatibleUnitError: attr_quantity=3.0 dimensionless should have units of angstrom

    You can attach a custom converter to an attribute.

    >>> class MyParameter:
    ...     # Both strings and integers convert nicely to floats with float().
    ...     attr_all_to_float = ParameterAttribute(converter=float)
    ...     attr_int_to_float = ParameterAttribute()
    ...     @attr_int_to_float.converter
    ...     def attr_int_to_float(self, attr, value):
    ...         # This converter converts only integers to float
    ...         # and raise an exception for the other types.
    ...         if isinstance(value, int):
    ...             return float(value)
    ...         elif not isinstance(value, float):
    ...             raise TypeError(f"Cannot convert '{value}' to float")
    ...         return value
    ...
    >>> my_par = MyParameter()

    attr_all_to_float accepts and convert to float both strings and integers

    >>> my_par.attr_all_to_float = 1
    >>> my_par.attr_all_to_float
    1.0
    >>> my_par.attr_all_to_float = '2.0'
    >>> my_par.attr_all_to_float
    2.0

    The custom converter associated to attr_int_to_float converts only integers instead.

    >>> my_par.attr_int_to_float = 3
    >>> my_par.attr_int_to_float
    3.0
    >>> my_par.attr_int_to_float = '4.0'
    Traceback (most recent call last):
    ...
    TypeError: Cannot convert '4.0' to float

    """

    class UNDEFINED:
        """Custom type used by ``ParameterAttribute`` to differentiate between ``None`` and undeclared default."""

        pass

    def __init__(self, default=UNDEFINED, unit=None, converter=None, docstring=""):
        self.default = default
        self._unit = unit
        self._converter = converter
        self.__doc__ = docstring

    def __set_name__(self, owner, name):
        self._name = "_" + name

    @property
    def name(self):
        # Get rid of the initial underscore.
        return self._name[1:]

    def __get__(self, instance, owner):
        if instance is None:
            # This is called from the class. Return the descriptor object.
            return self

        try:
            return getattr(instance, self._name)
        except AttributeError:
            # The attribute has not initialized. Check if there's a default.
            if self.default is ParameterAttribute.UNDEFINED:
                raise
            return self.default

    def __set__(self, instance, value):
        # Convert and validate the value.
        value = self._convert_and_validate(instance, value)
        setattr(instance, self._name, value)

    def converter(self, converter):
        """Create a new ParameterAttribute with an associated converter.

        This is meant to be used as a decorator (see main examples).
        """
        return self.__class__(default=self.default, converter=converter)

    def _convert_and_validate(self, instance, value):
        """Convert to Quantity, validate units, and call custom converter."""
        # The default value is always allowed.
        if self._is_valid_default(value):
            return value
        # Convert and validate units.
        value = self._validate_units(value)
        # Call the custom converter before setting the value.
        value = self._call_converter(value, instance)
        return value

    def _is_valid_default(self, value):
        """Return True if this is a defined default value."""
        return (
            self.default is not ParameterAttribute.UNDEFINED and value == self.default
        )

    def _validate_units(self, value):
        """Convert strings expressions to Quantity and validate the units if requested."""
        if self._unit is not None:
            # Convert eventual strings to Quantity objects.
            value = object_to_quantity(value)

            # Check if units are compatible.
            try:
                if not self._unit.is_compatible(value.unit):
                    raise IncompatibleUnitError(
                        f"{self.name}={value} should have units of {self._unit}"
                    )
            except AttributeError:
                # This is not a Quantity object.
                raise IncompatibleUnitError(
                    f"{self.name}={value} should have units of {self._unit}"
                )
        return value

    def _call_converter(self, value, instance):
        """Correctly calls static and instance converters."""
        if self._converter is not None:
            try:
                # Static function.
                return self._converter(value)
            except TypeError:
                # Instance method.
                return self._converter(instance, self, value)
        return value


class IndexedParameterAttribute(ParameterAttribute):
    """The attribute of a parameter with an unspecified number of terms.

    Some parameters can be associated to multiple terms, For example,
    torsions have parameters such as k1, k2, ..., and ``IndexedParameterAttribute``
    can be used to encapsulate the sequence of terms.

    The only substantial difference with ``ParameterAttribute`` is that
    only sequences are supported as values and converters and units are
    checked on each element of the sequence.

    Currently, the descriptor makes the sequence immutable. This is to
    avoid that an element of the sequence could be set without being
    properly validated. In the future, the data could be wrapped in a
    safe list that would safely allow mutability.

    Parameters
    ----------
    default : object, optional
        When specified, the descriptor makes this attribute optional by
        attaching a default value to it.
    unit : openmm.unit.Quantity, optional
        When specified, only sequences of quantities with compatible units
        are allowed to be set.
    converter : callable, optional
        An optional function that can be used to validate and cast each
        element of the sequence before setting the attribute.

    See Also
    --------
    ParameterAttribute
        A simple parameter attribute.
    MappedParameterAttribute
        A parameter attribute representing a mapping.
    IndexedMappedParameterAttribute
        A parameter attribute representing a sequence, each term of which is a mapping.

    Examples
    --------

    Create an optional indexed attribute with unit of angstrom.

    >>> from openmm import unit
    >>> class MyParameter:
    ...     length = IndexedParameterAttribute(default=None, unit=unit.angstrom)
    ...
    >>> my_par = MyParameter()
    >>> my_par.length is None
    True

    Strings are parsed into Quantity objects.

    >>> my_par.length = ['1 * angstrom', 0.5 * unit.nanometer]
    >>> my_par.length[0]
    Quantity(value=1, unit=angstrom)

    Similarly, custom converters work as with ``ParameterAttribute``, but
    they are used to validate each value in the sequence.

    >>> class MyParameter:
    ...     attr_indexed = IndexedParameterAttribute(converter=float)
    ...
    >>> my_par = MyParameter()
    >>> my_par.attr_indexed = [1, '1.0', '1e-2', 4.0]
    >>> my_par.attr_indexed
    [1.0, 1.0, 0.01, 4.0]

    """

    def _convert_and_validate(self, instance, value):
        """Overwrite ParameterAttribute._convert_and_validate to make the value a ValidatedList."""
        # The default value is always allowed.
        if self._is_valid_default(value):
            return value

        # We push the converters into a ValidatedList so that we can make
        # sure that elements are validated correctly when they are modified
        # after their initialization.
        # ValidatedList expects converters that take the value as a single
        # argument so we create a partial function with the instance assigned.
        static_converter = functools.partial(self._call_converter, instance=instance)
        value = ValidatedList(value, converter=[self._validate_units, static_converter])

        return value


class MappedParameterAttribute(ParameterAttribute):
    """The attribute of a parameter in which each term is a mapping.

    The substantial difference with ``IndexedParameterAttribute`` is that, unlike
    indexing, the mapping can be based on artbitrary references, like indices but
    can starting at non-zero values and include non-adjacent keys.

    Parameters
    ----------
    default : object, optional
        When specified, the descriptor makes this attribute optional by
        attaching a default value to it.
    unit : openmm.unit.Quantity, optional
        When specified, only sequences of mappings where values are quantities with
        compatible units are allowed to be set.
    converter : callable, optional
        An optional function that can be used to validate and cast each
        component of each element of the sequence before setting the attribute.

    See Also
    --------
    IndexedParameterAttribute
        A parameter attribute representing a sequence.
    IndexedMappedParameterAttribute
        A parameter attribute representing a sequence, each term of which is a mapping.

    Examples
    --------

    Create an optional indexed attribute with unit of angstrom.

    >>> from openmm import unit
    >>> class MyParameter:
    ...     length = MappedParameterAttribute(default=None, unit=unit.angstrom)
    ...
    >>> my_par = MyParameter()
    >>> my_par.length is None
    True

    Like other ParameterAttribute objects, strings are parsed into Quantity objects.

    >>> my_par.length = {1:'1.5 * angstrom', 2: '1.4 * angstrom'}
    >>> my_par.length[1]
    Quantity(value=1.5, unit=angstrom)

    Unlike other ParameterAttribute objects, the reference points can do not need ot be
    zero-indexed, non-adjancent, such as interpolating defining a bond parameter for
    interpolation by defining references values and bond orders 2 and 3:

    >>> my_par.length = {2:'1.42 * angstrom', 3: '1.35 * angstrom'}
    >>> my_par.length[2]
    Quantity(value=1.42, unit=angstrom)

    """

    def _convert_and_validate(self, instance, value):
        if self._is_valid_default(value):
            return value
        static_converter = functools.partial(self._call_converter, instance=instance)
        value = ValidatedDict(value, converter=[self._validate_units, static_converter])
        return value


class IndexedMappedParameterAttribute(ParameterAttribute):
    """The attribute of a parameter with an unspecified number of terms, where
    each term is a mapping.

    Some parameters can be associated to multiple terms,
    where those terms have multiple components.
    For example, torsions with fractional bond orders have parameters such as
    k1_bondorder1, k1_bondorder2, k2_bondorder1, k2_bondorder2, ..., and
    ``IndexedMappedParameterAttribute`` can be used to encapsulate the sequence of
    terms as mappings (typically, ``dict``\ s) of their components.

    The only substantial difference with ``IndexedParameterAttribute`` is that
    only sequences of mappings are supported as values and converters and units are
    checked on each component of each element in the sequence.

    Currently, the descriptor makes the sequence immutable. This is to
    avoid that an element of the sequence could be set without being
    properly validated. In the future, the data could be wrapped in a
    safe list that would safely allow mutability.

    Parameters
    ----------
    default : object, optional
        When specified, the descriptor makes this attribute optional by
        attaching a default value to it.
    unit : openmm.unit.Quantity, optional
        When specified, only sequences of mappings where values are quantities with
        compatible units are allowed to be set.
    converter : callable, optional
        An optional function that can be used to validate and cast each
        component of each element of the sequence before setting the attribute.

    See Also
    --------
    IndexedParameterAttribute
        A parameter attribute representing a sequence.
    MappedParameterAttribute
        A parameter attribute representing a mapping.

    Examples
    --------

    Create an optional indexed attribute with unit of angstrom.

    >>> from openmm import unit
    >>> class MyParameter:
    ...     length = IndexedMappedParameterAttribute(default=None, unit=unit.angstrom)
    ...
    >>> my_par = MyParameter()
    >>> my_par.length is None
    True

    Strings are parsed into Quantity objects.

    >>> my_par.length = [{1:'1 * angstrom'}, {1: 0.5 * unit.nanometer}]
    >>> my_par.length[0]
    {1: Quantity(value=1, unit=angstrom)}

    Similarly, custom converters work as with ``ParameterAttribute``, but
    they are used to validate each value in the sequence.

    >>> class MyParameter:
    ...     attr_indexed = IndexedMappedParameterAttribute(converter=float)
    ...
    >>> my_par = MyParameter()
    >>> my_par.attr_indexed = [{1: 1}, {2: '1.0', 3: '1e-2'}, {4: 4.0}]
    >>> my_par.attr_indexed
    [{1: 1.0}, {2: 1.0, 3: 0.01}, {4: 4.0}]

    """

    def _convert_and_validate(self, instance, value):
        """Overwrite ParameterAttribute._convert_and_validate to make the value a ValidatedList."""
        # The default value is always allowed.
        if self._is_valid_default(value):
            return value

        # We push the converters into a ValidatedListMapping so that we can make
        # sure that elements are validated correctly when they are modified
        # after their initialization.
        # ValidatedListMapping expects converters that take the value as a single
        # argument so we create a partial function with the instance assigned.
        static_converter = functools.partial(self._call_converter, instance=instance)

        value = ValidatedList(
            [
                ValidatedDict(
                    element, converter=[self._validate_units, static_converter]
                )
                for element in value
            ],
            converter=self._index_converter,
        )

        return value

    @staticmethod
    def _index_converter(x):
        return ValidatedDict(x)


class _ParameterAttributeHandler:
    """A base class for ``ParameterType`` and ``ParameterHandler`` objects.

    Encapsulate shared code of ``ParameterType`` and ``ParameterHandler``.
    In particular, this base class provides an ``__init__`` method that
    automatically initialize the attributes defined through the ``ParameterAttribute``
    and ``IndexedParameterAttribute`` descriptors, as well as handling
    cosmetic attributes.

    See Also
    --------
    ParameterAttribute
        A simple parameter attribute.
    IndexedParameterAttribute
        A parameter attribute with multiple terms.

    Examples
    --------

    This base class was design to encapsulate shared code between ``ParameterType``
    and ``ParameterHandler``, which both need to deal with parameter and cosmetic
    attributes.

    To create a new type/handler, you can use the ``ParameterAttribute`` descriptors.

    >>> class ParameterTypeOrHandler(_ParameterAttributeHandler):
    ...     length = ParameterAttribute(unit=unit.angstrom)
    ...     k = ParameterAttribute(unit=unit.kilocalorie_per_mole / unit.angstrom**2)
    ...

    ``_ParameterAttributeHandler`` and the descriptors take care of performing
    sanity checks on initialization and assignment of the single attributes. Because
    we attached units to the parameters, we need to pass them with compatible units.

    >>> my_par = ParameterTypeOrHandler(
    ...     length='1.01 * angstrom',
    ...     k=5 * unit.kilocalorie_per_mole / unit.angstrom**2
    ... )

    Note that ``_ParameterAttributeHandler`` took care of implementing
    a constructor, and that unit parameters support string assignments.
    These are automatically converted to ``Quantity`` objects.

    >>> my_par.length
    Quantity(value=1.01, unit=angstrom)

    While assigning incompatible units is forbidden.

    >>> my_par.k = 3.0 * unit.gram
    Traceback (most recent call last):
    ...
    openff.toolkit.utils.utils.IncompatibleUnitError: k=3.0 g should have units of kilocalorie/(angstrom**2*mole)

    On top of type checking, the constructor implemented in ``_ParameterAttributeHandler``
    checks if some required parameters are not given.

    >>> ParameterTypeOrHandler(length=3.0*unit.nanometer)
    Traceback (most recent call last):
    ...
    openff.toolkit.typing.engines.smirnoff.parameters.SMIRNOFFSpecError: <class 'openff.toolkit.typing.engines.smirnoff.parameters.ParameterTypeOrHandler'> require the following missing parameters: ['k']. Defined kwargs are ['length']

    Each attribute can be made optional by specifying a default value,
    and you can attach a converter function by passing a callable as an
    argument or through the decorator syntax.

    >>> class ParameterTypeOrHandler(_ParameterAttributeHandler):
    ...     attr_optional = ParameterAttribute(default=2)
    ...     attr_all_to_float = ParameterAttribute(converter=float)
    ...     attr_int_to_float = ParameterAttribute()
    ...
    ...     @attr_int_to_float.converter
    ...     def attr_int_to_float(self, attr, value):
    ...         # This converter converts only integers to floats
    ...         # and raise an exception for the other types.
    ...         if isinstance(value, int):
    ...             return float(value)
    ...         elif not isinstance(value, float):
    ...             raise TypeError(f"Cannot convert '{value}' to float")
    ...         return value
    ...
    >>> my_par = ParameterTypeOrHandler(attr_all_to_float='3.0', attr_int_to_float=1)
    >>> my_par.attr_optional
    2
    >>> my_par.attr_all_to_float
    3.0
    >>> my_par.attr_int_to_float
    1.0

    The float() function can convert strings to integers, but our custom
    converter forbids it

    >>> my_par.attr_all_to_float = '2.0'
    >>> my_par.attr_int_to_float = '4.0'
    Traceback (most recent call last):
    ...
    TypeError: Cannot convert '4.0' to float

    Parameter attributes that can be indexed can be handled with the
    ``IndexedParameterAttribute``. These support unit validation and
    converters exactly as ``ParameterAttribute``s, but the validation/conversion
    is performed for each indexed attribute.

    >>> class MyTorsionType(_ParameterAttributeHandler):
    ...     periodicity = IndexedParameterAttribute(converter=int)
    ...     k = IndexedParameterAttribute(unit=unit.kilocalorie_per_mole)
    ...
    >>> my_par = MyTorsionType(
    ...     periodicity1=2,
    ...     k1=5 * unit.kilocalorie_per_mole,
    ...     periodicity2='3',
    ...     k2=6 * unit.kilocalorie_per_mole,
    ... )
    >>> my_par.periodicity
    [2, 3]

    Indexed attributes, can be accessed both as a list or as their indexed
    parameter name.

    >>> my_par.periodicity2 = 6
    >>> my_par.periodicity[0] = 1
    >>> my_par.periodicity
    [1, 6]

    """

    def __init__(self, allow_cosmetic_attributes=False, **kwargs):
        """
        Initialize parameter and cosmetic attributes.

        Parameters
        ----------
        allow_cosmetic_attributes : bool optional. Default = False
            Whether to permit non-spec kwargs ("cosmetic attributes").
            If True, non-spec kwargs will be stored as an attribute of
            this parameter which can be accessed and written out. Otherwise,
            an exception will be raised.

        """
        # A list that may be populated to record the cosmetic attributes
        # read from a SMIRNOFF data source.
        self._cosmetic_attribs = []

        # Do not modify the original data.
        smirnoff_data = copy.deepcopy(kwargs)

        (
            smirnoff_data,
            indexed_mapped_attr_lengths,
        ) = self._process_indexed_mapped_attributes(smirnoff_data)
        smirnoff_data = self._process_indexed_attributes(
            smirnoff_data, indexed_mapped_attr_lengths
        )

        smirnoff_data = self._process_mapped_attributes(smirnoff_data)

        # Check for missing required arguments.
        given_attributes = set(smirnoff_data.keys())
        required_attributes = set(self._get_required_parameter_attributes().keys())
        missing_attributes = required_attributes.difference(given_attributes)
        if len(missing_attributes) != 0:
            msg = (
                f"{self.__class__} require the following missing parameters: {sorted(missing_attributes)}."
                f" Defined kwargs are {sorted(smirnoff_data.keys())}"
            )
            raise SMIRNOFFSpecError(msg)

        # Finally, set attributes of this ParameterType and handle cosmetic attributes.
        allowed_attributes = set(self._get_parameter_attributes().keys())
        for key, val in smirnoff_data.items():
            if key in allowed_attributes:
                setattr(self, key, val)
            # Handle all unknown kwargs as cosmetic so we can write them back out
            elif allow_cosmetic_attributes:
                self.add_cosmetic_attribute(key, val)
            else:
                msg = (
                    f"Unexpected kwarg ({key}: {val})  passed to {self.__class__} constructor. "
                    "If this is a desired cosmetic attribute, consider setting "
                    "'allow_cosmetic_attributes=True'"
                )
                raise SMIRNOFFSpecError(msg)

    def _process_mapped_attributes(self, smirnoff_data):
        kwargs = list(smirnoff_data.keys())
        for kwarg in kwargs:
            attr_name, key = self._split_attribute_mapping(kwarg)

            # Check if this is a mapped attribute
            if key is not None and attr_name in self._get_mapped_parameter_attributes():
                if attr_name not in smirnoff_data:
                    smirnoff_data[attr_name] = dict()

                smirnoff_data[attr_name][key] = smirnoff_data[kwarg]
                del smirnoff_data[kwarg]

        return smirnoff_data

    def _process_indexed_mapped_attributes(self, smirnoff_data):
        # TODO: construct data structure for holding indexed_mapped attrs, which
        # will get fed into setattr
        indexed_mapped_attr_lengths = {}
        reindex = set()
        reverse = defaultdict(dict)

        kwargs = list(smirnoff_data.keys())
        for kwarg in kwargs:
            attr_name, index, key = self._split_attribute_index_mapping(kwarg)

            # Check if this is an indexed_mapped attribute.
            if (
                (key is not None)
                and (index is not None)
                and attr_name in self._get_indexed_mapped_parameter_attributes()
            ):

                # we start with a dict because have no guarantee of order
                # in which we will see each kwarg
                # we'll switch this to a list later
                if attr_name not in smirnoff_data:
                    smirnoff_data[attr_name] = dict()
                    reindex.add(attr_name)

                if index not in smirnoff_data[attr_name]:
                    smirnoff_data[attr_name][index] = dict()

                smirnoff_data[attr_name][index][key] = smirnoff_data[kwarg]
                del smirnoff_data[kwarg]

                # build reverse mapping; needed for contiguity check below
                if index not in reverse[attr_name]:
                    reverse[attr_name][index] = dict()
                reverse[attr_name][index][key] = kwarg

        # turn all our top-level dicts into lists
        # catch cases where we skip an index,
        # e.g. k1_bondorder*, k3_bondorder* defined, but not k2_bondorder*
        for attr_name in reindex:
            indexed_mapping = []
            j = 0
            for i in sorted(smirnoff_data[attr_name].keys()):
                if int(i) == j:
                    indexed_mapping.append(smirnoff_data[attr_name][i])
                    j += 1
                else:
                    # any key will do; we are sensitive only to top-level index
                    key = sorted(reverse[attr_name][i].keys())[0]
                    kwarg = reverse[attr_name][i][key]
                    val = smirnoff_data[attr_name][i][key]

                    msg = (
                        f"Unexpected kwarg ({kwarg}: {val})  passed to {self.__class__} constructor. "
                        "If this is a desired cosmetic attribute, consider setting "
                        "'allow_cosmetic_attributes=True'"
                    )
                    raise SMIRNOFFSpecError(msg)

            smirnoff_data[attr_name] = indexed_mapping

            # keep track of lengths; used downstream for checking against other
            # indexed attributes
            indexed_mapped_attr_lengths[attr_name] = len(smirnoff_data[attr_name])

        return smirnoff_data, indexed_mapped_attr_lengths

    def _process_indexed_attributes(self, smirnoff_data, indexed_attr_lengths=None):
        # Check for indexed attributes and stack them into a list.
        # Keep track of how many indexed attribute we find to make sure they all have the same length.

        # TODO: REFACTOR ME; try looping over contents of `smirnoff_data`, using
        # `split_attribute_index` to extract values

        if indexed_attr_lengths is None:
            indexed_attr_lengths = {}

        for attrib_basename in self._get_indexed_parameter_attributes().keys():
            index = 1
            while True:
                attrib_w_index = "{}{}".format(attrib_basename, index)

                # Exit the while loop if the indexed attribute is not given.
                # this is the stop condition
                try:
                    attrib_w_index_value = smirnoff_data[attrib_w_index]
                except KeyError:
                    break

                # Check if this is the first iteration.
                if index == 1:
                    # Check if this attribute has been specified with and without index.
                    if attrib_basename in smirnoff_data:
                        err_msg = (
                            f"The attribute '{attrib_basename}' has been specified "
                            f"with and without index: '{attrib_w_index}'"
                        )
                        raise TypeError(err_msg)

                    # Otherwise create the list object.
                    smirnoff_data[attrib_basename] = list()

                # Append the new value to the list.
                smirnoff_data[attrib_basename].append(attrib_w_index_value)

                # Remove the indexed attribute from the kwargs as it will
                # be exposed only as an element of the list.
                del smirnoff_data[attrib_w_index]
                index += 1

            # Update the lengths with this attribute (if it was found).
            if index > 1:
                indexed_attr_lengths[attrib_basename] = len(
                    smirnoff_data[attrib_basename]
                )

        # Raise an error if we there are different indexed
        # attributes with a different number of terms.
        if len(set(indexed_attr_lengths.values())) > 1:
            raise TypeError(
                "The following indexed attributes have "
                f"different lengths: {indexed_attr_lengths}"
            )

        return smirnoff_data

    def to_dict(self, discard_cosmetic_attributes=False, duplicate_attributes=None):
        """
        Convert this object to dict format.

        The returning dictionary contains all the ``ParameterAttribute``
        and ``IndexedParameterAttribute`` as well as cosmetic attributes
        if ``discard_cosmetic_attributes`` is ``False``.

        Parameters
        ----------
        discard_cosmetic_attributes : bool, optional. Default = False
            Whether to discard non-spec attributes of this object
        duplicate_attributes : list of string, optional. Default = None
            A list of names of attributes that redundantly decsribe
            data and should be discarded during serializaiton

        Returns
        -------
        smirnoff_dict : dict
            The SMIRNOFF-compliant dict representation of this object.

        """
        # Make a list of all attribs that should be included in the
        # returned dict (call list() to make a copy). We discard
        # optional attributes that are set to None defaults.
        attribs_to_return = list(self._get_defined_parameter_attributes().keys())

        if duplicate_attributes is not None:
            for duplicate in duplicate_attributes:
                attribs_to_return.pop(attribs_to_return.index(duplicate))

        # Start populating a dict of the attribs.
        indexed_attribs = set(self._get_indexed_parameter_attributes().keys())
        mapped_attribs = set(self._get_mapped_parameter_attributes().keys())
        indexed_mapped_attribs = set(
            self._get_indexed_mapped_parameter_attributes().keys()
        )
        smirnoff_dict = OrderedDict()

        # If attribs_to_return is ordered here, that will effectively be an informal output ordering
        for attrib_name in attribs_to_return:
            attrib_value = getattr(self, attrib_name)

            if attrib_name in indexed_mapped_attribs:
                for idx, mapping in enumerate(attrib_value):
                    for key, val in mapping.items():
                        attrib_name_indexed, attrib_name_mapped = attrib_name.split("_")
                        smirnoff_dict[
                            f"{attrib_name_indexed}{str(idx+1)}_{attrib_name_mapped}{key}"
                        ] = val
            elif attrib_name in indexed_attribs:
                for idx, val in enumerate(attrib_value):
                    smirnoff_dict[attrib_name + str(idx + 1)] = val
            elif attrib_name in mapped_attribs:
                for key, val in attrib_value.items():
                    smirnoff_dict[f"{attrib_name}{str(key)}"] = val
            else:
                smirnoff_dict[attrib_name] = attrib_value

        # Serialize cosmetic attributes.
        if not (discard_cosmetic_attributes):
            for cosmetic_attrib in self._cosmetic_attribs:
                smirnoff_dict[cosmetic_attrib] = getattr(self, "_" + cosmetic_attrib)

        return smirnoff_dict

    def __getattr__(self, item):
        """Take care of mapping indexed attributes to their respective list elements."""

        # Try matching the case where there are two indices
        # this indicates a index_mapped parameter
        attr_name, index, key = self._split_attribute_index_mapping(item)

        # Check if this is an indexed_mapped attribute.
        if (
            (key is not None)
            and (index is not None)
            and attr_name in self._get_indexed_mapped_parameter_attributes()
        ):
            indexed_mapped_attr_value = getattr(self, attr_name)
            try:
                return indexed_mapped_attr_value[index][key]
            except (IndexError, KeyError) as err:
                if not err.args:
                    err.args = ("",)
                err.args = err.args + (
                    f"'{item}' is out of bound for indexed attribute '{attr_name}'",
                )
                raise

        # Otherwise, try indexed attribute
        # Separate the indexed attribute name from the list index.
        attr_name, index = self._split_attribute_index(item)

        # Check if this is an indexed attribute.
        if (
            index is not None
        ) and attr_name in self._get_indexed_parameter_attributes():
            indexed_attr_value = getattr(self, attr_name)
            try:
                return indexed_attr_value[index]
            except IndexError:
                raise IndexError(
                    f"'{item}' is out of bound for indexed attribute '{attr_name}'"
                )

        # Otherwise, forward the search to the next class in the MRO.
        try:
            return super().__getattr__(item)
        except AttributeError as e:
            # If this fails because the next classes in the MRO do not
            # implement __getattr__(), then raise the standard Attribute error.
            if "__getattr__" in str(e):
                raise AttributeError(
                    f"{self.__class__} object has no attribute '{item}'"
                )
            # Otherwise, re-raise the error from the class in the MRO.
            raise

    def __setattr__(self, key, value):
        """Take care of mapping indexed attributes to their respective list elements."""

        # Try matching the case where there are two indices
        # this indicates a index_mapped parameter
        attr_name, index, mapkey = self._split_attribute_index_mapping(key)

        # Check if this is an index_mapped attribute. avoiding an infinite
        # recursion by calling getattr() with non-existing keys.
        if (
            (mapkey is not None)
            and (index is not None)
            and attr_name in self._get_indexed_mapped_parameter_attributes()
        ):
            indexed_mapped_attr_value = getattr(self, attr_name)
            try:
                indexed_mapped_attr_value[index][mapkey] = value
                return
            except (IndexError, KeyError) as err:
                if not err.args:
                    err.args = ("",)
                err.args = err.args + (
                    f"'{key}' is out of bound for indexed attribute '{attr_name}'",
                )
                raise

        # Otherwise, try indexed attribute
        # Separate the indexed attribute name from the list index.
        attr_name, index = self._split_attribute_index(key)

        # Check if this is an indexed attribute. avoiding an infinite
        # recursion by calling getattr() with non-existing keys.
        if (index is not None) and (
            attr_name in self._get_indexed_parameter_attributes()
        ):
            indexed_attr_value = getattr(self, attr_name)
            try:
                indexed_attr_value[index] = value
                return
            except IndexError:
                raise IndexError(
                    f"'{key}' is out of bound for indexed attribute '{attr_name}'"
                )

        # Forward the request to the next class in the MRO.
        super().__setattr__(key, value)

    def add_cosmetic_attribute(self, attr_name, attr_value):
        """
        Add a cosmetic attribute to this object.

        This attribute will not have a functional effect on the object
        in the Open Force Field Toolkit, but can be written out during
        output.

        .. warning :: The API for modifying cosmetic attributes is experimental
           and may change in the future (see issue #338).

        Parameters
        ----------
        attr_name : str
            Name of the attribute to define for this object.
        attr_value : str
            The value of the attribute to define for this object.

        """
        setattr(self, "_" + attr_name, attr_value)
        self._cosmetic_attribs.append(attr_name)

    def delete_cosmetic_attribute(self, attr_name):
        """
        Delete a cosmetic attribute from this object.

        .. warning :: The API for modifying cosmetic attributes is experimental
           and may change in the future (see issue #338).

        Parameters
        ----------
        attr_name : str
            Name of the cosmetic attribute to delete.
        """
        # TODO: Can we handle this by overriding __delattr__ instead?
        #  Would we also need to override __del__ as well to cover both deletation methods?
        delattr(self, "_" + attr_name)
        self._cosmetic_attribs.remove(attr_name)

    def attribute_is_cosmetic(self, attr_name):
        """
        Determine whether an attribute of this object is cosmetic.

        .. warning :: The API for modifying cosmetic attributes is experimental
           and may change in the future (see issue #338).

        Parameters
        ----------
        attr_name : str
            The attribute name to check

        Returns
        -------
        is_cosmetic : bool
            Returns True if the attribute is defined and is cosmetic. Returns False otherwise.
        """
        return attr_name in self._cosmetic_attribs

    @staticmethod
    def _split_attribute_index(item):
        """Split the attribute name from the final index.

        For example, the method takes 'k2' and returns the tuple ('k', 1).
        If attribute_name doesn't end with an integer, it returns (item, None).
        """

        # Match any number (\d+) at the end of the string ($).
        match = re.search(r"\d+$", item)
        if match is None:
            return item, None

        index = match.group()  # This is a str.
        attr_name = item[: -len(index)]
        index = int(match.group()) - 1
        return attr_name, index

    @staticmethod
    def _split_attribute_index_mapping(item):
        """Split the attribute name from the final index.

        For example, the method takes 'k1_bondorder2' and returns the tuple ('k_bondorder', 0, 2).
        If attribute_name doesn't end with an integer, it returns (item, None, None).
        """
        # Match items of the form <item><index>_<mapping><key>
        # where <index> and <key> always integers
        match = re.search(r"\d+_[A-z]+\d+$", item)
        if match is None:
            return item, None, None

        # Match any number (\d+) at the end of the string ($).
        i_match = r"\d+$"

        indexed, mapped = item.split("_")

        # process indexed component
        match_indexed = re.search(i_match, indexed)
        index = match_indexed.group()  # This is a str.
        attr_name = indexed[: -len(index)]
        index = int(index) - 1

        # process mapped component
        match_mapping = re.search(i_match, mapped)
        key = match_mapping.group()  # This is a str.
        attr_name = f"{attr_name}_{mapped[:-len(key)]}"
        key = int(key)  # we don't subtract 1 here, because these are keys, not indices

        return attr_name, index, key

    @staticmethod
    def _split_attribute_mapping(item):
        """Split the attribute name from the and its mapping.

        For example, the method takes 'k_foo2' and returns the tuple ('k_foo', 2).
        If attribute_name doesn't end with an integer, it returns (item, None).

        """
        # TODO: Can these three splitting functions be collapsed down into one?
        # Match any number (\d+) at the end of the string ($).
        map_match = r"\d+$"

        match_mapping = re.search(map_match, item)
        if match_mapping is None:
            return item, None

        key = match_mapping.group()
        attr_name = item[: -len(key)]
        key = int(key)

        return attr_name, key

    @classmethod
    def _get_parameter_attributes(cls, filter=None):
        """Return all the attributes of the parameters.

        This is constructed dynamically by introspection gathering all
        the descriptors that are instances of the ParameterAttribute class.
        Parent classes of the parameter types are inspected as well.

        Note that since Python 3.6 the order of the class attribute definition
        is preserved (see PEP 520) so this function will return the attribute
        in their declaration order.

        Parameters
        ----------
        filter : Callable, optional
            An optional function with signature filter(ParameterAttribute) -> bool.
            If specified, only attributes for which this functions returns
            True are returned.

        Returns
        -------
        parameter_attributes : Dict[str, ParameterAttribute]
            A map from the name of the controlled parameter to the
            ParameterAttribute descriptor handling it.

        Examples
        --------
        >>> parameter_attributes = ParameterType._get_parameter_attributes()
        >>> sorted(parameter_attributes.keys())
        ['id', 'parent_id', 'smirks']
        >>> isinstance(parameter_attributes['id'], ParameterAttribute)
        True

        """
        # If no filter is specified, get all the parameters.
        if filter is None:
            filter = lambda x: True  # noqa

        # Go through MRO and retrieve also parents descriptors. The function
        # inspect.getmembers() automatically resolves the MRO, but it also
        # sorts the attribute alphabetically by name. Here we want the order
        # to be the same as the declaration order, which is guaranteed by PEP 520,
        # starting from the parent class.
        parameter_attributes = OrderedDict(
            (name, descriptor)
            for c in reversed(inspect.getmro(cls))
            for name, descriptor in c.__dict__.items()
            if isinstance(descriptor, ParameterAttribute) and filter(descriptor)
        )
        return parameter_attributes

    @classmethod
    def _get_indexed_mapped_parameter_attributes(cls):
        """Shortcut to retrieve only IndexedMappedParameterAttributes."""
        return cls._get_parameter_attributes(
            filter=lambda x: isinstance(x, IndexedMappedParameterAttribute)
        )

    @classmethod
    def _get_indexed_parameter_attributes(cls):
        """Shortcut to retrieve only IndexedParameterAttributes."""
        return cls._get_parameter_attributes(
            filter=lambda x: isinstance(x, IndexedParameterAttribute)
        )

    @classmethod
    def _get_mapped_parameter_attributes(cls):
        """Shortcut to retrieve only IndexedParameterAttributes."""
        return cls._get_parameter_attributes(
            filter=lambda x: isinstance(x, MappedParameterAttribute)
        )

    @classmethod
    def _get_required_parameter_attributes(cls):
        """Shortcut to retrieve only required ParameterAttributes."""
        return cls._get_parameter_attributes(filter=lambda x: x.default is x.UNDEFINED)

    @classmethod
    def _get_optional_parameter_attributes(cls):
        """Shortcut to retrieve only required ParameterAttributes."""
        return cls._get_parameter_attributes(
            filter=lambda x: x.default is not x.UNDEFINED
        )

    def _get_defined_parameter_attributes(self):
        """Returns all the attributes except for the optional attributes that have None default value.

        This returns first the required attributes and then the defined optional
        attribute in their respective declaration order.
        """
        required = self._get_required_parameter_attributes()
        optional = self._get_optional_parameter_attributes()
        # Filter the optional parameters that are set to their default.
        optional = OrderedDict(
            (name, descriptor)
            for name, descriptor in optional.items()
            if not (
                descriptor.default is None and getattr(self, name) == descriptor.default
            )
        )
        required.update(optional)
        return required


# We can't actually make this derive from dict, because it's possible for the user to change SMIRKS
# of parameters already in the list, which would cause the ParameterType object's SMIRKS and
# the dictionary key's SMIRKS to be out of sync.
class ParameterList(list):
    """
    Parameter list that also supports accessing items by SMARTS string.

    .. warning :: This API is experimental and subject to change.

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

        input_parameter_list = input_parameter_list or []
        # TODO: Should a ParameterList only contain a single kind of ParameterType?
        for input_parameter in input_parameter_list:
            self.append(input_parameter)

    def append(self, parameter):
        """
        Add a ParameterType object to the end of the ParameterList

        Parameters
        ----------
        parameter : a ParameterType object

        """
        # TODO: Ensure that newly added parameter is the same type as existing?
        super().append(parameter)

    def extend(self, other):
        """
        Add a ParameterList object to the end of the ParameterList

        Parameters
        ----------
        other : a ParameterList

        """
        if not isinstance(other, ParameterList):
            msg = (
                "ParameterList.extend(other) expected instance of ParameterList, "
                "but received {} (type {}) instead".format(other, type(other))
            )
            raise TypeError(msg)
        # TODO: Check if other ParameterList contains the same ParameterTypes?
        super().extend(other)

    def index(self, item):
        """
        Get the numerical index of a ParameterType object or SMIRKS in this ParameterList.
        Raises ParameterLookupError if the item is not found.

        Parameters
        ----------
        item : ParameterType object or str
            The parameter or SMIRKS to look up in this ParameterList

        Returns
        -------
        index : int
            The index of the found item

        Raises
        ------
        ParameterLookupError if SMIRKS pattern is passed in but not found

        """
        if isinstance(item, ParameterType):
            return super().index(item)
        else:
            for parameter in self:
                if parameter.smirks == item:
                    return self.index(parameter)
            raise ParameterLookupError(
                "SMIRKS {item} not found in ParameterList".format(item=item)
            )

    def insert(self, index, parameter):
        """
        Add a ParameterType object as if this were a list

        Parameters
        ----------
        index : int
            The numerical position to insert the parameter at
        parameter : a ParameterType object
            The parameter to insert
        """
        # TODO: Ensure that newly added parameter is the same type as existing?
        super().insert(index, parameter)

    def __delitem__(self, item):
        """
        Delete item by index or SMIRKS.

        Parameters
        ----------
        item : str or int
            SMIRKS or numerical index of item in this ParameterList
        """
        if type(item) is int:
            index = item
        else:
            # Try to find by SMIRKS
            index = self.index(item)
        super().__delitem__(index)

    def __getitem__(self, item):
        """
        Retrieve item by index or SMIRKS

        Parameters
        ----------
        item : str or int
            SMIRKS or numerical index of item in this ParameterList
        """
        if type(item) is int:
            index = item
        elif type(item) is slice:
            index = item
        elif isinstance(item, str):
            index = self.index(item)
        elif isinstance(item, ParameterType) or issubclass(item, ParameterType):
            raise ParameterLookupError("Lookup by instance is not supported")
        return super().__getitem__(index)

    # TODO: Override __setitem__ and __del__ to ensure we can slice by SMIRKS as well
    # This is needed for pickling. See https://github.com/openforcefield/openff-toolkit/issues/411
    # for more details.
    # TODO: Is there a cleaner way (getstate/setstate perhaps?) to allow FFs to be
    #       pickled?
    def __reduce__(self):
        return (__class__, (list(self),), self.__dict__)

    def __contains__(self, item):
        """Check to see if either Parameter or SMIRKS is contained in parameter list.

        Parameters
        ----------
        item : str
            SMIRKS of item in this ParameterList
        """
        if isinstance(item, str):
            # Special case for SMIRKS strings
            if item in [result.smirks for result in self]:
                return True
        # Fall back to traditional access
        return list.__contains__(self, item)

    def to_list(self, discard_cosmetic_attributes=True):
        """
        Render this ParameterList to a normal list, serializing each ParameterType object in it to dict.

        Parameters
        ----------

        discard_cosmetic_attributes : bool, optional. Default = True
            Whether to discard non-spec attributes of each ParameterType object.

        Returns
        -------
        parameter_list : List[dict]
            A serialized representation of a ParameterList, with each ParameterType it contains converted to dict.
        """
        parameter_list = list()

        for parameter in self:
            parameter_dict = parameter.to_dict(
                discard_cosmetic_attributes=discard_cosmetic_attributes
            )
            parameter_list.append(parameter_dict)

        return parameter_list


# TODO: Rename to better reflect role as parameter base class?
class ParameterType(_ParameterAttributeHandler):
    """
    Base class for SMIRNOFF parameter types.

    This base class provides utilities to create new parameter types. See
    the below for examples of how to do this.

    .. warning :: This API is experimental and subject to change.

    Attributes
    ----------
    smirks : str
        The SMIRKS pattern that this parameter matches.
    id : str or None
        An optional identifier for the parameter.
    parent_id : str or None
        Optionally, the identifier of the parameter of which this parameter
        is a specialization.

    See Also
    --------
    ParameterAttribute
    IndexedParameterAttribute

    Examples
    --------

    This class allows to define new parameter types by just listing its
    attributes. In the example below, ``_VALENCE_TYPE`` AND ``_ELEMENT_NAME``
    are used for the validation of the SMIRKS pattern associated to the
    parameter and the automatic serialization/deserialization into a ``dict``.

    >>> class MyBondParameter(ParameterType):
    ...     _VALENCE_TYPE = 'Bond'
    ...     _ELEMENT_NAME = 'Bond'
    ...     length = ParameterAttribute(unit=unit.angstrom)
    ...     k = ParameterAttribute(unit=unit.kilocalorie_per_mole / unit.angstrom**2)
    ...

    The parameter automatically inherits the required smirks attribute
    from ``ParameterType``. Associating a ``unit`` to a ``ParameterAttribute``
    cause the attribute to accept only values in compatible units and to
    parse string expressions.

    >>> my_par = MyBondParameter(
    ...     smirks='[*:1]-[*:2]',
    ...     length='1.01 * angstrom',
    ...     k=5 * unit.kilocalorie_per_mole / unit.angstrom**2
    ... )
    >>> my_par.length
    Quantity(value=1.01, unit=angstrom)
    >>> my_par.k = 3.0 * unit.gram
    Traceback (most recent call last):
    ...
    openff.toolkit.utils.utils.IncompatibleUnitError: k=3.0 g should have units of kilocalorie/(angstrom**2*mole)

    Each attribute can be made optional by specifying a default value,
    and you can attach a converter function by passing a callable as an
    argument or through the decorator syntax.

    >>> class MyParameterType(ParameterType):
    ...     _VALENCE_TYPE = 'Atom'
    ...     _ELEMENT_NAME = 'Atom'
    ...
    ...     attr_optional = ParameterAttribute(default=2)
    ...     attr_all_to_float = ParameterAttribute(converter=float)
    ...     attr_int_to_float = ParameterAttribute()
    ...
    ...     @attr_int_to_float.converter
    ...     def attr_int_to_float(self, attr, value):
    ...         # This converter converts only integers to floats
    ...         # and raise an exception for the other types.
    ...         if isinstance(value, int):
    ...             return float(value)
    ...         elif not isinstance(value, float):
    ...             raise TypeError(f"Cannot convert '{value}' to float")
    ...         return value
    ...
    >>> my_par = MyParameterType(smirks='[*:1]', attr_all_to_float='3.0', attr_int_to_float=1)
    >>> my_par.attr_optional
    2
    >>> my_par.attr_all_to_float
    3.0
    >>> my_par.attr_int_to_float
    1.0

    The float() function can convert strings to integers, but our custom
    converter forbids it

    >>> my_par.attr_all_to_float = '2.0'
    >>> my_par.attr_int_to_float = '4.0'
    Traceback (most recent call last):
    ...
    TypeError: Cannot convert '4.0' to float

    Parameter attributes that can be indexed can be handled with the
    ``IndexedParameterAttribute``. These support unit validation and
    converters exactly as ``ParameterAttribute``\ s, but the validation/conversion
    is performed for each indexed attribute.

    >>> class MyTorsionType(ParameterType):
    ...     _VALENCE_TYPE = 'ProperTorsion'
    ...     _ELEMENT_NAME = 'Proper'
    ...     periodicity = IndexedParameterAttribute(converter=int)
    ...     k = IndexedParameterAttribute(unit=unit.kilocalorie_per_mole)
    ...
    >>> my_par = MyTorsionType(
    ...     smirks='[*:1]-[*:2]-[*:3]-[*:4]',
    ...     periodicity1=2,
    ...     k1=5 * unit.kilocalorie_per_mole,
    ...     periodicity2='3',
    ...     k2=6 * unit.kilocalorie_per_mole,
    ... )
    >>> my_par.periodicity
    [2, 3]

    Indexed attributes, can be accessed both as a list or as their indexed
    parameter name.

    >>> my_par.periodicity2 = 6
    >>> my_par.periodicity[0] = 1
    >>> my_par.periodicity
    [1, 6]

    """

    # ChemicalEnvironment valence type string expected by SMARTS string for this Handler
    _VALENCE_TYPE = None
    # The string mapping to this ParameterType in a SMIRNOFF data source
    _ELEMENT_NAME = None

    # Parameter attributes shared among all parameter types.
    smirks = ParameterAttribute()
    id = ParameterAttribute(default=None)
    parent_id = ParameterAttribute(default=None)

    @smirks.converter
    def smirks(self, attr, smirks):
        # Validate the SMIRKS string to ensure it matches the expected
        # parameter type, raising an exception if it is invalid or doesn't
        # tag a valid set of atoms.

        # TODO: Add check to make sure we can't make tree non-hierarchical
        #       This would require parameter type knows which ParameterList it belongs to
        ChemicalEnvironment.validate_smirks(smirks, validate_valence_type=True)
        return smirks

    def __init__(self, smirks, allow_cosmetic_attributes=False, **kwargs):
        """
        Create a ParameterType.

        Parameters
        ----------
        smirks : str
            The SMIRKS match for the provided parameter type.
        allow_cosmetic_attributes : bool optional. Default = False
            Whether to permit non-spec kwargs ("cosmetic attributes"). If True, non-spec kwargs will be stored as
            an attribute of this parameter which can be accessed and written out. Otherwise an exception will
            be raised.

        """
        # This is just to make smirks a required positional argument.
        kwargs["smirks"] = smirks
        super().__init__(allow_cosmetic_attributes=allow_cosmetic_attributes, **kwargs)

    def __repr__(self):
        ret_str = "<{} with ".format(self.__class__.__name__)
        for attr, val in self.to_dict().items():
            ret_str += f"{attr}: {val}  "
        ret_str += ">"
        return ret_str


# ======================================================================
# PARAMETER HANDLERS
#
# The following classes are Handlers that know how to create Force
# subclasses and add them to an OpenMM System that is being created. Each
# Handler class must define three methods:
# 1) a constructor which takes as input hierarchical dictionaries of data
#    conformant to the SMIRNOFF spec;
# 2) a create_force() method that constructs the Force object and adds it
#    to the System; and
# 3) a labelForce() method that provides access to which terms are applied
#    to which atoms in specified mols.
# ======================================================================

# TODO: Should we have a parameter handler registry?


class ParameterHandler(_ParameterAttributeHandler):
    """Base class for parameter handlers.

    Parameter handlers are configured with some global parameters for a
    given section. They may also contain a :class:`ParameterList` populated
    with :class:`ParameterType` objects if they are responsible for assigning
    SMIRKS-based parameters.

    .. warning

       Parameter handler objects can only belong to a single :class:`ForceField` object.
       If you need to create a copy to attach to a different :class:`ForceField` object,
       use ``create_copy()``.

    .. warning :: This API is experimental and subject to change.

    """

    _TAGNAME = None  # str of section type handled by this ParameterHandler (XML element name for SMIRNOFF XML representation)
    _INFOTYPE = None  # container class with type information that will be stored in self._parameters
    _OPENMMTYPE = None  # OpenMM Force class (or None if no equivalent)
    _DEPENDENCIES = (
        None  # list of ParameterHandler classes that must precede this, or None
    )

    _KWARGS = []  # Kwargs to catch when create_force is called
    _SMIRNOFF_VERSION_INTRODUCED = (
        0.0  # the earliest version of SMIRNOFF spec that supports this ParameterHandler
    )
    _SMIRNOFF_VERSION_DEPRECATED = (
        None  # if deprecated, the first SMIRNOFF version number it is no longer used
    )
    _MIN_SUPPORTED_SECTION_VERSION = 0.3
    _MAX_SUPPORTED_SECTION_VERSION = 0.3

    version = ParameterAttribute()

    @version.converter
    def version(self, attr, new_version):
        """
        Raise a parsing exception if the given section version is unsupported.

        Raises
        ------
        SMIRNOFFVersionError if an incompatible version is passed in.

        """
        import packaging.version

        from openff.toolkit.typing.engines.smirnoff import SMIRNOFFVersionError

        # Use PEP-440 compliant version number comparison, if requested
        if (
            packaging.version.parse(str(new_version))
            > packaging.version.parse(str(self._MAX_SUPPORTED_SECTION_VERSION))
        ) or (
            packaging.version.parse(str(new_version))
            < packaging.version.parse(str(self._MIN_SUPPORTED_SECTION_VERSION))
        ):
            raise SMIRNOFFVersionError(
                f"SMIRNOFF offxml file was written with version {new_version}, but this version "
                f"of ForceField only supports version {self._MIN_SUPPORTED_SECTION_VERSION} "
                f"to version {self._MAX_SUPPORTED_SECTION_VERSION}"
            )
        return new_version

    def __init__(
        self, allow_cosmetic_attributes=False, skip_version_check=False, **kwargs
    ):
        """
        Initialize a ParameterHandler, optionally with a list of parameters and other kwargs.

        Parameters
        ----------
        allow_cosmetic_attributes : bool, optional. Default = False
            Whether to permit non-spec kwargs. If True, non-spec kwargs will be stored as attributes of this object
            and can be accessed and modified. Otherwise an exception will be raised if a non-spec kwarg is encountered.
        skip_version_check: bool, optional. Default = False
            If False, the SMIRNOFF section version will not be checked, and the ParameterHandler will be initialized
            with version set to _MAX_SUPPORTED_SECTION_VERSION.
        **kwargs : dict
            The dict representation of the SMIRNOFF data source

        """
        # Skip version check if requested.
        if "version" not in kwargs:
            if skip_version_check:
                kwargs["version"] = self._MAX_SUPPORTED_SECTION_VERSION
            else:
                raise SMIRNOFFSpecError(
                    f"Missing version while trying to construct {self.__class__}. "
                    f"0.3 SMIRNOFF spec requires each parameter section to have its own version."
                )

        # List of ParameterType objects (also behaves like an OrderedDict where keys are SMARTS).
        self._parameters = ParameterList()

        # Initialize ParameterAttributes and cosmetic attributes.
        super().__init__(allow_cosmetic_attributes=allow_cosmetic_attributes, **kwargs)

    def _add_parameters(self, section_dict, allow_cosmetic_attributes=False):
        """
        Extend the ParameterList in this ParameterHandler using a SMIRNOFF data source.

        Parameters
        ----------
        section_dict : dict
            The dict representation of a SMIRNOFF data source containing parameters to att to this ParameterHandler
        allow_cosmetic_attributes : bool, optional. Default = False
            Whether to allow non-spec fields in section_dict. If True, non-spec kwargs will be stored as an
            attribute of the parameter. If False, non-spec kwargs will raise an exception.

        """
        unitless_kwargs, attached_units = extract_serialized_units_from_dict(
            section_dict
        )
        smirnoff_data = attach_units(unitless_kwargs, attached_units)

        for key, val in smirnoff_data.items():
            if self._INFOTYPE is not None:
                element_name = self._INFOTYPE._ELEMENT_NAME
                # Skip sections that aren't the parameter list
                if key != element_name:
                    break
            # If there are multiple parameters, this will be a list. If there's just one, make it a list
            if not (isinstance(val, list)):
                val = [val]

            # If we're reading the parameter list, iterate through and attach units to
            # each parameter_dict, then use it to initialize a ParameterType
            for unitless_param_dict in val:

                param_dict = attach_units(unitless_param_dict, attached_units)
                new_parameter = self._INFOTYPE(
                    **param_dict, allow_cosmetic_attributes=allow_cosmetic_attributes
                )
                self._parameters.append(new_parameter)

    @property
    def parameters(self):
        """The ParameterList that holds this ParameterHandler's parameter objects"""
        return self._parameters

    @property
    def TAGNAME(self):
        """
        The name of this ParameterHandler corresponding to the SMIRNOFF tag name

        Returns
        -------
        handler_name : str
            The name of this parameter handler

        """
        return self._TAGNAME

    # TODO: Do we need to return these, or can we handle this internally
    @property
    def known_kwargs(self):
        """List of kwargs that can be parsed by the function."""
        # TODO: Should we use introspection to inspect the function signature instead?
        return set(self._KWARGS)

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
    def add_parameter(
        self, parameter_kwargs=None, parameter=None, after=None, before=None
    ):
        """Add a parameter to the force field, ensuring all parameters are valid.

        Parameters
        ----------
        parameter_kwargs: dict, optional
            The kwargs to pass to the ParameterHandler.INFOTYPE (a ParameterType) constructor
        parameter: ParameterType, optional
            A ParameterType to add to the ParameterHandler
        after : str or int, optional
            The SMIRKS pattern (if str) or index (if int) of the parameter directly before where
            the new parameter will be added
        before : str, optional
            The SMIRKS pattern (if str) or index (if int) of the parameter directly after where
            the new parameter will be added

        Note the following behavior:
          * Either `parameter_kwargs` or `parameter` must be specified.
          * When `before` and `after` are both `None`, the new parameter will be appended
            to the **END** of the parameter list.
          * When `before` and `after` are both specified, the new parameter will be added immediately
            after the parameter matching the `after` pattern or index.
          * The order of parameters in a parameter list can have significant impacts on parameter assignment. For details,
            see the [SMIRNOFF](https://open-forcefield-toolkit.readthedocs.io/en/latest/smirnoff.html#smirnoff-parameter-specification-is-hierarchical)
            specification.

        Examples
        --------

        Add a ParameterType to an existing ParameterList at a specified position.

        Given an existing parameter handler and a new parameter to add to it:

        >>> from openmm import unit
        >>> bh = BondHandler(skip_version_check=True)
        >>> length = 1.5 * unit.angstrom
        >>> k = 100 * unit.kilocalorie_per_mole / unit.angstrom ** 2
        >>> bh.add_parameter({'smirks': '[*:1]-[*:2]', 'length': length, 'k': k, 'id': 'b1'})
        >>> bh.add_parameter({'smirks': '[*:1]=[*:2]', 'length': length, 'k': k, 'id': 'b2'})
        >>> bh.add_parameter({'smirks': '[*:1]#[*:2]', 'length': length, 'k': k, 'id': 'b3'})
        >>> [p.id for p in bh.parameters]
        ['b1', 'b2', 'b3']

        >>> param = {'smirks': '[#1:1]-[#6:2]', 'length': length, 'k': k, 'id': 'b4'}

        Add a new parameter immediately after the parameter with the smirks '[*:1]=[*:2]'

        >>> bh.add_parameter(param, after='[*:1]=[*:2]')
        >>> [p.id for p in bh.parameters]
        ['b1', 'b2', 'b4', 'b3']
        """
        for val in [before, after]:
            if val and not isinstance(val, (str, int)):
                raise TypeError

        # If a dict was passed, construct it; if a ParameterType was passed, do nothing
        if parameter_kwargs:
            new_parameter = self._INFOTYPE(**parameter_kwargs)
        elif parameter:
            new_parameter = parameter
        else:
            raise ValueError("One of (parameter, parameter_kwargs) must be specified")

        if new_parameter.smirks in [p.smirks for p in self._parameters]:
            msg = f"A parameter SMIRKS pattern {new_parameter.smirks} already exists."
            raise DuplicateParameterError(msg)

        if before is not None:
            if isinstance(before, str):
                before_index = self._parameters.index(before)
            elif isinstance(before, int):
                before_index = before

        if after is not None:
            if isinstance(after, str):
                after_index = self._parameters.index(after)
            elif isinstance(after, int):
                after_index = after

        if None not in (before, after):
            if after_index > before_index:
                raise ValueError("before arg must be before after arg")

        if after is not None:
            self._parameters.insert(after_index + 1, new_parameter)
        elif before is not None:
            self._parameters.insert(before_index, new_parameter)
        else:
            self._parameters.append(new_parameter)

    def get_parameter(self, parameter_attrs):
        """
        Return the parameters in this ParameterHandler that match the parameter_attrs argument.
        When multiple attrs are passed, parameters that have any (not all) matching attributes
        are returned.

        Parameters
        ----------
        parameter_attrs : dict of {attr: value}
            The attrs mapped to desired values (for example {"smirks": "[*:1]~[#16:2]=,:[#6:3]~[*:4]", "id": "t105"} )

        Returns
        -------
        params : list of ParameterType objects
            A list of matching ParameterType objects

        Examples
        --------

        Create a parameter handler and populate it with some data.

        >>> from openmm import unit
        >>> handler = BondHandler(skip_version_check=True)
        >>> handler.add_parameter(
        ...     {
        ...         'smirks': '[*:1]-[*:2]',
        ...         'length': 1*unit.angstrom,
        ...         'k': 10*unit.kilocalorie_per_mole/unit.angstrom**2,
        ...     }
        ... )

        Look up, from this handler, all parameters matching some SMIRKS pattern

        >>> handler.get_parameter({'smirks': '[*:1]-[*:2]'})
        [<BondType with smirks: [*:1]-[*:2]  length: 1 A  k: 10 kcal/(A**2 mol)  >]

        """
        params = list()
        for attr, value in parameter_attrs.items():
            for param in self.parameters:
                if param in params:
                    continue
                # TODO: Cleaner accessing of cosmetic attributes
                # See issue #338
                if param.attribute_is_cosmetic(attr):
                    attr = "_" + attr
                if hasattr(param, attr):
                    if getattr(param, attr) == value:
                        params.append(param)
        return params

    class _Match:
        """Represents a ParameterType which has been matched to
        a given chemical environment.
        """

        @property
        def parameter_type(self):
            """ParameterType: The matched parameter type."""
            return self._parameter_type

        @property
        def environment_match(self):
            """Topology._ChemicalEnvironmentMatch: The environment which matched the type."""
            return self._environment_match

        def __init__(self, parameter_type, environment_match):
            """Constructs a new ParameterHandlerMatch object.

            Parameters
            ----------
            parameter_type: ParameterType
                The matched parameter type.
            environment_match: Topology._ChemicalEnvironmentMatch
                The environment which matched the type.
            """
            self._parameter_type = parameter_type
            self._environment_match = environment_match

    def find_matches(self, entity, unique=False):
        """Find the elements of the topology/molecule matched by a parameter type.

        Parameters
        ----------
        entity : openff.toolkit.topology.Topology
            Topology to search.
        unique : bool, default=False
            If False, SMARTS matching will enumerate every valid permutation of matching atoms.
            If True, only one order of each unique match will be returned.

        Returns
        ---------
        matches : ValenceDict[Tuple[int], ParameterHandler._Match]
            ``matches[particle_indices]`` is the ``ParameterType`` object
            matching the tuple of particle indices in ``entity``.
        """

        # TODO: Right now, this method is only ever called with an entity that is a Topology.
        #  Should we reduce its scope and have a check here to make sure entity is a Topology?
        return self._find_matches(entity, unique=unique)

    def _find_matches(
        self,
        entity,
        transformed_dict_cls=ValenceDict,
        unique=False,
    ):
        """Implement find_matches() and allow using a difference valence dictionary.
        Parameters
        ----------
        entity : openff.toolkit.topology.Topology
            Topology to search.
        transformed_dict_cls: class
            The type of dictionary to store the matches in. This
            will determine how groups of atom indices are stored
            and accessed (e.g for angles indices should be 0-1-2
            and not 2-1-0).
        unique : bool, default=False
            If False, SMARTS matching will enumerate every valid permutation of matching atoms.
            If True, only one order of each unique match will be returned.

        Returns
        ---------
        matches : `transformed_dict_cls` of ParameterHandlerMatch
            ``matches[particle_indices]`` is the ``ParameterType`` object
            matching the tuple of particle indices in ``entity``.
        """
        logger.debug("Finding matches for {}".format(self.__class__.__name__))

        matches = transformed_dict_cls()

        # TODO: There are probably performance gains to be had here
        #       by performing this loop in reverse order, and breaking early once
        #       all environments have been matched.
        for parameter_type in self._parameters:
            matches_for_this_type = {}

            for environment_match in entity.chemical_environment_matches(
                parameter_type.smirks,
                unique=unique,
            ):
                # Update the matches for this parameter type.
                handler_match = self._Match(parameter_type, environment_match)
                matches_for_this_type[
                    environment_match.topology_atom_indices
                ] = handler_match

            # Update matches of all parameter types.
            matches.update(matches_for_this_type)

            logger.debug(
                "{:64} : {:8} matches".format(
                    parameter_type.smirks, len(matches_for_this_type)
                )
            )

        logger.debug("{} matches identified".format(len(matches)))
        return matches

    @staticmethod
    def _assert_correct_connectivity(match, expected_connectivity=None):
        """A more performant version of the `topology.assert_bonded` method
        to ensure that the results of `_find_matches` are valid.

        Raises
        ------
        ValueError
            Raise an exception when the atoms in the match don't have
            the correct connectivity.

        Parameters
        ----------
        match: ParameterHandler._Match
            The match found by `_find_matches`
        connectivity: list of tuple of int, optional
            The expected connectivity of the match (e.g. for a torsion
            expected_connectivity=[(0, 1), (1, 2), (2, 3)]). If `None`,
            a connectivity of [(0, 1), ... (n - 1, n)] is assumed.
        """

        # I'm not 100% sure this is really necessary... but this should do
        # the same checks as the more costly assert_bonded method in the
        # ParameterHandler.create_force methods.
        if expected_connectivity is None:
            return

        reference_molecule = match.environment_match.reference_molecule

        for connectivity in expected_connectivity:

            atom_i = match.environment_match.reference_atom_indices[connectivity[0]]
            atom_j = match.environment_match.reference_atom_indices[connectivity[1]]

            reference_molecule.get_bond_between(atom_i, atom_j)

    def assign_parameters(self, topology, system):
        """Assign parameters for the given Topology to the specified OpenMM ``System`` object.

        Parameters
        ----------
        topology : openff.toolkit.topology.Topology
            The Topology for which parameters are to be assigned.
            Either a new Force will be created or parameters will be appended to an existing Force.
        system : openmm.System
            The OpenMM System object to add the Force (or append new parameters) to.
        """
        pass

    def postprocess_system(self, topology, system, **kwargs):
        """Allow the force to perform a a final post-processing pass on the OpenMM ``System`` following parameter assignment, if needed.

        Parameters
        ----------
        topology : openff.toolkit.topology.Topology
            The Topology for which parameters are to be assigned.
            Either a new Force will be created or parameters will be appended to an existing Force.
        system : openmm.System
            The OpenMM System object to add the Force (or append new parameters) to.
        """
        pass

    def to_dict(self, discard_cosmetic_attributes=False):
        """
        Convert this ParameterHandler to an OrderedDict, compliant with the SMIRNOFF data spec.

        Parameters
        ----------
        discard_cosmetic_attributes : bool, optional. Default = False.
            Whether to discard non-spec parameter and header attributes in this ParameterHandler.

        Returns
        -------
        smirnoff_data : OrderedDict
            SMIRNOFF-spec compliant representation of this ParameterHandler and its internal ParameterList.

        """
        smirnoff_data = OrderedDict()

        # Populate parameter list
        parameter_list = self._parameters.to_list(
            discard_cosmetic_attributes=discard_cosmetic_attributes
        )

        # NOTE: This assumes that a ParameterHandler will have just one homogenous ParameterList under it
        if self._INFOTYPE is not None:
            # smirnoff_data[self._INFOTYPE._ELEMENT_NAME] = unitless_parameter_list
            smirnoff_data[self._INFOTYPE._ELEMENT_NAME] = parameter_list

        # Collect parameter and cosmetic attributes.
        header_attribute_dict = super().to_dict(
            discard_cosmetic_attributes=discard_cosmetic_attributes
        )
        smirnoff_data.update(header_attribute_dict)

        return smirnoff_data

    # -------------------------------
    # Utilities for children classes.
    # -------------------------------

    @classmethod
    def _check_all_valence_terms_assigned(
        cls,
        assigned_terms,
        valence_terms,
        exception_cls=UnassignedValenceParameterException,
    ):
        """Check that all valence terms have been assigned and print a user-friendly error message.

        Parameters
        ----------
        assigned_terms : ValenceDict
            Atom index tuples defining added valence terms.
        valence_terms : Iterable[TopologyAtom] or Iterable[Iterable[TopologyAtom]]
            Atom or atom tuples defining topological valence terms.
        exception_cls : UnassignedValenceParameterException
            A specific exception class to raise to allow catching only specific
            types of errors.

        """
        from openff.toolkit.topology import TopologyAtom

        # Provided there are no duplicates in either list,
        # or something weird like a bond has been added to
        # a torsions list - this should work just fine I think.
        # If we expect either of those assumptions to be incorrect,
        # (i.e len(not_found_terms) > 0) we have bigger issues
        # in the code and should be catching those cases elsewhere!
        # The fact that we graph match all topol molecules to ref
        # molecules should avoid the len(not_found_terms) > 0 case.

        if len(assigned_terms) == len(valence_terms):
            return

        # Convert the valence term to a valence dictionary to make sure
        # the order of atom indices doesn't matter for comparison.
        valence_terms_dict = assigned_terms.__class__()
        for atoms in valence_terms:
            try:
                # valence_terms is a list of TopologyAtom tuples.
                atom_indices = (a.topology_particle_index for a in atoms)
            except TypeError:
                # valence_terms is a list of TopologyAtom.
                atom_indices = (atoms.topology_particle_index,)
            valence_terms_dict[atom_indices] = atoms

        # Check that both valence dictionaries have the same keys (i.e. terms).
        assigned_terms_set = set(assigned_terms.keys())
        valence_terms_set = set(valence_terms_dict.keys())
        unassigned_terms = valence_terms_set.difference(assigned_terms_set)
        not_found_terms = assigned_terms_set.difference(valence_terms_set)

        # Raise an error if there are unassigned terms.
        err_msg = ""

        if len(unassigned_terms) > 0:

            unassigned_topology_atom_tuples = []

            # Gain access to the relevant topology
            if type(valence_terms[0]) is TopologyAtom:
                topology = valence_terms[0].topology_molecule.topology
            else:
                topology = valence_terms[0][0].topology_molecule.topology
            unassigned_str = ""
            for unassigned_tuple in unassigned_terms:
                unassigned_str += "\n- Topology indices " + str(unassigned_tuple)
                unassigned_str += ": names and elements "

                unassigned_topology_atoms = []

                # Pull and add additional helpful info on missing terms
                for atom_idx in unassigned_tuple:
                    topology_atom = topology.atom(atom_idx)
                    unassigned_topology_atoms.append(topology_atom)
                    unassigned_str += f"({topology_atom.atom.name} {topology_atom.atom.element.symbol}), "
                unassigned_topology_atom_tuples.append(tuple(unassigned_topology_atoms))
            err_msg += (
                "{parameter_handler} was not able to find parameters for the following valence terms:\n"
                "{unassigned_str}"
            ).format(parameter_handler=cls.__name__, unassigned_str=unassigned_str)
        if len(not_found_terms) > 0:
            if err_msg != "":
                err_msg += "\n"
            not_found_str = "\n- ".join([str(x) for x in not_found_terms])
            err_msg += (
                "{parameter_handler} assigned terms that were not found in the topology:\n"
                "- {not_found_str}"
            ).format(parameter_handler=cls.__name__, not_found_str=not_found_str)
        if err_msg != "":
            err_msg += "\n"
            exception = exception_cls(err_msg)
            exception.unassigned_topology_atom_tuples = unassigned_topology_atom_tuples
            exception.handler_class = cls
            raise exception

    def _check_attributes_are_equal(
        self, other, identical_attrs=(), tolerance_attrs=(), tolerance=1e-6
    ):
        """Utility function to check that the given attributes of the two handlers are equal.

        Parameters
        ----------
        identical_attrs : List[str]
            Names of the parameters that must be checked with the equality operator.
        tolerance_attrs : List[str]
            Names of the parameters that must be equal up to a tolerance.
        tolerance : float
            The absolute tolerance used to compare the parameters.
        """

        def get_unitless_values(attr):
            this_val = getattr(self, attr)
            other_val = getattr(other, attr)
            # Strip quantities of their units before comparison.
            try:
                u = this_val.unit
            except AttributeError:
                return this_val, other_val
            return this_val / u, other_val / u

        for attr in identical_attrs:
            this_val, other_val = get_unitless_values(attr)

            if this_val != other_val:
                raise IncompatibleParameterError(
                    "{} values are not identical. "
                    "(handler value: {}, incompatible value: {}".format(
                        attr, this_val, other_val
                    )
                )

        for attr in tolerance_attrs:
            this_val, other_val = get_unitless_values(attr)
            if abs(this_val - other_val) > tolerance:
                raise IncompatibleParameterError(
                    "Difference between '{}' values is beyond allowed tolerance {}. "
                    "(handler value: {}, incompatible value: {}".format(
                        attr, tolerance, this_val, other_val
                    )
                )

    @staticmethod
    def check_partial_bond_orders_from_molecules_duplicates(pb_mols):
        if len(set(map(Molecule.to_smiles, pb_mols))) < len(pb_mols):
            raise ValueError(
                "At least two user-provided fractional bond order "
                "molecules are isomorphic"
            )

    @staticmethod
    def assign_partial_bond_orders_from_molecules(topology, pbo_mols):
        # for each reference molecule in our topology, we'll walk through the provided partial bond order molecules
        # if we find a match, we'll apply the partial bond orders and skip to the next molecule
        for ref_mol in topology.reference_molecules:
            for pbo_mol in pbo_mols:
                # we are as stringent as we are in the ElectrostaticsHandler
                # TODO: figure out whether bond order matching is redundant with aromatic matching
                isomorphic, topology_atom_map = Molecule.are_isomorphic(
                    ref_mol,
                    pbo_mol,
                    return_atom_map=True,
                    aromatic_matching=True,
                    formal_charge_matching=True,
                    bond_order_matching=True,
                    atom_stereochemistry_matching=True,
                    bond_stereochemistry_matching=True,
                )

                # if matching, assign bond orders and skip to next molecule
                # first match wins
                if isomorphic:
                    # walk through bonds on reference molecule
                    for bond in ref_mol.bonds:
                        # use atom mapping to translate to pbo_molecule bond
                        pbo_bond = pbo_mol.get_bond_between(
                            topology_atom_map[bond.atom1_index],
                            topology_atom_map[bond.atom2_index],
                        )
                        # extract fractional bond order
                        # assign fractional bond order to reference molecule bond
                        if pbo_bond.fractional_bond_order is None:
                            raise ValueError(
                                f"Molecule '{ref_mol}' was requested to be parameterized "
                                f"with user-provided fractional bond orders from '{pbo_mol}', but not "
                                "all bonds were provided with `fractional_bond_order` specified"
                            )

                        bond.fractional_bond_order = pbo_bond.fractional_bond_order

                    break
                # not necessary, but explicit
                else:
                    continue

    def __getitem__(self, val):
        """
        Syntax sugar for lookikng up a ParameterType in a ParameterHandler
        based on its SMIRKS.
        """
        return self.parameters[val]


class ConstraintHandler(ParameterHandler):
    """Handle SMIRNOFF ``<Constraints>`` tags

    ``ConstraintHandler`` must be applied before ``BondHandler`` and ``AngleHandler``,
    since those classes add constraints for which equilibrium geometries are needed from those tags.

    .. warning :: This API is experimental and subject to change.
    """

    class ConstraintType(ParameterType):
        """A SMIRNOFF constraint type

        .. warning :: This API is experimental and subject to change.
        """

        _VALENCE_TYPE = "Bond"
        _ELEMENT_NAME = "Constraint"

        distance = ParameterAttribute(default=None, unit=unit.angstrom)

    _TAGNAME = "Constraints"
    _INFOTYPE = ConstraintType
    _OPENMMTYPE = None  # don't create a corresponding OpenMM Force class

    def create_force(self, system, topology, **kwargs):
        constraint_matches = self.find_matches(topology)
        for (atoms, constraint_match) in constraint_matches.items():
            # Update constrained atom pairs in topology
            # topology.add_constraint(*atoms, constraint.distance)
            # If a distance is specified (constraint.distance != True), add the constraint here.
            # Otherwise, the equilibrium bond length will be used to constrain the atoms in HarmonicBondHandler
            constraint = constraint_match.parameter_type

            if constraint.distance is None:
                topology.add_constraint(*atoms, True)
            else:
                system.addConstraint(*atoms, constraint.distance)
                topology.add_constraint(*atoms, constraint.distance)


class BondHandler(ParameterHandler):
    """Handle SMIRNOFF ``<Bonds>`` tags

    .. warning :: This API is experimental and subject to change.
    """

    class BondType(ParameterType):
        """A SMIRNOFF bond type

        .. warning :: This API is experimental and subject to change.
        """

        # ChemicalEnvironment valence type string expected by SMARTS string for this Handler
        _VALENCE_TYPE = "Bond"
        _ELEMENT_NAME = "Bond"

        length = ParameterAttribute(default=None, unit=unit.angstrom)
        k = ParameterAttribute(
            default=None, unit=unit.kilocalorie_per_mole / unit.angstrom ** 2
        )

        # fractional bond order params
        length_bondorder = MappedParameterAttribute(default=None, unit=unit.angstrom)
        k_bondorder = MappedParameterAttribute(
            default=None, unit=unit.kilocalorie_per_mole / unit.angstrom ** 2
        )

        def __init__(self, **kwargs):
            # these checks enforce mutually-exclusive parameterattribute specifications
            has_k = "k" in kwargs.keys()
            has_k_bondorder = any(["k_bondorder" in key for key in kwargs.keys()])
            has_length = "length" in kwargs.keys()
            has_length_bondorder = any(
                ["length_bondorder" in key for key in kwargs.keys()]
            )

            # Are these errors too general? What about ParametersMissingError/ParametersOverspecifiedError?
            if has_k:
                if has_k_bondorder:
                    raise SMIRNOFFSpecError(
                        "BOTH k and k_bondorder* cannot be specified simultaneously."
                    )
            else:
                if not has_k_bondorder:
                    raise SMIRNOFFSpecError(
                        "Either k or k_bondorder* must be specified."
                    )
            if has_length:
                if has_length_bondorder:
                    raise SMIRNOFFSpecError(
                        "BOTH length and length_bondorder* cannot be specified simultaneously."
                    )
            else:
                if not has_length_bondorder:
                    raise SMIRNOFFSpecError(
                        "Either length or length_bondorder* must be specified."
                    )

            super().__init__(**kwargs)

    _TAGNAME = "Bonds"  # SMIRNOFF tag name to process
    _INFOTYPE = BondType  # class to hold force type info
    _OPENMMTYPE = openmm.HarmonicBondForce  # OpenMM force class to create
    _DEPENDENCIES = [ConstraintHandler]  # ConstraintHandler must be executed first
    _MAX_SUPPORTED_SECTION_VERSION = 0.4

    # Use the _allow_only filter here because this class's implementation contains all the information about supported
    # potentials for this handler.
    potential = ParameterAttribute(
        default="overridden in init",
        converter=_allow_only(["harmonic", "(k/2)*(r-length)^2"]),
    )
    # The default value for fractional_bondorder_method depends on the section version and is overwritten in __init__.
    # Do not use the allow_only filter here since ToolkitWrappers may be imported that support additional fractional
    # bondorder methods.
    fractional_bondorder_method = ParameterAttribute(default="overridden in init")
    # Use the _allow_only filter here because this class's implementation contains all the information about supported
    # interpolation types.
    fractional_bondorder_interpolation = ParameterAttribute(
        default="linear", converter=_allow_only(["linear"])
    )

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # Default value for fractional_bondorder_interpolation depends on section version
        if self.version == 0.3 and "fractional_bondorder_interpolation" not in kwargs:
            self.fractional_bondorder_method = "none"
        elif self.version == 0.4 and "fractional_bondorder_interpolation" not in kwargs:
            self.fractional_bondorder_method = "AM1-Wiberg"

        # Default value for potential depends on section version
        if self.version == 0.3 and "potential" not in kwargs:
            self.potential = "harmonic"
        elif self.version == 0.4 and "potential" not in kwargs:
            self.potential = "(k/2)*(r-length)^2"

    def check_handler_compatibility(self, other_handler):
        """
        Checks whether this ParameterHandler encodes compatible physics as another ParameterHandler. This is
        called if a second handler is attempted to be initialized for the same tag.

        Parameters
        ----------
        other_handler : a ParameterHandler object
            The handler to compare to.

        Raises
        ------
        IncompatibleParameterError if handler_kwargs are incompatible with existing parameters.
        """
        string_attrs_to_compare = [
            "fractional_bondorder_method",
            "fractional_bondorder_interpolation",
        ]
        self._check_attributes_are_equal(
            other_handler, identical_attrs=string_attrs_to_compare
        )

        # potential="harmonic" and potential="(k/2)*(r-length)^2" should be considered identical
        self_has_harmonic_potential = (
            self.potential == "harmonic" or self.potential == "(k/2)*(r-length)^2"
        )
        other_has_harmonic_potential = (
            other_handler.potential == "harmonic"
            or other_handler.potential == "(k/2)*(r-length)^2"
        )
        if not (self_has_harmonic_potential and other_has_harmonic_potential):
            if self.potential != other_handler.potential:
                raise IncompatibleParameterError(
                    f"potential values are not identical. "
                    f"(handler value: {self.potential}, incompatible value: {other_handler.potential}"
                )

    def create_force(self, system, topology, **kwargs):
        # Create or retrieve existing OpenMM Force object
        # TODO: The commented line below should replace the system.getForce search
        # force = super(BondHandler, self).create_force(system, topology, **kwargs)
        existing = [system.getForce(i) for i in range(system.getNumForces())]
        existing = [f for f in existing if type(f) == self._OPENMMTYPE]
        if len(existing) == 0:
            force = self._OPENMMTYPE()
            system.addForce(force)
        else:
            force = existing[0]

        # Do not trust previously-calculated partial bond orders, since we don't know
        # what method was used to assign them
        # TODO: Jeff tried implementing a way to mark how bond orders were assigned on the
        # topology, but realized that there's already a hierarchy of assignment
        # methods. That is, if a molecule was assigned using PBOs_from_mols, then
        # a different fractional bondorder method SHOULD NOT attempt
        # recalculation, whereas if the previous method was simply DIFFERENT,
        # then the old results should be erased/cached and overwritten with the
        # new ones. It will be easier to handle this at the level of caching
        # the results of molecule.assign_fractional_bond_orders
        for top_bond in topology.topology_bonds:
            top_bond.bond.fractional_bond_order = None

        # check whether any of the reference molecules in the topology
        # are in the partial_bond_orders_from_molecules list
        if "partial_bond_orders_from_molecules" in kwargs:
            # check whether molecules in the partial_bond_orders_from_molecules
            # list have any duplicates
            self.check_partial_bond_orders_from_molecules_duplicates(
                kwargs["partial_bond_orders_from_molecules"]
            )

            self.assign_partial_bond_orders_from_molecules(
                topology, kwargs["partial_bond_orders_from_molecules"]
            )

        # Add all bonds to the system.
        bond_matches = self.find_matches(topology)

        skipped_constrained_bonds = (
            0  # keep track of how many bonds were constrained (and hence skipped)
        )
        for (topology_atom_indices, bond_match) in bond_matches.items():
            # Get corresponding particle indices in Topology
            # particle_indices = tuple([ atom.particle_index for atom in atoms ])

            # Ensure atoms are actually bonded correct pattern in Topology
            try:
                self._assert_correct_connectivity(bond_match)
            except NotBondedError as e:
                smirks = bond_match.parameter_type.smirks
                raise NotBondedError(
                    f"While processing bond with SMIRKS {smirks}: " + e.msg
                )

            # topology.assert_bonded(atoms[0], atoms[1])
            bond_params = bond_match.parameter_type
            match = bond_match.environment_match

            # Compute equilibrium bond length and spring constant.
            bond = match.reference_molecule.get_bond_between(
                *match.reference_atom_indices
            )

            length_requires_interpolation = (
                getattr(bond_params, "length_bondorder", None) is not None
            )
            k_requires_interpolation = (
                getattr(bond_params, "k_bondorder", None) is not None
            )

            # Calculate fractional bond orders for this molecule only if needed.
            if (
                length_requires_interpolation or k_requires_interpolation
            ) and bond.fractional_bond_order is None:
                toolkit_registry = kwargs.get(
                    "toolkit_registry", GLOBAL_TOOLKIT_REGISTRY
                )
                match.reference_molecule.assign_fractional_bond_orders(
                    toolkit_registry=toolkit_registry,
                    bond_order_model=self.fractional_bondorder_method.lower(),
                )

            if not length_requires_interpolation:
                length = bond_params.length
            else:
                # Interpolate length using fractional bond orders
                bond_order = bond.fractional_bond_order
                if self.fractional_bondorder_interpolation == "linear":
                    if len(bond_params.length_bondorder) < 2:
                        raise SMIRNOFFSpecError(
                            "In order to use bond order interpolation, 2 or more parameters "
                            f"must be present. Found {len(bond_params.length_bondorder)} parameters."
                        )
                    length = _linear_inter_or_extrapolate(
                        points_dict=bond_params.length_bondorder,
                        x_query=bond_order,
                    )
                else:
                    # TODO: This code is effectively unreachable due to the the _allow_only converter used in this
                    #       ParameterAttribute's definition, which only allows "linear". Remove?
                    raise FractionalBondOrderInterpolationMethodUnsupportedError(
                        "Fractional bondorder interpolation method {} is not implemented.".format(
                            self.fractional_bondorder_interpolation
                        )
                    )
            if not k_requires_interpolation:
                k = bond_params.k
            else:
                # Interpolate k using fractional bond orders
                bond_order = bond.fractional_bond_order
                if self.fractional_bondorder_interpolation == "linear":
                    if len(bond_params.k_bondorder) < 2:
                        raise SMIRNOFFSpecError(
                            "In order to use bond order interpolation, 2 or more parameters "
                            f"must be present. Found {len(bond_params.k_bondorder)} parameters."
                        )
                    k = _linear_inter_or_extrapolate(
                        points_dict=bond_params.k_bondorder,
                        x_query=bond_order,
                    )
                else:
                    # TODO: This code is effectively unreachable due to the the _allow_only converter used in this
                    #       ParameterAttribute's definition, which only allows "linear". Remove?
                    raise FractionalBondOrderInterpolationMethodUnsupportedError(
                        "Fractional bondorder interpolation method {} is not implemented.".format(
                            self.fractional_bondorder_interpolation
                        )
                    )

            # If this pair of atoms is subject to a constraint, only use the length
            is_constrained = topology.is_constrained(*topology_atom_indices)
            if not is_constrained:
                # Add harmonic bond to HarmonicBondForce
                force.addBond(*topology_atom_indices, length, k)
            else:
                # Handle constraints.
                # Atom pair is constrained; we don't need to add a bond term.
                skipped_constrained_bonds += 1
                # Check if we need to add the constraint here to the equilibrium bond length.
                if is_constrained is True:
                    # Mark that we have now assigned a specific constraint distance to this constraint.
                    topology.add_constraint(*topology_atom_indices, length)
                    # Add the constraint to the System.
                    system.addConstraint(*topology_atom_indices, length)
                    # system.addConstraint(*particle_indices, length)

        logger.info(
            "{} bonds added ({} skipped due to constraints)".format(
                len(bond_matches) - skipped_constrained_bonds, skipped_constrained_bonds
            )
        )

        # Check that no topological bonds are missing force parameters.
        valence_terms = [list(b.atoms) for b in topology.topology_bonds]
        self._check_all_valence_terms_assigned(
            assigned_terms=bond_matches,
            valence_terms=valence_terms,
            exception_cls=UnassignedBondParameterException,
        )


class AngleHandler(ParameterHandler):
    """Handle SMIRNOFF ``<AngleForce>`` tags

    .. warning :: This API is experimental and subject to change.
    """

    class AngleType(ParameterType):
        """A SMIRNOFF angle type.

        .. warning :: This API is experimental and subject to change.
        """

        _VALENCE_TYPE = "Angle"  # ChemicalEnvironment valence type string expected by SMARTS string for this Handler
        _ELEMENT_NAME = "Angle"

        angle = ParameterAttribute(unit=unit.degree)
        k = ParameterAttribute(unit=unit.kilocalorie_per_mole / unit.degree ** 2)

    _TAGNAME = "Angles"  # SMIRNOFF tag name to process
    _INFOTYPE = AngleType  # class to hold force type info
    _OPENMMTYPE = openmm.HarmonicAngleForce  # OpenMM force class to create
    _DEPENDENCIES = [ConstraintHandler]  # ConstraintHandler must be executed first

    potential = ParameterAttribute(default="harmonic")

    def check_handler_compatibility(self, other_handler):
        """
        Checks whether this ParameterHandler encodes compatible physics as another ParameterHandler. This is
        called if a second handler is attempted to be initialized for the same tag.

        Parameters
        ----------
        other_handler : a ParameterHandler object
            The handler to compare to.

        Raises
        ------
        IncompatibleParameterError if handler_kwargs are incompatible with existing parameters.
        """
        string_attrs_to_compare = ["potential"]
        self._check_attributes_are_equal(
            other_handler, identical_attrs=string_attrs_to_compare
        )

    def create_force(self, system, topology, **kwargs):
        # force = super(AngleHandler, self).create_force(system, topology, **kwargs)
        existing = [system.getForce(i) for i in range(system.getNumForces())]
        existing = [f for f in existing if type(f) == self._OPENMMTYPE]
        if len(existing) == 0:
            force = self._OPENMMTYPE()
            system.addForce(force)
        else:
            force = existing[0]

        # Add all angles to the system.
        angle_matches = self.find_matches(topology)
        skipped_constrained_angles = (
            0  # keep track of how many angles were constrained (and hence skipped)
        )
        for (atoms, angle_match) in angle_matches.items():

            # Ensure atoms are actually bonded correct pattern in Topology
            # for (i, j) in [(0, 1), (1, 2)]:
            #     topology.assert_bonded(atoms[i], atoms[j])
            try:
                self._assert_correct_connectivity(angle_match)
            except NotBondedError as e:
                smirks = angle_match.parameter_type.smirks
                raise NotBondedError(
                    f"While processing angle with SMIRKS {smirks}: " + e.msg
                )

            if (
                topology.is_constrained(atoms[0], atoms[1])
                and topology.is_constrained(atoms[1], atoms[2])
                and topology.is_constrained(atoms[0], atoms[2])
            ):
                # Angle is constrained; we don't need to add an angle term.
                skipped_constrained_angles += 1
                continue

            angle = angle_match.parameter_type
            force.addAngle(*atoms, angle.angle, angle.k)

        logger.info(
            "{} angles added ({} skipped due to constraints)".format(
                len(angle_matches) - skipped_constrained_angles,
                skipped_constrained_angles,
            )
        )

        # Check that no topological angles are missing force parameters
        self._check_all_valence_terms_assigned(
            assigned_terms=angle_matches,
            valence_terms=list(topology.angles),
            exception_cls=UnassignedAngleParameterException,
        )


# TODO: There's a lot of duplicated code in ProperTorsionHandler and ImproperTorsionHandler
class ProperTorsionHandler(ParameterHandler):
    """Handle SMIRNOFF ``<ProperTorsionForce>`` tags

    .. warning :: This API is experimental and subject to change.
    """

    class ProperTorsionType(ParameterType):
        """A SMIRNOFF torsion type for proper torsions.

        .. warning :: This API is experimental and subject to change.
        """

        _VALENCE_TYPE = "ProperTorsion"
        _ELEMENT_NAME = "Proper"

        periodicity = IndexedParameterAttribute(converter=int)
        phase = IndexedParameterAttribute(unit=unit.degree)
        k = IndexedParameterAttribute(default=None, unit=unit.kilocalorie_per_mole)
        idivf = IndexedParameterAttribute(default=None, converter=float)

        # fractional bond order params
        k_bondorder = IndexedMappedParameterAttribute(
            default=None, unit=unit.kilocalorie_per_mole
        )

    _TAGNAME = "ProperTorsions"  # SMIRNOFF tag name to process
    _KWARGS = ["partial_bond_orders_from_molecules"]
    _INFOTYPE = ProperTorsionType  # info type to store
    _OPENMMTYPE = openmm.PeriodicTorsionForce  # OpenMM force class to create
    _MAX_SUPPORTED_SECTION_VERSION = 0.4

    potential = ParameterAttribute(
        default="k*(1+cos(periodicity*theta-phase))",
        converter=_allow_only(["k*(1+cos(periodicity*theta-phase))"]),
    )
    default_idivf = ParameterAttribute(default="auto")
    fractional_bondorder_method = ParameterAttribute(default="AM1-Wiberg")
    fractional_bondorder_interpolation = ParameterAttribute(
        default="linear", converter=_allow_only(["linear"])
    )

    def check_handler_compatibility(self, other_handler):
        """
        Checks whether this ParameterHandler encodes compatible physics as another ParameterHandler. This is
        called if a second handler is attempted to be initialized for the same tag.

        Parameters
        ----------
        other_handler : a ParameterHandler object
            The handler to compare to.

        Raises
        ------
        IncompatibleParameterError if handler_kwargs are incompatible with existing parameters.
        """
        float_attrs_to_compare = []
        string_attrs_to_compare = [
            "potential",
            "fractional_bondorder_method",
            "fractional_bondorder_interpolation",
        ]

        if self.default_idivf == "auto":
            string_attrs_to_compare.append("default_idivf")
        else:
            float_attrs_to_compare.append("default_idivf")

        self._check_attributes_are_equal(
            other_handler,
            identical_attrs=string_attrs_to_compare,
            tolerance_attrs=float_attrs_to_compare,
        )

    def create_force(self, system, topology, **kwargs):
        # force = super(ProperTorsionHandler, self).create_force(system, topology, **kwargs)
        existing = [system.getForce(i) for i in range(system.getNumForces())]
        existing = [f for f in existing if type(f) == self._OPENMMTYPE]

        if len(existing) == 0:
            force = self._OPENMMTYPE()
            system.addForce(force)
        else:
            force = existing[0]

        # Do not trust previously-calculated partial bond orders, since we don't know
        # what method was used to assign them
        # TODO: Jeff tried implementing a way to mark how bond orders were assigned on the
        # topology, but realized that there's already a hierarchy of assignment
        # methods. That is, if a molecule was assigned using PBOs_from_mols, then
        # a different fractional bondorder method SHOULD NOT attempt
        # recalculation, whereas if the previous method was simply DIFFERENT,
        # then the old results should be erased/cached and overwritten with the
        # new ones. It will be easier to handle this at the level of caching
        # the results of molecule.assign_fractional_bond_orders
        for top_bond in topology.topology_bonds:
            top_bond.bond.fractional_bond_order = None

        # check whether any of the reference molecules in the topology
        # are in the partial_bond_orders_from_molecules list
        if "partial_bond_orders_from_molecules" in kwargs:
            # check whether molecules in the partial_bond_orders_from_molecules
            # list have any duplicates
            self.check_partial_bond_orders_from_molecules_duplicates(
                kwargs["partial_bond_orders_from_molecules"]
            )

            self.assign_partial_bond_orders_from_molecules(
                topology, kwargs["partial_bond_orders_from_molecules"]
            )

        # find all proper torsions for which we have parameters
        # operates on reference molecules in topology
        # but gives back matches for atoms for instance molecules
        torsion_matches = self.find_matches(topology)

        for (atom_indices, torsion_match) in torsion_matches.items():
            # Ensure atoms are actually bonded correct pattern in Topology
            # Currently does nothing
            try:
                self._assert_correct_connectivity(torsion_match)
            except NotBondedError as e:
                smirks = torsion_match.parameter_type.smirks
                raise NotBondedError(
                    f"While processing torsion with SMIRKS {smirks}: " + e.msg
                )

            if torsion_match.parameter_type.k_bondorder is None:
                # TODO: add a check here that we have same number of terms for
                # `kX_bondorder*`, `periodicityX`, `phaseX`
                # only count a given `kX_bondorder*` once

                # assign torsion with no interpolation
                self._assign_torsion(atom_indices, torsion_match, force)
            else:
                # TODO: add a check here that we have same number of terms for
                # `kX_bondorder*`, `periodicityX`, `phaseX`
                # only count a given `kX_bondorder*` once

                # assign torsion with interpolation
                self._assign_fractional_bond_orders(
                    atom_indices, torsion_match, force, **kwargs
                )

        logger.info("{} torsions added".format(len(torsion_matches)))

        # Check that no topological torsions are missing force parameters

        # I can see the appeal of these kind of methods as an 'absolute' check
        # that things have gone well, but I think just making sure that the
        # reference molecule has been fully parametrised should have the same
        # effect! It would be good to eventually refactor things so that everything
        # is focused on the single unique molecules, and then simply just cloned
        # onto the system. It seems like John's proposed System object would do
        # exactly this.
        self._check_all_valence_terms_assigned(
            assigned_terms=torsion_matches,
            valence_terms=list(topology.propers),
            exception_cls=UnassignedProperTorsionParameterException,
        )

    def _assign_torsion(self, atom_indices, torsion_match, force):

        torsion_params = torsion_match.parameter_type

        for (periodicity, phase, k, idivf) in zip(
            torsion_params.periodicity,
            torsion_params.phase,
            torsion_params.k,
            torsion_params.idivf,
        ):

            if idivf == "auto":
                # TODO: Implement correct "auto" behavior
                raise NotImplementedError(
                    "The OpenForceField toolkit hasn't implemented "
                    "support for the torsion `idivf` value of 'auto'"
                )

            force.addTorsion(
                atom_indices[0],
                atom_indices[1],
                atom_indices[2],
                atom_indices[3],
                periodicity,
                phase,
                k / idivf,
            )

    def _assign_fractional_bond_orders(
        self, atom_indices, torsion_match, force, **kwargs
    ):
        from openff.toolkit.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY

        torsion_params = torsion_match.parameter_type
        match = torsion_match.environment_match

        for (periodicity, phase, k_bondorder, idivf) in zip(
            torsion_params.periodicity,
            torsion_params.phase,
            torsion_params.k_bondorder,
            torsion_params.idivf,
        ):

            if len(k_bondorder) < 2:
                raise ValueError(
                    "At least 2 bond order values required for `k_bondorder`; "
                    "got {}".format(len(k_bondorder))
                )

            if idivf == "auto":
                # TODO: Implement correct "auto" behavior
                raise NotImplementedError(
                    "The OpenForceField toolkit hasn't implemented "
                    "support for the torsion `idivf` value of 'auto'"
                )

            # get central bond for reference molecule
            central_bond = match.reference_molecule.get_bond_between(
                match.reference_atom_indices[1], match.reference_atom_indices[2]
            )

            # if fractional bond order not calculated yet, we calculate it
            # should only happen once per reference molecule for which we care
            # about fractional bond interpolation
            # and not at all for reference molecules we don't
            if central_bond.fractional_bond_order is None:
                toolkit_registry = kwargs.get(
                    "toolkit_registry", GLOBAL_TOOLKIT_REGISTRY
                )
                match.reference_molecule.assign_fractional_bond_orders(
                    toolkit_registry=toolkit_registry,
                    bond_order_model=self.fractional_bondorder_method.lower(),
                )

            # scale k based on the bondorder of the central bond
            if self.fractional_bondorder_interpolation == "linear":
                # we only interpolate on k
                k = _linear_inter_or_extrapolate(
                    k_bondorder, central_bond.fractional_bond_order
                )
            else:
                # TODO: This code is effectively unreachable due to the the _allow_only converter used in this
                #       ParameterAttribute's definition, which only allows "linear". Remove?
                raise FractionalBondOrderInterpolationMethodUnsupportedError(
                    "Fractional bondorder interpolation method {} is not implemented.".format(
                        self.fractional_bondorder_interpolation
                    )
                )

            # add a torsion with given parameters for topology atoms
            force.addTorsion(
                atom_indices[0],
                atom_indices[1],
                atom_indices[2],
                atom_indices[3],
                periodicity,
                phase,
                k / idivf,
            )


# TODO: There's a lot of duplicated code in ProperTorsionHandler and ImproperTorsionHandler
class ImproperTorsionHandler(ParameterHandler):
    """Handle SMIRNOFF ``<ImproperTorsionForce>`` tags

    .. warning :: This API is experimental and subject to change.
    """

    class ImproperTorsionType(ParameterType):
        """A SMIRNOFF torsion type for improper torsions.

        .. warning :: This API is experimental and subject to change.
        """

        _VALENCE_TYPE = "ImproperTorsion"
        _ELEMENT_NAME = "Improper"

        periodicity = IndexedParameterAttribute(converter=int)
        phase = IndexedParameterAttribute(unit=unit.degree)
        k = IndexedParameterAttribute(unit=unit.kilocalorie_per_mole)
        idivf = IndexedParameterAttribute(default=None, converter=float)

    _TAGNAME = "ImproperTorsions"  # SMIRNOFF tag name to process
    _INFOTYPE = ImproperTorsionType  # info type to store
    _OPENMMTYPE = openmm.PeriodicTorsionForce  # OpenMM force class to create

    potential = ParameterAttribute(
        default="k*(1+cos(periodicity*theta-phase))",
        converter=_allow_only(["k*(1+cos(periodicity*theta-phase))"]),
    )
    default_idivf = ParameterAttribute(default="auto")

    def check_handler_compatibility(self, other_handler):
        """
        Checks whether this ParameterHandler encodes compatible physics as another ParameterHandler. This is
        called if a second handler is attempted to be initialized for the same tag.

        Parameters
        ----------
        other_handler : a ParameterHandler object
            The handler to compare to.

        Raises
        ------
        IncompatibleParameterError if handler_kwargs are incompatible with existing parameters.
        """
        float_attrs_to_compare = []
        string_attrs_to_compare = ["potential"]

        if self.default_idivf == "auto":
            string_attrs_to_compare.append("default_idivf")
        else:
            float_attrs_to_compare.append("default_idivf")

        self._check_attributes_are_equal(
            other_handler,
            identical_attrs=string_attrs_to_compare,
            tolerance_attrs=float_attrs_to_compare,
        )

    def find_matches(self, entity, unique=False):
        """Find the improper torsions in the topology/molecule matched by a parameter type.

        Parameters
        ----------
        entity : openff.toolkit.topology.Topology
            Topology to search.

        Returns
        ---------
        matches : ImproperDict[Tuple[int], ParameterHandler._Match]
            ``matches[atom_indices]`` is the ``ParameterType`` object
            matching the 4-tuple of atom indices in ``entity``.

        """
        return self._find_matches(
            entity, transformed_dict_cls=ImproperDict, unique=unique
        )

    def create_force(self, system, topology, **kwargs):
        # force = super(ImproperTorsionHandler, self).create_force(system, topology, **kwargs)
        # force = super().create_force(system, topology, **kwargs)
        existing = [system.getForce(i) for i in range(system.getNumForces())]
        existing = [f for f in existing if type(f) == openmm.PeriodicTorsionForce]
        if len(existing) == 0:
            force = openmm.PeriodicTorsionForce()
            system.addForce(force)
        else:
            force = existing[0]

        # Add all improper torsions to the system
        improper_matches = self.find_matches(topology)
        for (atom_indices, improper_match) in improper_matches.items():
            # Ensure atoms are actually bonded correct pattern in Topology
            # For impropers, central atom is atom 1
            # for (i, j) in [(0, 1), (1, 2), (1, 3)]:
            #     topology.assert_bonded(atom_indices[i], atom_indices[j])
            try:
                self._assert_correct_connectivity(
                    improper_match, [(0, 1), (1, 2), (1, 3)]
                )
            except NotBondedError as e:
                smirks = improper_match.parameter_type.smirks
                raise NotBondedError(
                    f"While processing improper with SMIRKS {smirks}: " + e.msg
                )

            improper = improper_match.parameter_type

            # TODO: This is a lazy hack. idivf should be set according to the ParameterHandler's default_idivf attrib
            if improper.idivf is None:
                improper.idivf = [3 for item in improper.k]
            # Impropers are applied in three paths around the trefoil having the same handedness
            for (
                improper_periodicity,
                improper_phase,
                improper_k,
                improper_idivf,
            ) in zip(improper.periodicity, improper.phase, improper.k, improper.idivf):
                # TODO: Implement correct "auto" behavior
                if improper_idivf == "auto":
                    improper_idivf = 3
                    logger.warning(
                        "The OpenForceField toolkit hasn't implemented "
                        "support for the torsion `idivf` value of 'auto'."
                        "Currently assuming a value of '3' for impropers."
                    )
                # Permute non-central atoms
                others = [atom_indices[0], atom_indices[2], atom_indices[3]]
                # ((0, 1, 2), (1, 2, 0), and (2, 0, 1)) are the three paths around the trefoil
                for p in [
                    (others[i], others[j], others[k])
                    for (i, j, k) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]
                ]:
                    # The torsion force gets added three times, since the k is divided by three
                    force.addTorsion(
                        atom_indices[1],
                        p[0],
                        p[1],
                        p[2],
                        improper_periodicity,
                        improper_phase,
                        improper_k / improper_idivf,
                    )
        logger.info(
            "{} impropers added, each applied in a six-fold trefoil".format(
                len(improper_matches)
            )
        )


class _NonbondedHandler(ParameterHandler):
    """Base class for ParameterHandlers that deal with OpenMM NonbondedForce objects."""

    _OPENMMTYPE = openmm.NonbondedForce

    def create_force(self, system, topology, **kwargs):
        # If we aren't yet keeping track of which molecules' charges have been assigned by which charge methods,
        # initialize a dict for that here.
        # TODO: This should be an attribute of the _system_, not the _topology_. However, since we're still using
        #  OpenMM's System class, I am storing this data on the OFF Topology until we make an OFF System class.
        if not hasattr(topology, "_ref_mol_to_charge_method"):
            topology._ref_mol_to_charge_method = {
                ref_mol: None for ref_mol in topology.reference_molecules
            }

        # Retrieve the system's OpenMM NonbondedForce
        existing = [system.getForce(i) for i in range(system.getNumForces())]
        existing = [f for f in existing if type(f) == self._OPENMMTYPE]

        # If there isn't yet one, initialize it and populate it with particles
        if len(existing) == 0:
            force = self._OPENMMTYPE()
            system.addForce(force)
            # Create all atom particles. Virtual site particles are handled in
            # in its own handler
            for _ in topology.topology_atoms:
                force.addParticle(0.0, 1.0, 0.0)
        else:
            force = existing[0]

        return force

    def mark_charges_assigned(self, ref_mol, topology):
        """
        Record that charges have been assigned for a reference molecule.

        Parameters
        ----------
        ref_mol : openff.toolkit.topology.Molecule
            The molecule to mark as having charges assigned
        topology : openff.toolkit.topology.Topology
            The topology to record this information on.
        """
        # TODO: Change this to interface with system object instead of topology once we move away from OMM's System
        topology._ref_mol_to_charge_method[ref_mol] = self.__class__

    @staticmethod
    def check_charges_assigned(ref_mol, topology):
        """
        Check whether charges have been assigned for a reference molecule.

        Parameters
        ----------
        ref_mol : openff.toolkit.topology.Molecule
            The molecule to check for having charges assigned
        topology : openff.toolkit.topology.Topology
            The topology to query for this information

        Returns
        -------
        charges_assigned : bool
            Whether charges have already been assigned to this molecule

        """
        # TODO: Change this to interface with system object instead of topology once we move away from OMM's System
        return topology._ref_mol_to_charge_method[ref_mol] is not None


class vdWHandler(_NonbondedHandler):
    """Handle SMIRNOFF ``<vdW>`` tags

    .. warning :: This API is experimental and subject to change.
    """

    class vdWType(ParameterType):
        """A SMIRNOFF vdWForce type.

        .. warning :: This API is experimental and subject to change.
        """

        _VALENCE_TYPE = "Atom"  # ChemicalEnvironment valence type expected for SMARTS
        _ELEMENT_NAME = "Atom"

        epsilon = ParameterAttribute(unit=unit.kilocalorie_per_mole)
        sigma = ParameterAttribute(default=None, unit=unit.angstrom)
        rmin_half = ParameterAttribute(default=None, unit=unit.angstrom)

        def __init__(self, **kwargs):
            sigma = kwargs.get("sigma", None)
            rmin_half = kwargs.get("rmin_half", None)
            if (sigma is None) and (rmin_half is None):
                raise SMIRNOFFSpecError("Either sigma or rmin_half must be specified.")
            if (sigma is not None) and (rmin_half is not None):
                raise SMIRNOFFSpecError(
                    "BOTH sigma and rmin_half cannot be specified simultaneously."
                )

            super().__init__(**kwargs)

            if sigma:
                self._extra_nb_var = "rmin_half"
            if rmin_half:
                self._extra_nb_var = "sigma"

        def __setattr__(self, name, value):
            super().__setattr__(key=name, value=value)
            if name == "rmin_half":
                if type(value) == str:
                    value = object_to_quantity(value)
                super().__setattr__("sigma", 2.0 * value / 2 ** (1 / 6))
                self._extra_nb_var = "sigma"

            if name == "sigma":
                if type(value) == str:
                    value = object_to_quantity(value)
                super().__setattr__("rmin_half", value * 2 ** (1 / 6) / 2.0)
                self._extra_nb_var = "rmin_half"

        def to_dict(
            self,
            discard_cosmetic_attributes=False,
            duplicate_attributes=None,
        ):
            return super().to_dict(
                discard_cosmetic_attributes=discard_cosmetic_attributes,
                duplicate_attributes=[self._extra_nb_var],
            )

    _TAGNAME = "vdW"  # SMIRNOFF tag name to process
    _INFOTYPE = vdWType  # info type to store
    # _KWARGS = ['ewaldErrorTolerance',
    #            'useDispersionCorrection',
    #            'usePbc'] # Kwargs to catch when create_force is called

    potential = ParameterAttribute(
        default="Lennard-Jones-12-6", converter=_allow_only(["Lennard-Jones-12-6"])
    )
    combining_rules = ParameterAttribute(
        default="Lorentz-Berthelot", converter=_allow_only(["Lorentz-Berthelot"])
    )

    scale12 = ParameterAttribute(default=0.0, converter=float)
    scale13 = ParameterAttribute(default=0.0, converter=float)
    scale14 = ParameterAttribute(default=0.5, converter=float)
    scale15 = ParameterAttribute(default=1.0, converter=float)

    cutoff = ParameterAttribute(default=9.0 * unit.angstroms, unit=unit.angstrom)
    switch_width = ParameterAttribute(default=1.0 * unit.angstroms, unit=unit.angstrom)
    method = ParameterAttribute(
        default="cutoff", converter=_allow_only(["cutoff", "PME"])
    )

    # TODO: Use _allow_only when ParameterAttribute will support multiple converters (it'll be easy when we switch to use the attrs library)
    @scale12.converter
    def scale12(self, attrs, new_scale12):
        if new_scale12 != 0.0:
            raise SMIRNOFFSpecError(
                "Current OFF toolkit is unable to handle scale12 values other than 0.0. "
                "Specified 1-2 scaling was {}".format(self.scale12)
            )
        return new_scale12

    @scale13.converter
    def scale13(self, attrs, new_scale13):
        if new_scale13 != 0.0:
            raise SMIRNOFFSpecError(
                "Current OFF toolkit is unable to handle scale13 values other than 0.0. "
                "Specified 1-3 scaling was {}".format(self.scale13)
            )
        return new_scale13

    @scale15.converter
    def scale15(self, attrs, new_scale15):
        if new_scale15 != 1.0:
            raise SMIRNOFFSpecError(
                "Current OFF toolkit is unable to handle scale15 values other than 1.0. "
                "Specified 1-5 scaling was {}".format(self.scale15)
            )
        return new_scale15

    # Tolerance when comparing float attributes for handler compatibility.
    _SCALETOL = 1e-5

    def check_handler_compatibility(self, other_handler):
        """
        Checks whether this ParameterHandler encodes compatible physics as another ParameterHandler. This is
        called if a second handler is attempted to be initialized for the same tag.

        Parameters
        ----------
        other_handler : a ParameterHandler object
            The handler to compare to.

        Raises
        ------
        IncompatibleParameterError if handler_kwargs are incompatible with existing parameters.
        """
        float_attrs_to_compare = ["scale12", "scale13", "scale14", "scale15"]
        string_attrs_to_compare = ["potential", "combining_rules", "method"]
        unit_attrs_to_compare = ["cutoff"]

        self._check_attributes_are_equal(
            other_handler,
            identical_attrs=string_attrs_to_compare,
            tolerance_attrs=float_attrs_to_compare + unit_attrs_to_compare,
            tolerance=self._SCALETOL,
        )

    def create_force(self, system, topology, **kwargs):
        force = super().create_force(system, topology, **kwargs)

        # If we're using PME, then the only possible openMM Nonbonded type is LJPME
        if self.method == "PME":
            # If we're given a nonperiodic box, we always set NoCutoff. Later we'll add support for CutoffNonPeriodic
            if topology.box_vectors is None:
                force.setNonbondedMethod(openmm.NonbondedForce.NoCutoff)
                # if (topology.box_vectors is None):
                #     raise SMIRNOFFSpecError("If vdW method is  PME, a periodic Topology "
                #                             "must be provided")
            else:
                force.setNonbondedMethod(openmm.NonbondedForce.LJPME)
                force.setCutoffDistance(self.cutoff)
                force.setEwaldErrorTolerance(1.0e-4)

        # If method is cutoff, then we currently support openMM's PME for periodic system and NoCutoff for nonperiodic
        elif self.method == "cutoff":
            # If we're given a nonperiodic box, we always set NoCutoff. Later we'll add support for CutoffNonPeriodic
            if topology.box_vectors is None:
                force.setNonbondedMethod(openmm.NonbondedForce.NoCutoff)
            else:
                force.setNonbondedMethod(openmm.NonbondedForce.PME)
                force.setUseDispersionCorrection(True)
                force.setCutoffDistance(self.cutoff)

        # Iterate over all defined Lennard-Jones types, allowing later matches to override earlier ones.
        atom_matches = self.find_matches(topology)

        # Set the particle Lennard-Jones terms.
        for atom_key, atom_match in atom_matches.items():
            atom_idx = atom_key[0]
            ljtype = atom_match.parameter_type
            if ljtype.sigma is None:
                sigma = 2.0 * ljtype.rmin_half / (2.0 ** (1.0 / 6.0))
            else:
                sigma = ljtype.sigma
            force.setParticleParameters(atom_idx, 0.0, sigma, ljtype.epsilon)

        # Check that no atoms (n.b. not particles) are missing force parameters.
        self._check_all_valence_terms_assigned(
            assigned_terms=atom_matches, valence_terms=list(topology.topology_atoms)
        )


class ElectrostaticsHandler(_NonbondedHandler):
    """Handles SMIRNOFF ``<Electrostatics>`` tags.

    .. warning :: This API is experimental and subject to change.
    """

    _TAGNAME = "Electrostatics"
    _DEPENDENCIES = [vdWHandler]
    _KWARGS = ["charge_from_molecules", "allow_nonintegral_charges"]

    scale12 = ParameterAttribute(default=0.0, converter=float)
    scale13 = ParameterAttribute(default=0.0, converter=float)
    scale14 = ParameterAttribute(default=0.833333, converter=float)
    scale15 = ParameterAttribute(default=1.0, converter=float)
    cutoff = ParameterAttribute(default=9.0 * unit.angstrom, unit=unit.angstrom)
    switch_width = ParameterAttribute(default=0.0 * unit.angstrom, unit=unit.angstrom)
    method = ParameterAttribute(
        default="PME", converter=_allow_only(["Coulomb", "PME", "reaction-field"])
    )

    # TODO: Use _allow_only when ParameterAttribute will support multiple converters (it'll be easy when we switch to use the attrs library)
    @scale12.converter
    def scale12(self, attrs, new_scale12):
        if new_scale12 != 0.0:
            raise SMIRNOFFSpecError(
                "Current OFF toolkit is unable to handle scale12 values other than 0.0. "
                "Specified 1-2 scaling was {}".format(self.scale12)
            )
        return new_scale12

    @scale13.converter
    def scale13(self, attrs, new_scale13):
        if new_scale13 != 0.0:
            raise SMIRNOFFSpecError(
                "Current OFF toolkit is unable to handle scale13 values other than 0.0. "
                "Specified 1-3 scaling was {}".format(self.scale13)
            )
        return new_scale13

    @scale15.converter
    def scale15(self, attrs, new_scale15):
        if new_scale15 != 1.0:
            raise SMIRNOFFSpecError(
                "Current OFF toolkit is unable to handle scale15 values other than 1.0. "
                "Specified 1-5 scaling was {}".format(self.scale15)
            )
        return new_scale15

    @switch_width.converter
    def switch_width(self, attr, new_switch_width):
        if self.switch_width != 0.0 * unit.angstrom:
            raise IncompatibleParameterError(
                "The current implementation of the Open Force Field Toolkit can not "
                "support an electrostatic switching width. Currently only `0.0 angstroms` "
                f"is supported (SMIRNOFF data specified {new_switch_width})"
            )

    # Tolerance when comparing float attributes for handler compatibility.
    _SCALETOL = 1e-5

    def check_handler_compatibility(self, other_handler):
        """
        Checks whether this ParameterHandler encodes compatible physics as another ParameterHandler. This is
        called if a second handler is attempted to be initialized for the same tag.

        Parameters
        ----------
        other_handler : a ParameterHandler object
            The handler to compare to.

        Raises
        ------
        IncompatibleParameterError if handler_kwargs are incompatible with existing parameters.
        """
        float_attrs_to_compare = ["scale12", "scale13", "scale14", "scale15"]
        string_attrs_to_compare = ["method"]
        unit_attrs_to_compare = ["cutoff", "switch_width"]

        self._check_attributes_are_equal(
            other_handler,
            identical_attrs=string_attrs_to_compare,
            tolerance_attrs=float_attrs_to_compare + unit_attrs_to_compare,
            tolerance=self._SCALETOL,
        )

    def assign_charge_from_molecules(self, molecule, charge_mols):
        """
        Given an input molecule, checks against a list of molecules for an isomorphic match. If found, assigns
        partial charges from the match to the input molecule.

        Parameters
        ----------
        molecule : an openff.toolkit.topology.FrozenMolecule
            The molecule to have partial charges assigned if a match is found.
        charge_mols : list of [openff.toolkit.topology.FrozenMolecule]
            A list of molecules with charges already assigned.

        Returns
        -------
        match_found : bool
            Whether a match was found. If True, the input molecule will have been modified in-place.
        """
        # Check each charge_mol for whether it's isomorphic to the input molecule
        for charge_mol in charge_mols:
            ismorphic, topology_atom_map = Molecule.are_isomorphic(
                molecule,
                charge_mol,
                return_atom_map=True,
                aromatic_matching=True,
                formal_charge_matching=True,
                bond_order_matching=True,
                atom_stereochemistry_matching=True,
                bond_stereochemistry_matching=True,
            )
            # if they are isomorphic then use the mapping
            if ismorphic:
                # Take the first valid atom indexing map
                # Set the partial charges
                # Make a copy of the charge molecule's charges array (this way it's the right shape)
                temp_mol_charges = copy.deepcopy(
                    openmm.unit.Quantity(charge_mol.partial_charges)
                )
                for charge_idx, ref_idx in topology_atom_map.items():
                    temp_mol_charges[ref_idx] = charge_mol.partial_charges[charge_idx]
                molecule.partial_charges = temp_mol_charges
                return True

        # If no match was found, return False
        return False

    def create_force(self, system, topology, **kwargs):
        from openff.toolkit.topology import TopologyAtom, TopologyVirtualSite

        force = super().create_force(system, topology, **kwargs)

        # See if each molecule should have charges assigned by the charge_from_molecules kwarg
        for ref_mol in topology.reference_molecules:

            # If charges were already assigned, skip this molecule
            if self.check_charges_assigned(ref_mol, topology):
                continue

            # First, check whether any of the reference molecules in the topology are in the charge_from_mol list
            charges_from_charge_mol = False
            if "charge_from_molecules" in kwargs:
                charges_from_charge_mol = self.assign_charge_from_molecules(
                    ref_mol, kwargs["charge_from_molecules"]
                )

            # If this reference molecule wasn't in the charge_from_molecules list, end this iteration
            if not (charges_from_charge_mol):
                continue

            # Otherwise, the molecule is in the charge_from_molecules list, and we should assign charges to all
            # instances of it in this topology.
            for topology_molecule in topology._reference_molecule_to_topology_molecules[
                ref_mol
            ]:

                for topology_particle in topology_molecule.particles:

                    if type(topology_particle) is TopologyAtom:
                        ref_mol_particle_index = (
                            topology_particle.atom.molecule_particle_index
                        )
                    elif type(topology_particle) is TopologyVirtualSite:
                        ref_mol_particle_index = (
                            topology_particle.virtual_site.molecule_particle_index
                        )
                    else:
                        raise ValueError(
                            f"Particles of type {type(topology_particle)} are not supported"
                        )

                    topology_particle_index = topology_particle.topology_particle_index

                    particle_charge = ref_mol._partial_charges[ref_mol_particle_index]

                    # Retrieve nonbonded parameters for reference atom (charge not set yet)
                    _, sigma, epsilon = force.getParticleParameters(
                        topology_particle_index
                    )
                    # Set the nonbonded force with the partial charge
                    force.setParticleParameters(
                        topology_particle_index, particle_charge, sigma, epsilon
                    )

            # Finally, mark that charges were assigned for this reference molecule
            self.mark_charges_assigned(ref_mol, topology)

        # Set the nonbonded method
        current_nb_method = force.getNonbondedMethod()

        # First, check whether the vdWHandler set the nonbonded method to LJPME, because that means
        # that electrostatics also has to be PME
        if (current_nb_method == openmm.NonbondedForce.LJPME) and (
            self.method != "PME"
        ):
            raise IncompatibleParameterError(
                "In current Open Force Field Toolkit implementation, if vdW "
                "treatment is set to LJPME, electrostatics must also be PME "
                "(electrostatics treatment currently set to {}".format(self.method)
            )

        # Then, set nonbonded methods based on method keyword
        if self.method == "PME":
            # Check whether the topology is nonperiodic, in which case we always switch to NoCutoff
            # (vdWHandler will have already set this to NoCutoff)
            # TODO: This is an assumption right now, and a bad one. See issue #219
            if topology.box_vectors is None:
                assert current_nb_method == openmm.NonbondedForce.NoCutoff
                force.setCutoffDistance(self.cutoff)
                # raise IncompatibleParameterError("Electrostatics handler received PME method keyword, but a nonperiodic"
                #                                  " topology. Use of PME electrostatics requires a periodic topology.")
            else:
                if current_nb_method == openmm.NonbondedForce.LJPME:
                    pass
                    # There's no need to check for matching cutoff/tolerance here since both are hard-coded defaults
                else:
                    force.setNonbondedMethod(openmm.NonbondedForce.PME)
                    force.setCutoffDistance(self.cutoff)
                    force.setEwaldErrorTolerance(1.0e-4)

        # If vdWHandler set the nonbonded method to NoCutoff, then we don't need to change anything
        elif self.method == "Coulomb":
            if topology.box_vectors is None:
                # (vdWHandler will have already set this to NoCutoff)
                assert current_nb_method == openmm.NonbondedForce.NoCutoff
            else:
                raise IncompatibleParameterError(
                    "Electrostatics method set to Coulomb, and topology is periodic. "
                    "In the future, this will lead to use of OpenMM's CutoffPeriodic "
                    "Nonbonded force method, however this is not supported in the "
                    "current Open Force Field Toolkit."
                )

        # If the vdWHandler set the nonbonded method to PME, then ensure that it has the same cutoff
        elif self.method == "reaction-field":
            if topology.box_vectors is None:
                raise SMIRNOFFSpecError(
                    "Electrostatics method reaction-field can only be applied to a periodic system."
                )

            else:
                raise SMIRNOFFSpecUnimplementedError(
                    "Electrostatics method reaction-field is supported in the SMIRNOFF specification "
                    "but not yet implemented in the OpenFF Toolkit."
                )

    def postprocess_system(self, system, topology, **kwargs):
        force = super().create_force(system, topology, **kwargs)
        # Check to ensure all molecules have had charges assigned
        uncharged_mols = []
        for ref_mol in topology.reference_molecules:
            if not self.check_charges_assigned(ref_mol, topology):
                uncharged_mols.append(ref_mol)

        if len(uncharged_mols) != 0:
            msg = "The following molecules did not have charges assigned by any ParameterHandler in the ForceField:\n"
            for ref_mol in uncharged_mols:
                msg += f"{ref_mol.to_smiles()}\n"
            raise UnassignedMoleculeChargeException(msg)

        # Unless check is disabled, ensure that the sum of partial charges on a molecule
        # add up to approximately its formal charge
        allow_nonintegral_charges = kwargs.get("allow_nonintegral_charges", False)
        for top_mol in topology.topology_molecules:
            # Skip check if user manually disables it.
            if allow_nonintegral_charges:
                continue
            formal_charge_sum = top_mol.reference_molecule.total_charge
            partial_charge_sum = 0.0 * unit.elementary_charge
            for top_particle in top_mol.particles:
                q, _, _ = force.getParticleParameters(
                    top_particle.topology_particle_index
                )
                partial_charge_sum += q
            if (
                abs(formal_charge_sum - partial_charge_sum)
                > 0.01 * unit.elementary_charge
            ):
                msg = (
                    f"Partial charge sum ({partial_charge_sum}) "
                    f"for molecule '{top_mol.reference_molecule.name}' (SMILES "
                    f"{top_mol.reference_molecule.to_smiles()} does not equal formal charge sum "
                    f"({formal_charge_sum}). To override this error, provide the "
                    f"'allow_nonintegral_charges=True' keyword to ForceField.create_openmm_system"
                )
                raise NonintegralMoleculeChargeException(msg)


class LibraryChargeHandler(_NonbondedHandler):
    """Handle SMIRNOFF ``<LibraryCharges>`` tags

    .. warning :: This API is experimental and subject to change.
    """

    class LibraryChargeType(ParameterType):
        """A SMIRNOFF Library Charge type.

        .. warning :: This API is experimental and subject to change.
        """

        _VALENCE_TYPE = None  # This disables the connectivity check when parsing LibraryChargeType objects
        _ELEMENT_NAME = "LibraryCharge"

        name = ParameterAttribute(default=None)
        charge = IndexedParameterAttribute(unit=unit.elementary_charge)

        def __init__(self, **kwargs):
            super().__init__(**kwargs)
            unique_tags, connectivity = GLOBAL_TOOLKIT_REGISTRY.call(
                "get_tagged_smarts_connectivity", self.smirks
            )
            if len(self.charge) != len(unique_tags):
                raise SMIRNOFFSpecError(
                    f"LibraryCharge {self} was initialized with unequal number of "
                    f"tagged atoms and charges"
                )

        @classmethod
        def from_molecule(cls, molecule):
            """Construct a LibraryChargeType from a molecule with existing partial charges."""
            if molecule.partial_charges is None:
                raise ValueError("Input molecule is missing partial charges.")

            smirks = molecule.to_smiles(mapped=True)
            charges = molecule.partial_charges

            library_charge_type = cls(smirks=smirks, charge=charges)

            return library_charge_type

    _TAGNAME = "LibraryCharges"  # SMIRNOFF tag name to process
    _INFOTYPE = LibraryChargeType  # info type to store
    _DEPENDENCIES = [vdWHandler, ElectrostaticsHandler]

    def find_matches(self, entity, unique=True):
        """Find the elements of the topology/molecule matched by a parameter type.

        Parameters
        ----------
        entity : openff.toolkit.topology.Topology
            Topology to search.

        Returns
        ---------
        matches : ValenceDict[Tuple[int], ParameterHandler._Match]
            ``matches[particle_indices]`` is the ``ParameterType`` object
            matching the tuple of particle indices in ``entity``.
        """

        # TODO: Right now, this method is only ever called with an entity that is a Topology.
        #  Should we reduce its scope and have a check here to make sure entity is a Topology?
        return self._find_matches(
            entity,
            transformed_dict_cls=dict,
            unique=unique,
        )

    def create_force(self, system, topology, **kwargs):
        force = super().create_force(system, topology, **kwargs)

        # Iterate over all defined library charge parameters, allowing later matches to override earlier ones.
        atom_matches = self.find_matches(topology)

        # Create a set of all the topology atom indices for which library charges can be applied
        assignable_atoms = set()
        atom_assignments = dict()
        # TODO: This assumes that later matches should always override earlier ones. This may require more
        #       thought, since matches can be partially overlapping
        for topology_indices, library_charge in atom_matches.items():
            for charge_idx, top_idx in enumerate(topology_indices):
                if top_idx in assignable_atoms:
                    logger.debug(
                        f"Multiple library charge assignments found for atom {top_idx}"
                    )
                assignable_atoms.add(top_idx)
                atom_assignments[top_idx] = library_charge.parameter_type.charge[
                    charge_idx
                ]
        # TODO: Should header include a residue separator delimiter? Maybe not, since it's not clear how having
        #       multiple LibraryChargeHandlers could return a single set of matches while respecting different
        #       separators.

        # Keep track of the reference molecules that this successfully assigns charges to, so we can
        # mark them and subsequent charge generation handlers won't override the values
        ref_mols_assigned = set()

        # Check to see whether the set contains any complete molecules, and remove the matches if not.
        for top_mol in topology.topology_molecules:

            # Make a temporary copy of ref_mol to assign charges from charge_mol

            # If charges were already assigned, skip this molecule
            if self.check_charges_assigned(top_mol.reference_molecule, topology):
                continue

            # Ensure all of the atoms in this mol are covered, otherwise skip it
            top_particle_idxs = [atom.topology_particle_index for atom in top_mol.atoms]
            if (
                len(set(top_particle_idxs).intersection(assignable_atoms))
                != top_mol.n_atoms
            ):
                logger.debug(
                    "Entire molecule is not covered. Skipping library charge assignment."
                )
                continue

            # If we pass both tests above, go ahead and assign charges
            # TODO: We could probably save a little time by looking up this TopologyMolecule's _reference molecule_
            #       and assigning charges to all other instances of it in this topology
            for top_particle_idx in top_particle_idxs:
                _, sigma, epsilon = force.getParticleParameters(top_particle_idx)
                force.setParticleParameters(
                    top_particle_idx, atom_assignments[top_particle_idx], sigma, epsilon
                )

            ref_mols_assigned.add(top_mol.reference_molecule)

        # Finally, mark that charges were assigned for this reference molecule
        for assigned_mol in ref_mols_assigned:
            self.mark_charges_assigned(assigned_mol, topology)


class ToolkitAM1BCCHandler(_NonbondedHandler):
    """Handle SMIRNOFF ``<ToolkitAM1BCC>`` tags

    .. warning :: This API is experimental and subject to change.
    """

    _TAGNAME = "ToolkitAM1BCC"  # SMIRNOFF tag name to process
    _DEPENDENCIES = [vdWHandler, ElectrostaticsHandler, LibraryChargeHandler]
    _KWARGS = ["toolkit_registry"]  # Kwargs to catch when create_force is called

    def check_handler_compatibility(
        self, other_handler, assume_missing_is_default=True
    ):
        """
        Checks whether this ParameterHandler encodes compatible physics as another ParameterHandler. This is
        called if a second handler is attempted to be initialized for the same tag.

        Parameters
        ----------
        other_handler : a ParameterHandler object
            The handler to compare to.

        Raises
        ------
        IncompatibleParameterError if handler_kwargs are incompatible with existing parameters.
        """
        pass

    def create_force(self, system, topology, **kwargs):
        import warnings

        from openff.toolkit.topology import TopologyAtom, TopologyVirtualSite
        from openff.toolkit.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY

        force = super().create_force(system, topology, **kwargs)

        for ref_mol in topology.reference_molecules:

            # If charges were already assigned, skip this molecule
            if self.check_charges_assigned(ref_mol, topology):
                continue

            # If the molecule wasn't already assigned charge values, calculate them here
            toolkit_registry = kwargs.get("toolkit_registry", GLOBAL_TOOLKIT_REGISTRY)
            try:
                # We don't need to generate conformers here, since that will be done by default in
                # compute_partial_charges with am1bcc if the use_conformers kwarg isn't defined
                ref_mol.assign_partial_charges(
                    partial_charge_method="am1bcc", toolkit_registry=toolkit_registry
                )
            except Exception as e:
                warnings.warn(str(e), Warning)
                continue

            # Assign charges to relevant atoms
            for topology_molecule in topology._reference_molecule_to_topology_molecules[
                ref_mol
            ]:
                for topology_particle in topology_molecule.atoms:
                    if type(topology_particle) is TopologyAtom:
                        ref_mol_particle_index = (
                            topology_particle.atom.molecule_particle_index
                        )
                    elif type(topology_particle) is TopologyVirtualSite:
                        ref_mol_particle_index = (
                            topology_particle.virtual_site.molecule_particle_index
                        )
                    else:
                        raise ValueError(
                            f"Particles of type {type(topology_particle)} are not supported"
                        )

                    topology_particle_index = topology_particle.topology_particle_index

                    particle_charge = ref_mol._partial_charges[ref_mol_particle_index]

                    # Retrieve nonbonded parameters for reference atom (charge not set yet)
                    _, sigma, epsilon = force.getParticleParameters(
                        topology_particle_index
                    )
                    # Set the nonbonded force with the partial charge
                    force.setParticleParameters(
                        topology_particle_index, particle_charge, sigma, epsilon
                    )
            # Finally, mark that charges were assigned for this reference molecule
            self.mark_charges_assigned(ref_mol, topology)

    # TODO: Move chargeModel and library residue charges to SMIRNOFF spec
    def postprocess_system(self, system, topology, **kwargs):

        bond_matches = self.find_matches(topology)

        # Apply bond charge increments to all appropriate force groups
        # QUESTION: Should we instead apply this to the Topology in a preprocessing step, prior to spreading out charge onto virtual sites?
        for force in system.getForces():
            if force.__class__.__name__ in [
                "NonbondedForce"
            ]:  # TODO: We need to apply this to all Force types that involve charges, such as (Custom)GBSA forces and CustomNonbondedForce
                for (atoms, bond_match) in bond_matches.items():
                    # Get corresponding particle indices in Topology
                    bond = bond_match.parameter_type

                    particle_indices = tuple([atom.particle_index for atom in atoms])
                    # Retrieve parameters
                    [charge0, sigma0, epsilon0] = force.getParticleParameters(
                        particle_indices[0]
                    )
                    [charge1, sigma1, epsilon1] = force.getParticleParameters(
                        particle_indices[1]
                    )
                    # Apply bond charge increment
                    charge0 -= bond.increment
                    charge1 += bond.increment
                    # Update charges
                    force.setParticleParameters(
                        particle_indices[0], charge0, sigma0, epsilon0
                    )
                    force.setParticleParameters(
                        particle_indices[1], charge1, sigma1, epsilon1
                    )
                    # TODO: Calculate exceptions


class ChargeIncrementModelHandler(_NonbondedHandler):
    """Handle SMIRNOFF ``<ChargeIncrementModel>`` tags

    .. warning :: This API is experimental and subject to change.
    """

    class ChargeIncrementType(ParameterType):
        """A SMIRNOFF bond charge correction type.

        .. warning :: This API is experimental and subject to change.
        """

        _VALENCE_TYPE = None  # This disables the connectivity check when parsing LibraryChargeType objects
        _ELEMENT_NAME = "ChargeIncrement"

        charge_increment = IndexedParameterAttribute(unit=unit.elementary_charge)

        def __init__(self, **kwargs):
            super().__init__(**kwargs)
            unique_tags, connectivity = GLOBAL_TOOLKIT_REGISTRY.call(
                "get_tagged_smarts_connectivity", self.smirks
            )

            n_tags = len(unique_tags)
            n_increments = len(self.charge_increment)
            diff = n_tags - n_increments

            if diff < 0 or diff > 1:
                # TODO: Consider dealing with diff > 2 by smearing charges across
                # all un-specified increments
                raise SMIRNOFFSpecError(
                    f"ChargeIncrement {self} was initialized with an invalid combination "
                    f"of tagged atoms and charge increments"
                )

    _TAGNAME = "ChargeIncrementModel"  # SMIRNOFF tag name to process
    _INFOTYPE = ChargeIncrementType  # info type to store
    _DEPENDENCIES = [
        vdWHandler,
        ElectrostaticsHandler,
        LibraryChargeHandler,
        ToolkitAM1BCCHandler,
    ]
    _MAX_SUPPORTED_SECTION_VERSION = 0.4

    number_of_conformers = ParameterAttribute(default=1, converter=int)

    partial_charge_method = ParameterAttribute(default="AM1-Mulliken", converter=str)

    def check_handler_compatibility(
        self, other_handler, assume_missing_is_default=True
    ):
        """
        Checks whether this ParameterHandler encodes compatible physics as another ParameterHandler. This is
        called if a second handler is attempted to be initialized for the same tag.

        Parameters
        ----------
        other_handler : a ParameterHandler object
            The handler to compare to.

        Raises
        ------
        IncompatibleParameterError if handler_kwargs are incompatible with existing parameters.
        """

        int_attrs_to_compare = ["number_of_conformers"]
        string_attrs_to_compare = ["partial_charge_method"]

        self._check_attributes_are_equal(
            other_handler,
            identical_attrs=string_attrs_to_compare + int_attrs_to_compare,
        )

    def find_matches(self, entity, unique=False):
        """Find the elements of the topology/molecule matched by a parameter type.

        Parameters
        ----------
        entity : openff.toolkit.topology.Topology
            Topology to search.

        Returns
        ---------
        matches : ValenceDict[Tuple[int], ParameterHandler._Match]
            ``matches[particle_indices]`` is the ``ParameterType`` object
            matching the tuple of particle indices in ``entity``.
        """
        matches = self._find_matches(
            entity, transformed_dict_cls=TagSortedDict, unique=unique
        )
        return matches

    def create_force(self, system, topology, **kwargs):
        import warnings

        from openff.toolkit.topology import TopologyAtom, TopologyVirtualSite

        # We only want one instance of this force type
        existing = [system.getForce(i) for i in range(system.getNumForces())]
        existing = [f for f in existing if type(f) == self._OPENMMTYPE]
        if len(existing) == 0:
            force = self._OPENMMTYPE()
            system.addForce(force)
        else:
            force = existing[0]

        for ref_mol in topology.reference_molecules:

            # If charges were already assigned, skip this molecule
            if self.check_charges_assigned(ref_mol, topology):
                continue

            toolkit_registry = kwargs.get("toolkit_registry", GLOBAL_TOOLKIT_REGISTRY)
            try:
                # If the molecule wasn't assigned parameters from a manually-input charge_mol, calculate them here
                ref_mol.generate_conformers(n_conformers=self.number_of_conformers)
                ref_mol.assign_partial_charges(
                    partial_charge_method=self.partial_charge_method,
                    toolkit_registry=toolkit_registry,
                )
            except Exception as e:
                warnings.warn(str(e), Warning)
                continue

            charges_to_assign = {}

            # Assign initial, un-incremented charges to relevant atoms
            for topology_molecule in topology._reference_molecule_to_topology_molecules[
                ref_mol
            ]:
                for topology_particle in topology_molecule.particles:
                    topology_particle_index = topology_particle.topology_particle_index
                    if type(topology_particle) is TopologyAtom:
                        ref_mol_particle_index = (
                            topology_particle.atom.molecule_particle_index
                        )
                    if type(topology_particle) is TopologyVirtualSite:
                        ref_mol_particle_index = (
                            topology_particle.virtual_site.molecule_particle_index
                        )
                    particle_charge = ref_mol._partial_charges[ref_mol_particle_index]
                    charges_to_assign[topology_particle_index] = particle_charge

            # Find SMARTS-based matches for charge increments
            charge_increment_matches = self.find_matches(topology)

            # We ignore the atom index order in the keys here, since they have been
            # sorted in order to deduplicate matches and let us identify when one parameter overwrites another
            # in the SMIRNOFF parameter hierarchy. Since they are sorted, the position of the atom index
            # in the key tuple DOES NOT correspond to the appropriate charge_incrementX value.
            # Instead, the correct ordering of the match indices is found in
            # charge_increment_match.environment_match.topology_atom_indices
            for (_, charge_increment_match) in charge_increment_matches.items():
                # Adjust the values in the charges_to_assign dict by adding any
                # charge increments onto the existing values
                atom_indices = (
                    charge_increment_match.environment_match.topology_atom_indices
                )
                charge_increments = copy.deepcopy(
                    charge_increment_match.parameter_type.charge_increment
                )

                # If we've been provided with one less charge increment value than tagged atoms, assume the last
                # tagged atom offsets the charge of the others to make the chargeincrement net-neutral
                if len(atom_indices) - len(charge_increments) == 1:
                    charge_increment_sum = 0.0 * unit.elementary_charge
                    for ci in charge_increments:
                        charge_increment_sum += ci
                    charge_increments.append(-charge_increment_sum)
                elif len(atom_indices) - len(charge_increments) == 0:
                    pass
                else:
                    raise SMIRNOFFSpecError(
                        f"Trying to apply chargeincrements {charge_increment_match.parameter_type} "
                        f"to tagged atoms {atom_indices}, but the number of chargeincrements "
                        f"must be either the same as- or one less than the number of tagged atoms."
                    )

                for top_particle_idx, charge_increment in zip(
                    atom_indices, charge_increments
                ):
                    if top_particle_idx in charges_to_assign:
                        charges_to_assign[top_particle_idx] += charge_increment

            # Set the incremented charges on the System particles
            for topology_particle_index, charge_to_assign in charges_to_assign.items():
                _, sigma, epsilon = force.getParticleParameters(topology_particle_index)
                force.setParticleParameters(
                    topology_particle_index, charge_to_assign, sigma, epsilon
                )

            # Finally, mark that charges were assigned for this reference molecule
            self.mark_charges_assigned(ref_mol, topology)


class GBSAHandler(ParameterHandler):
    """Handle SMIRNOFF ``<GBSA>`` tags

    .. warning :: This API is experimental and subject to change.
    """

    class GBSAType(ParameterType):
        """A SMIRNOFF GBSA type.

        .. warning :: This API is experimental and subject to change.
        """

        _VALENCE_TYPE = "Atom"
        _ELEMENT_NAME = "Atom"

        radius = ParameterAttribute(unit=unit.angstrom)
        scale = ParameterAttribute(converter=float)

    _TAGNAME = "GBSA"
    _INFOTYPE = GBSAType
    _OPENMMTYPE = openmm.GBSAOBCForce
    # It's important that this runs AFTER partial charges are assigned to all particles, since this will need to
    # collect and assign them to the GBSA particles
    _DEPENDENCIES = [
        vdWHandler,
        ElectrostaticsHandler,
        ToolkitAM1BCCHandler,
        ChargeIncrementModelHandler,
        LibraryChargeHandler,
    ]

    gb_model = ParameterAttribute(
        default="OBC1", converter=_allow_only(["HCT", "OBC1", "OBC2"])
    )
    solvent_dielectric = ParameterAttribute(default=78.5, converter=float)
    solute_dielectric = ParameterAttribute(default=1, converter=float)
    sa_model = ParameterAttribute(default="ACE", converter=_allow_only(["ACE", None]))
    surface_area_penalty = ParameterAttribute(
        default=5.4 * unit.calorie / unit.mole / unit.angstrom ** 2,
        unit=unit.calorie / unit.mole / unit.angstrom ** 2,
    )
    solvent_radius = ParameterAttribute(default=1.4 * unit.angstrom, unit=unit.angstrom)

    def _validate_parameters(self):
        """
        Checks internal attributes, raising an exception if they are configured in an invalid way.
        """
        # If we're using HCT via GBSAHCTForce(CustomAmberGBForceBase):, then we need to ensure that:
        #   surface_area_energy is 5.4 cal/mol/A^2
        #   solvent_radius is 1.4 A
        # Justification at https://github.com/openforcefield/openff-toolkit/pull/363
        if self.gb_model == "HCT":
            if (
                self.surface_area_penalty
                != 5.4 * unit.calorie / unit.mole / unit.angstrom ** 2
            ) and (self.sa_model is not None):
                raise IncompatibleParameterError(
                    f"The current implementation of HCT GBSA does not "
                    f"support surface_area_penalty values other than 5.4 "
                    f"cal/mol A^2 (data source specified value of "
                    f"{self.surface_area_penalty})"
                )

            if (self.solvent_radius != 1.4 * unit.angstrom) and (
                self.sa_model is not None
            ):
                raise IncompatibleParameterError(
                    f"The current implementation of HCT GBSA does not "
                    f"support solvent_radius values other than 1.4 "
                    f"A (data source specified value of "
                    f"{self.solvent_radius})"
                )

        # If we're using OBC1 via GBSAOBC1Force(CustomAmberGBForceBase), then we need to ensure that:
        #   surface_area_energy is 5.4 cal/mol/A^2
        #   solvent_radius is 1.4 A
        # Justification at https://github.com/openforcefield/openff-toolkit/pull/363
        if self.gb_model == "OBC1":
            if (
                self.surface_area_penalty
                != 5.4 * unit.calorie / unit.mole / unit.angstrom ** 2
            ) and (self.sa_model is not None):
                raise IncompatibleParameterError(
                    f"The current implementation of OBC1 GBSA does not "
                    f"support surface_area_penalty values other than 5.4 "
                    f"cal/mol A^2 (data source specified value of "
                    f"{self.surface_area_penalty})"
                )

            if (self.solvent_radius != 1.4 * unit.angstrom) and (
                self.sa_model is not None
            ):
                raise IncompatibleParameterError(
                    f"The current implementation of OBC1 GBSA does not "
                    f"support solvent_radius values other than 1.4 "
                    f"A (data source specified value of "
                    f"{self.solvent_radius})"
                )

        # If we're using OBC2 via GBSAOBCForce, then we need to ensure that
        #   solvent_radius is 1.4 A
        # Justification at https://github.com/openforcefield/openff-toolkit/pull/363
        if self.gb_model == "OBC2":

            if (self.solvent_radius != 1.4 * unit.angstrom) and (
                self.sa_model is not None
            ):
                raise IncompatibleParameterError(
                    f"The current implementation of OBC1 GBSA does not "
                    f"support solvent_radius values other than 1.4 "
                    f"A (data source specified value of "
                    f"{self.solvent_radius})"
                )

    # Tolerance when comparing float attributes for handler compatibility.
    _SCALETOL = 1e-5

    def check_handler_compatibility(self, other_handler):
        """
        Checks whether this ParameterHandler encodes compatible physics as another ParameterHandler. This is
        called if a second handler is attempted to be initialized for the same tag.

        Parameters
        ----------
        other_handler : a ParameterHandler object
            The handler to compare to.

        Raises
        ------
        IncompatibleParameterError if handler_kwargs are incompatible with existing parameters.
        """
        float_attrs_to_compare = ["solvent_dielectric", "solute_dielectric"]
        string_attrs_to_compare = ["gb_model", "sa_model"]
        unit_attrs_to_compare = ["surface_area_penalty", "solvent_radius"]

        self._check_attributes_are_equal(
            other_handler,
            identical_attrs=string_attrs_to_compare,
            tolerance_attrs=float_attrs_to_compare + unit_attrs_to_compare,
            tolerance=self._SCALETOL,
        )

    def create_force(self, system, topology, **kwargs):
        import simtk

        self._validate_parameters()

        # Grab the existing nonbonded force (which will have particle charges)
        existing = [system.getForce(i) for i in range(system.getNumForces())]
        existing = [f for f in existing if type(f) == openmm.NonbondedForce]
        assert len(existing) == 1

        nonbonded_force = existing[0]

        # No previous GBSAForce should exist, so we're safe just making one here.
        force_map = {
            "HCT": simtk.openmm.app.internal.customgbforces.GBSAHCTForce,
            "OBC1": simtk.openmm.app.internal.customgbforces.GBSAOBC1Force,
            "OBC2": simtk.openmm.GBSAOBCForce,
            # It's tempting to do use the class below, but the customgbforce
            # version of OBC2 doesn't provide setSolventRadius()
            # 'OBC2': simtk.openmm.app.internal.customgbforces.GBSAOBC2Force,
        }
        openmm_force_type = force_map[self.gb_model]

        if nonbonded_force.getNonbondedMethod() == openmm.NonbondedForce.NoCutoff:
            amber_cutoff = None
        else:
            amber_cutoff = nonbonded_force.getCutoffDistance().value_in_unit(
                unit.nanometer
            )

        if self.gb_model == "OBC2":
            gbsa_force = openmm_force_type()

        else:
            # We set these values in the constructor if we use the internal AMBER GBSA type wrapper
            gbsa_force = openmm_force_type(
                solventDielectric=self.solvent_dielectric,
                soluteDielectric=self.solute_dielectric,
                SA=self.sa_model,
                cutoff=amber_cutoff,
                kappa=0,
            )
            # WARNING: If using a CustomAmberGBForce, the functional form is affected by whether
            # the cutoff kwarg is None *during initialization*. So, if you initialize it with a
            # non-None value, and then try to change it to None, you're likely to get unphysical results.

        # Set the GBSAForce to have the same cutoff as NonbondedForce
        # gbsa_force.setCutoffDistance(nonbonded_force.getCutoffDistance())
        if amber_cutoff is not None:
            gbsa_force.setCutoffDistance(amber_cutoff)

        if nonbonded_force.usesPeriodicBoundaryConditions():
            # WARNING: The lines below aren't equivalent. The NonbondedForce and
            # CustomGBForce NonbondedMethod enums have different meanings.
            # More details:
            # http://docs.openmm.org/latest/api-python/generated/openmm.openmm.NonbondedForce.html
            # http://docs.openmm.org/latest/api-python/generated/openmm.openmm.GBSAOBCForce.html
            # http://docs.openmm.org/latest/api-python/generated/openmm.openmm.CustomGBForce.html

            # gbsa_force.setNonbondedMethod(simtk.openmm.NonbondedForce.CutoffPeriodic)
            gbsa_force.setNonbondedMethod(simtk.openmm.CustomGBForce.CutoffPeriodic)
        else:
            # gbsa_force.setNonbondedMethod(simtk.openmm.NonbondedForce.NoCutoff)
            gbsa_force.setNonbondedMethod(simtk.openmm.CustomGBForce.NoCutoff)

        # Add all GBSA terms to the system. Note that this will have been done above
        if self.gb_model == "OBC2":
            gbsa_force.setSolventDielectric(self.solvent_dielectric)
            gbsa_force.setSoluteDielectric(self.solute_dielectric)
            if self.sa_model is None:
                gbsa_force.setSurfaceAreaEnergy(0)
            else:
                gbsa_force.setSurfaceAreaEnergy(self.surface_area_penalty)

        # Iterate over all defined GBSA types, allowing later matches to override earlier ones.
        atom_matches = self.find_matches(topology)

        # Create all particles.

        # !!! WARNING: CustomAmberGBForceBase expects different per-particle parameters
        # depending on whether you use addParticle or setParticleParameters. In
        # setParticleParameters, we have to apply the offset and scale BEFORE setting
        # parameters, whereas in addParticle, the offset is applied automatically, and the particle
        # parameters are not set until an auxillary finalize() method is called. !!!

        # To keep it simple, we DO NOT pre-populate the particles in the GBSA force here.
        # We call addParticle further below instead.
        # These lines are commented out intentionally as an example of what NOT to do.
        # for topology_particle in topology.topology_particles:
        # gbsa_force.addParticle([0.0, 1.0, 0.0])

        params_to_add = [[] for _ in topology.topology_particles]
        for atom_key, atom_match in atom_matches.items():
            atom_idx = atom_key[0]
            gbsatype = atom_match.parameter_type
            charge, _, _2 = nonbonded_force.getParticleParameters(atom_idx)
            params_to_add[atom_idx] = [charge, gbsatype.radius, gbsatype.scale]

        if self.gb_model == "OBC2":
            for particle_param in params_to_add:
                gbsa_force.addParticle(*particle_param)
        else:
            for particle_param in params_to_add:
                gbsa_force.addParticle(particle_param)
            # We have to call finalize() for models that inherit from CustomAmberGBForceBase,
            # otherwise the added particles aren't actually passed to the underlying CustomGBForce
            gbsa_force.finalize()

        # Check that no atoms (n.b. not particles) are missing force parameters.
        self._check_all_valence_terms_assigned(
            assigned_terms=atom_matches, valence_terms=list(topology.topology_atoms)
        )

        system.addForce(gbsa_force)


class VirtualSiteHandler(_NonbondedHandler):
    """Handle SMIRNOFF ``<VirtualSites>`` tags

    TODO: Add example usage/documentation

    .. warning :: This API is experimental and subject to change.
    """

    # Virtual Site exclusions policies
    ############################################################################
    # none: do not add any exclusions

    # minimal: only add exclusions between vsite particles and their "single"
    # parent atom. This is the atom that the vsite's origin is defined as

    # parents: only add exclusions between vsite particles and all of the
    # associated parent atoms

    # local: add exclusions between vsites that share exactly the same atoms.

    # neighbors: add exclusions between vsites and atoms that share the same
    # "clique" of virtual sites. For example, if 1-2-3 and 3-4-5 each have a
    # vsite, then the vsite on 1-2-3 will be excluded from atoms 4 and 5
    # since they share atom 3.

    # connected: add exclusions between the vsite and all atoms connected to
    # the parents, e.g the entire molecule making it an interaction only
    # between two nonbonded fragments.

    # all: exclude all interactions, effectively turning vsites off.

    # Note: TODO: only up to parents is implemented!

    class _ExclusionPolicy(Enum):
        NONE = 1
        MINIMAL = 2
        PARENTS = 3
        LOCAL = 4
        NEIGHBORS = 5
        CONNECTED = 6
        ALL = 7

    _parameter_to_policy = {
        "none": _ExclusionPolicy.NONE,
        "minimal": _ExclusionPolicy.MINIMAL,
        "parents": _ExclusionPolicy.PARENTS,
        "local": _ExclusionPolicy.LOCAL,
        "neighbors": _ExclusionPolicy.NEIGHBORS,
        "connected": _ExclusionPolicy.CONNECTED,
        "all": _ExclusionPolicy.ALL,
    }

    exclusion_policy = ParameterAttribute(default="parents")  # has custom converter

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self._virtual_site_types = dict()
        for vsite_cls in all_subclasses(self.__class__.VirtualSiteType):
            # catch classes which are not actual implementations, which should return None
            vtype = vsite_cls.vsite_type()
            if vtype:
                self.register_virtual_site_type(vtype, vsite_cls)

    def add_parameter(
        self, parameter_kwargs=None, parameter=None, after=None, before=None
    ):
        """Add a parameter to the force field, ensuring all parameters are valid.
        This method differs from other handlers in that it uses a plugin-style
        enable/disable type system


        Parameters
        ----------
        parameter_kwargs: dict, optional
            The kwargs to pass to the ParameterHandler.INFOTYPE (a ParameterType) constructor
        parameter: ParameterType, optional
            A ParameterType to add to the ParameterHandler
        after : str or int, optional
            The SMIRKS pattern (if str) or index (if int) of the parameter directly before where
            the new parameter will be added
        before : str, optional
            The SMIRKS pattern (if str) or index (if int) of the parameter directly after where
            the new parameter will be added

        Note that one of (parameter_kwargs, parameter) must be specified
        Note that when `before` and `after` are both None, the new parameter will be appended
            to the END of the parameter list.
        Note that when `before` and `after` are both specified, the new parameter
            will be added immediately after the parameter matching the `after` pattern or index.

        """

        # TODO: This function need unit tests

        for val in [before, after]:
            if val and not isinstance(val, (str, int)):
                raise TypeError

        # If a dict was passed, construct it; if a ParameterType was passed, do nothing
        if parameter_kwargs:
            vsite_type = parameter_kwargs["type"]
            if (
                vsite_type in self._virtual_site_types
                and self._virtual_site_types[vsite_type] is not None
            ):
                new_parameter = self._virtual_site_types[vsite_type](**parameter_kwargs)
            else:
                raise ValueError(
                    f"Virtual site type {vsite_type} not enabled or implemented in this handler {self.__class__}"
                )
        elif parameter:
            new_parameter = parameter
            # As a convenience, if parameter type not present, register it
            if parameter.type not in self._virtual_site_types:
                self.register_virtual_site_type(parameter.type, type(parameter))
            # additionally, if the type was previously disabled (set to None),
            # reenable it with this new type
            elif self._virtual_site_types.get(parameter.type, None) is None:
                self.register_virtual_site_type(
                    parameter.type, type(parameter), replace=True
                )

        else:
            raise ValueError("One of (parameter, parameter_kwargs) must be specified")

        if (
            (new_parameter.smirks in [p.smirks for p in self._parameters])
            and (new_parameter.type in [p.type for p in self._parameters])
            and (new_parameter.name in [p.name for p in self._parameters])
        ):
            msg = f"A parameter SMIRKS pattern {new_parameter.smirks} already exists for type {new_parameter.type} and name {new_parameter.name}."
            raise DuplicateParameterError(msg)

        if before is not None:
            if isinstance(before, str):
                before_index = self._parameters.index(before)
            elif isinstance(before, int):
                before_index = before

        if after is not None:
            if isinstance(after, str):
                after_index = self._parameters.index(after)
            elif isinstance(after, int):
                after_index = after

        if None not in (before, after):
            if after_index > before_index:
                raise ValueError("before arg must be before after arg")

        if after is not None:
            self._parameters.insert(after_index + 1, new_parameter)
        elif before is not None:
            self._parameters.insert(before_index, new_parameter)
        else:
            self._parameters.append(new_parameter)

    def _add_parameters(self, section_dict, allow_cosmetic_attributes=False):
        """
        Extend the ParameterList in this VirtualSiteHandler using a SMIRNOFF data source.

        Parameters
        ----------
        section_dict : dict
            The dict representation of a SMIRNOFF data source containing parameters to att to this VirtualSiteHandler
        allow_cosmetic_attributes : bool, optional. Default = False
            Whether to allow non-spec fields in section_dict. If True, non-spec kwargs will be stored as an
            attribute of the parameter. If False, non-spec kwargs will raise an exception.

        """

        # Most of this is exactly the same as the base _add_parameters. The only
        # difference is how INFOTYPE is implemented, see the comment below

        unitless_kwargs, attached_units = extract_serialized_units_from_dict(
            section_dict
        )
        smirnoff_data = attach_units(unitless_kwargs, attached_units)

        for key, val in smirnoff_data.items():
            if self._INFOTYPE is not None:
                element_name = self._INFOTYPE._ELEMENT_NAME
                # Skip sections that aren't the parameter list
                if key != element_name:
                    break
            # If there are multiple parameters, this will be a list. If there's just one, make it a list
            if not (isinstance(val, list)):
                val = [val]

            # If we're reading the parameter list, iterate through and attach units to
            # each parameter_dict, then use it to initialize a ParameterType
            for unitless_param_dict in val:

                param_dict = attach_units(unitless_param_dict, attached_units)

                # This differs from other handlers in that we use both a
                # dynamic version of INFOTYPE, and also allow a plugin-style
                # system where we allow visibility of virtual site types to
                # control which ones are allowed to be activated
                vsite_type = param_dict["type"]
                if (
                    vsite_type in self._virtual_site_types
                    and self._virtual_site_types[vsite_type] is not None
                ):
                    new_parameter = self._virtual_site_types[vsite_type](
                        **param_dict,
                        allow_cosmetic_attributes=allow_cosmetic_attributes,
                    )
                    self._parameters.append(new_parameter)
                else:
                    raise ValueError(
                        f"Virtual site type {vsite_type} not enabled or implemented in this handler {self.__class__}"
                    )

    def register_virtual_site_type(self, vsite_name, vsite_cls, replace=False):
        """
        Register an implemented virtual site type. Doing this must be done to
        pass the validation and option checking. To disable a type, register the
        name with None and replace=True

        Parameters
        ----------
        vsite_name : str
            The name of the type. This name must be what is found in the "type"
            attribute in the OFFXML format

        vsite_cls : Any
            The class to register the name with that implements the type.

        Returns
        -------
        None
        """

        if vsite_name in self._virtual_site_types and not replace:
            raise DuplicateVirtualSiteTypeException(
                "VirtualSite type {} already registered for handler {} and replace=False".format(
                    vsite_name, self.__class__
                )
            )
        self._virtual_site_types[vsite_name] = vsite_cls
        self._parameters = ParameterList(
            [param for param in self._parameters if param.type != vsite_name]
        )

    @property
    def virtual_site_types(self):
        """
        Return the dictionary of registered virtual site types

        Parameters
        ----------
        None

        Returns
        -------
        virtual_site_types : dict
            A list of virtual site types already registered, paired with their
            class that implements them
        """

        return self._virtual_site_types

    @exclusion_policy.converter
    def exclusion_policy(self, attr, policy):
        """
        Convert and validate the exclusion policy specified in the VirtualSiteHandler

        Parameters
        ----------
        attr : openff.toolkit.typing.engines.smirnoff.parameters.ParameterAttribute
            The underlying ParameterAttribute
        policy : Any
            The policy name to validate

        Returns
        -------
        policy : str
            The policy name if it is valid

        Raises
        ------
        SMIRNOFFSpecError if the value of policy did not match the SMIRNOFF Specification
        ValueError if policy cannot be converted to a string

        .. warning :: This API is experimental and subject to change.
        """

        try:
            policy = str(policy)
        except ValueError:
            raise

        _exclusion_policies_implemented = ["none", "minimal", "parents"]
        if policy in _exclusion_policies_implemented:
            return policy
        else:

            raise SMIRNOFFSpecError(
                'VirtualSiteHandler exclusion policy not understood. Set "exclusion_policy" to one of {}'.format(
                    _exclusion_policies_implemented
                )
            )

    class _VirtualSiteTypeSelector:
        """A SMIRNOFF virtual site base selector

        This is a placeholder class that will dynamically choose the correct
        virtual site to create based on the type specified in attributes.
        Normally, the OFFXML element explicity defines the type, but here it
        depends on the type attribute as well, which needs the introspection
        implemented here.

        """

        _VALENCE_TYPE = None
        # This is needed so VirtualSite elements are parsed correctly
        # using this generic selector as a type
        _ELEMENT_NAME = "VirtualSite"

        _enable_types = {}

        def __new__(cls, **attrs):

            VSH = VirtualSiteHandler

            vsite_type = attrs["type"]

            if vsite_type == "BondCharge":
                cls = VSH.VirtualSiteBondChargeType

            elif vsite_type == "MonovalentLonePair":
                cls = VSH.VirtualSiteMonovalentLonePairType

            elif vsite_type == "DivalentLonePair":
                cls = VSH.VirtualSiteDivalentLonePairType

            elif vsite_type == "TrivalentLonePair":
                cls = VSH.VirtualSiteTrivalentLonePairType
            else:
                raise SMIRNOFFSpecError(
                    'VirtualSite type not understood. Choose one of "BondCharge", "MonovalentLonePair", "DivalentLonePair", "TrivalentLonePair"'
                )

            return cls(**attrs)

    class VirtualSiteType(vdWHandler.vdWType, abc.ABC):
        """A SMIRNOFF virtual site base type

        .. warning :: This API is experimental and subject to change.
        """

        # The attributes that we expect in the OFFXML
        name = ParameterAttribute(default="EP", converter=str)
        distance = ParameterAttribute(unit=unit.angstrom)
        charge_increment = IndexedParameterAttribute(unit=unit.elementary_charge)
        # Type has a delayed converter/validator to support a plugin-style enable/disable system
        type = ParameterAttribute(converter=str)
        match = ParameterAttribute(default="all_permutations")  # has custom converter
        epsilon = ParameterAttribute(
            default=0.0 * unit.kilocalorie_per_mole, unit=unit.kilocalorie_per_mole
        )
        sigma = ParameterAttribute(default=0.0 * unit.angstrom, unit=unit.angstrom)
        rmin_half = ParameterAttribute(default=None, unit=unit.angstrom)

        # Here we define the default sorting behavior if we need to sort the
        # atom key into a canonical ordering
        transformed_dict_cls = ValenceDict

        # Value of None indicates "not a valid type" or "not an actual implemented type".
        # To enable/register a new virtual site type, make a subclass and set its
        # _vsite_type to what would need to be provided in the OFFXML "type" attr,
        # e.g. type="BondCharge" would mean _vsite_type="BondCharge"
        _vsite_type = None

        @classmethod
        def vsite_type(cls):
            """
            The type of this virtual site as represented in the SMIRNOFF specification

            .. warning :: This API is experimental and subject to change.
            """
            return cls._vsite_type

        def __init__(self, **kwargs):
            """
            Create a virtual site parameter type
            """

            # Need to create default vdW parameters if not specified, since they are optional
            sigma = kwargs.get("sigma", None)
            rmin_half = kwargs.get("rmin_half", None)

            if (sigma is None) and (rmin_half is None):
                kwargs["sigma"] = 0.0 * unit.angstrom

            if kwargs.get("epsilon", None) is None:
                kwargs["epsilon"] = 0.0 * unit.kilocalorie_per_mole

            super().__init__(**kwargs)

            if sigma:
                self._extra_nb_var = "rmin_half"
            if rmin_half:
                self._extra_nb_var = "sigma"

        # @type.converter
        # def type(self, attr, vsite_type):
        #     """
        #     Convert and validate the virtual site type specified in the VirtualSite element

        #     Parameters
        #     ----------
        #     attr : openff.toolkit.typing.engines.smirnoff.parameters.ParameterAttribute
        #         The underlying ParameterAttribute
        #     vsite_type : Any
        #         The virtual site type to validate

        #     Returns
        #     -------
        #     vsite_type : str
        #         The virtual site type if it is valid

        #     Raises
        #     ------
        #     SMIRNOFFSpecError if the value of policy did not match the SMIRNOFF Specification
        #     ValueError if policy cannot be converted to a string

        #     .. warning :: This API is experimental and subject to change.
        #     """

        #     try:
        #         vsite_type = str(vsite_type)
        #     except ValueError:
        #         raise

        #     if vsite_type in VirtualSiteHandler.virtual_site_types():
        #         return vsite_type
        #     else:
        #         valid_types = ", ".join(
        #             [str('"' + vtype + '"') for vtype in self._virtual_site_types]
        #         )
        #         raise SMIRNOFFSpecError(
        #             "VirtualSite not given a type. Set type to one of:\n" + valid_types
        #         )

        @match.converter
        def match(self, attr, match):
            """
            Convert and validate the virtual site type specified in the VirtualSite element

            Parameters
            ----------
            match : Any
                The virtual site type to validate

            Returns
            -------
            match : str
                The virtual site type if it is valid

            Raises
            ------
            SMIRNOFFSpecError
            ValueError

            .. warning :: This API is experimental and subject to change.
            """

            try:
                match = str(match)
            except ValueError:
                raise

            if match == "once" or match == "all_permutations":
                return match
            else:
                raise SMIRNOFFSpecError(
                    'VirtualSite type must specify "match" as either "once" or "all_permutations"'
                )

        def __eq__(self, obj):
            if type(self) != type(obj):
                return False
            A = ["name"]
            are_equal = [getattr(self, a) == getattr(obj, a) for a in A]
            return all(are_equal)

        def _add_virtual_site(self, fn, atoms, orientations, *args, **kwargs):
            """

            Parameters
            ----------
            fn : callable
                The underlying OpenFF function that should be called to create
                the virtual site in the toolkit. Currently, these are:
                    * `Molecule._add_bond_charge_virtual_site`
                    * `Molecule._add_monovalent_lone_pair_virtual_site`
                    * `Molecule._add_divalent_lone_pair_virtual_site`
                    * `Molecule._add_trivalent_lone_pair_virtual_site`

            orientations : list of int tuples
                The permutations corresponding to each virtual particle in the
                virtual site.
            Returns
            -------
                The index of the created virtual site
            """

            args = [atoms, self.distance] + list(args)

            # This needs to be dealt with better
            # Since we cannot save state with this object during find_matches,
            # we have no idea here which permutations actually matched.
            # Since we are past the point where we can examine chemical
            # environment matches to determine the possible orientations, we
            # must default and take the explicit interpretation:
            # "all_permutations" will try to make a virtual particle for every
            # permutation, and "once" is the canonical sorting of the atom
            # indices.

            # This means that, using the "match" option in the spec, it is not
            # possible to choose specific permutations. For the current cases,
            # this should be fine and works well.

            # The API above takes all given orientations, but the OFFXML
            # has the match setting, which ultimately decides which orientations
            # to include.
            if self.match == "once":
                key = self.transformed_dict_cls.key_transform(orientations[0])
                orientations = [key]
                # else all matches wanted, so keep whatever was matched.

            base_args = {
                "name": self.name,
                "charge_increments": self.charge_increment,
                "epsilon": self.epsilon,
                "sigma": self.sigma,
                "rmin_half": self.rmin_half,
                "orientations": orientations,
                "replace": kwargs.pop("replace", False),
            }
            kwargs.update(base_args)
            kwargs.pop(self._extra_nb_var)

            return fn(*args, **kwargs)

    class VirtualSiteBondChargeType(VirtualSiteType):
        """A SMIRNOFF virtual site bond charge type

        .. warning :: This API is experimental and subject to change.
        """

        _vsite_type = "BondCharge"

        def add_virtual_site(self, molecule, orientations, replace=False):
            """
            Add a virtual site to the molecule

            Parameters
            ----------
            molecule : openff.toolkit.topology.molecule.Molecule
                The molecule to add the virtual site to
            orientations : List[Tuple[int]]
                A list of orientation tuples which define the permuations used
                to contruct the geometry of the virtual site particles
            replace : bool, default=False
                Replace this virtual site if it already exists in the molecule

            Returns
            -------
            off_idx : int
                The index of the first particle added due to this virtual site

            .. warning :: This API is experimental and subject to change.
            """

            fn = molecule._add_bond_charge_virtual_site
            ref_key = self.transformed_dict_cls.key_transform(orientations[0])
            atoms = list([molecule.atoms[i] for i in ref_key])
            args = (atoms, orientations)
            off_idx = super()._add_virtual_site(fn, *args, replace=replace)
            return off_idx

    class VirtualSiteMonovalentLonePairType(VirtualSiteType):
        """A SMIRNOFF monovalent lone pair virtual site type

        .. warning :: This API is experimental and subject to change.
        """

        outOfPlaneAngle = ParameterAttribute(unit=unit.degree)
        inPlaneAngle = ParameterAttribute(unit=unit.degree)

        _vsite_type = "MonovalentLonePair"

        def add_virtual_site(self, molecule, orientations, replace=False):
            """
            Add a virtual site to the molecule

            Parameters
            ----------
            molecule : openff.toolkit.topology.molecule.Molecule
                The molecule to add the virtual site to
            orientations : List[Tuple[int]]
                A list of orientation tuples which define the permuations used
                to contruct the geometry of the virtual site particles
            replace : bool, default=False
                Replace this virtual site if it already exists in the molecule

            Returns
            -------
            off_idx : int
                The index of the first particle added due to this virtual site

            .. warning :: This API is experimental and subject to change.
            """

            fn = molecule._add_monovalent_lone_pair_virtual_site
            ref_key = self.transformed_dict_cls.key_transform(orientations[0])
            atoms = list([molecule.atoms[i] for i in ref_key])
            args = (atoms, orientations, self.outOfPlaneAngle, self.inPlaneAngle)
            off_idx = super()._add_virtual_site(fn, *args, replace=replace)
            return off_idx

    class VirtualSiteDivalentLonePairType(VirtualSiteType):
        """A SMIRNOFF divalent lone pair virtual site type

        .. warning :: This API is experimental and subject to change.
        """

        outOfPlaneAngle = ParameterAttribute(unit=unit.degree)

        _vsite_type = "DivalentLonePair"

        def add_virtual_site(self, molecule, orientations, replace=False):
            """
            Add a virtual site to the molecule

            Parameters
            ----------
            molecule : openff.toolkit.topology.molecule.Molecule
                The molecule to add the virtual site to
            orientations : List[Tuple[int]]
                A list of orientation tuples which define the permuations used
                to contruct the geometry of the virtual site particles
            replace : bool, default=False
                Replace this virtual site if it already exists in the molecule

            Returns
            -------
            off_idx : int
                The index of the first particle added due to this virtual site

            .. warning :: This API is experimental and subject to change.
            """
            fn = molecule._add_divalent_lone_pair_virtual_site
            ref_key = self.transformed_dict_cls.key_transform(orientations[0])
            atoms = list([molecule.atoms[i] for i in ref_key])
            args = (atoms, orientations, self.outOfPlaneAngle)
            off_idx = super()._add_virtual_site(fn, *args, replace=replace)
            return off_idx

    class VirtualSiteTrivalentLonePairType(VirtualSiteType):
        """A SMIRNOFF trivalent lone pair virtual site type

        .. warning :: This API is experimental and subject to change.
        """

        transformed_dict_cls = ImproperDict

        _vsite_type = "TrivalentLonePair"

        def __init__(self, **kwargs):
            """
            Special init method for TrivalentLonePairSites that ensures that match="all_permutations"
            """
            super().__init__(**kwargs)
            if self.match != "once":
                raise SMIRNOFFSpecError(
                    f"TrivalentLonePair virtual site defined with match attribute set to {self.match}. "
                    f"Only supported value is 'once'."
                )

        def add_virtual_site(self, molecule, orientations, replace=False):
            """
            Add a virtual site to the molecule

            Parameters
            ----------
            molecule : openff.toolkit.topology.molecule.Molecule
                The molecule to add the virtual site to
            orientations : List[Tuple[int]]
                A list of orientation tuples which define the permuations used
                to contruct the geometry of the virtual site particles
            replace : bool, default=False
                Replace this virtual site if it already exists in the molecule

            Returns
            -------
            off_idx : int
                The index of the first particle added due to this virtual site

            .. warning :: This API is experimental and subject to change.
            """
            fn = molecule._add_trivalent_lone_pair_virtual_site
            ref_key = self.transformed_dict_cls.key_transform(orientations[0])
            atoms = list([molecule.atoms[i] for i in ref_key])

            # Trivalents should never need multiple orientations as long
            # as there are no angle parameters
            args = (atoms, orientations)
            off_idx = super()._add_virtual_site(fn, *args, replace=replace)
            return off_idx

    _DEPENDENCIES = [
        ElectrostaticsHandler,
        LibraryChargeHandler,
        ChargeIncrementModelHandler,
        ToolkitAM1BCCHandler,
    ]

    _TAGNAME = "VirtualSites"  # SMIRNOFF tag name to process

    # Trying to create an instance of this selector will cause
    # some introspection to be done on the type attr passed in, and
    # will dispatch the appropriate virtual site type.
    _INFOTYPE = _VirtualSiteTypeSelector  # class to hold force type info

    def _find_matches(
        self,
        entity,
        transformed_dict_cls=UnsortedDict,
        use_named_slots=False,
        expand_permutations=False,
    ):
        """Implement find_matches() and allow using a difference valence dictionary.

        Parameters
        ----------
        entity : openff.toolkit.topology.Topology
            Topology to search.
        transformed_dict_cls: Union[Dict, ValenceDict, ImproperDict]
            The type of dictionary to store the matches in. This
            will determine how groups of atom indices are stored
            and accessed (e.g for angles indices should be 0-1-2
            and not 2-1-0).

        Returns
        -------
        matches : `transformed_dict_cls` of ParameterHandlerMatch
            ``matches[particle_indices]`` is the ``ParameterType`` object
            matching the tuple of particle indices in ``entity``.
        """
        from collections import defaultdict

        logger.debug("Finding matches for {}".format(self.__class__.__name__))

        matches = transformed_dict_cls()

        for parameter_type in self._parameters:

            matches_for_this_type = defaultdict(list)

            ce_matches = entity.chemical_environment_matches(parameter_type.smirks)

            # Split the groups into unique sets i.e. 13,14 and 13,15
            # Needed for vsites, where a vsite could match C-H with for a CH2 group
            # Since these are distinct vsite definitions, we need to split them
            # up into separate groups (match_groups)
            match_groups_set = [m.topology_atom_indices for m in ce_matches]
            match_groups = []
            for key in set(match_groups_set):
                distinct_atom_pairs = [
                    x
                    for x in ce_matches
                    if sorted(x.topology_atom_indices) == sorted(key)
                ]
                match_groups.append(distinct_atom_pairs)

            for ce_matches in match_groups:
                for environment_match in ce_matches:
                    # Update the matches for this parameter type.
                    handler_match = self._Match(parameter_type, environment_match)
                    key = environment_match.topology_atom_indices

                    # only a match if orientation matches
                    if not hasattr(handler_match._parameter_type, "match"):
                        # Probably should never get here
                        raise SMIRNOFFSpecError(
                            "The match keyword not found in this parameter?!"
                        )
                    else:

                        # The possible orders of this match
                        # We must check that the tuple of atoms are the same
                        # as they can be different in e.g. formaldehyde
                        orders = [m.topology_atom_indices for m in ce_matches]
                        orientation_flag = handler_match._parameter_type.match

                        tdc = handler_match._parameter_type.transformed_dict_cls
                        index_of_key = tdc.index_of(key, possible=orders)

                        if orientation_flag == "once":
                            orientation = [0]
                        elif orientation_flag == "all_permutations":
                            orientation = [
                                tdc.index_of(k, possible=orders) for k in orders
                            ]
                        else:
                            # Probably will never reach here since validation
                            # happens elsewhere
                            raise Exception(
                                "VirtualSite match keyword not understood. Choose from 'once' or 'all_permutations'. This error should be impossible to reach; please submit an issue at https://github.com/openforcefield/openff-toolkit"
                            )

                        orders = [
                            order for order in orders if sorted(key) == sorted(order)
                        ]

                        # Find these matches is from the older implementation that allows
                        # specifying specific orientations. Leaving in for now..
                        if len(orientation) > len(orders):
                            error_msg = (
                                "For parameter of type\n{:s}\norientations {} "
                                + "exceeds length of possible orders "
                                + "({:d}):\n{:s}"
                            ).format(
                                str(parameter_type),
                                orientation,
                                len(orders),
                                str(orders),
                            )
                            raise IndexError(error_msg)

                        if not expand_permutations:
                            key = tdc.key_transform(key)

                        hit = sum([index_of_key == ornt for ornt in orientation])
                        assert (
                            hit < 2
                        ), "VirtualSite orientation for {:s} indices invalid: Has duplicates".format(
                            parameter_type.__repr__
                        )
                        if hit == 1:
                            matches_for_this_type[key].append(handler_match)

                # Resolve match overriding by the use of the name attribute
                # If two parameters match but have the same name, use most recent,
                # but if the names are different, keep and apply both parameters
                if use_named_slots:
                    for k in matches_for_this_type:
                        if k not in matches:
                            matches[k] = {}
                    for k, v in matches_for_this_type.items():
                        marginal_matches = []
                        for new_match in v:
                            unique = True
                            new_item = new_match._parameter_type
                            for idx, (name, existing_match) in enumerate(
                                matches[k].items()
                            ):
                                existing_item = existing_match._parameter_type
                                same_parameter = False

                                same_type = type(existing_item) == type(new_item)
                                if same_type:
                                    same_parameter = existing_item == new_item

                                # same, so replace it to have a FIFO priority
                                # and the last parameter matching wins
                                if same_parameter:
                                    matches[k][new_item.name] = new_match
                                    unique = False
                            if unique:
                                marginal_matches.append(new_match)
                        matches[k].update(
                            {p._parameter_type.name: p for p in marginal_matches}
                        )
                else:
                    matches.update(matches_for_this_type)

            logger.debug(
                "{:64} : {:8} matches".format(
                    parameter_type.smirks, len(matches_for_this_type)
                )
            )

        logger.debug("{} matches identified".format(len(matches)))

        if use_named_slots:
            for k, v in matches.items():
                matches[k] = list(v.values())

        return matches

    def create_force(self, system, topology, **kwargs):
        """

        Parameters
        ----------

        Returns
        -------
        """
        force = super().create_force(system, topology, **kwargs)

        # Separate the logic of adding vsites in the oFF state and the OpenMM
        # system. Operating on the topology is not ideal (a hack), so hopefully
        # this loop, which adds the oFF vsites to the topology, will live
        # somewhere else

        logger.debug("Creating OpenFF virtual site representations...")
        self.create_openff_virtual_sites(topology)

        # The toolkit now has a representation of the vsites in the topology,
        # and here we create the OpenMM parameters/objects/exclusions
        logger.debug("Creating OpenMM VSite particles...")
        for ref_mol in topology.reference_molecules:
            logger.debug("Adding vsites for reference mol: {}".format(str(ref_mol)))
            self._create_openmm_virtual_sites(system, force, topology, ref_mol)

    def check_handler_compatibility(self, other_handler):
        """
        Checks whether this ParameterHandler encodes compatible physics as
        another ParameterHandler. This is called if a second handler is
        attempted to be initialized for the same tag.

        Parameters
        ----------
        other_handler : a ParameterHandler object
            The handler to compare to.

        Raises
        ------
        IncompatibleParameterError if handler_kwargs are incompatible with
        existing parameters.
        """
        string_attrs_to_compare = [
            "exclusion_policy",
        ]
        self._check_attributes_are_equal(
            other_handler, identical_attrs=string_attrs_to_compare
        )

    def find_matches(self, entity, expand_permutations=True, unique=False):
        """Find the virtual sites in the topology/molecule matched by a
        parameter type.

        Parameters
        ----------
        entity : openff.toolkit.topology.Topology
            Topology to search.

        Returns
        ---------
        matches : Dict[Tuple[int], ParameterHandler._Match]
            ``matches[atom_indices]`` is the ``ParameterType`` object
            matching the n-tuple of atom indices in ``entity``.

        """
        if unique:
            raise NotImplementedError(
                "`unique=True` not implemented in VirtualSiteHandler"
            )
        return self._find_matches(
            entity,
            transformed_dict_cls=UnsortedDict,
            use_named_slots=True,
            expand_permutations=expand_permutations,
        )

    def _apply_charge_increment(self, force, atom_key, charge_increment):
        vsite_charge = charge_increment[0]
        vsite_charge *= 0.0
        for charge_i, atom in enumerate(atom_key):
            o_charge, o_sigma, o_epsilon = force.getParticleParameters(atom)
            vsite_charge -= charge_increment[charge_i]
            o_charge += charge_increment[charge_i]
            force.setParticleParameters(atom, o_charge, o_sigma, o_epsilon)
        return vsite_charge

    def _same_virtual_site_type(self, vs_i, vs_j):
        if type(vs_i) != type(vs_j):
            return False
        if vs_i.name != vs_j.name:
            return False
        return True

    def _reduce_virtual_particles_to_sites(self, atom_matches):
        combined_orientations = []

        # These are the indices representing the tuples (VSITE_TYPE, KEY_LIST).
        # Basically an ordered dictionary with ints as keys
        VSITE_TYPE = 0
        KEY_LIST = 1

        for key, atom_match_lst in atom_matches.items():
            for match in atom_match_lst:

                # For each match, loop through existing virtual sites found,
                # and determine if this match is a unique virtual site,
                # or a member of an existing virtual site (e.g. TIP5)
                vs_i = match.parameter_type
                found = False

                for i, vsite_struct in enumerate(combined_orientations):

                    vs_j = vsite_struct[VSITE_TYPE]

                    # The logic to determine if the particles should be
                    # combined into a single virtual site. If the atoms
                    # are the same, the vsite is the same, but the absolute
                    # ordering of the match is different, then we say
                    # this is part of the same virtual site.
                    # Note that the expand_permutations=True above is what
                    # returns the different orders for each match (normally,
                    # this is not the case for e.g. bonds where 1-2 is the
                    # same parameter as 2-1 and is always returned as 1-2.
                    same_atoms = all(
                        [sorted(key) == sorted(k) for k in vsite_struct[KEY_LIST]]
                    )
                    diff_keys = key not in vsite_struct[KEY_LIST]
                    same_vsite = self._same_virtual_site_type(vs_i, vs_j)

                    if same_atoms and same_vsite and diff_keys:
                        combined_orientations[i][KEY_LIST].append(key)
                        found = True

                    # Skip out early since there is no reason to keep
                    # searching since we will never add the same
                    # particle twice
                    if found:
                        break

                # If the entire loop did not produce a hit, then this is
                # a brand new virtual site
                if not found:
                    newsite = [None, None]
                    newsite[VSITE_TYPE] = vs_i
                    newsite[KEY_LIST] = [key]
                    combined_orientations.append(newsite)

        return combined_orientations

    def create_openff_virtual_sites(self, topology):
        """
        Modifies the input topology to contain VirtualSites assigned by this handler.

        Parameters
        ----------
        topology : openff.toolkit.topology.Topology
            Topology to add virtual sites to.
        """

        for molecule in topology.reference_molecules:

            """The following two lines below should be avoided but is left
            until a better solution is found (see #699). The issue is that a
            topology should not be required since `find_matches` works on
            FrozenMolecules. However, the signature is different, as they return
            different results.

            Also, we are using a topology to retrieve the indices for the
            matches, but then using those indices as a direct `Atom` object
            lookup in the molecule. This is unsafe because there is no reason to
            believe that the indices should be consistent. However, there is
            currently no `FrozenMolecule.atom(index)` functionality so using the
            topology match indices is the only clear way forward.  See the
            contents of `add_virtual_site` called below for the code that shows
            this."""

            top_mol = Topology.from_molecules([molecule])
            matches = self.find_matches(top_mol, expand_permutations=True)

            virtual_sites = self._reduce_virtual_particles_to_sites(matches)

            # Now handle the vsites for this molecule
            # This call batches the key tuples into a single list, in order
            # for the virtual site to represent multiple particles
            for vsite_type, orientations in virtual_sites:
                vsite_type.add_virtual_site(molecule, orientations, replace=True)

    def _create_openmm_virtual_sites(self, system, force, topology, ref_mol):

        """
        Here we must assume that
            1. All atoms in the topology are already present
            2. The order we iterate these virtual sites is the order they
                appear in the OpenMM system
        If 1 is not met, then 2 will fail, and it will be quite difficult to
        find the mapping since we currently do not keep track of any OpenMM
        state, and it is unlikely that we will ever do so.  If 1 is met, then 2
        should fall into place naturally.

        This means that we will not check that our index matches the OpenMM
        index, as there is no reason, from a purely technical and/or API
        standpoint, to require them to be.
        """

        for vsite in ref_mol.virtual_sites:
            ref_key = [atom.molecule_atom_index for atom in vsite.atoms]
            logger.debug("VSite ref_key: {}".format(ref_key))

            ms = topology._reference_molecule_to_topology_molecules[ref_mol]
            for top_mol in ms:
                logger.debug("top_mol: {}".format(top_mol))

                ids = self._create_openmm_virtual_particle(
                    system, force, top_mol, vsite, ref_key
                )

                # Go and exclude each of the vsite particles; this makes
                # sense because these particles cannot "feel" forces, only
                # exert them
                policy = self._parameter_to_policy[self.exclusion_policy]
                if policy.value != self._ExclusionPolicy.NONE.value:
                    # Default here is to always exclude vsites which are
                    # of the same virtual site. Their positions are rigid,
                    # and so any energy that would be added to the system
                    # due to their pairwise interaction would not make sense.
                    for i, j in combinations(ids, 2):
                        logger.debug("Excluding vsite {} vsite {}".format(i, j))
                        force.addException(i, j, 0.0, 0.0, 0.0, replace=True)

    def _create_openmm_virtual_particle(self, system, force, top_mol, vsite, ref_key):

        policy = self._parameter_to_policy[self.exclusion_policy]

        ids = []

        for vp in vsite.particles:
            orientation = vp.orientation
            sort_key = [orientation.index(i) for i in ref_key]
            atom_key = [ref_key[i] for i in sort_key]
            logger.debug("sort_key: {}".format(sort_key))
            atom_key = [top_mol.atom_start_topology_index + i for i in atom_key]

            omm_vsite = vsite.get_openmm_virtual_site(atom_key)
            vsite_q = self._apply_charge_increment(
                force, atom_key, vsite.charge_increments
            )

            ljtype = vsite
            if ljtype.sigma is None:
                sigma = 2.0 * ljtype.rmin_half / (2.0 ** (1.0 / 6.0))
            else:
                sigma = ljtype.sigma

            # create the vsite particle
            mass = 0.0
            vsite_idx = system.addParticle(mass)
            ids.append(vsite_idx)

            logger.debug(
                "vsite_id: {} orientation: {} atom_key: {}".format(
                    vsite_idx, orientation, atom_key
                )
            )

            system.setVirtualSite(vsite_idx, omm_vsite)
            force.addParticle(vsite_q, sigma, ljtype.epsilon)

            logger.debug(f"Added virtual site particle with charge {vsite_q}")
            logger.debug(f"  charge_increments: {vsite.charge_increments}")

            # add exclusion to the "parent" atom of the vsite
            if policy.value >= self._ExclusionPolicy.MINIMAL.value:
                keylen = len(atom_key)
                if keylen == 2:
                    owning_atom = atom_key[0]
                elif keylen == 3:
                    owning_atom = atom_key[1]
                else:
                    # The second atom of an improper is considered the
                    # "owning" atom
                    owning_atom = atom_key[1]
                logger.debug(f"Excluding vsite {vsite_idx} atom {owning_atom}")
                force.addException(owning_atom, vsite_idx, 0.0, 0.0, 0.0, replace=True)

            # add exclusions to all atoms in the vsite definition (the parents)
            if policy.value >= self._ExclusionPolicy.PARENTS.value:
                for i in atom_key:
                    if i == owning_atom:
                        continue
                    logger.debug(f"Excluding vsite {vsite_idx} atom {i}")
                    force.addException(i, vsite_idx, 0.0, 0.0, 0.0, replace=True)

            if policy.value > self._ExclusionPolicy.PARENTS.value:
                raise NotImplementedError(
                    "Only the 'parents', 'minimal', and 'none' exclusion_policies are implemented"
                )

        return ids


if __name__ == "__main__":
    import doctest

    doctest.testmod()
    # doctest.run_docstring_examples(_ParameterAttributeHandler, globals())
