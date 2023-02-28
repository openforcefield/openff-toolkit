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
    "NotEnoughPointsForInterpolationError",
    "ParameterLookupError",
    "SMIRNOFFSpecError",
    "SMIRNOFFSpecUnimplementedError",
    "UnassignedAngleParameterException",
    "UnassignedBondParameterException",
    "UnassignedMoleculeChargeException",
    "UnassignedProperTorsionParameterException",
    "UnassignedValenceParameterException",
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
    "ParameterType",
    "ConstraintType",
    "BondType",
    "AngleType",
    "ProperTorsionType",
    "ImproperTorsionType",
    "vdWType",
    "LibraryChargeType",
    "GBSAType",
    "ChargeIncrementType",
    "VirtualSiteType",
]
import copy
import functools
import inspect
import logging
import re
from collections import defaultdict
from typing import (
    Any,
    Callable,
    Dict,
    List,
    Literal,
    Optional,
    Tuple,
    Union,
    cast,
    get_args,
)

import numpy as np
from openff.units import unit
from packaging.version import Version

from openff.toolkit.topology import ImproperDict, TagSortedDict, Topology, ValenceDict
from openff.toolkit.topology.molecule import Molecule
from openff.toolkit.utils.collections import ValidatedDict, ValidatedList
from openff.toolkit.utils.exceptions import (
    DuplicateParameterError,
    DuplicateVirtualSiteTypeException,
    FractionalBondOrderInterpolationMethodUnsupportedError,
    IncompatibleParameterError,
    IncompatibleUnitError,
    MissingIndexedAttributeError,
    MissingPartialChargesError,
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
from openff.toolkit.utils.utils import object_to_quantity

logger = logging.getLogger(__name__)


_cal_mol_a2 = unit.calorie / unit.mole / unit.angstrom**2


def _linear_inter_or_extrapolate(points_dict, x_query):
    """
    Linearly interpolate or extrapolate based on a piecewise linear function
    defined by a set of points. This function is designed to work with
    key:value pairs where the value may be a ``openff.units.Quantity``.

    Parameters
    ----------
    points_dict : dict{float: float or float-valued openff.units.Quantity}
        A dictionary with each item representing a point, where the key is the X value and the value is the Y value.
    x_query : float
        The X value of the point to interpolate or extrapolate.

    Returns
    -------
    y_value : float or float-valued openff.units.Quantity
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
        k = points_dict[bond_orders[0]] - (
            (points_dict[bond_orders[1]] - points_dict[bond_orders[0]])
            / (bond_orders[1] - bond_orders[0])
        ) * (bond_orders[0] - x_query)
        return k

    # extrapolate for fractional bond orders above our highest defined bond order
    elif above is None:
        bond_orders = sorted(points_dict)
        k = points_dict[bond_orders[-1]] + (
            (points_dict[bond_orders[-1]] - points_dict[bond_orders[-2]])
            / (bond_orders[-1] - bond_orders[-2])
        ) * (x_query - bond_orders[-1])
        return k


# TODO: This is technically a validator, not a converter, but ParameterAttribute doesn't support them yet
#       (it'll be easy if we switch to use the attrs library).
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


def _validate_units(attr, value: Union[str, unit.Quantity], units: unit.Unit):
    value = object_to_quantity(value)

    try:
        if not units.is_compatible_with(value.units):
            raise IncompatibleUnitError(
                f"{attr.name}={value} should have units of {units}"
            )
    except AttributeError:
        raise IncompatibleUnitError(f"{attr.name}={value} should have units of {units}")
    return value


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
    unit : openff.units.Quantity, optional
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

    >>> from openff.units import unit
    >>> class MyParameter:
    ...     attr_quantity = ParameterAttribute(unit=unit.angstrom)
    ...
    >>> my_par = MyParameter()
    >>> my_par.attr_quantity = '1.0 * nanometer'
    >>> my_par.attr_quantity
    <Quantity(1.0, 'nanometer')>
    >>> my_par.attr_quantity = 3.0
    Traceback (most recent call last):
    ...
    openff.toolkit.utils.exceptions.IncompatibleUnitError:
    attr_quantity=3.0 dimensionless should have units of angstrom

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

    def __init__(
        self,
        default: Any = UNDEFINED,
        unit: Optional[unit.Unit] = None,
        converter: Optional[Callable] = None,
        docstring: str = "",
    ):
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
                if not self._unit.is_compatible_with(value.units):
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
    unit : openff.units.Quantity, optional
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

    >>> from openff.units import unit
    >>> class MyParameter:
    ...     length = IndexedParameterAttribute(default=None, unit=unit.angstrom)
    ...
    >>> my_par = MyParameter()
    >>> my_par.length is None
    True

    Strings are parsed into Quantity objects.

    >>> my_par.length = ['1 * angstrom', 0.5 * unit.nanometer]
    >>> my_par.length[0]
    <Quantity(1, 'angstrom')>

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
    unit : openff.units.Quantity, optional
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

    >>> from openff.units import unit
    >>> class MyParameter:
    ...     length = MappedParameterAttribute(default=None, unit=unit.angstrom)
    ...
    >>> my_par = MyParameter()
    >>> my_par.length is None
    True

    Like other ParameterAttribute objects, strings are parsed into Quantity objects.

    >>> my_par.length = {1:'1.5 * angstrom', 2: '1.4 * angstrom'}
    >>> my_par.length[1]
    <Quantity(1.5, 'angstrom')>

    Unlike other ParameterAttribute objects, the reference points can do not need ot be
    zero-indexed, non-adjancent, such as interpolating defining a bond parameter for
    interpolation by defining references values and bond orders 2 and 3:

    >>> my_par.length = {2:'1.42 * angstrom', 3: '1.35 * angstrom'}
    >>> my_par.length[2]
    <Quantity(1.42, 'angstrom')>

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
    unit : openff.units.Quantity, optional
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

    >>> from openff.units import unit
    >>> class MyParameter:
    ...     length = IndexedMappedParameterAttribute(default=None, unit=unit.angstrom)
    ...
    >>> my_par = MyParameter()
    >>> my_par.length is None
    True

    Strings are parsed into Quantity objects.

    >>> my_par.length = [{1:'1 * angstrom'}, {1: 0.5 * unit.nanometer}]
    >>> my_par.length[0]
    {1: <Quantity(1, 'angstrom')>}

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
    ...     k = ParameterAttribute(unit=unit.kilocalorie / unit.mole / unit.angstrom**2)
    ...

    ``_ParameterAttributeHandler`` and the descriptors take care of performing
    sanity checks on initialization and assignment of the single attributes. Because
    we attached units to the parameters, we need to pass them with compatible units.

    >>> my_par = ParameterTypeOrHandler(
    ...     length='1.01 * angstrom',
    ...     k=5 * unit.kilocalorie / unit.mole / unit.angstrom**2
    ... )

    Note that ``_ParameterAttributeHandler`` took care of implementing
    a constructor, and that unit parameters support string assignments.
    These are automatically converted to ``Quantity`` objects.

    >>> my_par.length
    <Quantity(1.01, 'angstrom')>

    While assigning incompatible units is forbidden.

    >>> my_par.k = 3.0 * unit.gram
    Traceback (most recent call last):
    ...
    openff.toolkit.utils.exceptions.IncompatibleUnitError:
    k=3.0 gram should have units of kilocalorie / angstrom ** 2 / mole

    On top of type checking, the constructor implemented in ``_ParameterAttributeHandler``
    checks if some required parameters are not given.

    >>> ParameterTypeOrHandler(length=3.0*unit.nanometer)
    Traceback (most recent call last):
    ...
    openff.toolkit.utils.exceptions.SMIRNOFFSpecError:
    <class '...ParameterTypeOrHandler'> require the following missing
    parameters: ['k']. Defined kwargs are ['length']

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
    ...     k = IndexedParameterAttribute(unit=unit.kilocalorie / unit.mole)
    ...
    >>> my_par = MyTorsionType(
    ...     periodicity1=2,
    ...     k1=5 * unit.kilocalorie / unit.mole,
    ...     periodicity2='3',
    ...     k2=6 * unit.kilocalorie / unit.mole,
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
                attrib_w_index = f"{attrib_basename}{index}"

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
                try:
                    attribs_to_return.pop(attribs_to_return.index(duplicate))
                except ValueError:
                    # The attribute was not in the list
                    continue

        # Start populating a dict of the attribs.
        indexed_attribs = set(self._get_indexed_parameter_attributes().keys())
        mapped_attribs = set(self._get_mapped_parameter_attributes().keys())
        indexed_mapped_attribs = set(
            self._get_indexed_mapped_parameter_attributes().keys()
        )
        smirnoff_dict = dict()

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
            elif attrib_name == "version":
                smirnoff_dict[attrib_name] = str(attrib_value)
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
            key is not None
            and index is not None
            and attr_name in self._get_indexed_mapped_parameter_attributes()
        ):
            indexed_mapped_attr_value = getattr(self, attr_name)
            try:
                return indexed_mapped_attr_value[index][key]
            except (IndexError, KeyError) as err:
                raise MissingIndexedAttributeError(
                    f"{str(err)} '{item}' is out of bounds for indexed attribute '{attr_name}'"
                )

        # Otherwise, try indexed attribute
        # Separate the indexed attribute name from the list index.
        attr_name, index = self._split_attribute_index(item)

        # Check if this is an indexed attribute.
        if index is not None and attr_name in self._get_indexed_parameter_attributes():
            indexed_attr_value = getattr(self, attr_name)
            try:
                return indexed_attr_value[index]
            except IndexError:
                raise MissingIndexedAttributeError(
                    f"'{item}' is out of bounds for indexed attribute '{attr_name}'"
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
                raise MissingIndexedAttributeError(
                    f"{str(err)} '{key}' is out of bounds for indexed attribute '{attr_name}'"
                )

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
                raise MissingIndexedAttributeError(
                    f"'{key}' is out of bounds for indexed attribute '{attr_name}'"
                )

        # Forward the request to the next class in the MRO.
        super().__setattr__(key, value)

    def add_cosmetic_attribute(self, attr_name, attr_value):
        """
        Add a cosmetic attribute to this object.

        This attribute will not have a functional effect on the object
        in the OpenFF Toolkit, but can be written out during
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

            def filter(x):
                return True

        # Go through MRO and retrieve also parents descriptors. The function
        # inspect.getmembers() automatically resolves the MRO, but it also
        # sorts the attribute alphabetically by name. Here we want the order
        # to be the same as the declaration order, which is guaranteed by PEP 520,
        # starting from the parent class.
        parameter_attributes = dict(
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
        optional = dict(
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
                f"but received {other} (type {type(other)}) instead"
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
            raise ParameterLookupError(f"SMIRKS {item} not found in ParameterList")

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
    attributes. In the example below, ``_ELEMENT_NAME`` is used to
    describe the SMIRNOFF parameter being defined, and is used during
    automatic serialization/deserialization into a ``dict``.

    >>> class MyBondParameter(ParameterType):
    ...     _ELEMENT_NAME = 'Bond'
    ...     length = ParameterAttribute(unit=unit.angstrom)
    ...     k = ParameterAttribute(unit=unit.kilocalorie / unit.mole / unit.angstrom**2)
    ...

    The parameter automatically inherits the required smirks attribute
    from ``ParameterType``. Associating a ``unit`` to a ``ParameterAttribute``
    cause the attribute to accept only values in compatible units and to
    parse string expressions.

    >>> my_par = MyBondParameter(
    ...     smirks='[*:1]-[*:2]',
    ...     length='1.01 * angstrom',
    ...     k=5 * unit.kilocalorie / unit.mole / unit.angstrom**2
    ... )
    >>> my_par.length
    <Quantity(1.01, 'angstrom')>
    >>> my_par.k = 3.0 * unit.gram
    Traceback (most recent call last):
    ...
    openff.toolkit.utils.exceptions.IncompatibleUnitError:
    k=3.0 gram should have units of kilocalorie / angstrom ** 2 / mole

    Each attribute can be made optional by specifying a default value,
    and you can attach a converter function by passing a callable as an
    argument or through the decorator syntax.

    >>> class MyParameterType(ParameterType):
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
    ...     _ELEMENT_NAME = 'Proper'
    ...     periodicity = IndexedParameterAttribute(converter=int)
    ...     k = IndexedParameterAttribute(unit=unit.kilocalorie / unit.mole)
    ...
    >>> my_par = MyTorsionType(
    ...     smirks='[*:1]-[*:2]-[*:3]-[*:4]',
    ...     periodicity1=2,
    ...     k1=5 * unit.kilocalorie / unit.mole,
    ...     periodicity2='3',
    ...     k2=6 * unit.kilocalorie / unit.mole,
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

    # The string mapping to this ParameterType in a SMIRNOFF data source
    _ELEMENT_NAME: Optional[str] = None

    # Parameter attributes shared among all parameter types.
    smirks = ParameterAttribute()
    id = ParameterAttribute(default=None)
    parent_id = ParameterAttribute(default=None)

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
        ret_str = f"<{self.__class__.__name__} with "
        for attr, val in self.to_dict().items():
            ret_str += f"{attr}: {val}  "
        ret_str += ">"
        return ret_str


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

    # str of section type handled by this ParameterHandler (XML element name for SMIRNOFF XML representation)
    _TAGNAME: Optional[str] = None
    # container class with type information that will be stored in self._parameters
    _INFOTYPE: Optional[Any] = None
    # OpenMM Force class (or None if no equivalent)
    _OPENMMTYPE: Optional[str] = None
    # list of ParameterHandler classes that must precede this, or None
    _DEPENDENCIES: Optional[Any] = None

    # Kwargs to catch when create_force is called
    _KWARGS: List[str] = []
    # the earliest version of SMIRNOFF spec that supports this ParameterHandler
    _SMIRNOFF_VERSION_INTRODUCED = 0.0
    _SMIRNOFF_VERSION_DEPRECATED = None
    # if deprecated, the first SMIRNOFF version number it is no longer used
    _MIN_SUPPORTED_SECTION_VERSION = Version("0.3")
    _MAX_SUPPORTED_SECTION_VERSION = Version("0.3")

    version = ParameterAttribute()

    @version.converter
    def version(self, attr, new_version):
        """
        Raise a parsing exception if the given section version is unsupported.

        Raises
        ------
        SMIRNOFFVersionError if an incompatible version is passed in.

        """
        from openff.toolkit.utils.exceptions import SMIRNOFFVersionError

        if isinstance(new_version, Version):
            pass
        elif isinstance(new_version, str):
            new_version = Version(new_version)
        elif isinstance(new_version, (float, int)):
            new_version = Version(str(new_version))
        else:
            raise Exception(f"Could not convert type {type(new_version)}")

        # Use PEP-440 compliant version number comparison, if requested
        if (new_version > self._MAX_SUPPORTED_SECTION_VERSION) or (
            new_version < self._MIN_SUPPORTED_SECTION_VERSION
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
        for key, val in section_dict.items():
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
            for param_dict in val:
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

    def _index_of_parameter(
        self, parameter: Optional[ParameterType] = None, key: Optional[Any] = None
    ) -> Optional[int]:
        """Attempts to find the index of a parameter in the parameters list.

        By default, two parameters are considered 'the same' if they have the same
        SMIRKS pattern.

        Parameters
        ----------
        parameter
            The parameter to find the index of. This argument is mutually exclusive with
            ``key``.
        key
            The SMIRKS pattern associated with the parameter to find the index
            of. This argument is mutually exclusive with ``parameter``.

        Returns
        -------
            The index of the parameter if found, otherwise ``None``.
        """

        if (key is None and parameter is None) or (
            key is not None and parameter is not None
        ):
            raise ValueError("`key` and `parameter` are mutually exclusive arguments")

        key = key if parameter is None else parameter.smirks

        for index, existing_parameter in enumerate(self._parameters):
            if existing_parameter.smirks != key:
                continue

            return index

        return None

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
          * The order of parameters in a parameter list can have significant impacts on parameter assignment. For
            details, see the SMIRNOFF specification:
            https://openforcefield.github.io/standards/standards/smirnoff/#smirnoff-parameter-specification-is-hierarchical

        Examples
        --------

        Add a ParameterType to an existing ParameterList at a specified position.

        Given an existing parameter handler and a new parameter to add to it:

        >>> from openff.units import unit
        >>> bh = BondHandler(skip_version_check=True)
        >>> length = 1.5 * unit.angstrom
        >>> k = 100 * unit.kilocalorie / unit.mole / unit.angstrom ** 2
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

        if self._index_of_parameter(new_parameter) is not None:
            msg = f"A parameter SMIRKS pattern {new_parameter.smirks} already exists."
            raise DuplicateParameterError(msg)

        before_index, after_index = None, None

        if before is not None:
            if isinstance(before, int):
                before_index = before
            else:
                before_index = self._index_of_parameter(key=before)

        if after is not None:
            if isinstance(after, int):
                after_index = after
            else:
                after_index = self._index_of_parameter(key=after)

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

        >>> from openff.units import unit
        >>> handler = BondHandler(skip_version_check=True)
        >>> handler.add_parameter(
        ...     {
        ...         'smirks': '[*:1]-[*:2]',
        ...         'length': 1*unit.angstrom,
        ...         'k': 10*unit.kilocalorie / unit.mole/unit.angstrom**2,
        ...     }
        ... )

        Look up, from this handler, all parameters matching some SMIRKS pattern

        >>> handler.get_parameter({'smirks': '[*:1]-[*:2]'})
        [<BondType with smirks: [*:1]-[*:2]  length: 1 angstrom  k: 10.0 kilocalorie / angstrom ** 2 / mole  >]

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
            ``matches[atom_indices]`` is the ``ParameterType`` object
            matching the tuple of atom indices in ``entity``.
        """

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
            ``matches[atom_indices]`` is the ``ParameterType`` object
            matching the tuple of atom indices in ``entity``.
        """
        logger.debug(f"Finding matches for {self.__class__.__name__}")

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

        logger.debug(f"{len(matches)} matches identified")
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

    def create_force(self, *args, **kwarsg):
        """
        .. deprecated:: 0.11.0

            This method was deprecated in v0.11.0, no longer has any
            functionality, and will soon be removed. Use the `OpenFF Interchange
            <https://docs.openforcefield.org/interchange>`_ package instead.
        """
        raise NotImplementedError(
            "`ParameterHandler`s no longer create OpenMM forces. Use `openff-interchange` instead."
        )

    def to_dict(self, discard_cosmetic_attributes=False):
        """
        Convert this ParameterHandler to a dict, compliant with the SMIRNOFF data spec.

        Parameters
        ----------
        discard_cosmetic_attributes : bool, optional. Default = False.
            Whether to discard non-spec parameter and header attributes in this ParameterHandler.

        Returns
        -------
        smirnoff_data : dict
            SMIRNOFF-spec compliant representation of this ParameterHandler and its internal ParameterList.

        """
        smirnoff_data = dict()

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
                this_val.units
            except AttributeError:
                return this_val, other_val
            assert this_val.units == other_val.units
            return this_val.m, other_val.m

        for attr in identical_attrs:
            this_val, other_val = get_unitless_values(attr)

            if this_val != other_val:
                msg = (
                    f"{attr} values are not identical. (handler value: '{this_val}', "
                    f"incompatible value: '{other_val}').\n"
                )
                if "AM1-Wiberg" in (this_val, other_val) and "none" in (
                    this_val,
                    other_val,
                ):
                    msg += (
                        "This likely results from mixing bond handlers with different versions "
                        "(0.3 and 0.4). Consider upgrading bond handlers to version 0.4 or "
                        "manually setting `fractional_bondorder_method='AM1-Wiberg'`."
                    )
                raise IncompatibleParameterError(msg)

        for attr in tolerance_attrs:
            try:
                this_val, other_val = get_unitless_values(attr)
            except AttributeError:
                raise AttributeError(
                    f"Mismatch found with attr={attr}, this_val={this_val}, "
                    f"other_val={other_val}"
                )
            if abs(this_val - other_val) > tolerance:
                raise IncompatibleParameterError(
                    "Difference between '{}' values is beyond allowed tolerance {}. "
                    "(handler value: {}, incompatible value: {}".format(
                        attr, tolerance, this_val, other_val
                    )
                )

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

        _ELEMENT_NAME = "Constraint"

        distance = ParameterAttribute(default=None, unit=unit.angstrom)

    _TAGNAME = "Constraints"
    _INFOTYPE = ConstraintType
    _OPENMMTYPE = None  # don't create a corresponding OpenMM Force class


class BondHandler(ParameterHandler):
    """Handle SMIRNOFF ``<Bonds>`` tags

    .. warning :: This API is experimental and subject to change.
    """

    class BondType(ParameterType):
        """A SMIRNOFF bond type

        .. warning :: This API is experimental and subject to change.
        """

        _ELEMENT_NAME = "Bond"

        length = ParameterAttribute(default=None, unit=unit.angstrom)
        k = ParameterAttribute(
            default=None, unit=unit.kilocalorie / unit.mole / unit.angstrom**2
        )

        # fractional bond order params
        length_bondorder = MappedParameterAttribute(default=None, unit=unit.angstrom)
        k_bondorder = MappedParameterAttribute(
            default=None, unit=unit.kilocalorie / unit.mole / unit.angstrom**2
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
    _OPENMMTYPE = "HarmonicBondForce"
    _DEPENDENCIES = [ConstraintHandler]  # ConstraintHandler must be executed first
    _MAX_SUPPORTED_SECTION_VERSION = Version("0.4")

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
        if (
            self.version == Version("0.3")
            and "fractional_bondorder_method" not in kwargs
        ):
            self.fractional_bondorder_method = "none"
        elif (
            self.version == Version("0.4")
            and "fractional_bondorder_method" not in kwargs
        ):
            self.fractional_bondorder_method = "AM1-Wiberg"

        # Default value for potential depends on section version
        if self.version == Version("0.3") and "potential" not in kwargs:
            self.potential = "harmonic"
        elif self.version == Version("0.4") and "potential" not in kwargs:
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


class AngleHandler(ParameterHandler):
    """Handle SMIRNOFF ``<AngleForce>`` tags

    .. warning :: This API is experimental and subject to change.
    """

    class AngleType(ParameterType):
        """A SMIRNOFF angle type.

        .. warning :: This API is experimental and subject to change.
        """

        _ELEMENT_NAME = "Angle"

        angle = ParameterAttribute(unit=unit.degree)
        k = ParameterAttribute(unit=unit.kilocalorie / unit.mole / unit.degree**2)

    _TAGNAME = "Angles"  # SMIRNOFF tag name to process
    _INFOTYPE = AngleType  # class to hold force type info
    _OPENMMTYPE = "HarmonicAngleForce"
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


# TODO: There's a lot of duplicated code in ProperTorsionHandler and ImproperTorsionHandler
class ProperTorsionHandler(ParameterHandler):
    """Handle SMIRNOFF ``<ProperTorsionForce>`` tags

    .. warning :: This API is experimental and subject to change.
    """

    class ProperTorsionType(ParameterType):
        """A SMIRNOFF torsion type for proper torsions.

        .. warning :: This API is experimental and subject to change.
        """

        _ELEMENT_NAME = "Proper"

        periodicity = IndexedParameterAttribute(converter=int)
        phase = IndexedParameterAttribute(unit=unit.degree)
        k = IndexedParameterAttribute(default=None, unit=unit.kilocalorie / unit.mole)
        idivf = IndexedParameterAttribute(default=None, converter=float)

        # fractional bond order params
        k_bondorder = IndexedMappedParameterAttribute(
            default=None, unit=unit.kilocalorie / unit.mole
        )

    _TAGNAME = "ProperTorsions"  # SMIRNOFF tag name to process
    _KWARGS = ["partial_bond_orders_from_molecules"]
    _INFOTYPE = ProperTorsionType  # info type to store
    _OPENMMTYPE = "PeriodicTorsionForce"
    _MAX_SUPPORTED_SECTION_VERSION = Version("0.4")

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


# TODO: There's a lot of duplicated code in ProperTorsionHandler and ImproperTorsionHandler
class ImproperTorsionHandler(ParameterHandler):
    """Handle SMIRNOFF ``<ImproperTorsionForce>`` tags

    .. warning :: This API is experimental and subject to change.
    """

    class ImproperTorsionType(ParameterType):
        """A SMIRNOFF torsion type for improper torsions.

        .. warning :: This API is experimental and subject to change.
        """

        _ELEMENT_NAME = "Improper"

        periodicity = IndexedParameterAttribute(converter=int)
        phase = IndexedParameterAttribute(unit=unit.degree)
        k = IndexedParameterAttribute(unit=unit.kilocalorie / unit.mole)
        idivf = IndexedParameterAttribute(default=None, converter=float)

    _TAGNAME = "ImproperTorsions"  # SMIRNOFF tag name to process
    _INFOTYPE = ImproperTorsionType  # info type to store
    _OPENMMTYPE = "PeriodicTorsionForce"

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


class _NonbondedHandler(ParameterHandler):
    """Base class for ParameterHandlers that deal with OpenMM NonbondedForce objects."""

    _OPENMMTYPE = "NonbondedForce"


class vdWHandler(_NonbondedHandler):
    """Handle SMIRNOFF ``<vdW>`` tags

    .. warning :: This API is experimental and subject to change.
    """

    class vdWType(ParameterType):
        """A SMIRNOFF vdWForce type.

        .. warning :: This API is experimental and subject to change.
        """

        _ELEMENT_NAME = "Atom"

        epsilon = ParameterAttribute(unit=unit.kilocalorie / unit.mole)
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
                duplicate_attributes=[
                    *([] if duplicate_attributes is None else duplicate_attributes),
                    self._extra_nb_var,
                ],
            )

    _TAGNAME = "vdW"  # SMIRNOFF tag name to process
    _INFOTYPE = vdWType  # info type to store

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

    # TODO: Use _allow_only when ParameterAttribute will support multiple converters
    #       (it'll be easy when we switch to use the attrs library)
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
        unit_attrs_to_compare = ["cutoff", "switch_width"]

        self._check_attributes_are_equal(
            other_handler,
            identical_attrs=string_attrs_to_compare,
            tolerance_attrs=float_attrs_to_compare + unit_attrs_to_compare,
            tolerance=self._SCALETOL,
        )


class ElectrostaticsHandler(_NonbondedHandler):
    """Handles SMIRNOFF ``<Electrostatics>`` tags.

    .. warning :: This API is experimental and subject to change.
    """

    _TAGNAME = "Electrostatics"
    _DEPENDENCIES = [vdWHandler]
    _KWARGS = ["charge_from_molecules"]
    _MAX_SUPPORTED_SECTION_VERSION = Version("0.4")

    # Tolerance when comparing float attributes for handler compatibility.
    _SCALETOL = 1e-5
    _DEFAULT_REACTION_FIELD_EXPRESSION = (
        "charge1*charge2/(4*pi*epsilon0)*(1/r + k_rf*r^2 - c_rf);"
        "k_rf=(cutoff^(-3))*(solvent_dielectric-1)/(2*solvent_dielectric+1);"
        "c_rf=cutoff^(-1)*(3*solvent_dielectric)/(2*solvent_dielectric+1)"
    )

    scale12 = ParameterAttribute(default=0.0, converter=float)
    scale13 = ParameterAttribute(default=0.0, converter=float)
    scale14 = ParameterAttribute(default=0.833333, converter=float)
    scale15 = ParameterAttribute(default=1.0, converter=float)
    cutoff = ParameterAttribute(default=9.0 * unit.angstrom, unit=unit.angstrom)
    switch_width = ParameterAttribute(default=0.0 * unit.angstrom, unit=unit.angstrom)
    solvent_dielectric = ParameterAttribute(default=None)

    # TODO: How to validate arbitrary algebra in a converter?
    periodic_potential = ParameterAttribute(
        default="Ewald3D-ConductingBoundary",
        converter=_allow_only(
            [
                "Ewald3D-ConductingBoundary",
                "Coulomb",
                _DEFAULT_REACTION_FIELD_EXPRESSION,
            ]
        ),
    )
    nonperiodic_potential = ParameterAttribute(
        default="Coulomb", converter=_allow_only(["Coulomb"])
    )
    exception_potential = ParameterAttribute(
        default="Coulomb", converter=_allow_only(["Coulomb"])
    )

    # TODO: Use _allow_only when ParameterAttribute will support multiple converters
    #       (it'll be easy when we switch to use the attrs library)
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
        if new_switch_width not in [0.0 * unit.angstrom, None, "None", "none"]:
            raise SMIRNOFFSpecUnimplementedError(
                "The current implementation of the OpenFF Toolkit does not support an electrostatic "
                f"switch width (passed a value of {new_switch_width}). Currently only `0.0 angstroms` is supported "
                "and no switching function will be applied to the resulting `NonbondedForce`. If this behavior is "
                "important to you, please raise an issue at https://github.com/openforcefield/openff-toolkit/issues."
            )

    @periodic_potential.converter
    def periodic_potential(self, attr, new_value):
        if new_value in ["PME", "Ewald3D-ConductingBoundary"]:
            return "Ewald3D-ConductingBoundary"
        elif new_value in ["reaction-field", self._DEFAULT_REACTION_FIELD_EXPRESSION]:
            return self._DEFAULT_REACTION_FIELD_EXPRESSION
        elif new_value.lower() == "coulomb":
            return "Coulomb"
        else:
            raise NotImplementedError(
                "Failed to process unexpected periodic potential value: {new_value}"
            )

    @solvent_dielectric.converter
    def solvent_dielectric(self, attr, new_value):
        if new_value is not None:
            raise SMIRNOFFSpecUnimplementedError(
                "The current implementation of the OpenFF Toolkit does not support any electrostatic "
                "functions that make use of `solvent_dielectric`. If this behavior is important to you, please raise"
                "raise an issue at https://github.com/openforcefield/openff-toolkit/issues."
            )

    def __init__(self, **kwargs):
        if kwargs.get("version") == 0.4:
            if "method" in kwargs:
                raise SMIRNOFFSpecError(
                    "`method` attribute has been removed in version 0.4 of the Electrostatics tag. Use "
                    "`periodic_potential`, `nonperiodic_potenetial`, and `exception_potential` instead. "
                    "See https://openforcefield.github.io/standards/standards/smirnoff/#electrostatics"
                )
        if kwargs.get("version") == 0.3:
            logger.info(
                "Attempting to up-convert Electrostatics section from 0.3 to 0.4"
            )
            # Default value in 0.3 is "PME", so we have to handle these cases identically
            if kwargs.get("method") in ["PME", None]:
                kwargs["periodic_potential"] = "Ewald3D-ConductingBoundary"
                kwargs["nonperiodic_potential"] = "Coulomb"
                kwargs["exception_potential"] = "Coulomb"
                kwargs["version"] = 0.4
                kwargs.pop("method", None)
                logger.info(
                    'Successfully up-converted Electrostatics section from 0.3 to 0.4. `method="PME"` '
                    'is now split into `periodic_potential="Ewald3D-ConductingBoundary"`, '
                    '`nonperiodic_potential="Coulomb"`, and `exception_potential="Coulomb"`.'
                )
            elif kwargs["method"] == "Coulomb":
                kwargs["periodic_potential"] = "Coulomb"
                kwargs["nonperiodic_potential"] = "Coulomb"
                kwargs["exception_potential"] = "Coulomb"
                kwargs["version"] = 0.4
                kwargs.pop("method", None)
                logger.info(
                    'Successfully up-converted Electrostatics section from 0.3 to 0.4. `method="Coulomb"` '
                    'is now split into `periodic_potential="Coulob"`, '
                    '`nonperiodic_potential="Coulomb"`, and `exception_potential="Coulomb"`.'
                )
            elif kwargs["method"] == "reaction-field":
                kwargs["periodic_potential"] = self._DEFAULT_REACTION_FIELD_EXPRESSION
                kwargs["nonperiodic_potential"] = "Coulomb"
                kwargs["exception_potential"] = "Coulomb"
                kwargs["version"] = 0.4
                kwargs.pop("method", None)
                logger.info(
                    'Successfully up-converted Electrostatics section from 0.3 to 0.4. `method="Coulomb"` '
                    f'is now split into `periodic_potential="{self._DEFAULT_REACTION_FIELD_EXPRESSION}"` '
                    '`nonperiodic_potential="Coulomb"`, and `exception_potential="Coulomb"`.'
                )
            else:
                raise NotImplementedError(
                    "Failed to up-convert Electrostatics section from 0.3 to 0.4. Did not know "
                    f"how to up-convert `method={kwargs['method']}`."
                )
        super().__init__(**kwargs)

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
        string_attrs_to_compare = [
            "periodic_potential",
            "nonperiodic_potential",
            "exception_potential",
        ]
        unit_attrs_to_compare = ["cutoff", "switch_width"]

        self._check_attributes_are_equal(
            other_handler,
            identical_attrs=string_attrs_to_compare,
            tolerance_attrs=float_attrs_to_compare + unit_attrs_to_compare,
            tolerance=self._SCALETOL,
        )


class LibraryChargeHandler(_NonbondedHandler):
    """Handle SMIRNOFF ``<LibraryCharges>`` tags

    .. warning :: This API is experimental and subject to change.
    """

    class LibraryChargeType(ParameterType):
        """A SMIRNOFF Library Charge type.

        .. warning :: This API is experimental and subject to change.
        """

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
        def from_molecule(cls, molecule: Molecule):
            """
            Construct a LibraryChargeType from a molecule with existing partial charges.

            Parameters
            ----------
            molecule : openff.toolkit.topology.molecule.Molecule
                The molecule to create the LibraryChargeType from. The molecule must have partial charges.

            Returns
            -------
            library_charge_type : LibraryChargeType
                A LibraryChargeType that is expected to match this molecule and its partial charges.

            Raises
            ------

            MissingPartialChargesError : If the input molecule does not have partial charges.
            """
            if molecule.partial_charges is None:
                raise MissingPartialChargesError(
                    "Input molecule is missing partial charges."
                )

            smirks = molecule.to_smiles(mapped=True)
            charges = molecule.partial_charges

            library_charge_type = cls(smirks=smirks, charge=charges)

            return library_charge_type

    _TAGNAME = "LibraryCharges"  # SMIRNOFF tag name to process
    _INFOTYPE = LibraryChargeType  # info type to store
    _DEPENDENCIES = [vdWHandler, ElectrostaticsHandler]

    def find_matches(self, entity, unique=False):
        """Find the elements of the topology/molecule matched by a parameter type.

        Parameters
        ----------
        entity : openff.toolkit.topology.Topology
            Topology to search.

        Returns
        ---------
        matches : ValenceDict[Tuple[int], ParameterHandler._Match]
            ``matches[atom_indices]`` is the ``ParameterType`` object
            matching the tuple of atom indices in ``entity``.
        """

        return self._find_matches(
            entity,
            transformed_dict_cls=dict,
            unique=unique,
        )


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


class ChargeIncrementModelHandler(_NonbondedHandler):
    """Handle SMIRNOFF ``<ChargeIncrementModel>`` tags

    .. warning :: This API is experimental and subject to change.
    """

    class ChargeIncrementType(ParameterType):
        """A SMIRNOFF bond charge correction type.

        .. warning :: This API is experimental and subject to change.
        """

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
    _MAX_SUPPORTED_SECTION_VERSION = Version("0.4")

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
            ``matches[atom_indices]`` is the ``ParameterType`` object
            matching the tuple of atom indices in ``entity``.
        """
        matches = self._find_matches(
            entity, transformed_dict_cls=TagSortedDict, unique=unique
        )
        return matches


class GBSAHandler(ParameterHandler):
    """Handle SMIRNOFF ``<GBSA>`` tags

    .. warning :: This API is experimental and subject to change.
    """

    class GBSAType(ParameterType):
        """A SMIRNOFF GBSA type.

        .. warning :: This API is experimental and subject to change.
        """

        _ELEMENT_NAME = "Atom"

        radius = ParameterAttribute(unit=unit.angstrom)
        scale = ParameterAttribute(converter=float)

    _TAGNAME = "GBSA"
    _INFOTYPE = GBSAType
    _OPENMMTYPE = "GBSAOBCForce"
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
        default=unit.Quantity(5.4, _cal_mol_a2),
        unit=_cal_mol_a2,
    )
    solvent_radius = ParameterAttribute(default=1.4 * unit.angstrom, unit=unit.angstrom)

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


_VirtualSiteType = Literal[
    "BondCharge",
    "MonovalentLonePair",
    "DivalentLonePair",
    "TrivalentLonePair",
]


class VirtualSiteHandler(_NonbondedHandler):
    """Handle SMIRNOFF ``<VirtualSites>`` tags
    TODO: Add example usage/documentation
    .. warning :: This API is experimental and subject to change.
    """

    class VirtualSiteType(vdWHandler.vdWType):
        _ELEMENT_NAME = "VirtualSite"

        name = ParameterAttribute(default="EP", converter=str)
        type = ParameterAttribute(converter=str)

        match = ParameterAttribute(converter=str)

        distance = ParameterAttribute(unit=unit.angstrom)
        outOfPlaneAngle = ParameterAttribute(unit=unit.degree)
        inPlaneAngle = ParameterAttribute(unit=unit.degree)

        epsilon = ParameterAttribute(
            default=0.0 * unit.kilocalorie_per_mole, unit=unit.kilocalorie_per_mole
        )
        sigma = ParameterAttribute(default=1.0 * unit.angstrom, unit=unit.angstrom)
        rmin_half = ParameterAttribute(default=None, unit=unit.angstrom)

        charge_increment = IndexedParameterAttribute(unit=unit.elementary_charge)

        @property
        def parent_index(self) -> int:
            """Returns the index of the atom matched by the SMIRKS pattern that should
            be considered the 'parent' to the virtual site.
            A value of ``0`` corresponds to the atom matched by the ``:1`` selector in
            the SMIRKS pattern, a value ``2`` the atom matched by ``:2`` and so on.
            """
            return self.type_to_parent_index(self.type)

        @classmethod
        def type_to_parent_index(cls, type_: _VirtualSiteType) -> int:
            """Returns the index of the atom matched by the SMIRKS pattern that should
            be considered the 'parent' to a given type of virtual site.
            A value of ``0`` corresponds to the atom matched by the ``:1`` selector in
            the SMIRKS pattern, a value ``2`` the atom matched by ``:2`` and so on.
            """

            if type_.replace("VirtualSite", "") in get_args(_VirtualSiteType):
                return 0

            raise NotImplementedError()

        @outOfPlaneAngle.converter
        def outOfPlaneAngle(self, attr, value):
            if value == "None":
                return

            supports_out_of_plane_angle = self._supports_out_of_plane_angle(self.type)

            if not supports_out_of_plane_angle and value is not None:
                raise SMIRNOFFSpecError(
                    f"'{self.type}' sites do not support `outOfPlaneAngle`"
                )
            elif supports_out_of_plane_angle:
                return _validate_units(attr, value, unit.degrees)

            return value

        @inPlaneAngle.converter
        def inPlaneAngle(self, attr, value):
            if value == "None":
                return

            supports_in_plane_angle = self._supports_in_plane_angle(self.type)

            if not supports_in_plane_angle and value is not None:
                raise SMIRNOFFSpecError(
                    f"'{self.type}' sites do not support `inPlaneAngle`"
                )
            elif supports_in_plane_angle:
                return _validate_units(attr, value, unit.degrees)

            return value

        def __init__(self, **kwargs):
            self._add_default_init_kwargs(kwargs)
            super().__init__(**kwargs)

        @classmethod
        def _add_default_init_kwargs(cls, kwargs):
            """Adds any missing default values to the ``kwargs`` dictionary, and
            partially validates any provided values that aren't easily validated with
            converters.
            """

            type_ = kwargs.get("type", None)

            if type_ is None:
                raise SMIRNOFFSpecError("the `type` keyword is missing")
            if type_ not in get_args(_VirtualSiteType):
                raise SMIRNOFFSpecError(
                    f"'{type_}' is not a supported virtual site type"
                )

            if "charge_increment" in kwargs:
                expected_num_charge_increments = cls._expected_num_charge_increments(
                    type_
                )
                num_charge_increments = len(kwargs["charge_increment"])
                if num_charge_increments != expected_num_charge_increments:
                    raise SMIRNOFFSpecError(
                        f"'{type_}' virtual sites expect exactly {expected_num_charge_increments} "
                        f"charge increments, but got {kwargs['charge_increment']} "
                        f"(length {num_charge_increments}) instead."
                    )

            supports_in_plane_angle = cls._supports_in_plane_angle(type_)
            supports_out_of_plane_angle = cls._supports_out_of_plane_angle(type_)

            if not supports_out_of_plane_angle:
                kwargs["outOfPlaneAngle"] = kwargs.get("outOfPlaneAngle", None)
            if not supports_in_plane_angle:
                kwargs["inPlaneAngle"] = kwargs.get("inPlaneAngle", None)

            match = kwargs.get("match", None)

            if match is None:
                raise SMIRNOFFSpecError("the `match` keyword is missing")

            out_of_plane_angle = kwargs.get("outOfPlaneAngle", 0.0 * unit.degree)
            is_in_plane = (
                None
                if not supports_out_of_plane_angle
                else np.isclose(out_of_plane_angle.m_as(unit.degree), 0.0)
            )

            if not cls._supports_match(type_, match, is_in_plane):
                raise SMIRNOFFSpecError(
                    f"match='{match}' not supported with type='{type_}'"
                    + ("" if is_in_plane is None else f" and is_in_plane={is_in_plane}")
                )

            if "rmin_half" not in kwargs:
                kwargs["sigma"] = kwargs.get("sigma", 0.0 * unit.angstrom)

            kwargs["epsilon"] = kwargs.get("epsilon", 0.0 * unit.kilocalorie_per_mole)

        @classmethod
        def _supports_in_plane_angle(cls, type_: _VirtualSiteType) -> bool:
            return type_ in {"MonovalentLonePair"}

        @classmethod
        def _supports_out_of_plane_angle(cls, type_: _VirtualSiteType) -> bool:
            return type_ in {"MonovalentLonePair", "DivalentLonePair"}

        @classmethod
        def _expected_num_charge_increments(cls, type_: _VirtualSiteType) -> int:
            if type_ == "BondCharge":
                return 2
            elif (type_ == "MonovalentLonePair") or (type_ == "DivalentLonePair"):
                return 3
            elif type_ == "TrivalentLonePair":
                return 4
            raise NotImplementedError()

        @classmethod
        def _supports_match(
            cls, type_: _VirtualSiteType, match: str, is_in_plane: Optional[bool] = None
        ) -> bool:
            is_in_plane = True if is_in_plane is None else is_in_plane

            if match == "once":
                return type_ == "TrivalentLonePair" or (
                    type_ == "DivalentLonePair" and is_in_plane
                )
            elif match == "all_permutations":
                return type_ in {"BondCharge", "MonovalentLonePair", "DivalentLonePair"}

            raise NotImplementedError()

    _TAGNAME = "VirtualSites"
    _INFOTYPE = VirtualSiteType
    _OPENMMTYPE = "NonbondedForce"
    _DEPENDENCIES = [
        ElectrostaticsHandler,
        LibraryChargeHandler,
        ChargeIncrementModelHandler,
        ToolkitAM1BCCHandler,
        vdWHandler,
    ]

    exclusion_policy = ParameterAttribute(default="parents")

    @classmethod
    def _validate_found_match(
        cls,
        atoms_by_index: Dict,
        matched_indices: Tuple[int, ...],
        parameter: VirtualSiteType,
    ):
        """
        We place limitations on the chemical environments that certain types of v-site
        can be applied to, e.g. we enforce that divalent lone pairs can only be applied
        to environments that look like a carboxyl group.
        These somewhat artificial restrictions limit the number of potential edge
        cases that need to be thought through, and significantly reduces the number
        of test cases / problematic choices that need to be made. If users meet a not
        supported exception, they should open an issue on GitHub explaining their exact
        use case so that we can ensure that exactly what they need is both supported
        and works as expected through expansion of the unit tests.
        """

        supported_connectivity = {
            # We currently expect monovalent lone pairs to be applied to something
            # like a carboxyl group, where the parent of the lone pair has a
            # connectivity of 1, while it neighbour has a connectivity of 3
            ("MonovalentLonePair", 0): 1,
            ("MonovalentLonePair", 1): 3,
            # We currently expect divalent lone pairs to be applied to something
            # like an sp2 nitrogen, or a hydroxyl oxygen
            ("DivalentLonePair", 0): 2,
            # We currently expect trivalent lone pairs to be applied to something
            # like an sp3 nitrogen
            ("TrivalentLonePair", 0): 3,
        }

        for smirks_index, atom_index in enumerate(matched_indices):
            if (parameter.type, smirks_index) not in supported_connectivity:
                # No restrictions placed on this matched atom.
                continue

            matched_atom = atoms_by_index[atom_index]
            connectivity = len(matched_atom.bonds)
            expected_connectivity = supported_connectivity[
                (parameter.type, smirks_index)
            ]
            if expected_connectivity == connectivity:
                continue

            raise NotImplementedError(
                f"{parameter.smirks} matched chemical environment that is currently "
                f"unsupported by virtual sites of type {parameter.type}. Atom with "
                f"smirks index={smirks_index} matched topology atom {atom_index} with "
                f"connectivity={connectivity}, but it was expected to have connectivity "
                f"{expected_connectivity}. If this is "
                f"a use case you would like supported, please describe what it is "
                f"you are trying to do in an issue on the OpenFF Toolkit GitHub: "
                f"https://github.com/openforcefield/openff-toolkit/issues"
            )

    def check_handler_compatibility(self, other_handler):
        self._check_attributes_are_equal(
            other_handler, identical_attrs=["exclusion_policy"]
        )

    def _index_of_parameter(
        self,
        parameter: Optional[ParameterType] = None,
        key: Optional[Any] = None,
    ) -> Optional[int]:
        """Attempts to find the index of a parameter in the parameters list.
        By default, two parameters are considered 'the same' if they have the same
        SMIRKS pattern, type, and name.
        Parameters
        ----------
        parameter
            The parameter to find the index of. This argument is mutually exclusive with
            ``key``.
        key
            A tuple of the type, SMIRKS, and name of the parameter to find the index
            of. This argument is mutually exclusive with ``parameter``.
        Returns
        -------
            The index of the parameter if found, otherwise ``None``.
        """

        if (key is None and parameter is None) or (
            key is not None and parameter is not None
        ):
            raise ValueError("`key` and `parameter` are mutually exclusive arguments")

        key = cast(
            Tuple[str, str, str],
            key
            if parameter is None
            else (parameter.type, parameter.smirks, parameter.name),
        )
        expected_type, expected_smirks, expected_name = key

        for i, existing_parameter in enumerate(self.parameters):
            if (
                existing_parameter.type != expected_type
                or existing_parameter.smirks != expected_smirks
                or existing_parameter.name != expected_name
            ):
                continue

            return i

        return None

    def _find_matches_by_parent(self, entity: Topology) -> Dict[int, list]:
        from collections import defaultdict

        topology_atoms = {
            i: topology_atom for i, topology_atom in enumerate(entity.atoms)
        }

        # We need to find all the parameters that would lead to a v-site being placed
        # onto a given 'parent atom'. We only allow each parent atom to be assigned one
        # v-site with a given 'name', whereby the last parameter to be matched wins.
        matches_by_parent: Dict = defaultdict(lambda: defaultdict(list))

        for parameter in self._parameters:
            for match in entity.chemical_environment_matches(parameter.smirks):
                parent_index = match.topology_atom_indices[parameter.parent_index]

                matches_by_parent[parent_index][parameter.name].append(
                    (parameter, match)
                )

        # we then need to find which parameter was the last one to be assigned to each
        # given 'parent' atom, and all the ways that that parameter matches the atoms
        # surrounding the parent. Whether we keep the different orientations or not
        # depends on the 'match' setting of the parameter.
        assigned_matches_by_parent = defaultdict(list)

        for parent_index, matches_by_name in matches_by_parent.items():
            for matches in matches_by_name.values():
                assigned_parameter, _ = matches[-1]  # last match wins

                match_orientations = [
                    match
                    for parameter_index, match in matches
                    if parameter_index == assigned_parameter
                ]

                if assigned_parameter.match == "once":
                    # the v-site types were designed such that we should be safe
                    # choosing an arbitrary ordering of the non-parent matched atoms.
                    match_orientations = [match_orientations[0]]
                elif assigned_parameter.match == "all_permutations":
                    pass
                else:
                    raise SMIRNOFFSpecError(
                        f"{assigned_parameter.match} match keyword is not supported"
                    )

                assigned_matches_by_parent[parent_index].append(
                    (assigned_parameter, match_orientations)
                )

                for match in match_orientations:
                    # make sure the match does not look like a weird edge case that we
                    # haven't tested to ensure 'sensible' behaviour in most cases.
                    self._validate_found_match(
                        topology_atoms, match.topology_atom_indices, assigned_parameter
                    )

        return assigned_matches_by_parent

    def _find_matches(
        self,
        entity: Topology,
        transformed_dict_cls=dict,
        unique=False,
    ) -> Dict[Tuple[int], List[ParameterHandler._Match]]:
        assigned_matches_by_parent = self._find_matches_by_parent(entity)
        return_dict = {}
        for parent_index, assigned_parameters in assigned_matches_by_parent.items():
            assigned_matches = []
            for assigned_parameter, match_orientations in assigned_parameters:
                for match in match_orientations:
                    assigned_matches.append(
                        ParameterHandler._Match(assigned_parameter, match)
                    )
            return_dict[(parent_index,)] = assigned_matches

        return return_dict


ConstraintType = ConstraintHandler.ConstraintType
BondType = BondHandler.BondType
AngleType = AngleHandler.AngleType
ProperTorsionType = ProperTorsionHandler.ProperTorsionType
ImproperTorsionType = ImproperTorsionHandler.ImproperTorsionType
vdWType = vdWHandler.vdWType
LibraryChargeType = LibraryChargeHandler.LibraryChargeType
GBSAType = GBSAHandler.GBSAType
ChargeIncrementType = ChargeIncrementModelHandler.ChargeIncrementType
VirtualSiteType = VirtualSiteHandler.VirtualSiteType


if __name__ == "__main__":
    import doctest

    doctest.testmod()
    # doctest.run_docstring_examples(_ParameterAttributeHandler, globals())
