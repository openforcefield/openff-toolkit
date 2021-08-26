import copy
import functools
import inspect
import re
from collections import OrderedDict, defaultdict

from openff.toolkit.typing.chemistry import ChemicalEnvironment
from openff.toolkit.utils.collections import ValidatedDict, ValidatedList
from openff.toolkit.utils.exceptions import ParameterLookupError, SMIRNOFFSpecError
from openff.toolkit.utils.utils import IncompatibleUnitError, object_to_quantity


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
    unit : simtk.unit.Quantity, optional
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

    >>> from simtk import unit
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
    unit : simtk.unit.Quantity, optional
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

    >>> from simtk import unit
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
    unit : simtk.unit.Quantity, optional
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

    >>> from simtk import unit
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
    unit : simtk.unit.Quantity, optional
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

    >>> from simtk import unit
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
            filter = lambda x: True

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
