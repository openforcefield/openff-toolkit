"""
Core classes for OpenFF Units
"""

from __future__ import annotations

import ast
import json
import operator as op
import uuid
import warnings
from typing import TYPE_CHECKING, Any, Literal, overload

import numpy  # possible to make this optional?
import pint
from openff.utilities import has_package, requires_package
from pint import Measurement as _Measurement
from pint import Quantity as _Quantity
from pint import Unit as _Unit
from pint.facets.plain.quantity import PlainQuantity as PintQuantity

from openff.toolkit.utils.exceptions import OpenFFToolkitException

if TYPE_CHECKING:
    import openmm.unit

try:
    from pydantic import GetCoreSchemaHandler
    from pydantic_core import core_schema

    has_pydantic = True
except ImportError:
    has_pydantic = False


__all__ = (
    "DEFAULT_UNIT_REGISTRY",
    "Measurement",
    "Quantity",
    "Unit",
    "ensure_quantity",
    "from_openmm",
    "openmm_unit_to_string",
    "string_to_openmm_unit",
    "to_openmm",
    "unit",
)


class MissingOpenMMUnitError(OpenFFToolkitException):
    """Raised when a unit cannot be converted to an equivalent OpenMM unit"""


class NoneQuantityError(OpenFFToolkitException):
    """Raised when attempting to convert `None` between unit packages as a quantity object"""


class NoneUnitError(OpenFFToolkitException):
    """Raised when attempting to convert `None` between unit packages as a unit object"""


"""
from openff.toolkit.utils.exceptions import (
    MissingOpenMMUnitError,
    NoneQuantityError,
    NoneUnitError,
)
"""


def get_defaults_path() -> str:
    """Get the full path to the ``defaults.txt`` file"""
    from openff.utilities import get_data_file_path

    return get_data_file_path("data/units/defaults.txt", "openff.toolkit")


class Unit(pint.UnitRegistry.Unit):
    """A unit of measure."""

    pass


if has_pydantic:

    class _QuantityMixin:
        @classmethod
        def serialize(
            cls,
            v: PintQuantity,
            info: core_schema.SerializationInfo | None = None,
        ) -> dict | str | PintQuantity:
            to_json = info is not None and info.mode_is_json()

            if to_json:
                magnitude = v.magnitude

                # storing numpy arrays natively works fine in memory in Python,
                # but must be list-ified when serializing to JSON.
                if isinstance(magnitude, numpy.ndarray):
                    magnitude = v.magnitude.tolist()

                # TODO: I think this is necessary for handling unit-wrapped arrays, but it is
                #       not so performant. Scalar quantities can be directly serialized to much
                #       shorter strings
                return json.dumps(
                    {
                        "magnitude": magnitude,
                        "units": str(v.units),
                    }
                )

            return {
                "magnitude": v.magnitude,
                "units": str(v.units),
            }

        @classmethod
        def validate(
            cls,
            v: dict | str | PintQuantity,
        ):
            if isinstance(v, Quantity):
                return v
            elif isinstance(v, str):
                # TODO: A significant wart is that we have to try to guess whether the string is a
                #       JSON-serialized quantity or a simple string representation of a quantity.
                #       For example: input of "0.9 nanometer" can be passed directly to the
                #       Quantity constructor, but "{"magnitude": 0.9, "units": "nanometer"}" cannot
                #       as it needs to be unwrapped. A better solution would require a better way
                #       of serializing unit-wrapped arrays to JSON
                if "{" in v:
                    deserialized = json.loads(v)
                    return Quantity(
                        deserialized["magnitude"],
                        deserialized["units"],
                    )
                else:
                    return Quantity(v)
            elif isinstance(v, dict):
                return Quantity(v["magnitude"], v["units"])
            else:
                # this cannot be accessed with the current core_schema definition - the types of
                # the `v` argument to this method **happen** to be identical to the supported types
                # in the core_schema. If **either** is changed, this clause may be hit
                raise ValueError(f"Invalid type {type(v)} for Quantity")

        @classmethod
        def __get_pydantic_core_schema__(
            cls,
            source_type: Any,
            handler: GetCoreSchemaHandler,
        ) -> core_schema.CoreSchema:

            validate_schema = core_schema.chain_schema(
                [
                    core_schema.union_schema(
                        [
                            core_schema.is_instance_schema(PintQuantity),
                            core_schema.str_schema(),
                            core_schema.dict_schema(),
                            # any other types that could be accepted by a quantity field?
                        ]
                    ),
                    core_schema.no_info_plain_validator_function(cls.validate),
                ]
            )

            validate_json_schema = core_schema.chain_schema(
                [
                    core_schema.union_schema(
                        [
                            core_schema.str_schema(coerce_numbers_to_str=True),
                            core_schema.dict_schema(),
                        ]
                    ),
                    core_schema.no_info_plain_validator_function(cls.validate),
                ]
            )

            serialize_schema = core_schema.plain_serializer_function_ser_schema(
                cls.serialize,
                info_arg=True,
            )

            return core_schema.json_or_python_schema(
                json_schema=validate_json_schema,
                python_schema=validate_schema,
                serialization=serialize_schema,
            )
else:

    class _QuantityMixin:  # type: ignore[no-redef]
        pass


class Quantity(_QuantityMixin, PintQuantity):
    """A value with associated units."""

    def __dask_tokenize__(self):
        return uuid.uuid4().hex

    @staticmethod
    def _dask_finalize(results, func, args, units):
        values = func(results, *args)
        return Quantity(values, units)


@requires_package("openmm")
def _to_openmm(self) -> openmm.unit.Quantity:
    """Convert the quantity to an ``openmm.unit.Quantity``.

    Returns
    -------
    openmm_quantity : openmm.unit.quantity.Quantity
        The OpenMM compatible quantity.
    """
    from openff.toolkit.utils.units import to_openmm

    return to_openmm(self)


class Measurement(pint.UnitRegistry.Measurement):  # type: ignore
    """A value with associated units and uncertainty."""

    def __dask_tokenize__(self):
        return uuid.uuid4().hex

    @staticmethod
    def _dask_finalize(results, func, args, units):
        values = func(results, *args)
        return Measurement(values, units)


class UnitRegistry(pint.UnitRegistry):
    _quantity_class = Quantity
    _unit_class = Unit
    _measurement_class = Measurement


DEFAULT_UNIT_REGISTRY = UnitRegistry(get_defaults_path())

unit = DEFAULT_UNIT_REGISTRY

Unit: type[_Unit] = DEFAULT_UNIT_REGISTRY.Unit  # type: ignore[no-redef]
Quantity: type[_Quantity] = DEFAULT_UNIT_REGISTRY.Quantity  # type: ignore[no-redef]
Measurement: type[_Measurement] = DEFAULT_UNIT_REGISTRY.Measurement  # type: ignore

Quantity.to_openmm = _to_openmm  # type: ignore[attr-defined]

if has_pydantic:
    # Re-attach the Pydantic magic to our new Quantity class, which itself was
    # dynamically created by Pint's magic (and lost these methods in the process).
    Quantity.__get_pydantic_core_schema__ = _QuantityMixin.__get_pydantic_core_schema__  # type: ignore
    Quantity.validate = _QuantityMixin.validate  # type: ignore
    Quantity.serialize = _QuantityMixin.serialize  # type: ignore


pint.set_application_registry(DEFAULT_UNIT_REGISTRY)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    Quantity([])


# below this line was openff/units/openmm.py

if has_package("openmm.unit"):
    import openmm.unit

    EitherQuantity: type[Quantity | openmm.unit.Quantity] = Quantity | openmm.unit.Quantity  # type: ignore[no-redef]
else:
    EitherQuantity: type[Quantity] = Quantity  # type: ignore[no-redef]


@requires_package("openmm.unit")
def openmm_unit_to_string(input_unit: openmm.unit.Unit) -> str:
    """
    Convert a openmm.unit.Unit to a string representation.

    Parameters
    ----------
    input_unit : A openmm.unit
        The unit to serialize

    Returns
    -------
    unit_string : str
        The serialized unit.
    """
    if input_unit is None:
        raise NoneUnitError("Input is None, expected an (OpenMM) Unit object.")

    if input_unit == openmm.unit.dimensionless:
        return "dimensionless"

    if input_unit == openmm.unit.dalton:
        return "g/mol"

    # Decompose output_unit into a tuples of (base_dimension_unit, exponent)
    unit_string = ""

    for unit_component in input_unit.iter_base_or_scaled_units():
        unit_component_name = unit_component[0].name
        # Convert, for example "elementary charge" --> "elementary_charge"
        unit_component_name = unit_component_name.replace(" ", "_")
        if unit_component[1] == 1:
            contribution = f"{unit_component_name}"
        else:
            contribution = f"{unit_component_name}**{int(unit_component[1])}"
        if unit_string == "":
            unit_string = contribution
        else:
            unit_string += f" * {contribution}"

    return unit_string


def _ast_eval(node):
    """
    Performs an algebraic syntax tree evaluation of a unit.

    Parameters
    ----------
    node : An ast parsing tree node

    Raises
    ------
    openff.units.exceptions.MissingOpenMMUnitError
        if the unit is unavailable in OpenMM.
    """

    operators = {
        ast.Add: op.add,
        ast.Sub: op.sub,
        ast.Mult: op.mul,
        ast.Div: op.truediv,
        ast.Pow: op.pow,
        ast.BitXor: op.xor,
        ast.USub: op.neg,
    }

    if isinstance(node, ast.Constant):  # <number>
        return node.value
    elif isinstance(node, ast.BinOp):  # <left> <operator> <right>
        return operators[type(node.op)](_ast_eval(node.left), _ast_eval(node.right))
    elif isinstance(node, ast.UnaryOp):  # <operator> <operand> e.g., -1
        return operators[type(node.op)](_ast_eval(node.operand))
    elif isinstance(node, ast.Name):
        # see if this is a openmm unit
        try:
            b = getattr(openmm.unit, node.id)
        except AttributeError:
            raise MissingOpenMMUnitError(node.id)
        return b
    # TODO: This toolkit code that had a hack to cover some edge behavior
    #       not clear which tests trigger it
    elif isinstance(node, ast.List):
        return ast.literal_eval(node)
    else:
        raise TypeError(node)


def string_to_openmm_unit(unit_string: str) -> openmm.unit.Unit:
    """
    Deserializes a openmm.unit.Quantity from a string representation, for
    example: "kilocalories_per_mole / angstrom ** 2"


    Parameters
    ----------
    unit_string : dict
        Serialized representation of a openmm.unit.Quantity.

    Returns
    -------
    output_unit: openmm.unit.Quantity
        The deserialized unit from the string

    Raises
    ------
    openff.units.exceptions.MissingOpenMMUnitError
        if the unit is unavailable in OpenMM.
    """
    if unit_string == "standard_atmosphere":
        return openmm.unit.atmosphere

    output_unit = _ast_eval(ast.parse(unit_string, mode="eval").body)
    return output_unit


@requires_package("openmm.unit")
def from_openmm(openmm_quantity: openmm.unit.Quantity) -> Quantity:
    """Convert an OpenMM ``Quantity`` to an OpenFF ``Quantity``

    :class:`openmm.unit.quantity.Quantity` from OpenMM and
    :class:`openff.units.Quantity` from this package both represent a numerical
    value with units.

    Examples
    --------

    >>> from openff.units import Quantity as OpenFFQuantity
    >>> from openff.units import from_openmm
    >>> from openmm import unit
    >>> length = unit.Quantity(9.0, unit.angstrom)
    >>> from_openmm(length)
    <Quantity(9.0, 'angstrom')>
    >>> assert isinstance(from_openmm(length), OpenFFQuantity)

    """
    if openmm_quantity is None:
        raise NoneQuantityError("Input is None, expected an (OpenMM) Quantity object.")

    if isinstance(openmm_quantity, list):
        openmm_quantity = openmm.unit.Quantity(openmm_quantity)
    openmm_unit_ = openmm_quantity.unit
    openmm_value = openmm_quantity.value_in_unit(openmm_unit_)

    target_unit = openmm_unit_to_string(openmm_unit_)
    target_unit = Unit(target_unit)

    return Quantity(openmm_value, target_unit)


@requires_package("openmm.unit")
def to_openmm(quantity: Quantity) -> openmm.unit.Quantity:
    """Convert an OpenFF ``Quantity`` to an OpenMM ``Quantity``

    :class:`openmm.unit.quantity.Quantity` from OpenMM and
    :class:`openff.units.Quantity` from this package both represent a numerical
    value with units. The units available in the two packages differ; when a
    unit is missing from the target package, the resulting quantity will be in
    base units (kg/m/s/A/K/mole), which are shared between both packages. This
    may cause the resulting value to be slightly different to the input due to
    the limited precision of floating point numbers.

    Examples
    --------

    >>> from openff.units import unit
    >>> from openff.units import to_openmm
    >>> from openmm import unit as openmm_unit
    >>> length = unit.Quantity(9.0, unit.angstrom)
    >>> to_openmm(length)
    9.0 A
    >>> assert isinstance(to_openmm(length), openmm_unit.Quantity)

    """
    if quantity is None:
        raise NoneQuantityError("Input is None, expected an (OpenFF) Quantity object.")

    def to_openmm_inner(quantity) -> openmm.unit.Quantity:
        value = quantity.m

        unit_string = str(quantity.units._units)
        openmm_unit_ = string_to_openmm_unit(unit_string)

        return value * openmm_unit_

    assert isinstance(quantity, Quantity)

    try:
        return to_openmm_inner(quantity)
    except MissingOpenMMUnitError:
        return to_openmm_inner(quantity.to_base_units())


@requires_package("openmm.unit")
def _ensure_openmm_quantity(
    unknown_quantity: EitherQuantity,
) -> openmm.unit.Quantity:
    if "openmm" in str(type(unknown_quantity)):
        from openmm import unit as openmm_unit

        if isinstance(unknown_quantity, openmm_unit.Quantity):
            return unknown_quantity
        else:
            raise ValueError(f"Failed to process input of type {type(unknown_quantity)}.")
    elif isinstance(unknown_quantity, Quantity):
        return to_openmm(unknown_quantity)
    else:
        import openmm.unit

        try:
            return openmm.unit.Quantity(
                unknown_quantity,
                openmm.unit.dimensionless,
            )
        except Exception as e:
            raise ValueError(f"Failed to process input of type {type(unknown_quantity)}.") from e


def _ensure_openff_quantity(
    unknown_quantity: EitherQuantity,
) -> Quantity:
    if isinstance(unknown_quantity, Quantity):
        return unknown_quantity
    elif "openmm" in str(type(unknown_quantity)):
        import openmm.unit

        if isinstance(unknown_quantity, openmm.unit.Quantity):
            return from_openmm(unknown_quantity)
        else:
            raise ValueError(f"Failed to process input of type {type(unknown_quantity)}.")
    else:
        try:
            return Quantity(
                unknown_quantity,
                "dimensionless",
            )
        except Exception as e:
            raise ValueError(f"Failed to process input of type {type(unknown_quantity)}.") from e


@overload
def ensure_quantity(
    unknown_quantity: EitherQuantity,
    type_to_ensure: Literal["openmm"],
) -> openmm.unit.Quantity: ...


@overload
def ensure_quantity(
    unknown_quantity: EitherQuantity,
    type_to_ensure: Literal["openff"],
) -> Quantity: ...


def ensure_quantity(
    unknown_quantity: EitherQuantity,
    type_to_ensure: Literal["openmm", "openff"],
) -> EitherQuantity:
    """
    Given a quantity that could be of a variety of types, attempt to coerce into a given type.

    Examples
    --------

    >>> import numpy
    >>> from openmm import unit as openmm_unit
    >>> from openff.units import unit
    >>> from openff.units import ensure_quantity
    >>> # Create a 9 Angstrom quantity with each registry
    >>> length1 = unit.Quantity(9.0, unit.angstrom)
    >>> length2 = openmm_unit.Quantity(9.0, openmm_unit.angstrom)
    >>> # Similar quantities are be coerced into requested type
    >>> assert type(ensure_quantity(length1, "openmm")) == openmm_unit.Quantity
    >>> assert type(ensure_quantity(length2, "openff")) == unit.Quantity
    >>> # Seemingly-redundant "conversions" short-circuit
    >>> assert ensure_quantity(length1, "openff") == ensure_quantity(length2, "openff")
    >>> assert ensure_quantity(length1, "openmm") == ensure_quantity(length2, "openmm")
    >>> # NumPy arrays and some primitives are automatically up-converted to `Quantity` objects
    >>> # Note that their units are set to "dimensionless"
    >>> ensure_quantity(numpy.array([1, 2]), "openff")
    <Quantity([1 2], 'dimensionless')>
    >>> ensure_quantity(4.0, "openmm")
    4.0 dimensionless

    """
    if type_to_ensure == "openmm":
        return _ensure_openmm_quantity(unknown_quantity)
    elif type_to_ensure == "openff":
        return _ensure_openff_quantity(unknown_quantity)
    else:
        raise ValueError(f"Unsupported `type_to_ensure` found. Given {type_to_ensure}, expected 'openff' or 'openmm'.")
