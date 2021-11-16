#!/usr/bin/env python
"""
Utility subroutines

"""

__all__ = [
    "requires_package",
    "inherit_docstrings",
    "all_subclasses",
    "temporary_cd",
    "get_data_file_path",
    "unit_to_string",
    "quantity_to_string",
    "string_to_unit",
    "string_to_quantity",
    "object_to_quantity",
    "check_units_are_compatible",
    "extract_serialized_units_from_dict",
    "attach_units",
    "detach_units",
    "serialize_numpy",
    "deserialize_numpy",
    "convert_all_quantities_to_string",
    "convert_all_strings_to_quantity",
    "convert_0_1_smirnoff_to_0_2",
    "convert_0_2_smirnoff_to_0_3",
    "get_molecule_parameterIDs",
]


import contextlib
import functools
import logging

from openff.units import unit
from openmm import unit as openmm_unit

from openff.toolkit.utils.exceptions import (
    IncompatibleUnitError,
    MissingDependencyError,
)

# =============================================================================================
# CONFIGURE LOGGER
# =============================================================================================

logger = logging.getLogger(__name__)


# =============================================================================================
# UTILITY SUBROUTINES
# =============================================================================================


def requires_package(package_name):
    """
    Helper function to denote that a funciton requires some optional
    dependency. A function decorated with this decorator will raise
    `MissingDependencyError` if the package is not found by
    `importlib.import_module()`.

    Parameters
    ----------
    package_name : str
        The directory path to enter within the context

    Raises
    ------
    MissingDependencyError

    """

    def inner_decorator(function):
        @functools.wraps(function)
        def wrapper(*args, **kwargs):
            import importlib

            try:
                importlib.import_module(package_name)
            except (ImportError, ModuleNotFoundError):
                raise MissingDependencyError(package_name)
            except Exception as e:
                raise e

            return function(*args, **kwargs)

        return wrapper

    return inner_decorator


def inherit_docstrings(cls):
    """Inherit docstrings from parent class"""
    from inspect import getmembers, isfunction

    for name, func in getmembers(cls, isfunction):
        if func.__doc__:
            continue
        for parent in cls.__mro__[1:]:
            if hasattr(parent, name):
                func.__doc__ = getattr(parent, name).__doc__
    return cls


def all_subclasses(cls):
    """Recursively retrieve all subclasses of the specified class"""
    return cls.__subclasses__() + [
        g for s in cls.__subclasses__() for g in all_subclasses(s)
    ]


@contextlib.contextmanager
def temporary_cd(dir_path):
    """Context to temporary change the working directory.

    Parameters
    ----------

    dir_path : str
        The directory path to enter within the context

    Examples
    --------

    >>> dir_path = '/tmp'
    >>> with temporary_cd(dir_path):
    ...     pass  # do something in dir_path

    """
    import os

    prev_dir = os.getcwd()
    os.chdir(os.path.abspath(dir_path))
    try:
        yield
    finally:
        os.chdir(prev_dir)


def get_data_file_path(relative_path):
    """Get the full path to one of the reference files in testsystems.
    In the source distribution, these files are in ``openff/toolkit/data/``,
    but on installation, they're moved to somewhere in the user's python
    site-packages directory.

    Parameters
    ----------

    name : str
        Name of the file to load (with respect to the repex folder).

    """

    import os

    from pkg_resources import resource_filename

    fn = resource_filename("openff.toolkit", os.path.join("data", relative_path))

    if not os.path.exists(fn):
        raise ValueError(
            f"Sorry! {fn} does not exist. If you just added it, you'll have to re-install"
        )

    return fn


def unit_to_string(input_unit):
    return str(input_unit)


def _unit_to_string(input_unit):
    """
    Serialize a openmm.unit.Unit and return it as a string.

    Parameters
    ----------
    input_unit : A openmm.unit.Unit
        The Unit object to serialize

    Returns
    -------
    unit_string : str
        The serialized unit.
    """

    if input_unit == openmm_unit.dimensionless:
        return "dimensionless"

    # Decompose output_unit into a tuples of (base_dimension_unit, exponent)
    unit_string = None

    for unit_component in input_unit.iter_base_or_scaled_units():
        unit_component_name = unit_component[0].name
        # Convert, for example "elementary charge" --> "elementary_charge"
        unit_component_name = unit_component_name.replace(" ", "_")
        if unit_component[1] == 1:
            contribution = "{}".format(unit_component_name)
        else:
            contribution = "{}**{}".format(unit_component_name, int(unit_component[1]))
        if unit_string is None:
            unit_string = contribution
        else:
            unit_string += " * {}".format(contribution)

    return unit_string


def quantity_to_dict(input_quantity):
    return {
        "value": input_quantity.magnitude,
        "unit": str(input_quantity.units),
    }


def dict_to_quantity(input_dict):
    return input_dict["value"] * unit.Unit(input_dict["unit"])


def quantity_to_string(input_quantity):
    return str(input_quantity)


def _quantity_to_string(input_quantity):
    """
    Serialize a openmm.unit.Quantity to a string.

    Parameters
    ----------
    input_quantity : openmm.unit.Quantity
        The quantity to serialize

    Returns
    -------
    output_string : str
        The serialized quantity

    """
    import numpy as np

    if input_quantity is None:
        return None
    unitless_value = input_quantity.value_in_unit(input_quantity.unit)
    # The string representation of a numpy array doesn't have commas and breaks the
    # parser, thus we convert any arrays to list here
    if isinstance(unitless_value, np.ndarray):
        unitless_value = list(unitless_value)
    unit_string = unit_to_string(input_quantity.unit)
    output_string = "{} * {}".format(unitless_value, unit_string)
    return output_string


def _ast_eval(node):
    """
    Performs an algebraic syntax tree evaluation of a unit.

    Parameters
    ----------
    node : An ast parsing tree node
    """
    import ast
    import operator as op

    operators = {
        ast.Add: op.add,
        ast.Sub: op.sub,
        ast.Mult: op.mul,
        ast.Div: op.truediv,
        ast.Pow: op.pow,
        ast.BitXor: op.xor,
        ast.USub: op.neg,
    }

    if isinstance(node, ast.Num):  # <number>
        return node.n
    elif isinstance(node, ast.BinOp):  # <left> <operator> <right>
        return operators[type(node.op)](_ast_eval(node.left), _ast_eval(node.right))
    elif isinstance(node, ast.UnaryOp):  # <operator> <operand> e.g., -1
        return operators[type(node.op)](_ast_eval(node.operand))
    elif isinstance(node, ast.Name):
        # see if this is a openmm unit
        b = getattr(unit, node.id)
        return b
    # TODO: This was a quick hack that surprisingly worked. We should validate this further.
    elif isinstance(node, ast.List):
        return ast.literal_eval(node)
    else:
        raise TypeError(node)


def string_to_unit(unit_string):
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
    """
    return unit.Unit(unit_string)
    import ast

    output_unit = _ast_eval(ast.parse(unit_string, mode="eval").body)
    return output_unit

    # if (serialized['unitless_value'] is None) and (serialized['unit'] is None):
    #    return None
    # quantity_unit = None
    # for unit_name, power in serialized['unit']:
    #    unit_name = unit_name.replace(
    #        ' ', '_')  # Convert eg. 'elementary charge' to 'elementary_charge'
    #    if quantity_unit is None:
    #        quantity_unit = (getattr(unit, unit_name)**power)
    #    else:
    #        quantity_unit *= (getattr(unit, unit_name)**power)
    # quantity = unit.Quantity(serialized['unitless_value'], quantity_unit)
    # return quantity


def string_to_quantity(quantity_string):
    from tokenize import TokenError

    import pint

    try:
        return unit.Quantity(quantity_string)
    except (pint.errors.UndefinedUnitError, TokenError):
        return quantity_string


def _string_to_quantity(quantity_string):
    """
    Takes a string representation of a quantity and returns a openmm.unit.Quantity

    Parameters
    ----------
    quantity_string : str
        The quantity to deserialize

    Returns
    -------
    output_quantity : openmm.unit.Quantity
        The deserialized quantity
    """
    if quantity_string is None:
        return None
    # This can be the exact same as string_to_unit
    import ast

    output_quantity = _ast_eval(ast.parse(quantity_string, mode="eval").body)
    return output_quantity


def convert_all_strings_to_quantity(smirnoff_data):
    """
    Traverses a SMIRNOFF data structure, attempting to convert all
    quantity-defining strings into openmm.unit.Quantity objects.

    Integers and floats are ignored and not converted into a dimensionless
    ``openmm.unit.Quantity`` object.

    Parameters
    ----------
    smirnoff_data : dict
        A hierarchical dict structured in compliance with the SMIRNOFF spec

    Returns
    -------
    converted_smirnoff_data : dict
        A hierarchical dict structured in compliance with the SMIRNOFF spec,
        with quantity-defining strings converted to openmm.unit.Quantity objects
    """
    if isinstance(smirnoff_data, dict):
        for key, value in smirnoff_data.items():
            smirnoff_data[key] = convert_all_strings_to_quantity(value)
        obj_to_return = smirnoff_data

    elif isinstance(smirnoff_data, list):
        for index, item in enumerate(smirnoff_data):
            smirnoff_data[index] = convert_all_strings_to_quantity(item)
        obj_to_return = smirnoff_data

    elif isinstance(smirnoff_data, int) or isinstance(smirnoff_data, float):
        obj_to_return = smirnoff_data

    else:
        try:
            obj_to_return = object_to_quantity(smirnoff_data)
        except (AttributeError, TypeError, SyntaxError):
            obj_to_return = smirnoff_data

    return obj_to_return


def convert_all_quantities_to_string(smirnoff_data):
    """
    Traverses a SMIRNOFF data structure, attempting to convert all
    quantities into strings.

    Parameters
    ----------
    smirnoff_data : dict
        A hierarchical dict structured in compliance with the SMIRNOFF spec

    Returns
    -------
    converted_smirnoff_data : dict
        A hierarchical dict structured in compliance with the SMIRNOFF spec,
        with openmm.unit.Quantitys converted to string
    """

    if isinstance(smirnoff_data, dict):
        for key, value in smirnoff_data.items():
            smirnoff_data[key] = convert_all_quantities_to_string(value)
        obj_to_return = smirnoff_data
    elif isinstance(smirnoff_data, list):
        for index, item in enumerate(smirnoff_data):
            smirnoff_data[index] = convert_all_quantities_to_string(item)
        obj_to_return = smirnoff_data
    elif isinstance(smirnoff_data, openmm_unit.Quantity):
        obj_to_return = quantity_to_string(smirnoff_data)
    else:
        obj_to_return = smirnoff_data

    return obj_to_return


@functools.singledispatch
def object_to_quantity(object):
    """
    Attempts to turn the provided object into openmm.unit.Quantity(s).

    Can handle float, int, strings, quantities, or iterators over
    the same. Raises an exception if unable to convert all inputs.

    Parameters
    ----------
    object : int, float, string, quantity, or iterator of strings of quantities
        The object to convert to a ``openmm.unit.Quantity`` object.

    Returns
    -------
    converted_object : openmm.unit.Quantity or List[openmm.unit.Quantity]

    """
    # If we can't find a custom type, we treat this as a generic iterator.
    return [object_to_quantity(sub_obj) for sub_obj in object]


@object_to_quantity.register(openmm_unit.Quantity)
def _(obj):
    from openff.units.openmm import from_openmm

    return from_openmm(obj)


@object_to_quantity.register(unit.Quantity)
def _(obj):
    return obj


@object_to_quantity.register(str)
def _(obj):
    import pint

    try:
        return string_to_quantity(obj)
    except pint.errors.UndefinedUnitError:
        raise ValueError


@object_to_quantity.register(int)
@object_to_quantity.register(float)
def _(obj):
    return unit.Quantity(obj)


def check_units_are_compatible(object_name, object, unit_to_check, context=None):
    """
    Checks whether a openmm.unit.Quantity or list of openmm.unit.Quantitys is compatible with given unit.

    Parameters
    ----------
    object_name : string
        Name of object, used in printing exception.
    object : A openmm.unit.Quantity or list of openmm.unit.Quantitys
    unit_to_check : A openmm.unit.Unit
    context : string, optional. Default=None
        Additional information to provide at the beginning of the exception message if raised

    Raises
    ------
    IncompatibleUnitError
    """

    # If context is not provided, explicitly make it a blank string
    if context is None:
        context = ""
    # Otherwise add a space after the end of it to correct message printing
    else:
        context += " "

    if isinstance(object, list):
        for sub_object in object:
            check_units_are_compatible(
                object_name, sub_object, unit_to_check, context=context
            )
    elif isinstance(object, openmm_unit.Quantity):
        if not object.unit.is_compatible(unit_to_check):
            msg = (
                f"{context}{object_name} with "
                f"value {object} is incompatible with expected unit {unit_to_check}"
            )
            raise IncompatibleUnitError(msg)
    else:
        msg = (
            f"{context}{object_name} with "
            f"value {object} is incompatible with expected unit {unit_to_check}"
        )
        raise IncompatibleUnitError(msg)


def extract_serialized_units_from_dict(input_dict):
    """
    Create a mapping of (potentially unit-bearing) quantities from a dictionary, where some keys exist in pairs like
    {'length': 8, 'length_unit':'angstrom'}.

    Parameters
    ----------
    input_dict : dict
       Dictionary where some keys are paired like {'X': 1.0, 'X_unit': angstrom}.

    Returns
    -------
    unitless_dict : dict
       input_dict, but with keys ending in ``_unit`` removed.
    attached_units : dict str : openmm.unit.Unit
       ``attached_units[parameter_name]`` is the openmm.unit.Unit combination that should be attached to corresponding
       parameter ``parameter_name``. For example ``attached_units['X'] = openmm.unit.angstrom.

    """

    # TODO: Should this scheme also convert "1" to int(1) and "8.0" to float(8.0)?
    from collections import OrderedDict

    attached_units = OrderedDict()
    unitless_dict = input_dict.copy()
    keys_to_delete = []
    for key in input_dict.keys():
        if key.endswith("_unit"):
            parameter_name = key[:-5]
            parameter_units_string = input_dict[key]
            try:
                parameter_units = string_to_unit(parameter_units_string)
            except Exception as e:
                e.msg = (
                    "Could not parse units {}\n".format(parameter_units_string) + e.msg
                )
                raise e
            attached_units[parameter_name] = parameter_units
            # Remember this key and delete it later (we break the dict if we delete a key in the loop)
            keys_to_delete.append(key)
    # Clean out the '*_unit' keys that we processed
    for key in keys_to_delete:
        del unitless_dict[key]

    return unitless_dict, attached_units


def attach_units(unitless_dict, attached_units):
    """
    Attach units to dict entries for which units are specified.

    Parameters
    ----------
    unitless_dict : dict
       Dictionary, where some items are to have units applied.
    attached_units : dict [str : openmm.unit.Unit]
       ``attached_units[parameter_name]`` is the openmm.unit.Unit combination that should be attached to corresponding
       parameter ``parameter_name``

    Returns
    -------
    unit_bearing_dict : dict
       Updated dict with openmm.unit.Unit units attached to values for which units were specified for their keys

    """
    temp_dict = unitless_dict.copy()
    for parameter_name, units_to_attach in attached_units.items():
        if parameter_name in temp_dict.keys():
            parameter_attrib_string = temp_dict[parameter_name]
            try:
                temp_dict[parameter_name] = (
                    float(parameter_attrib_string) * units_to_attach
                )
            except ValueError as e:
                e.msg = (
                    "Expected numeric value for parameter '{}',"
                    "instead found '{}' when trying to attach units '{}'\n"
                ).format(parameter_name, parameter_attrib_string, units_to_attach)
                raise e

        # Now check for matches like "phase1", "phase2"
        c = 1
        while (parameter_name + str(c)) in temp_dict.keys():
            indexed_parameter_name = parameter_name + str(c)
            parameter_attrib_string = temp_dict[indexed_parameter_name]
            try:
                temp_dict[indexed_parameter_name] = (
                    float(parameter_attrib_string) * units_to_attach
                )
            except ValueError as e:
                e.msg = "Expected numeric value for parameter '{}', instead found '{}' when trying to attach units '{}'\n".format(
                    indexed_parameter_name, parameter_attrib_string, units_to_attach
                )
                raise e
            c += 1
    return temp_dict


def detach_units(unit_bearing_dict, output_units=None):
    """
    Given a dict which may contain some openmm.unit.Quantity objects, return the same dict with the Quantities
    replaced with unitless values, and a new dict containing entries with the suffix "_unit" added, containing
    the units.

    Parameters
    ----------
    unit_bearing_dict : dict
        A dictionary potentially containing openmm.unit.Quantity objects as values.
    output_units : dict[str : openmm.unit.Unit], optional. Default = None
        A mapping from parameter fields to the output unit its value should be converted to.
        For example, {'length_unit': unit.angstrom}. If no output_unit is defined for a key:value pair in which
        the value is a openmm.unit.Quantity, the output unit will be the Quantity's unit, and this information
        will be included in the unit_dict return value.

    Returns
    -------
    unitless_dict : dict
        The input smirnoff_dict object, with all openmm.unit.Quantity values converted to unitless values.
    unit_dict : dict
        A dictionary in which keys are keys of openmm.unit.Quantity values in unit_bearing_dict,
        but suffixed with "_unit". Values are openmm.unit.Unit .
    """

    if output_units is None:
        output_units = {}

    # initialize dictionaries for outputs
    unit_dict = {}
    unitless_dict = unit_bearing_dict.copy()

    for key, value in unit_bearing_dict.items():
        # If no conversion is needed, skip this item
        if not isinstance(value, unit.Quantity):
            continue

        # If conversion is needed, see if the user has requested an output unit
        unit_key = key + "_unit"

        if unit_key in output_units:
            output_unit = output_units[unit_key]
        else:
            output_unit = value.units
        if not (output_unit.is_compatible_with(value.units)):
            raise ValueError(
                "Requested output unit {} is not compatible with "
                "quantity unit {}.".format(output_unit, value.units)
            )
        unitless_dict[key] = value.m_as(output_unit)
        unit_dict[unit_key] = output_unit

    return unitless_dict, unit_dict


def serialize_numpy(np_array):
    """
    Serializes a numpy array into a JSON-compatible string. Leverages the numpy.save function,
    thereby preserving the shape of the input array

    from https://stackoverflow.com/questions/30698004/how-can-i-serialize-a-numpy-array-while-preserving-matrix-dimensions#30699208

    Parameters
    ----------
    np_array : A numpy array
        Input numpy array

    Returns
    -------
    serialized : str
        A serialized representation of the numpy array.
    shape : tuple of ints
        The shape of the serialized array
    """

    bigendian_array = np_array.newbyteorder(">")
    serialized = bigendian_array.tobytes()
    shape = np_array.shape
    return serialized, shape


def deserialize_numpy(serialized_np, shape):
    """
    Deserializes a numpy array from a JSON-compatible string.

    from https://stackoverflow.com/questions/30698004/how-can-i-serialize-a-numpy-array-while-preserving-matrix-dimensions#30699208

    Parameters
    ----------
    serialized_np : bytes or list
        A byte or list serialized representation of a numpy array
    shape : tuple of ints
        The shape of the serialized array
    Returns
    -------
    np_array : numpy.ndarray
        The deserialized numpy array
    """

    import numpy as np

    if isinstance(serialized_np, list):
        np_array = np.array(serialized_np)
    if isinstance(serialized_np, bytes):
        dt = np.dtype("float")
        dt.newbyteorder(">")  # set to big-endian
        np_array = np.frombuffer(serialized_np, dtype=dt)
    np_array = np_array.reshape(shape)
    return np_array


def convert_0_2_smirnoff_to_0_3(smirnoff_data_0_2):
    """
    Convert an 0.2-compliant SMIRNOFF dict to an 0.3-compliant one.
    This involves removing units from header tags and adding them
    to attributes of child elements.
    It also requires converting ProperTorsions and ImproperTorsions
    potentials from "charmm" to "fourier".

    Parameters
    ----------
    smirnoff_data_0_2 : dict
        Hierarchical dict representing a SMIRNOFF data structure according the the 0.2 spec

    Returns
    -------
    smirnoff_data_0_3
        Hierarchical dict representing a SMIRNOFF data structure according the the 0.3 spec
    """
    # Legacy force fields sometimes specify the NonbondedForce's sigma_unit value, but then provide
    # atom size as rmin_half. Here we correct for this behavior by explicitly defining both as
    # the same unit if either one is defined.
    if "vdW" in smirnoff_data_0_2["SMIRNOFF"].keys():
        rmh_unit = smirnoff_data_0_2["SMIRNOFF"]["vdW"].get("rmin_half_unit", None)
        sig_unit = smirnoff_data_0_2["SMIRNOFF"]["vdW"].get("sigma_unit", None)
        if (rmh_unit is not None) and (sig_unit is None):
            smirnoff_data_0_2["SMIRNOFF"]["vdW"]["sigma_unit"] = rmh_unit
        elif (sig_unit is not None) and (rmh_unit is None):
            smirnoff_data_0_2["SMIRNOFF"]["vdW"]["rmin_half_unit"] = sig_unit
        # If both are None, or both are defined, don't overwrite anything
        else:
            pass

    # Recursively attach unit strings
    smirnoff_data = recursive_attach_unit_strings(smirnoff_data_0_2, {})

    # Change TorsionHandler potential from "charmm" to "k*(1+cos(periodicity*theta-phase))". Note that, scientifically,
    # we should have used "k*(1+cos(periodicity*theta-phase))" all along, since "charmm" technically
    # implies that we would support a harmonic potential for torsion terms with periodicity 0
    # More at: https://github.com/openforcefield/openff-toolkit/issues/303#issuecomment-490156779
    if "ProperTorsions" in smirnoff_data["SMIRNOFF"]:
        if "potential" in smirnoff_data["SMIRNOFF"]["ProperTorsions"]:
            if smirnoff_data["SMIRNOFF"]["ProperTorsions"]["potential"] == "charmm":
                smirnoff_data["SMIRNOFF"]["ProperTorsions"][
                    "potential"
                ] = "k*(1+cos(periodicity*theta-phase))"
    if "ImproperTorsions" in smirnoff_data["SMIRNOFF"]:
        if "potential" in smirnoff_data["SMIRNOFF"]["ImproperTorsions"]:
            if smirnoff_data["SMIRNOFF"]["ImproperTorsions"]["potential"] == "charmm":
                smirnoff_data["SMIRNOFF"]["ImproperTorsions"][
                    "potential"
                ] = "k*(1+cos(periodicity*theta-phase))"

    # Add per-section tag
    sections_not_to_version_0_3 = ["Author", "Date", "version", "aromaticity_model"]
    for l1_tag in smirnoff_data["SMIRNOFF"].keys():
        if l1_tag not in sections_not_to_version_0_3:

            if smirnoff_data["SMIRNOFF"][l1_tag] is None:
                # Handle empty entries, such as the ToolkitAM1BCC handler.
                smirnoff_data["SMIRNOFF"][l1_tag] = {}

            smirnoff_data["SMIRNOFF"][l1_tag]["version"] = 0.3

    # Update top-level tag
    smirnoff_data["SMIRNOFF"]["version"] = 0.3

    return smirnoff_data


def convert_0_1_smirnoff_to_0_2(smirnoff_data_0_1):
    """
    Convert an 0.1-compliant SMIRNOFF dict to an 0.2-compliant one.
    This involves renaming several tags, adding Electrostatics and ToolkitAM1BCC tags, and
    separating improper torsions into their own section.

    Parameters
    ----------
    smirnoff_data_0_1 : dict
        Hierarchical dict representing a SMIRNOFF data structure according the the 0.1 spec

    Returns
    -------
    smirnoff_data_0_2
        Hierarchical dict representing a SMIRNOFF data structure according the the 0.2 spec
    """
    smirnoff_data = smirnoff_data_0_1.copy()

    l0_replacement_dict = {"SMIRFF": "SMIRNOFF"}
    l1_replacement_dict = {
        "HarmonicBondForce": "Bonds",
        "HarmonicAngleForce": "Angles",
        "PeriodicTorsionForce": "ProperTorsions",
        "NonbondedForce": "vdW",
    }
    for old_l0_tag, new_l0_tag in l0_replacement_dict.items():
        # Convert first-level smirnoff_data tags.
        # Right now this just changes the SMIRFF tag to SMIRNOFF
        if old_l0_tag in smirnoff_data.keys():
            smirnoff_data[new_l0_tag] = smirnoff_data[old_l0_tag]
            del smirnoff_data[old_l0_tag]

    # SMIRFF tag will have been converted to SMIRNOFF here
    # Convert second-level tags here
    for old_l1_tag, new_l1_tag in l1_replacement_dict.items():
        if old_l1_tag in smirnoff_data["SMIRNOFF"].keys():
            smirnoff_data["SMIRNOFF"][new_l1_tag] = smirnoff_data["SMIRNOFF"][
                old_l1_tag
            ]
            del smirnoff_data["SMIRNOFF"][old_l1_tag]

    # Add 'potential' field to each l1 tag
    default_potential = {
        "Bonds": "harmonic",
        "Angles": "harmonic",
        "ProperTorsions": "charmm",
        # Note that "charmm" isn't actually correct, and was later changed
        # in the 0.3 spec. More info at
        # https://github.com/openforcefield/openff-toolkit/pull/311#commitcomment-33494506
        "vdW": "Lennard-Jones-12-6",
    }
    for l1_tag in smirnoff_data["SMIRNOFF"].keys():
        if l1_tag in default_potential.keys():
            # Ensure that it isn't there already (shouldn't happen, but better to be safe)
            if "potential" in smirnoff_data["SMIRNOFF"][l1_tag].keys():
                assert smirnoff_data[l1_tag].keys == default_potential[l1_tag]
                continue
            # Issue an informative warning about assumptions made during conversion.
            logger.warning(
                f"0.1 SMIRNOFF spec file does not contain 'potential' attribute for '{l1_tag}' tag. "
                f"The SMIRNOFF spec converter is assuming it has a value of '{default_potential[l1_tag]}'"
            )
            smirnoff_data["SMIRNOFF"][l1_tag]["potential"] = default_potential[l1_tag]

    # Separate improper torsions from propers
    if "ProperTorsions" in smirnoff_data["SMIRNOFF"]:
        if "Improper" in smirnoff_data["SMIRNOFF"]["ProperTorsions"]:
            # First generate an ImproperTorsions header, taking the relevant values from the ProperTorsions header
            improper_section = {
                "k_unit": smirnoff_data["SMIRNOFF"]["ProperTorsions"]["k_unit"],
                "phase_unit": smirnoff_data["SMIRNOFF"]["ProperTorsions"]["phase_unit"],
                "potential": smirnoff_data["SMIRNOFF"]["ProperTorsions"]["potential"],
                "Improper": smirnoff_data["SMIRNOFF"]["ProperTorsions"]["Improper"],
            }

            # Then, attach the newly-made ImproperTorsions section
            smirnoff_data["SMIRNOFF"]["ImproperTorsions"] = improper_section
            del smirnoff_data["SMIRNOFF"]["ProperTorsions"]["Improper"]

    # Add Electrostatics tag, setting several values to their defaults and
    # warning about assumptions that are being made
    electrostatics_section = {
        "method": "PME",
        "scale12": 0.0,
        "scale13": 0.0,
        "scale15": 1.0,
        "cutoff": 9.0,
        "cutoff_unit": "angstrom",
    }
    logger.warning(
        "0.1 SMIRNOFF spec did not allow the 'Electrostatics' tag. Adding it in 0.2 spec conversion, and "
        "assuming the following values:"
    )
    for key, val in electrostatics_section.items():
        logger.warning(f"\t{key}: {val}")

    # Take electrostatics 1-4 scaling term from 0.1 spec's NonBondedForce tag
    electrostatics_section["scale14"] = smirnoff_data["SMIRNOFF"]["vdW"][
        "coulomb14scale"
    ]
    del smirnoff_data["SMIRNOFF"]["vdW"]["coulomb14scale"]
    smirnoff_data["SMIRNOFF"]["Electrostatics"] = electrostatics_section

    # Change vdW's lj14scale to 14scale, add other scaling terms
    vdw_section_additions = {
        "method": "cutoff",
        "combining_rules": "Lorentz-Berthelot",
        "scale12": "0.0",
        "scale13": "0.0",
        "scale15": "1",
        "switch_width": "1.0",
        "switch_width_unit": "angstrom",
        "cutoff": "9.0",
        "cutoff_unit": "angstrom",
    }
    for key, val in vdw_section_additions.items():
        if key not in smirnoff_data["SMIRNOFF"]["vdW"].keys():
            logger.warning(
                f"0.1 SMIRNOFF spec file does not contain '{key}' attribute for 'NonBondedMethod/vdW'' tag. "
                f"The SMIRNOFF spec converter is assuming it has a value of '{val}'"
            )
            smirnoff_data["SMIRNOFF"]["vdW"][key] = val

    # Rename L-J 1-4 scaling term from 0.1 spec's NonBondedForce tag to vdW's scale14
    smirnoff_data["SMIRNOFF"]["vdW"]["scale14"] = smirnoff_data["SMIRNOFF"]["vdW"][
        "lj14scale"
    ]
    del smirnoff_data["SMIRNOFF"]["vdW"]["lj14scale"]

    # Add <ToolkitAM1BCC/> tag
    smirnoff_data["SMIRNOFF"]["ToolkitAM1BCC"] = {}

    # Update top-level tag
    smirnoff_data["SMIRNOFF"]["version"] = 0.2

    return smirnoff_data


def recursive_attach_unit_strings(smirnoff_data, units_to_attach):
    """
    Recursively traverse a SMIRNOFF data structure, appending "* {unit}" to values in key:value pairs
    where "key_unit":"unit_string" is present at a higher level in the hierarchy.
    This function expects all items in smirnoff_data to be formatted as strings.

    Parameters
    ----------
    smirnoff_data : dict
        Any level of hierarchy that is part of a SMIRNOFF dict, with all data members
        formatted as string.
    units_to_attach : dict
        Dict of the form {key:unit_string}

    Returns
    -------
    unit_appended_smirnoff_data: dict
    """
    import re

    # Make a copy of units_to_attach so we don't modify the original (otherwise things like k_unit could
    # leak between sections)
    units_to_attach = units_to_attach.copy()
    # smirnoff_data = smirnoff_data.copy()

    # If we're working with a dict, see if there are any new unit entries and store them,
    # then operate recursively on the values in the dict.
    if isinstance(smirnoff_data, dict):

        # Go over all key:value pairs once to see if there are new units to attach.
        # Note that units to be attached can be defined in the same dict as the
        # key:value pair they will be attached to, so we need to complete this check
        # before we are able to check other items in the dict.
        for key, value in list(smirnoff_data.items()):
            if key[-5:] == "_unit":
                units_to_attach[key[:-5]] = value
                del smirnoff_data[key]

        # Go through once more to attach units as appropriate
        for key in smirnoff_data.keys():

            # We use regular expressions to catch possible indexed attributes
            attach_unit = None
            for unit_key, unit_string in units_to_attach.items():
                if re.match(f"{unit_key}[0-9]*", key):
                    attach_unit = unit_string

            if attach_unit is not None:
                smirnoff_data[key] = str(smirnoff_data[key]) + " * " + attach_unit

            # And recursively act on value, in case it's a deeper level of hierarchy
            smirnoff_data[key] = recursive_attach_unit_strings(
                smirnoff_data[key], units_to_attach
            )

    # If it's a list, operate on each member of the list
    elif isinstance(smirnoff_data, list):
        for index, value in enumerate(smirnoff_data):
            smirnoff_data[index] = recursive_attach_unit_strings(value, units_to_attach)

    # Otherwise, just return smirnoff_data unchanged
    else:
        pass

    return smirnoff_data


def get_molecule_parameterIDs(molecules, forcefield):
    """Process a list of molecules with a specified SMIRNOFF ffxml file and determine which parameters are used by
    which molecules, returning collated results.

    Parameters
    ----------
    molecules : list of openff.toolkit.topology.Molecule
        List of molecules (with explicit hydrogens) to parse
    forcefield : openff.toolkit.typing.engines.smirnoff.ForceField
        The ForceField to apply

    Returns
    -------
    parameters_by_molecule : dict
        Parameter IDs used in each molecule, keyed by isomeric SMILES
        generated from provided OEMols. Each entry in the dict is a list
        which does not necessarily have unique entries; i.e. parameter IDs
        which are used more than once will occur multiple times.

    parameters_by_ID : dict
        Molecules in which each parameter ID occur, keyed by parameter ID.
        Each entry in the dict is a set of isomeric SMILES for molecules
        in which that parameter occurs. No frequency information is stored.

    """

    from openff.toolkit.topology import Topology

    # Create storage
    parameters_by_molecule = dict()
    parameters_by_ID = dict()

    # Generate isomeric SMILES for each molecule, ensuring all molecules are unique
    isosmiles = [molecule.to_smiles() for molecule in molecules]
    already_seen = set()
    duplicates = set(
        smiles
        for smiles in isosmiles
        if smiles in already_seen or already_seen.add(smiles)
    )
    if len(duplicates) > 0:
        raise ValueError(
            "Error: get_molecule_parameterIDs has been provided a list of oemols which contains some duplicates: {}".format(
                duplicates
            )
        )

    # Assemble molecules into a Topology
    topology = Topology()
    for molecule in molecules:
        topology.add_molecule(molecule)

    # Label molecules
    labels = forcefield.label_molecules(topology)

    # Organize labels into output dictionary by looping over all molecules/smiles
    for idx in range(len(isosmiles)):
        # Pull smiles, initialize storage
        smi = isosmiles[idx]
        parameters_by_molecule[smi] = []

        # Organize data for this molecule
        data = labels[idx]
        for force_type in data.keys():
            for atom_indices, parameter_type in data[force_type].items():

                pid = parameter_type.id
                # Store pid to molecule
                parameters_by_molecule[smi].append(pid)

                # Store which molecule this pid occurred in
                if pid not in parameters_by_ID:
                    parameters_by_ID[pid] = set()
                    parameters_by_ID[pid].add(smi)
                else:
                    parameters_by_ID[pid].add(smi)

    return parameters_by_molecule, parameters_by_ID


def sort_smirnoff_dict(data):
    """
    Recursively sort the keys in a dict of SMIRNOFF data.

    Adapted from https://stackoverflow.com/a/47882384/4248961

    TODO: Should this live elsewhere?
    """
    sorted_dict = dict()
    for key, val in sorted(data.items()):
        if isinstance(val, dict):
            # This should hit each ParameterHandler and dicts within them
            sorted_dict[key] = sort_smirnoff_dict(val)
        elif isinstance(val, list):
            # Handle case of ParameterLists, which show up in
            # the smirnoff dicts as lists of OrderedDicts
            new_parameter_list = list()
            for param in val:
                new_parameter_list.append(sort_smirnoff_dict(param))
            sorted_dict[key] = new_parameter_list
        else:
            # Handle metadata or the bottom of a recursive dict
            sorted_dict[key] = val
    return sorted_dict
