#!/usr/bin/env python
"""
Utility subroutines

"""

__all__ = [
    'MessageException',
    'inherit_docstrings',
    'all_subclasses',
    'temporary_cd',
    'temporary_directory',
    'get_data_file_path',
    'unit_to_string',
    'quantity_to_string',
    'string_to_unit',
    'string_to_quantity',
    'extract_serialized_units_from_dict',
    'attach_units',
    'detach_units',
    'serialize_numpy',
    'deserialize_numpy',
]

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import contextlib
from simtk import unit



#=============================================================================================
# COMMON EXCEPTION TYPES
#=============================================================================================


class MessageException(Exception):
    """A base class for exceptions that print out a string given in their constructor"""
    def __init__(self, msg):
        super().__init__(self, msg)
        self.msg = msg

# =============================================================================================
# UTILITY SUBROUTINES
# =============================================================================================


def inherit_docstrings(cls):
    """Inherit docstrings from parent class"""
    from inspect import getmembers, isfunction
    for name, func in getmembers(cls, isfunction):
        if func.__doc__: continue
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


@contextlib.contextmanager
def temporary_directory():
    """Context for safe creation of temporary directories."""

    import tempfile
    tmp_dir = tempfile.mkdtemp()
    try:
        yield tmp_dir
    finally:
        import shutil
        shutil.rmtree(tmp_dir)


def get_data_file_path(relative_path):
    """Get the full path to one of the reference files in testsystems.
    In the source distribution, these files are in ``openforcefield/data/``,
    but on installation, they're moved to somewhere in the user's python
    site-packages directory.
    Parameters
    ----------
    name : str
        Name of the file to load (with respect to the repex folder).
    """

    from pkg_resources import resource_filename
    import os
    fn = resource_filename('openforcefield', os.path.join(
        'data', relative_path))

    if not os.path.exists(fn):
        raise ValueError(
            "Sorry! %s does not exist. If you just added it, you'll have to re-install"
            % fn)

    return fn


def unit_to_string(input_unit):
    """
    Serialize a simtk.unit.Unit and return it as a string.

    Parameters
    ----------
    input_unit : A simtk.unit
        The unit to serialize

    Returns
    -------
    unit_string : str
        The serialized unit.
    """


    # Decompose output_unit into a tuples of (base_dimension_unit, exponent)
    unit_string = None

    for unit_component in input_unit.iter_base_or_scaled_units():
        unit_component_name = unit_component[0].name
        # Convert, for example "elementary charge" --> "elementary_charge"
        unit_component_name = unit_component_name.replace(' ', '_')
        if unit_component[1] == 1:
            contribution = '{}'.format(unit_component_name)
        else:
            contribution = '{}**{}'.format(unit_component_name, int(unit_component[1]))
        if unit_string is None:
            unit_string = contribution
        else:
            unit_string += ' * {}'.format(contribution)

    return unit_string

def quantity_to_string(input_quantity):
    """
    Serialize a simtk.unit.Quantity to a string.

    Parameters
    ----------
    input_quantity : simtk.unit.Quantity
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
    # The string representaiton of a numpy array doesn't have commas and breaks the
    # parser, thus we convert any arrays to list here
    if isinstance(unitless_value, np.ndarray):
        unitless_value = list(unitless_value)
    unit_string = unit_to_string(input_quantity.unit)
    output_string = '{} * {}'.format(unitless_value, unit_string)
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

    operators = {ast.Add: op.add, ast.Sub: op.sub, ast.Mult: op.mul,
        ast.Div: op.truediv, ast.Pow: op.pow, ast.BitXor: op.xor,
        ast.USub: op.neg}

    if isinstance(node, ast.Num): # <number>
        return node.n
    elif isinstance(node, ast.BinOp): # <left> <operator> <right>
        return operators[type(node.op)](_ast_eval(node.left), _ast_eval(node.right))
    elif isinstance(node, ast.UnaryOp): # <operator> <operand> e.g., -1
        return operators[type(node.op)](_ast_eval(node.operand))
    elif isinstance(node, ast.Name):
        # see if this is a simtk unit
        b = getattr(unit, node.id)
        return b
    # TODO: This was a quick hack that surprisingly worked. We should validate this further.
    elif isinstance(node, ast.List):
        return ast.literal_eval(node)
    else:
        raise TypeError(node)


def string_to_unit(unit_string):
    """
    Deserializes a simtk.unit.Quantity from a string representation, for
    example: "kilocalories_per_mole / angstrom ** 2"


    Parameters
    ----------
    unit_string : dict
        Serialized representation of a simtk.unit.Quantity.

    Returns
    -------
    output_unit: simtk.unit.Quantity
        The deserialized unit from the string
    """
    import ast
    output_unit = _ast_eval(ast.parse(unit_string, mode='eval').body)
    return output_unit

    #if (serialized['unitless_value'] is None) and (serialized['unit'] is None):
    #    return None
    #quantity_unit = None
    #for unit_name, power in serialized['unit']:
    #    unit_name = unit_name.replace(
    #        ' ', '_')  # Convert eg. 'elementary charge' to 'elementary_charge'
    #    if quantity_unit is None:
    #        quantity_unit = (getattr(unit, unit_name)**power)
    #    else:
    #        quantity_unit *= (getattr(unit, unit_name)**power)
    #quantity = unit.Quantity(serialized['unitless_value'], quantity_unit)
    #return quantity

def string_to_quantity(quantity_string):
    """
    Takes a string representation of a quantity and returns a simtk.unit.Quantity

    Parameters
    ----------
    quantity_string : str
        The quantity to deserialize

    Returns
    -------
    output_quantity : simtk.unit.Quantity
        The deserialized quantity
    """
    if quantity_string is None:
        return None
    # This can be the exact same as string_to_unit
    import ast
    output_quantity = _ast_eval(ast.parse(quantity_string, mode='eval').body)
    return output_quantity


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
    attached_units : dict str : simtk.unit.Unit
       ``attached_units[parameter_name]`` is the simtk.unit.Unit combination that should be attached to corresponding
       parameter ``parameter_name``. For example ``attached_units['X'] = simtk.unit.angstrom.

    """

    # TODO: Should this scheme also convert "1" to int(1) and "8.0" to float(8.0)?
    from collections import OrderedDict
    attached_units = OrderedDict()
    unitless_dict = input_dict.copy()
    keys_to_delete = []
    for key in input_dict.keys():
        if key.endswith('_unit'):
            parameter_name = key[:-5]
            parameter_units_string = input_dict[key]
            try:
                parameter_units = string_to_unit(parameter_units_string)
            except Exception as e:
                e.msg = "Could not parse units {}\n".format(
                    parameter_units_string) + e.msg
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
    attached_units : dict [str : simtk.unit.Unit]
       ``attached_units[parameter_name]`` is the simtk.unit.Unit combination that should be attached to corresponding
       parameter ``parameter_name``

    Returns
    -------
    unit_bearing_dict : dict
       Updated dict with simtk.unit.Unit units attached to values for which units were specified for their keys

    """
    temp_dict = unitless_dict.copy()
    for parameter_name, units_to_attach in attached_units.items():
        if parameter_name in temp_dict.keys():
            parameter_attrib_string = temp_dict[parameter_name]
            try:
                temp_dict[parameter_name] = float(parameter_attrib_string) * units_to_attach
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
                temp_dict[indexed_parameter_name] = float(
                          parameter_attrib_string) * units_to_attach
            except ValueError as e:
                e.msg = "Expected numeric value for parameter '{}', instead found '{}' when trying to attach units '{}'\n".format(
                    indexed_parameter_name, parameter_attrib_string,
                    units_to_attach)
                raise e
            c += 1
    return temp_dict



def detach_units(unit_bearing_dict, output_units=None):
    """
    Given a dict which may contain some simtk.unit.Quantity objects, return the same dict with the Quantities
    replaced with unitless values, and a new dict containing entries with the suffix "_unit" added, containing
    the units.

    Parameters
    ----------
    unit_bearing_dict : dict
        A dictionary potentially containing simtk.unit.Quantity objects as values.
    output_units : dict[str : simtk.unit.Unit], optional. Default = None
        A mapping from parameter fields to the output unit its value should be converted to.
        For example, {'length_unit': unit.angstrom}. If no output_unit is defined for a key:value pair in which
        the value is a simtk.unit.Quantity, the output unit will be the Quantity's unit, and this information
        will be included in the unit_dict return value.

    Returns
    -------
    unitless_dict : dict
        The input smirnoff_dict object, with all simtk.unit.Quantity values converted to unitless values.
    unit_dict : dict
        A dictionary in which keys are keys of simtk.unit.Quantity values in unit_bearing_dict,
        but suffixed with "_unit". Values are simtk.unit.Unit .
    """
    from simtk import unit

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
        unit_key = key + '_unit'

        if unit_key in output_units:
            output_unit = output_units[unit_key]
        else:
            output_unit = value.unit
        if not (output_unit.is_compatible(value.unit)):
            raise ValueError("Requested output unit {} is not compatible with "
                             "quantity unit {} .".format(output_unit, value.unit))
        unitless_dict[key] = value.value_in_unit(output_unit)
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

    bigendian_array = np_array.newbyteorder('>')
    serialized = bigendian_array.tobytes()
    shape = np_array.shape
    return serialized, shape


def deserialize_numpy(serialized_np, shape):
    """
    Deserializes a numpy array from a JSON-compatible string.

    from https://stackoverflow.com/questions/30698004/how-can-i-serialize-a-numpy-array-while-preserving-matrix-dimensions#30699208

    Parameters
    ----------
    serialized_np : str
        A serialized numpy array
    shape : tuple of ints
        The shape of the serialized array
    Returns
    -------
    np_array : numpy.ndarray
        The deserialized numpy array
    """

    import numpy as np
    dt = np.dtype('float')
    dt.newbyteorder('>')  # set to big-endian
    np_array = np.frombuffer(serialized_np, dtype=dt)
    np_array = np_array.reshape(shape)
    return np_array
