#!/usr/bin/env python
"""
Utility subroutines

"""
#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import contextlib
from simtk import unit

#=============================================================================================
# UTILITY SUBROUTINES
#=============================================================================================


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
    ...     # do something in dir_path

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


def get_data_filename(relative_path):
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


def serialize_quantity(quantity):
    """
    Serialized a simtk.unit.Quantity into a dict of {'unitless_value': X, 'unit': Y}

    Parameters
    ----------
    quantity : A simtk.unit.Quantity-wrapped value or iterator over values
        The object to serialize

    Returns
    -------
    serialzied : dict
        The serialized object
    """

    serialized = dict()

    # If it's None, just return None in all fields
    if quantity is None:
        serialized['unitless_value'] = None
        serialized['unit'] = None
        return serialized

    # If it's not None, make sure it's a simtk.unit.Quantity
    assert (hasattr(quantity, 'unit'))

    quantity_unit = list()
    for base_unit in quantity.unit.iter_all_base_units():
        quantity_unit.append((base_unit[0].name, base_unit[1]))

    unitless_value = quantity / quantity.unit
    serialized['unitless_value'] = unitless_value
    serialized['unit'] = quantity_unit
    return serialized


def deserialize_quantity(serialized):
    """
    Deserializes a simtk.unit.Quantity.

    Parameters
    ----------
    serialized : dict
        Serialized representation of a simtk.unit.Quantity. Must have keys ["unitless_value", "unit"]

    Returns
    -------
    simtk.unit.Quantity
    """
    if (serialized['unitless_value'] is None) and (serialized['unit'] is None):
        return None
    quantity_unit = None
    for unit_name, power in serialized['unit']:
        unit_name = unit_name.replace(
            ' ', '_')  # Convert eg. 'elementary charge' to 'elementary_charge'
        if quantity_unit is None:
            quantity_unit = (getattr(unit, unit_name)**power)
        else:
            quantity_unit *= (getattr(unit, unit_name)**power)
    quantity = unit.Quantity(serialized['unitless_value'], quantity_unit)
    return quantity


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
    # This version ties us to numpy -- Making a more general solution
    #import io
    #import numpy
    #import json
    #memfile = io.BytesIO()
    #numpy.save(memfile, np_array)
    #memfile.seek(0)
    #serialized = json.dumps(memfile.read().decode('latin-1'))
    #dt = np.dtype('float')
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
    #import io
    #import numpy
    #import json
    #memfile = io.BytesIO()
    #memfile.write(json.loads(serialized_np).encode('latin-1'))
    #memfile.seek(0)
    #np_array = numpy.load(memfile)
    #return np_array
    import numpy as np
    dt = np.dtype('float')
    dt.newbyteorder('>')  # set to big-endian
    np_array = np.frombuffer(serialized_np, dtype=dt)
    np_array = np_array.reshape(shape)
    return np_array
