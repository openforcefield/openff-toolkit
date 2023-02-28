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
from typing import Dict, Iterable, List, Tuple, Union

import numpy as np
import pint
from openff.units import unit
from openff.utilities import requires_package

logger = logging.getLogger(__name__)


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


@pint.register_unit_format("simple")
def format_unit_simple(unit, registry, **options):
    return " * ".join(f"{u} ** {p}" for u, p in unit.items())


def unit_to_string(input_unit: unit.Unit) -> str:
    return f"{input_unit:simple}"


def quantity_to_dict(input_quantity):
    value = input_quantity.magnitude
    if isinstance(value, np.ndarray):
        value = value.tolist()

    return {
        "value": value,
        "unit": str(input_quantity.units),
    }


def dict_to_quantity(input_dict):
    return input_dict["value"] * unit.Unit(input_dict["unit"])


def quantity_to_string(input_quantity: unit.Quantity) -> str:
    """
    Serialize a openff.units.unit.Quantity to a string representation that is backwards-compatible
    with older versions of the OpenFF Toolkit. This includes a " * " between numerical values and
    their units and "A" being used in place of the unicode Å ("\N{ANGSTROM SIGN}").

    Parameters
    ----------
    input_quantity : openff.units.unit.Quantity
        The quantity to serialize

    Returns
    -------
    output_string : str
        The serialized quantity

    """
    unitless_value = input_quantity.m_as(input_quantity.units)
    # The string representation of a numpy array doesn't have commas and breaks the
    # parser, thus we convert any arrays to list here
    if isinstance(unitless_value, np.ndarray):
        unitless_value = list(unitless_value)
    unit_string = unit_to_string(input_quantity.units)
    output_string = "{} * {}".format(unitless_value, unit_string)
    return output_string

    return str(input_quantity)


def string_to_unit(unit_string):
    """
    Deserializes a openff.units.unit.Quantity from a string representation, for
    example: "kilocalories_per_mole / angstrom ** 2"


    Parameters
    ----------
    unit_string : dict
        Serialized representation of a openff.units.unit.Quantity.

    Returns
    -------
    output_unit: openff.units.unit.Quantity
        The deserialized unit from the string
    """
    return unit.Unit(unit_string)


def string_to_quantity(quantity_string) -> Union[str, int, float, unit.Quantity]:
    """Attempt to parse a string into a unit.Quantity.

    Note that dimensionless floats and ints are returns as floats or ints, not Quantity objects.
    """

    from tokenize import TokenError

    from pint import UndefinedUnitError

    try:
        quantity = unit.Quantity(quantity_string)
    except (TokenError, UndefinedUnitError):
        return quantity_string

    # TODO: Should intentionally unitless array-likes be Quantity objects
    #       or their raw representation?
    if (quantity.units == unit.dimensionless) and isinstance(quantity.m, (int, float)):
        return quantity.m
    else:
        return quantity


def convert_all_strings_to_quantity(
    smirnoff_data: Dict,
    ignore_keys: Iterable[str] = tuple(),
):
    """
    Traverses a SMIRNOFF data structure, attempting to convert all
    quantity-defining strings into openff.units.unit.Quantity objects.

    Integers and floats are ignored and not converted into a dimensionless
    ``openff.units.unit.Quantity`` object.

    Parameters
    ----------
    smirnoff_data : dict
        A hierarchical dict structured in compliance with the SMIRNOFF spec
    ignore_keys : iterable of str, optional, default=tuple()
        A list of keys to skip when converting strings to quantities

    Returns
    -------
    converted_smirnoff_data : dict
        A hierarchical dict structured in compliance with the SMIRNOFF spec,
        with quantity-defining strings converted to openff.units.unit.Quantity objects
    """
    from pint import DefinitionSyntaxError

    if isinstance(smirnoff_data, dict):
        for key, value in smirnoff_data.items():
            if key in ignore_keys:
                smirnoff_data[key] = value
            else:
                smirnoff_data[key] = convert_all_strings_to_quantity(
                    value,
                    ignore_keys=ignore_keys,
                )
        obj_to_return = smirnoff_data

    elif isinstance(smirnoff_data, list):
        for index, item in enumerate(smirnoff_data):
            smirnoff_data[index] = convert_all_strings_to_quantity(
                item,
                ignore_keys=ignore_keys,
            )
        obj_to_return = smirnoff_data

    elif isinstance(smirnoff_data, int) or isinstance(smirnoff_data, float):
        obj_to_return = smirnoff_data

    else:
        try:
            obj_to_return = object_to_quantity(smirnoff_data)
        except (TypeError, DefinitionSyntaxError):
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
        with openff.units.unit.Quantitys converted to string
    """

    if isinstance(smirnoff_data, dict):
        for key, value in smirnoff_data.items():
            smirnoff_data[key] = convert_all_quantities_to_string(value)
        obj_to_return = smirnoff_data
    elif isinstance(smirnoff_data, list):
        for index, item in enumerate(smirnoff_data):
            smirnoff_data[index] = convert_all_quantities_to_string(item)
        obj_to_return = smirnoff_data
    elif isinstance(smirnoff_data, unit.Quantity):
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


try:
    import openmm
    from openff.units.openmm import from_openmm

    @object_to_quantity.register(openmm.unit.Quantity)
    def _(obj):
        return from_openmm(obj)

except ImportError:
    pass  # pragma: nocover


def serialize_numpy(np_array) -> Tuple[bytes, Tuple[int]]:
    """
    Serializes a numpy array into a big-endian bytestring and tuple representing its shape.

    Parameters
    ----------
    np_array : A numpy array
        Input numpy array

    Returns
    -------
    serialized : bytes
        A big-endian bytestring of the NumPy array.
    shape : tuple of ints
        The shape of the serialized array
    """
    import numpy as np

    bigendian_float = np.dtype(float).newbyteorder(">")
    bigendian_array = np_array.astype(bigendian_float)
    serialized = bigendian_array.tobytes()
    shape = np_array.shape
    return serialized, shape


def deserialize_numpy(serialized_np: Union[bytes, List], shape: Tuple[int]):
    """
    Deserializes a numpy array from a bytestring or list. The input, if a bytestring, is
    assumed to be in big-endian byte order.

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
        dt = np.dtype("float").newbyteorder(">")
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
            "Error: get_molecule_parameterIDs has been provided a list of oemols which contains some duplicates: "
            f"{duplicates}"
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
            for parameter_type in data[force_type].values():
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
            # the smirnoff dicts as lists of OrderedDicts or dicts
            new_parameter_list = list()
            for param in val:
                new_parameter_list.append(sort_smirnoff_dict(param))
            sorted_dict[key] = new_parameter_list
        else:
            # Handle metadata or the bottom of a recursive dict
            sorted_dict[key] = val
    return sorted_dict
