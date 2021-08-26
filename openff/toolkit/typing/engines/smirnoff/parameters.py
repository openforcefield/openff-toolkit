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
from openff.toolkit.typing.engines.smirnoff.parameter_handlers import (
    AngleHandler,
    BondHandler,
    ChargeIncrementModelHandler,
    ConstraintHandler,
    ElectrostaticsHandler,
    GBSAHandler,
    ImproperTorsionHandler,
    LibraryChargeHandler,
    NonbondedMethod,
    ParameterHandler,
    ProperTorsionHandler,
    ToolkitAM1BCCHandler,
    VirtualSiteHandler,
    _allow_only,
    _linear_inter_or_extrapolate,
    _NonbondedHandler,
    vdWHandler,
)
from openff.toolkit.typing.engines.smirnoff.parameter_types import (
    IndexedMappedParameterAttribute,
    IndexedParameterAttribute,
    MappedParameterAttribute,
    ParameterAttribute,
    ParameterList,
    ParameterType,
    _ParameterAttributeHandler,
)
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
