"""
openff-toolkit
A modern, extensible library for molecular mechanics force field science from the Open Force Field Consortium.
"""

from ._version import get_versions  # type: ignore

__version__ = get_versions()["version"]

import openff.toolkit.topology as topology
import openff.toolkit.typing as typing
import openff.toolkit.utils as utils
from openff.toolkit.topology import Molecule, Topology
from openff.toolkit.typing.engines.smirnoff.forcefield import (
    ForceField,
    get_available_force_fields,
)

__all__ = [
    "__version__",
    "topology",
    "typing",
    "utils",
    "Molecule",
    "Topology",
    "ForceField",
    "get_available_force_fields",
]
