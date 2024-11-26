"""
openff-toolkit
A modern, extensible library for molecular mechanics force field science from the Open Force Field Consortium.
"""

import importlib
from importlib.metadata import version
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    # These types are imported lazily at runtime, but we need to tell type
    # checkers what they are
    from openff.units import Quantity, unit

    from openff.toolkit.topology import Molecule, Topology
    from openff.toolkit.typing.engines.smirnoff import (
        ForceField,
        get_available_force_fields,
    )
    from openff.toolkit.utils.toolkits import (
        GLOBAL_TOOLKIT_REGISTRY,
        AmberToolsToolkitWrapper,
        BuiltInToolkitWrapper,
        OpenEyeToolkitWrapper,
        RDKitToolkitWrapper,
        ToolkitRegistry,
    )

__version__ = version("openff.toolkit")

__all__ = [
    "__version__",
    "Molecule",
    "Topology",
    "ForceField",
    "get_available_force_fields",
    "GLOBAL_TOOLKIT_REGISTRY",
    "AmberToolsToolkitWrapper",
    "BuiltInToolkitWrapper",
    "OpenEyeToolkitWrapper",
    "RDKitToolkitWrapper",
    "ToolkitRegistry",
    "Quantity",
    "unit",
]

# Dictionary of objects to lazily import; maps the object's name to its module path
_lazy_imports_obj = {
    "ForceField": "openff.toolkit.typing.engines.smirnoff",
    "get_available_force_fields": "openff.toolkit.typing.engines.smirnoff",
    "Molecule": "openff.toolkit.topology",
    "Topology": "openff.toolkit.topology",
    "GLOBAL_TOOLKIT_REGISTRY": "openff.toolkit.utils.toolkits",
    "AmberToolsToolkitWrapper": "openff.toolkit.utils.toolkits",
    "BuiltInToolkitWrapper": "openff.toolkit.utils.toolkits",
    "OpenEyeToolkitWrapper": "openff.toolkit.utils.toolkits",
    "RDKitToolkitWrapper": "openff.toolkit.utils.toolkits",
    "ToolkitRegistry": "openff.toolkit.utils.toolkits",
    "Quantity": "openff.units",
    "unit": "openff.units",
    # Remember to add new lazy imports to __all__ and the if TYPE_CHECKING imports
}

# Dictionary of modules to lazily import; maps the modules's name to its path
_lazy_imports_mod = {
    "topology": "openff.toolkit.topology",
    "typing": "openff.toolkit.typing",
    "utils": "openff.toolkit.utils",
    "unit": "openff.units.units",
}


def __getattr__(name):
    """Lazily import objects from _lazy_imports_obj or _lazy_imports_mod

    Note that this method is only called by Python if the name cannot be found
    in the current module."""
    obj_mod = _lazy_imports_obj.get(name)
    if obj_mod is not None:
        mod = importlib.import_module(obj_mod)
        try:
            return mod.__dict__[name]
        except KeyError:  # account for lazy loaders
            return getattr(mod, name)

    lazy_mod = _lazy_imports_mod.get(name)
    if lazy_mod is not None:
        return importlib.import_module(lazy_mod)

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    """Add _lazy_imports_obj and _lazy_imports_mod to dir(<module>)"""
    keys = (*globals().keys(), *_lazy_imports_obj.keys(), *_lazy_imports_mod.keys())
    return sorted(keys)
