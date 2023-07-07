"""A project (and object) for storing, manipulating, and converting molecular mechanics data."""
import importlib
from types import ModuleType
from typing import TYPE_CHECKING

from openff.interchange._version import get_versions

if TYPE_CHECKING:
    # Type checkers can't see lazy-imported objects
    from openff.interchange.components.interchange import Interchange

# Handle versioneer
versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions

__all__: list[str] = [
    "Interchange",
]

_objects: dict[str, str] = {
    "Interchange": "openff.interchange.components.interchange",
}


def __getattr__(name) -> ModuleType:
    """
    Lazily import objects from submodules.

    Taken from openff/toolkit/__init__.py
    """
    module = _objects.get(name)
    if module is not None:
        return importlib.import_module(module).__dict__[name]

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    """Add _objects to dir()."""
    keys = (*globals().keys(), *_objects.keys())
    return sorted(keys)
