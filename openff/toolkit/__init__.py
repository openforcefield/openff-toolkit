"""
openff-toolkit
A modern, extensible library for molecular mechanics force field science from the Open Force Field Consortium.
"""

from ._version import get_versions  # type: ignore

__version__ = get_versions()["version"]
del get_versions
