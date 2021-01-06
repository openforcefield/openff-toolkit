import warnings

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions


warnings.warn(
    "Importing this package as `import openforcefield.XXX` and "
    "`from openforcefield import XXX` was marked for deprecation in version `0.8.3`. From version "
    "`0.9.0` onwards this package will need to be imported as "
    "`import openff.toolkit.XXX` and `from openff.toolkit import XXX`. See the `0.8.3` "
    "release notes for more information.",
    DeprecationWarning,
    stacklevel=2,
)
