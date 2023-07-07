"""
Collections for custom SMIRNOFF handlers.
"""

from openff.interchange.smirnoff._base import SMIRNOFFCollection

__all__ = [
    "SMIRNOFFCollection",
]


def load_smirnoff_plugins() -> list:
    """Load external potential handlers as plugins."""
    from importlib_metadata import entry_points

    plugin_handlers = list()

    for entry_point in entry_points().select(
        group="openff.interchange.plugins.collections",
    ):
        plugin_handlers.append(entry_point.load())

    return plugin_handlers
