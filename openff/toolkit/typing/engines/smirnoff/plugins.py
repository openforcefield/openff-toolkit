"""This module defines functions for loading parameter handler and parser classes which
have been registered through the `entrypoint plugin system <https://packaging.python.
org/guides/creating-and-discovering-plugins/#using-package-metadata>`_.

.. warning ::

   This feature is experimental and may be removed / altered in future versions.

Currently only ``ParameterHandler`` classes can be registered via the plugin
system. This is possible by registering the handler class through an entry
point in your projects ``setup.py`` file::

    setup(
        ...
        entry_points={
            'openff.toolkit.plugins.handlers': ['CustomHandler = myapp:CustomHandler']
        },
        ...
    )

where in this example your package is named ``myapp`` and contains a class which
inherits from ``ParameterHandler`` named ``CustomHandler``.
"""

import logging
from importlib.metadata import entry_points

from openff.toolkit.typing.engines.smirnoff.parameters import ParameterHandler

logger = logging.getLogger(__name__)

SUPPORTED_PLUGIN_NAMES = ["handlers"]  # io_handlers could be supported in future.


def _load_handler_plugins(handler_name: str, expected_type: type) -> list[type]:
    """Loads parameter handler plugins of a specified type which have been registered
    with the ``entrypoint`` plugin system.

    Parameters
    ----------
    handler_name
        The name of the hander plugin. This can currently be any of the names
        listed in ``SUPPORTED_PLUGIN_NAMES``.
    expected_type
        The expected class type of the plugin. E.g. when loading parameter io
        handler plugins the expected class type is ``ParameterIOHandler``. Any
        classes not matching the expected type will be skipped.
    """

    if handler_name not in SUPPORTED_PLUGIN_NAMES:
        raise NotImplementedError(
            f"Plugins named {handler_name} not supported. The list of currently "
            f"supported plugin names is: {SUPPORTED_PLUGIN_NAMES}."
        )

    PluginType = type[ParameterHandler]

    discovered_plugins: list[PluginType] = list()

    for entry_point in entry_points().select(
        group=f"openff.toolkit.plugins.{handler_name}"
    ):
        try:
            discovered_plugins.append(entry_point.load())
        except (ImportError, AttributeError):
            logger.exception(f"Could not load the {entry_point} plugin")
            continue

    valid_plugins: list[PluginType] = list()

    for discovered_plugin in discovered_plugins:
        if not issubclass(discovered_plugin, expected_type):
            logger.info(
                f"The {discovered_plugin.__name__} object has been registered as a "
                f"{handler_name} plugin, but does not inherit from "
                f"{expected_type.__name__}. This plugin will be skipped."
            )
            continue

        valid_plugins.append(discovered_plugin)

    return valid_plugins


def load_handler_plugins() -> list[type[ParameterHandler]]:
    """Loads any ``ParameterHandler`` class plugins which have been registered through
    the ``entrypoint`` plugin system.

    Returns
    -------
        The registered ``ParameterHandler`` plugins.
    """
    return _load_handler_plugins("handlers", ParameterHandler)
