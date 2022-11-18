"""
Test classes and function in module openff.toolkit.typing.engines.smirnoff.plugins
"""
import pytest

from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.typing.engines.smirnoff.plugins import load_handler_plugins


def test_force_field_custom_handler():
    """Tests a force field can make use of a custom parameter handler registered
    through the entrypoint plugin system.
    """

    # Construct a simple FF which only uses the custom handler.
    force_field_contents = "\n".join(
        [
            "<?xml version='1.0' encoding='ASCII'?>",
            "<SMIRNOFF version='0.3' aromaticity_model='OEAroModel_MDL'>",
            "  <CustomHandler version='0.3'></CustomHandler>",
            "</SMIRNOFF>",
        ]
    )

    # An exception should be raised when plugins aren't allowed.
    with pytest.raises(KeyError) as error_info:
        ForceField(force_field_contents)

    assert (
        "Cannot find a registered parameter handler class for tag 'CustomHandler'"
        in error_info.value.args[0]
    )

    # Otherwise the FF should be created as expected.
    force_field = ForceField(force_field_contents, load_plugins=True)

    parameter_handler = force_field.get_parameter_handler("CustomHandler")
    assert parameter_handler is not None
    assert parameter_handler.__class__.__name__ == "CustomHandler"

    assert parameter_handler.__class__ in force_field._plugin_parameter_handler_classes


def test_load_handler_plugins():
    """Tests that parameter handlers can be registered as plugins."""

    registered_plugins = load_handler_plugins()

    assert len(registered_plugins) == 1
    assert registered_plugins[0].__name__ == "CustomHandler"
