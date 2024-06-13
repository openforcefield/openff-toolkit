"""
Test classes and function in module openff.toolkit.typing.engines.smirnoff.plugins
"""

import pytest

from openff.toolkit import ForceField, Quantity
from openff.toolkit.typing.engines.smirnoff.plugins import (
    _load_handler_plugins,
    load_handler_plugins,
)


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

    assert len(registered_plugins) == 2

    registered_plugin_names = [plugin.__name__ for plugin in registered_plugins]

    assert "CustomHandler" in registered_plugin_names
    assert "BuckinghamHandler" in registered_plugin_names


def test_do_not_load_other_type():
    with pytest.raises(NotImplementedError, match="foobar"):
        _load_handler_plugins(handler_name="foobar", expected_type=type(None))


@pytest.mark.skip(reason="Broken around 8/29/23, unclear why")
def test_skip_wrong_subclass(caplog):
    import logging

    caplog.set_level(logging.INFO)
    load_handler_plugins()

    assert "does not inherit from ParameterHandler" in caplog.text, caplog.text


def test_buckingham_type():
    """Reproduce, in part, issue #1888."""
    from custom_plugins.handler_plugins import BuckinghamHandler

    parameter = BuckinghamHandler.BuckinghamType(
        smirks="[*:1]",
        a="2 kilojoule_per_mole",
        b="1/nanometer",
        c="-0.5 kilojoule_per_mole * nanometer**6",
    )

    for param in ["a", "b", "c"]:
        assert isinstance(getattr(parameter, param), Quantity)

    assert str(parameter.a.units) == "kilojoule_per_mole"
    assert str(parameter.b.units) == "1 / nanometer"
    assert str(parameter.c.units) == "kilojoule_per_mole * nanometer ** 6"
