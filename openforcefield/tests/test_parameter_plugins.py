"""
Test classes and function in module openforcefield.typing.engines.smirnoff.plugins
"""
import pkg_resources
import pytest

from openforcefield.typing.engines.smirnoff import ForceField
from openforcefield.typing.engines.smirnoff.plugins import load_handler_plugins


@pytest.yield_fixture()
def mock_entry_point_plugins():
    """Registers a fake parameter handler and io handler with the
    entry point plugin system.

    Notes
    -----
    This function is based on `this stack overflow answer
    <https://stackoverflow.com/a/48666503/11808960>`_
    """

    previous_entries = pkg_resources.working_set.entries
    previous_entry_keys = pkg_resources.working_set.entry_keys
    previous_by_key = pkg_resources.working_set.by_key

    # Create a fake distribution to insert into the global working_set
    distribution = pkg_resources.Distribution(__file__)

    # Create the fake entry point definitions. These include a parameter handler
    # which is supported, and an io parameter handler which should be skipped.
    handler_entry_point = pkg_resources.EntryPoint.parse(
        "CustomHandler = openforcefield.tests.plugins:CustomHandler", dist=distribution
    )
    io_handler_entry_point = pkg_resources.EntryPoint.parse(
        "CustomIOHandler = openforcefield.tests.plugins:CustomIOHandler",
        dist=distribution,
    )

    # Add the mapping to the fake EntryPoint
    distribution._ep_map = {
        "openff.toolkit.plugins.handlers": {
            "CustomHandler": handler_entry_point,
            "CustomIOHandler": io_handler_entry_point,
        }
    }

    # Add the fake distribution to the global working_set
    pkg_resources.working_set.add(distribution, "CustomHandler")
    pkg_resources.working_set.add(distribution, "CustomIOHandler")

    yield

    pkg_resources.working_set.entries = previous_entries
    pkg_resources.working_set.entry_keys = previous_entry_keys
    pkg_resources.working_set.by_key = previous_by_key


def test_force_field_custom_handler(mock_entry_point_plugins):  # doctest: +SKIP
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


def test_load_handler_plugins(mock_entry_point_plugins):
    """Tests that parameter handlers can be registered as plugins."""

    registered_plugins = load_handler_plugins()

    assert len(registered_plugins) == 1
    assert registered_plugins[0].__name__ == "CustomHandler"
