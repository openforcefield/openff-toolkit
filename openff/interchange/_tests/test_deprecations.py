import pytest

from openff.interchange import Interchange
from openff.interchange._tests import _BaseTest
from openff.interchange.warnings import InterchangeDeprecationWarning


class TestDeprecation(_BaseTest):
    def test_potential_handler_deprecation(self):
        from openff.interchange.components.potentials import Collection

        with pytest.warns(
            InterchangeDeprecationWarning,
            match="`PotentialHandler` has been renamed to `Collection`.",
        ):
            from openff.interchange.components.potentials import PotentialHandler

        assert PotentialHandler is Collection

    @pytest.fixture()
    def prepared_system(self, sage, water):
        return Interchange.from_smirnoff(sage, [water])

    def test_slot_map_deprecation(self, prepared_system):
        with pytest.warns(
            InterchangeDeprecationWarning,
            match="The `slot_map` attribute is deprecated. Use `key_map` instead.",
        ):
            prepared_system["vdW"].slot_map

    def test_handlers_deprecation(self, prepared_system):
        with pytest.warns(
            InterchangeDeprecationWarning,
            match="The `handlers` attribute is deprecated. Use `collections` instead.",
        ):
            prepared_system.handlers
