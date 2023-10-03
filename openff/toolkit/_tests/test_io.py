"""
Test classes and function in module openff.toolkit.typing.engines.smirnoff.io.

"""

import pytest

from openff.toolkit.typing.engines.smirnoff.io import XMLParameterIOHandler


class TestXMLParameterIOHandler:
    def test_from_string(self):
        pass

    def test_raise_file_not_found(self):
        """Raise FileNotFoundError when the file doesn't exist."""
        io_handler = XMLParameterIOHandler()
        with pytest.raises(FileNotFoundError):
            io_handler.parse_file("nonexisting_file.xml")
