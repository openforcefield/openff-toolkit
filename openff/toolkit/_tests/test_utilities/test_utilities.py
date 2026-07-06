import os

import pytest

from openff.toolkit._utilities.exceptions import MissingOptionalDependencyError
from openff.toolkit._utilities.testing import skip_if_missing
from openff.toolkit._utilities.utilities import (
    get_data_dir_path,
    get_data_file_path,
    has_executable,
    has_package,
    requires_oe_module,
    requires_package,
    temporary_cd,
)


def compare_paths(path_1: str, path_2: str) -> bool:
    """Checks whether two paths are the same.

    Parameters
    ----------
    path_1
        The first path.
    path_2
        The second path.

    Returns
    -------
    True if the paths are equivalent.
    """
    return os.path.normpath(path_1) == os.path.normpath(path_2)


def test_get_data_dir_path():
    """Tests that the `get_data_dir_path` can correctly find
    data directories.
    """

    # Test a path which should exist.
    data_file_path = get_data_dir_path("data/more", package_name="openff.toolkit._utilities")
    assert os.path.isdir(data_file_path)

    # Ensure a double-checking through data/ takes place
    data_file_path = get_data_dir_path("more", package_name="openff.toolkit._utilities")
    assert os.path.isdir(data_file_path)

    # Test a path which should not exist.
    with pytest.raises(NotADirectoryError):
        get_data_dir_path("missing", package_name="openff.toolkit._utilities")

    # Test that a file directory raises NotaDirectoryError
    with pytest.raises(NotADirectoryError):
        get_data_dir_path("data/data.dat", package_name="openff.toolkit._utilities")


def test_get_data_file_path():
    """Tests that the `get_data_file_path` can correctly find
    data files.
    """

    # Test a path which should exist.
    data_file_path = get_data_file_path("data/data.dat", package_name="openff.toolkit._utilities")
    assert os.path.isfile(data_file_path)

    # Ensure a double-checking through data/ takes place
    data_file_path = get_data_file_path("data.dat", package_name="openff.toolkit._utilities")
    assert os.path.isfile(data_file_path)

    # Test a path which should not exist.
    with pytest.raises(FileNotFoundError):
        get_data_file_path("missing.file", package_name="openff.toolkit._utilities")

    # Test that a directory raises FileNotFoundError
    with pytest.raises(FileNotFoundError):
        get_data_file_path("data/", package_name="openff.toolkit._utilities")


def test_temporary_cd():
    """Tests that temporary cd works as expected"""

    original_directory = os.getcwd()

    # Move to the parent directory
    with temporary_cd(os.pardir):
        current_directory = os.getcwd()
        expected_directory = os.path.abspath(os.path.join(original_directory, os.pardir))

        assert compare_paths(current_directory, expected_directory)

    assert compare_paths(os.getcwd(), original_directory)

    # Move to a temporary directory
    with temporary_cd():
        assert not compare_paths(os.getcwd(), original_directory)

    assert compare_paths(os.getcwd(), original_directory)

    # Move to the same directory
    with temporary_cd(""):
        assert compare_paths(os.getcwd(), original_directory)

    assert compare_paths(os.getcwd(), original_directory)

    with temporary_cd(os.curdir):
        assert compare_paths(os.getcwd(), original_directory)

    assert compare_paths(os.getcwd(), original_directory)


def test_has_package():
    assert has_package("os")
    assert has_package("pytest")
    assert not has_package("nummmmmmpy")


def test_has_executable():
    assert has_executable("pwd")
    assert has_executable("pytest")

    assert has_executable("/usr/bin/whoami")

    assert not has_package("pyyyyython")


def test_requires_package():
    """Tests that the ``requires_package`` utility behaves as expected."""

    def dummy_function():
        pass

    # sys should always be found so this should not raise an exception.
    requires_package("sys")(dummy_function)()

    with pytest.raises(MissingOptionalDependencyError) as error_info:
        requires_package("fake-lib")(dummy_function)()

    assert error_info.value.library_name == "fake-lib"


@skip_if_missing("openeye.oechem")
@pytest.mark.skipif("OE_LICENSE" not in os.environ, reason="Requires an OpenEye license is NOT set up")
def test_requires_oe_module():
    """Tests that the ``requires_package`` utility behaves as expected when an OpenEye license is set up."""

    def dummy_function():
        pass

    requires_oe_module("oechem")(dummy_function)()


@pytest.mark.parametrize("oe_module", ["oechem", "oeiupac", "oeomega", "oequacpac", "oedepict"])
@skip_if_missing("openeye")
@pytest.mark.skipif("OE_LICENSE" in os.environ, reason="Requires an OpenEye license is NOT set up")
def test_requires_oe_module_installed_missing_license(oe_module):
    """Tests that the ``requires_package`` utility behaves as expected while OpenEye toolkits are
    installed but no OpenEye license is set up."""

    def dummy_function():
        pass

    with pytest.raises(MissingOptionalDependencyError) as error_info:
        requires_oe_module(oe_module)(dummy_function)()

    assert oe_module in str(error_info.value)
    assert "conda-forge" not in str(error_info.value)


@pytest.mark.skipif(has_package("openeye.oechem"), reason="Requires OpenEye toolkits are NOT installed")
def test_requires_oe_module_not_installed():
    """Tests that the ``requires_package`` utility behaves as expected while OpenEye toolkits are
    installed."""

    def dummy_function():
        pass

    with pytest.raises(MissingOptionalDependencyError) as error_info:
        requires_oe_module("oechem")(dummy_function)()

    assert "oechem" in str(error_info.value)
    assert "conda-forge" not in str(error_info.value)
