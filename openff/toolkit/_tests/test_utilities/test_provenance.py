import importlib

import pytest

from openff.toolkit._utilities.provenance import (
    _get_conda_list_package_versions,
    get_ambertools_version,
)
from openff.toolkit._utilities.warnings import CondaExecutableNotFoundWarning


def mock_has_executable(executable_name: str) -> bool:
    return False


@pytest.mark.leaky
def test_conda_unavailable_returns_empty_dict(monkeypatch):
    """
    Test that if `conda|etc. list` fails, with or without AmberTools (properly) installed,
    get_ambertools_version returns `None`.

    This test leaks! Do not run it alongside other tests! Select either it "pytest -m leaky"
    or avoid it "pytest -m 'not leaky'"
    """
    with (
        monkeypatch.context() as m,
        pytest.warns(CondaExecutableNotFoundWarning, match="No .* executable found"),
    ):
        m.setenv("PIXI_IN_SHELL", "0")
        m.setenv("CONDA_SHLVL", "0")

        assert _get_conda_list_package_versions() == dict()


def test_conda_available_get_ambertools_version_found():
    # Skip if ambertools is not installed - sander would a better import test
    # for AmberTools (since parmed could come from the standalone
    # parmed package), however it's currently broken on ARM macs.
    pytest.importorskip("parmed")

    assert get_ambertools_version() is not None, (
        f"note that len of package versions is {len(_get_conda_list_package_versions())}"
    )


def test_conda_available_get_ambertools_version_not_found():
    try:
        importlib.import_module("parmed")
    except ImportError:
        assert len(_get_conda_list_package_versions()) > 0
        assert get_ambertools_version() is None
        return

    pytest.skip("only run when ambertools is not installed.")
