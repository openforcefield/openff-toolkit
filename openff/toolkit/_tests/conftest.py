"""
Configuration file for pytest.

This adds the following command line options.
- runslow: Run tests marked as slow (default is False).

"""

import logging

import pytest

logger = logging.getLogger(__name__)


def pytest_configure(config):
    """
    Initialization hook to register custom markers without a pytest.ini
    More info: https://docs.pytest.org/en/latest/reference.html#initialization-hooks
    """
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with `-m 'not slow'`)"
    )


def untar_full_alkethoh_and_freesolv_set():
    """When running slow tests, we unpack the full AlkEthOH and FreeSolv test
    sets in advance to speed things up.

    See
        test_forcefield.py::test_alkethoh_parameters_assignment
        test_forcefield.py::test_freesolv_parameters_assignment
    """
    import os
    import tarfile

    from openff.toolkit.utils import get_data_file_path

    molecule_dir_path = get_data_file_path("molecules")
    for tarfile_name in ["AlkEthOH_tripos.tar.gz", "FreeSolv.tar.gz"]:
        tarfile_path = os.path.join(molecule_dir_path, tarfile_name)
        with tarfile.open(tarfile_path, "r:gz") as tar:
            tar.extractall(path=molecule_dir_path)


def pytest_addoption(parser):
    """Add the pytest command line option --runslow and --failwip.

    If --runslow is not given, tests marked with pytest.mark.slow are
    skipped.

    If --failwip is not given, tests marked with pytest.mark.wip are
    xfailed.

    Parameters
    ----------
    parser : argparsing.parser
        The parser used by pytest to process arguments

    Returns
    -------
    None

    """

    # Loaded pytest plugins define their own arguments, and in certain cases
    # our options use the same name. Although this can define two different
    # behaviors for the same argument, the two options below are fairly
    # unambiguous in their interpretation, so we skip any errors and allow
    # the testing to proceed

    try:
        parser.addoption(
            "--runslow", action="store_true", default=False, help="run slow tests"
        )
    except ValueError:
        logger.warning(
            "Option --runslow already added elsewhere (from a plugin possibly?). Skipping..."
        )

    try:
        parser.addoption(
            "--failwip",
            action="store_true",
            default=False,
            help="fail work in progress tests",
        )
    except ValueError:
        logger.warning(
            "Option --failwip already added elsewhere (from a plugin possibly?). Skipping..."
        )


def pytest_collection_modifyitems(config, items):
    if config.getoption("runslow"):
        # If --runslow is given, we don't have to mark items for skipping,
        # but we need to extract the whole AlkEthOH and FreeSolv sets (see
        # test_forcefield::test_alkethoh/freesolv_parameters_assignment).
        untar_full_alkethoh_and_freesolv_set()
    else:
        # Mark for skipping all items marked as slow.
        skip_slow = pytest.mark.skip(
            reason="specify --runslow pytest option to run this test."
        )
        for item in items:
            if "slow" in item.keywords:
                item.add_marker(skip_slow)

    # Mark work-in-progress tests for xfail.
    if not config.getoption("failwip"):
        xfail_wip_reason = (
            "This is a work in progress test. Specify "
            "--failwip pytest option to make this test fail."
        )
        for item in items:
            if "wip" in item.keywords:
                # Augment original reason.
                reason = xfail_wip_reason + item.get_closest_marker("wip").kwargs.get(
                    "reason", ""
                )
                item.add_marker(pytest.mark.xfail(reason=reason))
