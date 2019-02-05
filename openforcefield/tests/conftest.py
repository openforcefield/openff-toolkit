#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Configuration file for pytest.

This adds the following command line options.
- runslow: Run tests marked as slow (default is False).

"""


#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import pytest


#=============================================================================================
# UTILITY FUNCTIONS
#=============================================================================================

def untar_full_alkethoh_set():
    """When running slow tests, we unpack the full AlkEthOH test in advance to speed things up."""
    import os
    import tarfile
    from openforcefield.utils import get_data_filename
    tarfile_path = os.path.join(get_data_filename('molecules'), 'AlkEthOH_tripos.tar.gz')
    with tarfile.open(tarfile_path, 'r:gz') as tar:
        tar.extractall(path=os.path.dirname(tarfile_path))


#=============================================================================================
# CONFIGURATION
#=============================================================================================

def pytest_addoption(parser):
    """Add the pytest command line option --runslow.

    If --runslow is not given, tests marked with pytest.mark.slow are
    skipped.
    """
    parser.addoption(
        "--runslow", action="store_true", default=False, help="run slow tests"
    )


def pytest_collection_modifyitems(config, items):
    if not config.getoption("runslow"):
        # Mark for skipping all items marked as slow.
        skip_slow = pytest.mark.skip(reason="need to set --runslow to run this test")
        for item in items:
            if "slow" in item.keywords:
                item.add_marker(skip_slow)
    else:
        # If --runslow is given, we don't have to mark items for skipping,
        # but we need to extract the whole AlkEthOH set (see
        # test_forcefield::test_alkethoh_parameters_assignment).
        untar_full_alkethoh_set()
