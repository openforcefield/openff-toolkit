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
    # If --runslow is given, we don't have to mark items for skipping.
    if config.getoption("--runslow"):
        return
    # Mark for skipping all items marked as slow.
    skip_slow = pytest.mark.skip(reason="need to set --runslow to run this test")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)
