"""
This is a configuration file for PyTest.

This file is designed to be used with nbval: https://nbval.readthedocs.io/

The nbval-sanitize.cfg file is also for use with nbval.
"""


def pytest_collectstart(collector):
    """
    Ignore certain outputs when validating documentation notebooks

    See https://nbval.readthedocs.io/en/latest/#Skipping-certain-output-types
    """
    if collector.fspath and collector.fspath.ext == ".ipynb":
        collector.skip_compare += (
            # SVGs are often produced by Molecule.visualize()
            # The precise SVG output varies with environment and RDkit version,
            # and failures of visualize() produce exceptions or text output, so
            # we don't lose much by ignoring the actual content.
            "image/svg+xml",
        )
