def pytest_collectstart(collector):
    """
    Ignore the contents of stderr and svg image outputs in documentation notebooks

    See https://nbval.readthedocs.io/en/latest/#Skipping-certain-output-types
    """
    if collector.fspath and collector.fspath.ext == ".ipynb":
        collector.skip_compare += (
            "image/svg+xml",
            "stderr",
        )
