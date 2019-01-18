#!/usr/bin/env python

import shutil
from os import path

from openforcefield.properties.estimator import runner
from openforcefield.properties.estimator.backends.dask import DaskLocalClusterBackend
from openforcefield.utils import setup_timestamp_logging


def start_property_estimator_server():
    """An integrated test of the property estimator"""

    setup_timestamp_logging()

    # Remove any existing data.
    if path.isdir('property-data'):
        shutil.rmtree('property-data')

    backend = DaskLocalClusterBackend(1, 1)

    property_server = runner.PropertyCalculationRunner(backend)
    property_server.run_until_killed()


if __name__ == "__main__":
    start_property_estimator_server()
