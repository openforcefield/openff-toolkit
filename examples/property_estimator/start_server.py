#!/usr/bin/env python

import shutil
from os import path

from openforcefield.properties.estimator import runner
from openforcefield.properties.estimator.backends.dask import DaskLocalClusterBackend
from openforcefield.properties.estimator.storage import LocalFileStorage
from openforcefield.utils import setup_timestamp_logging


def start_property_estimator_server():
    """An integrated test of the property estimator"""

    setup_timestamp_logging()

    working_directory = 'working-directory'

    # Remove any existing data.
    if path.isdir(working_directory):
        shutil.rmtree(working_directory)

    calculation_backend = DaskLocalClusterBackend(1, 1)
    storage_backend = LocalFileStorage(root_key='stored_data')

    property_server = runner.PropertyCalculationRunner(calculation_backend,
                                                       storage_backend,
                                                       working_directory=working_directory)

    property_server.run_until_killed()


if __name__ == "__main__":
    start_property_estimator_server()
