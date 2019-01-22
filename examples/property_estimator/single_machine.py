#!/usr/bin/env python

import logging
import shutil
from os import path

from openforcefield.properties.datasets import ThermoMLDataSet
from openforcefield.properties.estimator import client
from openforcefield.properties.estimator import runner
from openforcefield.properties.estimator.backends.dask import DaskLocalClusterBackend
from openforcefield.properties.estimator.storage import LocalFileStorage
from openforcefield.typing.engines import smirnoff
from openforcefield.utils import get_data_filename, setup_timestamp_logging


def run_property_estimator():
    """An integrated test of the property estimator"""

    setup_timestamp_logging()

    working_directory = 'working-directory'

    # Remove any existing data.
    if path.isdir(working_directory):
        shutil.rmtree(working_directory)

    data_set = ThermoMLDataSet.from_file_list(get_data_filename('properties/single_density.xml'))
    # data_set = ThermoMLDataSet.from_file_list(get_data_filename('properties/single_dielectric.xml'))

    # data_set = ThermoMLDataSet.from_file_list(get_data_filename('properties/single_density.xml'),
    #                                           get_data_filename('properties/single_dielectric.xml'))

    # data_set = ThermoMLDataSet.from_file_list(get_data_filename('properties/density_dielectric.xml'))
    # data_set = ThermoMLDataSet.from_file_list(get_data_filename('properties/two_species.xml'))
    # data_set = ThermoMLDataSet.from_file_list(get_data_filename('properties/binary.xml'))
    # data_set = ThermoMLDataSet.from_file_list(get_data_filename('properties/fake_data.xml'))
    # data_set = ThermoMLDataSet.from_file_list(get_data_filename('properties/j.jct.2007.09.004.xml'))

    force_field = smirnoff.ForceField(get_data_filename('forcefield/smirnoff99Frosst.offxml'))

    calculation_backend = DaskLocalClusterBackend(1, 1)
    storage_backend = LocalFileStorage(root_key='stored_data')

    property_server = runner.PropertyCalculationRunner(calculation_backend,
                                                       storage_backend,
                                                       working_directory=working_directory)

    property_estimator = client.PropertyEstimator()
    property_estimator.submit_computations(data_set, force_field)

    property_server.run_until_complete()

    for property_id in property_server.finished_calculations:

        logging.info('Calculation {}:'.format(property_id))
        logging.info('{}'.format(property_server.finished_calculations[property_id].json()))


if __name__ == "__main__":
    run_property_estimator()
