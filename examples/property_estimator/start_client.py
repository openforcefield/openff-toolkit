#!/usr/bin/env python

import logging
import shutil
from os import path

from openforcefield.properties.datasets import ThermoMLDataSet
from openforcefield.properties.estimator import client
from openforcefield.typing.engines import smirnoff
from openforcefield.utils import get_data_filename, setup_timestamp_logging


def submit_calculation_to_server():
    """Submit calculations to a running server instance"""

    setup_timestamp_logging()

    data_set = ThermoMLDataSet.from_file_list(get_data_filename('properties/single_density.xml'))
    force_field = smirnoff.ForceField(get_data_filename('forcefield/smirnoff99Frosst.offxml'))

    property_estimator = client.PropertyEstimator()
    ticket_ids = property_estimator.submit_computations(data_set, force_field)

    logging.info('Ticket info: {}'.format(ticket_ids))
    result = property_estimator.wait_for_result(ticket_ids)

    logging.info('The server has returned a response: {}'.format(result))


if __name__ == "__main__":
    submit_calculation_to_server()
