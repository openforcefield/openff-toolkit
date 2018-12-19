import logging
import sys
import shutil

from os import path

from openforcefield.properties.datasets import ThermoMLDataSet

from openforcefield.typing.engines import smirnoff
from openforcefield.utils import get_data_filename

from openforcefield.properties.estimator import client


def run_property_estimator():
    """An integrated test of the property estimator"""

    # Remove any existing data.
    if path.isdir('property-data'):
        shutil.rmtree('property-data')

    # Set up time-based logging to help debug threading issues.
    formatter = logging.Formatter(fmt='%(asctime)s.%(msecs)03d %(levelname)-8s %(message)s',
                                  datefmt='%H:%M:%S')

    screen_handler = logging.StreamHandler(stream=sys.stdout)
    screen_handler.setFormatter(formatter)

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.addHandler(screen_handler)

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

    property_estimator = client.PropertyEstimator()

    results = property_estimator.compute_properties(data_set.properties, force_field, 1)
    client.PropertyEstimator.produce_calculation_report(data_set, results)


if __name__ == "__main__":
    run_property_estimator()
