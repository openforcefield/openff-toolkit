import logging
import sys
import shutil

from os import path

import simtk.unit as unit

from openforcefield import substances

from openforcefield.datasets import ThermoMLDataSet

from openforcefield.typing.engines import smirnoff
from openforcefield.utils import get_data_filename

from openforcefield.propertycalculator import protocols
from openforcefield.propertycalculator import runner


def run_property_calculator():

    if path.isdir('property-data'):
        shutil.rmtree('property-data')

    formatter = logging.Formatter(fmt='%(asctime)s.%(msecs)03d %(levelname)-8s %(message)s',
                                  datefmt='%H:%M:%S')

    screen_handler = logging.StreamHandler(stream=sys.stdout)
    screen_handler.setFormatter(formatter)

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.addHandler(screen_handler)

    # data_set = ThermoMLDataSet.from_file_list('../data/properties/fake_data.xml')
    data_set = ThermoMLDataSet.from_file_list('../data/properties/j.jct.2007.09.004.xml')
    force_field = smirnoff.ForceField(get_data_filename('forcefield/smirnoff99Frosst.offxml'))

    property_calculator = runner.PropertyCalculationRunner(2)
    results = property_calculator.run(data_set.measured_properties, force_field)


# test_mixture_protocol()