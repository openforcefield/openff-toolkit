import logging
import sys
import shutil
import pickle

from os import path

from openforcefield.properties import Density, DielectricConstant

from openforcefield.properties.datasets import ThermoMLDataSet
from openforcefield.properties.estimator.backends.dask import DaskLocalClusterBackend
from openforcefield.properties.estimator.layers import SimulationLayer

from openforcefield.typing.engines import smirnoff
from openforcefield.utils import get_data_filename

from openforcefield.properties.estimator import client, CalculationSchema
from openforcefield.properties.estimator import runner


def test_calculation_schema():
    """Tests serialisation and deserialization of a calculation schema."""
    density_schema = Density.get_default_calculation_schema()
    density_schema.validate_interfaces()

    density_json = density_schema.json()
    print(density_json)
    
    dielectric_schema = DielectricConstant.get_default_calculation_schema()
    dielectric_schema.validate_interfaces()

    dielectric_json = dielectric_schema.json()
    print(dielectric_json)

    density_schema_from_json = CalculationSchema.parse_raw(density_json)
    print(density_schema_from_json)

    dielectric_schema_from_json = CalculationSchema.parse_raw(dielectric_json)
    print(dielectric_schema_from_json)


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

    backend = DaskLocalClusterBackend(1, 1)
    property_server = runner.PropertyCalculationRunner(backend)

    property_estimator = client.PropertyEstimator()
    ticket_ids = property_estimator.submit_computations(data_set, force_field)

    property_server.run_until_complete()

    # results = property_estimator.compute_properties(data_set, force_field)
    # client.PropertyEstimator.produce_calculation_report(data_set, results)


def run_simulation_backend():

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

    dummy_pickle = b''
    data_model = pickle.loads(dummy_pickle)

    backend = DaskLocalClusterBackend(1, 1)
    backend.start()

    def dummy_callback(*args, **kwargs):
        print(args, kwargs)
        pass

    simulation_layer = SimulationLayer()
    simulation_layer.schedule_calculation(backend, data_model, {}, dummy_callback, True)

    backend._client.get_worker_logs()


if __name__ == "__main__":
    # test_calculation_schema()
    # run_simulation_backend()
    run_property_estimator()
