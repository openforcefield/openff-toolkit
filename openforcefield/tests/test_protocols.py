import logging
import sys
import shutil

from os import path

from openforcefield.properties import Density

from openforcefield.properties.datasets import ThermoMLDataSet

from openforcefield.typing.engines import smirnoff
from openforcefield.utils import get_data_filename

from openforcefield.properties.estimator.components import protocols, groups

from openforcefield.properties.estimator.components.protocols import available_protocols
from openforcefield.properties.estimator.components.groups import available_groups

from openforcefield.properties.estimator.components.protocols import ProtocolInputReference
from openforcefield.properties.estimator import client, CalculationSchema


def test_calculation_schema():
    """Tests serialisation and deserialization of a calculation schema."""
    schema = CalculationSchema()

    schema.property_type = str(type(Density))
    schema.id = 'DensitySchema'

    build_coordinates = protocols.BuildCoordinatesPackmol()
    build_coordinates.id = 'build_coordinates'

    build_coordinates.input_references = [
        # Globals
        ProtocolInputReference(input_property_name='substance',
                               output_protocol_id='global',
                               output_property_name='substance')
    ]

    schema.protocols[build_coordinates.id] = build_coordinates.schema

    assign_topology = protocols.BuildSmirnoffTopology()
    assign_topology.id = 'build_topology'

    assign_topology.input_references = [
        # Globals
        ProtocolInputReference(input_property_name='force_field',
                               output_protocol_id='global',
                               output_property_name='force_field'),
        # Locals
        ProtocolInputReference(input_property_name='topology',
                               output_protocol_id=build_coordinates.id,
                               output_property_name='topology'),

        ProtocolInputReference(input_property_name='molecules',
                               output_protocol_id=build_coordinates.id,
                               output_property_name='molecules')
    ]

    schema.protocols[assign_topology.id] = assign_topology.schema

    converge_uncertainty = groups.ConditionalGroup({
        build_coordinates.id: build_coordinates,
        assign_topology.id: assign_topology
    })

    converge_uncertainty.id = 'converge_uncertainty'

    schema.groups[converge_uncertainty.id] = converge_uncertainty.schema

    json_object = schema.json()
    print(json_object)

    schema_from_json = CalculationSchema.parse_raw(json_object)

    test_object = available_protocols[schema.protocols[build_coordinates.id].type]()
    test_object.schema = schema_from_json.protocols[build_coordinates.id]

    print(schema)


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
