import uuid

from simtk import unit

from openforcefield.properties import Density, DielectricConstant, PropertyPhase
from openforcefield.properties.estimator import CalculationSchema
from openforcefield.properties.estimator.client import PropertyEstimatorOptions
from openforcefield.properties.estimator.layers.simulation import DirectCalculation, DirectCalculationGraph
from openforcefield.properties.statistics import Statistics
from openforcefield.properties.substances import Mixture
from openforcefield.properties.thermodynamics import ThermodynamicState
from openforcefield.utils import get_data_filename, graph


def test_statistics_object():

    statistics_object = Statistics.from_openmm_csv(get_data_filename('properties/stats_openmm.csv'), 1*unit.atmosphere)
    print(statistics_object)

    statistics_object.save_as_pandas_csv('stats_pandas.csv')

    statistics_object = Statistics.from_pandas_csv('stats_pandas.csv')


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

    density_recreated_json = density_schema_from_json.json()
    assert density_json == density_recreated_json

    dielectric_schema_from_json = CalculationSchema.parse_raw(dielectric_json)
    print(dielectric_schema_from_json)

    dielectric_recreated_json = dielectric_schema_from_json.json()
    assert dielectric_json == dielectric_recreated_json


def test_density_merging():

    substance = Mixture()
    substance.add_component('C', 1.0)

    dummy_property = Density(thermodynamic_state=ThermodynamicState(temperature=298*unit.kelvin,
                                                                    pressure=1*unit.atmosphere),
                             phase=PropertyPhase.Liquid,
                             substance=substance,
                             value=10*unit.gram/unit.mole,
                             uncertainty=1*unit.gram/unit.mole,
                             id=str(uuid.uuid4()))

    calculation_schema = dummy_property.get_default_calculation_schema()

    density_calculation_a = DirectCalculation(dummy_property,
                                              get_data_filename('forcefield/smirnoff99Frosst.offxml'),
                                              calculation_schema,
                                              PropertyEstimatorOptions())

    density_calculation_b = DirectCalculation(dummy_property,
                                              get_data_filename('forcefield/smirnoff99Frosst.offxml'),
                                              calculation_schema,
                                              PropertyEstimatorOptions())

    # merge_order_a = graph.topological_sort(density_calculation_a.dependants_graph)
    # merge_order_b = graph.topological_sort(density_calculation_b.dependants_graph)
    #
    # for protocol_id_A, protocol_id_B in zip(merge_order_a, merge_order_b):
    #
    #     print('Trying to merge {} and {}.'.format(protocol_id_A, protocol_id_B))
    #
    #     protocol_a = density_calculation_a.protocols[protocol_id_A]
    #     protocol_b = density_calculation_b.protocols[protocol_id_B]
    #
    #     can_merge = protocol_a.can_merge(protocol_b)
    #     assert can_merge
    #
    #     density_calculation_b.replace_protocol(protocol_b, protocol_a)
    #     merged_ids = protocol_a.merge(protocol_b)
    #
    #     for protocol_id in density_calculation_b.protocols:
    #
    #         for old_id, new_id in merged_ids.items():
    #             density_calculation_b.protocols[protocol_id].replace_protocol(old_id, new_id)

    calculation_graph = DirectCalculationGraph('')

    calculation_graph.add_calculation(density_calculation_a)
    calculation_graph.add_calculation(density_calculation_b)

    merge_order_a = graph.topological_sort(density_calculation_a.dependants_graph)
    merge_order_b = graph.topological_sort(density_calculation_b.dependants_graph)

    assert len(calculation_graph._nodes_by_id) == len(density_calculation_a.protocols)

    for protocol_id in density_calculation_a.protocols:
        assert protocol_id in calculation_graph._nodes_by_id

    for protocol_id_A, protocol_id_B in zip(merge_order_a, merge_order_b):

        assert density_calculation_a.protocols[protocol_id_A].schema.json() == \
               density_calculation_b.protocols[protocol_id_B].schema.json()


def test_dielectric_merging():

    substance = Mixture()
    substance.add_component('C', 1.0)

    dummy_property = DielectricConstant(thermodynamic_state=ThermodynamicState(temperature=298*unit.kelvin,
                                                                               pressure=1*unit.atmosphere),
                                        phase=PropertyPhase.Liquid,
                                        substance=substance,
                                        value=10*unit.gram/unit.mole,
                                        uncertainty=1*unit.gram/unit.mole,
                                        id=str(uuid.uuid4()))

    calculation_schema = dummy_property.get_default_calculation_schema()

    dielectric_calculation_a = DirectCalculation(dummy_property,
                                                 get_data_filename('forcefield/smirnoff99Frosst.offxml'),
                                                 calculation_schema,
                                                 PropertyEstimatorOptions())

    dielectric_calculation_b = DirectCalculation(dummy_property,
                                                 get_data_filename('forcefield/smirnoff99Frosst.offxml'),
                                                 calculation_schema,
                                                 PropertyEstimatorOptions())

    calculation_graph = DirectCalculationGraph('')

    calculation_graph.add_calculation(dielectric_calculation_a)
    calculation_graph.add_calculation(dielectric_calculation_b)

    merge_order_a = graph.topological_sort(dielectric_calculation_a.dependants_graph)
    merge_order_b = graph.topological_sort(dielectric_calculation_b.dependants_graph)

    assert len(calculation_graph._nodes_by_id) == len(dielectric_calculation_a.protocols)

    for protocol_id in dielectric_calculation_a.protocols:
        assert protocol_id in calculation_graph._nodes_by_id

    for protocol_id_A, protocol_id_B in zip(merge_order_a, merge_order_b):
        assert dielectric_calculation_a.protocols[protocol_id_A].schema.json() == \
               dielectric_calculation_b.protocols[protocol_id_B].schema.json()


def test_density_dielectric_merging():

    substance = Mixture()
    substance.add_component('C', 1.0)

    density = Density(thermodynamic_state=ThermodynamicState(temperature=298*unit.kelvin,
                                                             pressure=1*unit.atmosphere),
                      phase=PropertyPhase.Liquid,
                      substance=substance,
                      value=10*unit.gram/unit.mole,
                      uncertainty=1*unit.gram/unit.mole,
                      id=str(uuid.uuid4()))

    dielectric = DielectricConstant(thermodynamic_state=ThermodynamicState(temperature=298*unit.kelvin,
                                                                           pressure=1*unit.atmosphere),
                                    phase=PropertyPhase.Liquid,
                                    substance=substance,
                                    value=10*unit.gram/unit.mole,
                                    uncertainty=1*unit.gram/unit.mole,
                                    id=str(uuid.uuid4()))

    density_schema = density.get_default_calculation_schema()
    dielectric_schema = dielectric.get_default_calculation_schema()

    density_calculation = DirectCalculation(density,
                                            get_data_filename('forcefield/smirnoff99Frosst.offxml'),
                                            density_schema,
                                            PropertyEstimatorOptions())

    dielectric_calculation = DirectCalculation(dielectric,
                                               get_data_filename('forcefield/smirnoff99Frosst.offxml'),
                                               dielectric_schema,
                                               PropertyEstimatorOptions())

    calculation_graph = DirectCalculationGraph('')

    calculation_graph.add_calculation(density_calculation)
    calculation_graph.add_calculation(dielectric_calculation)

    merge_order_a = graph.topological_sort(density_calculation.dependants_graph)
    merge_order_b = graph.topological_sort(dielectric_calculation.dependants_graph)

    for protocol_id_A, protocol_id_B in zip(merge_order_a, merge_order_b):

        if protocol_id_A.find('extract_traj') < 0:

            assert density_calculation.protocols[protocol_id_A].schema.json() == \
                   dielectric_calculation.protocols[protocol_id_B].schema.json()

        else:

            assert density_calculation.protocols[protocol_id_A].schema.json() != \
                   dielectric_calculation.protocols[protocol_id_B].schema.json()


def test_protocol_decorators():
    pass
    # build_coordinates = BuildCoordinatesPackmol('build_coordinates')
    #
    # build_coordinates.substance = ProtocolPath('substance', 'global')
    #
    # value = getattr(BuildCoordinatesPackmol, 'substance')
    #
    # assert value is None


def test_simulation_layer():
    """Manually test the simulation layer"""

    # if path.isdir('property-data'):
    #     shutil.rmtree('property-data')
    #
    # # Set up time-based logging to help debug threading issues.
    # formatter = logging.Formatter(fmt='%(asctime)s.%(msecs)03d %(levelname)-8s %(message)s',
    #                               datefmt='%H:%M:%S')
    #
    # screen_handler = logging.StreamHandler(stream=sys.stdout)
    # screen_handler.setFormatter(formatter)
    #
    # logger = logging.getLogger()
    # logger.setLevel(logging.INFO)
    # logger.addHandler(screen_handler)
    #
    # dummy_pickle = b''
    # data_model = pickle.loads(dummy_pickle)
    #
    # backend = DaskLocalClusterBackend(1, 1)
    # backend.start()
    #
    # def dummy_callback(*args, **kwargs):
    #     print(args, kwargs)
    #     pass
    #
    # simulation_layer = SimulationLayer()
    # simulation_layer.schedule_calculation(backend, data_model, {}, dummy_callback, True)
    pass
