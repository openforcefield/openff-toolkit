from openforcefield.properties import Density, DielectricConstant
from openforcefield.properties.estimator import CalculationSchema


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