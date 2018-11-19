import logging
import simtk.unit as unit

from openforcefield import substances

from openforcefield.datasets import ThermoMLDataSet

from openforcefield.typing.engines import smirnoff
from openforcefield.utils import get_data_filename

from openforcefield.propertycalculator import protocols


def test_mixture_protocol():

    logging.getLogger().setLevel(logging.INFO)

    data_set = ThermoMLDataSet.from_file_list('../data/properties/single_density.xml')
    force_field = smirnoff.ForceField(get_data_filename('forcefield/smirnoff99Frosst.offxml'))

    mixture_protocol = protocols.BuildLiquidCoordinates(512,
                                                        1.0 * unit.grams / unit.milliliters)

    topology_protocol = protocols.BuildSmirnoffTopology()
    energy_minimisation_protocol = protocols.RunEnergyMinimisation()
    nvt_protocol = protocols.RunNVTSimulation()
    npt_protocol = protocols.RunNPTSimulation()

    protocol_data = protocols.ProtocolData()

    protocol_data.root_directory = 'test_data/'

    protocol_data.substance = data_set.measured_properties[0].substance
    protocol_data.thermodynamic_state = data_set.measured_properties[0].thermodynamic_state
    protocol_data.force_field = force_field

    protocol_data = mixture_protocol.execute(protocol_data)

    if protocol_data is None:
        print('Failed to set-up the initial coordinates / topology.')

    protocol_data = topology_protocol.execute(protocol_data)

    if protocol_data is None:
        print('Failed to set-up the system using SMIRNOFF.')

    protocol_data = energy_minimisation_protocol.execute(protocol_data)

    if protocol_data is None:
        print('Failed to minimise the system energy.')

    protocol_data = nvt_protocol.execute(protocol_data)

    if protocol_data is None:
        print('Failed to run the NVT simulation.')

    protocol_data = npt_protocol.execute(protocol_data)

    if protocol_data is None:
        print('Failed to run the NPT simulation.')

    print(protocol_data)


test_mixture_protocol()