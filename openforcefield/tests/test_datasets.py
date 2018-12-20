from os import listdir
from os.path import isfile, join

import logging

from pydantic import ValidationError

from openforcefield.utils import get_data_filename

from openforcefield.properties import PhysicalProperty
from openforcefield.properties.datasets import ThermoMLDataSet

# TODO: Add tests for specific ThermoML data sets that give 100% coverage.
# These may need to be hand written.


def test_from_url():

    data_set = ThermoMLDataSet.from_url_list('https://trc.nist.gov/journals/jct/2005v37/i04/j.jct.2004.09.022.xml')
    assert data_set is not None

    assert len(data_set.properties) > 0

    data_set = ThermoMLDataSet.from_url_list('https://trc.nist.gov/journals/jct/2005v37/i04/j.jct.2004.09.022.xmld')
    assert data_set is None


def test_serialization():

    data_set = ThermoMLDataSet.from_doi_list('10.1016/j.jct.2016.10.001')
    assert data_set is not None

    assert len(data_set.properties) > 0

    for mixture_tag in data_set.properties:

        for physical_property in data_set.properties[mixture_tag]:
            physical_property_json = physical_property.json()
            print(physical_property_json)

            physical_property_recreated = PhysicalProperty.parse_raw(physical_property_json)
            print(physical_property_recreated)


def test_from_doi():

    data_set = ThermoMLDataSet.from_doi_list('10.1016/j.jct.2016.10.001')
    assert data_set is not None

    assert len(data_set.properties) > 0

    for mixture_tag in data_set.properties:

        for physical_property in data_set.properties[mixture_tag]:

            physical_property_json = physical_property.json()
            print(physical_property_json)

            physical_property_recreated = PhysicalProperty.parse_raw(physical_property_json)
            print(physical_property_recreated)

    data_set = ThermoMLDataSet.from_doi_list('10.1016/j.jct.2016.12.009')
    assert data_set is None

    data_set = ThermoMLDataSet.from_doi_list('10.1016/j.jct.2016.12.009x')
    assert data_set is None


def test_from_files():

    data_set = ThermoMLDataSet.from_file_list(get_data_filename('properties/j.jct.2004.09.014.xml'),
                                              get_data_filename('properties/j.jct.2004.09.022.xml'),
                                              get_data_filename('properties/j.jct.2007.09.004.xml'))
    assert data_set is not None

    assert len(data_set.properties) > 0

    data_set = ThermoMLDataSet.from_file_list('../data/properties/j.jct.2004.09.014.xmld')
    assert data_set is None


def parse_all_jct_files():

    logging.basicConfig(filename='data_sets.log', filemode='w', level=logging.INFO)

    data_path = '../data/properties/JCT'
    thermoml_files = []

    for file_path in listdir(data_path):

        full_path = join(data_path, file_path)

        if not isfile(full_path):
            continue

        thermoml_files.append(full_path)

    data_set = ThermoMLDataSet.from_file_list(*thermoml_files)


if __name__ == "__main__":
    test_from_doi()
