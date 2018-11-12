from openforcefield.datasets import ThermoMLDataSet

print('Loading real file:\n\n')
dataset = ThermoMLDataSet.from_url('https://trc.nist.gov/journals/jct/2005v37/i04/j.jct.2004.09.022.xml')

print('\n\nTrying to load non-existant file:\n\n')

dataset = ThermoMLDataSet.from_url('https://trc.nist.gov/journals/jct/2005v37/i04/j.jct.2004.09.022.xmld')

print('\n\nTrying to load as local file:\n\n')

with open('../data/properties/j.jct.2004.09.022.xml') as file:

    dataset = ThermoMLDataSet.from_xml(file.read())