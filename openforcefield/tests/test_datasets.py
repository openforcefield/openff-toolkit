from openforcefield.datasets import ThermoMLDataSet

print('Loading files from url:\n\n')

dataset = ThermoMLDataSet.from_url('https://trc.nist.gov/journals/jct/2005v37/i04/j.jct.2004.09.022.xml')

print('\n')

dataset = ThermoMLDataSet.from_url('https://trc.nist.gov/journals/jct/2005v37/i05/j.jct.2004.09.014.xml')

print('\n\nTrying to load from non-existent url:\n\n')

dataset = ThermoMLDataSet.from_url('https://trc.nist.gov/journals/jct/2005v37/i04/j.jct.2004.09.022.xmld')

print('\n\nTrying to load as local file:\n\n')

dataset = ThermoMLDataSet.from_file('../data/properties/binary.xml')

print('\n\nTrying to load from non-existent local file:\n\n')

dataset = ThermoMLDataSet.from_file('../data/properties/j.jct.2004.09.022.xmld')
