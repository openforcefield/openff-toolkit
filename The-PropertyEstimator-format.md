# Property calculation toolkit

Property calculation toolkit from the [Open Forcefield Consortium](http://openforcefield.org).

This toolkit provides an API for storing, manipulating, and computing measured physical properties from simulation data.

## API docs

### Physical property measurements

[**Physical property measurements**](https://en.wikipedia.org/wiki/Physical_property) are measured properties of a substance that provide some information about the physical parameters that define the interactions within the substance.
A physical property is defined by a combination of:
* A `Mixture` specifying the substance that the measurement was performed on
* A `ThermodynamicState` specifying the thermodynamic conditions under which the measurement was performed
* A `PhysicalProperty` is the physical property that was measured
* A `MeasurementMethod` specifying the kind of measurement that was performed

An example of each:
* `Mixture`: a 0.8 mole fraction mixture of ethanol and water
* `ThermodynamicState`: 298 kelvin, 1 atmosphere
* `PhysicalProperty`: mass density
* `MeasurementMethod`: vibrating tube method

#### Physical substances

We generally use the concept of a liquid or gas `Mixture`, which is a subclass of `Substance`.

A simple liquid has only one component:
```python
liquid = Mixture()
liquid.addComponent('water')
```

A binary mixture has two components:
```python
binary_mixture = Mixture()
binary_mixture.addComponent('water', mole_fraction=0.2)
binary_mixture.addComponent('methanol') # assumed to be rest of mixture if no mole_fraction specified
```

A ternary mixture has three components:
```python
ternary_mixture = Mixture()
ternary_mixture.addComponent('ethanol', mole_fraction=0.2)
ternary_mixture.addComponent('methanol', mole_fraction=0.2)
ternary_mixture.addComponent('water')
```

The infinite dilution of one solute within a solvent or mixture is also specified as a `Mixture`, where the solute has zero mole fraction:
```python
infinite_dilution = Mixture()
infinite_dilution.addComponent('phenol', impurity=True) # infinite dilution; one copy only of the impurity
infinite_dilution.addComponent('water')
```

For testing against synthetic data, we also make use of `isolatedMolecule` objects (also subclasses of `Substance`) that represent a single molecule:
```python
phenol = IsolatedMolecule(iupac='phenol')
ethane = IsolatedMolecule(smiles='CC')
```

You can iterate over the components in a mixture:
```python
for component in mixture.components:
    print (component.iupac_name, component.mole_fraction)
```
retrieve a component by name:
```python
component = mixture.components['ethanol']
```
or get the number of components in a mixture:
```python
ncomponents = mixture.ncomponents
```
or check if a component is an impurity:
```python
if component.impurity == True:
    ...
```

#### Thermodynamic states

A `ThermodynamicState` specifies a combination of thermodynamic parameters (e.g. temperature, pressure) at which a measurement is performed.
```python
from simtk import unit
thermodynamic_state = ThermodynamicState(pressure=500*unit.kilopascals, temperature=298.15*unit.kelvin)
```
We use the `simtk.unit` unit system from [OpenMM](http://openmm.org) for units (though we may later migrate to [`pint`](https://pint.readthedocs.io) for portability).

#### Measurement methods

A `MeasurementMethod` subclass has information specific to the particular method used to measure a property (such as experimental uncertainty guidance).

Some examples:
* `FlowCalorimetry` for `HeatCapacity` or `ExcessMolarEnthalpy`
* `VibratingTubeMethod` for `MassDensity`

#### Physical property measurements

A `MeasuredPhysicalProperty` is a combination of `Substance`, `ThermodynamicState`, and a unit-bearing measured property `value` and `uncertainty`:
```python
# Define mixture
mixture = Mixture()
mixture.addComponent('water', mole_fraction=0.2)
mixture.addComponent('methanol')
# Define thermodynamic state
thermodynamic_state = ThermodynamicState(pressure=500*unit.kilopascals, temperature=298.15*unit.kelvin)
# Define measurement
measurement = ExcessMolarEnthalpy(substance, thermodynamic_state, value=83.3863244*unit.kilojoules_per_mole, uncertainty=0.1220794866*unit.kilojoules_per_mole)
```
The various properties are all subclasses of `MeasuredPhysicalProperty` and generally follow the `<ePropName/>` ThermoML tag names.

Some examples of `MeasuredPhysicalProperty`:
* `MassDensity` - mass density
* `ExcessMolarEnthalpy` - excess partial apparent molar enthalpy
* `HeatCapacity` - molar heat capacity at constant pressure

There are also some examples of measured physical properties we use exclusively for testing against synthetic data:
* `MeanPotentialEnergy` - mean potential energy of a system
* `BondMoment` - the `n`th moment (mean `E[x]` if n==1, or central moment `E[(x-E[x])^n]` if n >= 2) of a specified bond of an isolated molecule
* `AngleMoment` - the `n`th moment (mean `E[x]` if n==1, or central moment `E[(x-E[x])^n]` if n >= 2) of a specified angle of an isolated molecule
* `TorsionMoment` - the `n`th circular moment (`E[e^{i*n*phi}]`) of a specified torsion of an isolated molecule

For example:
```python
molecule = IsolatedMolecule(smiles=smiles)
thermodynamic_state = ThermodynamicState(temperature=300*unit.kelvin)
mean_potential = MeanPotentialEnergy(molecule, thermodynamic_state, value=124.4*unit.kilojoules_per_mole, uncertainty=14.5*unit.kilojoules_per_mole)
bond_average = BondMoment(molecule, thermodynamic_state, value=1.12*unit.angstroms, uncertainty=0.02*unit.angstroms, moment=1, smirks='[#6:1]-[#6:2]')
bond_variance = BondMoment(molecule, thermodynamic_state, value=0.05*unit.angstroms**2, uncertainty=0.02*unit.angstroms**2, moment=2, smirks='[#6:1]-[#6:2]')
angle_average = AngleMoment(molecule, thermodynamic_state, value=20*unit.degrees, uncertainty=0.05*unit.degrees, moment=1, smirks='[#6:1]-[#6:2]-[#6:3]')
torsion_moment = TorsionMoment(molecule, thermodynamic_state, value=(0.5, 0.3), uncertainty=(0.05, 0.03), moment=1, smirks='[#6:1]-[#6:2]-[#6:3]-[#6:4]')
```
**QUESTIONS**
* Is it useful to average over any bonds that match the SMARTS strings? How would you even construct those SMARTS strings algorithmically for the molecules in the set, to avoid matching anything else?
* Note that we need the first noncentral moment (average), so I had an exception where we don't compute the central moments for `moment == 1`
* Do we really want the variance (second central moment) instead of the standard deviation?

A [roadmap of physical properties to be implemented](https://github.com/open-forcefield-group/open-forcefield-tools/wiki/Physical-Properties-for-Calculation) is available.
Please raise an issue if your physical property of interest is not listed!

Each `MeasuredPhysicalProperty` has several properties:
* `.substance` - the `Mixture` for which the measurement was made
* `.thermodynamic_state` - the `ThermodynamicState` at which the measurement was made
* `.measurement_method` - the `MeasurementMethod` used to measure the physical property
* `.value` - the unit-bearing measurement value
* `.uncertainty` - the standard uncertainty of the measurement
* `.reference` - the literature reference (if available) for the measurement
* `.DOI` - the literature reference DOI (if available) for the measurement

The value, uncertainty, reference, and DOI do not necessarily need to be defined for a dataset in order for property calculations to be performed.

### Physical property datasets

A `PhysicalPropertyDataset` is a collection of `MeasuredPhysicalProperty` objects that are related in some way.
```python
dataset = PhysicalPropertyDataset([measurement1, measurement2])
```
The dataset is iterable:
```python
dataset = PhysicalPropertyDataset([measurement1, measurement2])
for measurement in dataset:
    print measurement.value
```
and has accessors to retrieve DOIs and references associated with measurements in the dataset:
```python
# Print the DOIs associated with this dataset
print(dataset.DOIs)
# Print the references associated with this dataset
print(dataset.references)
```

For convenience, you can retrieve the dataset as a [pandas `DataFrame`](https://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.html):
```python
dataset.to_pandas()
```

#### ThermoML datasets

A `ThermoMLDataset` object represents a physical property dataset stored in the IUPAC-standard [ThermoML](http://trc.nist.gov/ThermoMLRecommendations.pdf) for specifying thermodynamic properties in XML format.
`ThermoMLDataset` is a subclass of `PhysicalPropertyDataset`, and provides the same API interface (in addition to some ThermoML-specfic methods).
Direct access to the [NIST ThermoML Archive](http://trc.nist.gov/ThermoML.html) is supported for obtaining physical property measurements in this format directly from the NIST TRC repository.

For example, to retrieve [the ThermoML dataset](http://trc.boulder.nist.gov/ThermoML/10.1016/j.jct.2005.03.012) that accompanies [this paper](http://www.sciencedirect.com/science/article/pii/S0021961405000741), we can simply use the DOI `10.1016/j.jct.2005.03.012` as a key for creating a `PhysicalPropertyDataset` subclassed object from the ThermoML Archive:
```python
dataset = ThermoMLDataset(keys='10.1016/j.jct.2005.03.012')
```
You can also specify multiple ThermoML Archive keys to create a dataset from multiple ThermoML files:
```python
thermoml_keys = ['10.1021/acs.jced.5b00365', '10.1021/acs.jced.5b00474']
dataset = ThermoMLDataset(keys=thermoml_keys)
```

It is also possible to specify ThermoML datasets housed at other locations, such as
```python
dataset = ThermoMLDataset(url='http://openforcefieldgroup.org/thermoml-datasets')
```
or
```python
dataset = ThermoMLDataset(url='file:///Users/choderaj/thermoml')
```
or
```python
dataset = ThermoMLDataset(keys=['10.1021/acs.jced.5b00365', '10.1021/acs.jced.5b00474'], url='http://openforcefieldgroup.org/thermoml-datasets')
```
or from ThermoML and a different URL:
```python
dataset = ThermoMLDataset(keys=thermoml_keys)
dataset.retrieve(keys=local_keys, url='http://openforcefieldgroup.org/thermoml-datasets')
```

You can see which DOIs contribute to the current `ThermoMLDataset` with the convenience functions:
```python
print(dataset.DOIs)
```

NIST has compiled a JSON frame of corrections to uncertainties.
These can be used to update or correct data uncertainties and discard outliers using `applyNISTUncertainties()`:
```python
# Modify uncertainties according to NIST evaluation
dataset.applyNISTUncertainties(nist_uncertainties, adjust_uncertainties=True, discard_outliers=True)
```

**TODO:**
* We should merge any other useful parts parts of the [ThermoPyL API](https://github.com/choderalab/thermopyl) in here.

#### Other datasets

In future, we will add interfaces to other online datasets, such as
* [BindingDB](https://www.bindingdb.org/bind/index.jsp) for retrieving [host-guest binding affinity](https://www.bindingdb.org/bind/HostGuest.jsp) datasets.

### Estimating properties

The `PropertyEstimator` class creates objects that handle property estimation of all of the properties in a dataset, given a set or sets of parameters.
The implementation will isolate the user from whatever backend (local machine, HPC cluster, [XSEDE resources](http://xsede.org), [Amazon EC2](https://aws.amazon.com/ec2)) is being used to compute the properties, as well as whether new simulations are being launched and analyzed or existing simulation data is being reweighted.
Different backends will take different optional arguments, but here is an example that will launch and use 10 worker processes on a cluster:
```python
estimator = PropertyEstimator(nworkers=10) # NOTE: multiple backends will be supported in the future
computed_properties = estimator.computeProperties(dataset, parameter_sets)
```
Here, `dataset` is a `PhysicalPropertyDataset` or subclass, and `parameter_sets` is a list containing `ForceField` objects used to parameterize the physical systems in the dataset.
This can be a single parameter set or multiple (usually closely related) parameter sets.

`PropertyEstimator.computeProperties(...)` returns a list of `ComputedPhysicalProperty` objects that provide access to several pieces of information:
* `property.value` - the computed property value, with appropriate units
* `property.uncertainty` - the statistical uncertainty in the computed property
* `property.parameters` - a reference to the parameter set used to compute this property
* `property.property` - a reference to the corresponding `MeasuredPhysicalProperty` this property was computed for

**TODO**:
* How should we instruct `computeProperties()` to provide gradients (or components of gradients)?

This API can be extended in the future to provide access to the simulation data used to estimate the property, such as
```python
# Attach to my compute and storage resources
estimator = PropertyEstimator(...)
# Estimate some properties
computed_properties = estimator.computeProperties(dataset, parameters)
# Get statistics about simulation data that went into each property
for property in computed_properties:
   # Get statistics about simulation data that was reweighted to estimate this property
   for simulation in property.simulations:
      print('The simulation was %.3f ns long' % (simulation.length / unit.nanoseconds))
      print('The simulation was run at %.1f K and %.1f atm' % (simulation.thermodynamic_state.temperature / unit.kelvin, simulation.thermodynamic_state.pressure / unit.atmospheres))
      # Get the ParameterSet that was used for this simulation
      parameters = simulation.parameters
      # what else do you want...?
```

In future, we will want to use a parallel key/value database like [cassandra](http://cassandra.apache.org) to store simulations, along with a distributed task management system like [celery](http://www.celeryproject.org) with [redis](https://www.google.com/search?client=safari&rls=en&q=redis&ie=UTF-8&oe=UTF-8).

## API Usage Examples

### Using the high-level API

In this example, datasets are retrieved from the ThermoML and filtered to retain certain properties.
The corresponding properties for a given parameter set filename are then computed for a SMIRFF parameter set and printed.
```python
# Define the input datasets from ThermoML
thermoml_keys = ['10.1016/j.jct.2005.03.012', ...]
dataset = ThermoMLDataset(thermoml_keys)
# Filter the dataset to include only molar heat capacities measured between 280-350 K
dataset.filter(ePropName='Excess molar enthalpy (molar enthalpy of mixing), kJ/mol') # filter to retain only this property name
dataset.filter(VariableType='eTemperature', min=280*unit.kelvin, max=350*kelvin) # retain only measurements with `eTemperature` in specified range
# Load an initial parameter set
parameter_set = [ SMIRFFParameterSet('smarty-initial.xml') ]
# Compute physical properties for these measurements
estimator = PropertyEstimator(nworkers=10) # NOTE: multiple backends will be supported in the future
computed_properties = estimator.computeProperties(dataset, parameter_set)
# Write out statistics about errors in computed properties
for (computed, measured) in (computed_properties, dataset):
   property_unit = measured.value.unit
   print('%24s : experiment %8.3f (%.3f) | calculated %8.3f (%.3f) %s' % (measured.value / property_unit, measured.uncertainty / property_unit, computed.value / property_unit, computed.uncertainty / property_unit, str(property_unit))
```

### Using the low-level API

TBD
