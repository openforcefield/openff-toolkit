# Property Estimator Data

This directory contains data for running the PropertyEstimator class. 

## Manifest

* **ThermoML.xsd**: ThermoML Schema downloaded from [NIST](http://trc.nist.gov/ThermoML.xsd). Current version from 
  Jan 26, 2018 fetched by `curl` (OSX) or `wget` (Linux).
* **thermoml_schema module**: Mini-module which tries to import a thermoml_schema compatible with your installed 
  pyxb version without th eneed to write `try...except` blocks in your own import statements
* **thermoml_schema/thermoml_schema12\*.py**: Importable pyxb generated schemas created from command 
  `pyxbgen -u ThermoML.xsd -m thermoml_schema12*.py` where `*` is replaced with pyxb minor version number.
   * e.g. `thermoml_schema125.py` was generated from pyxb version 1.2.5.
   * Pyxb requires is neither forward nor reverse compatible and requires specific versions