# Testing system conversion for SMIRNOFF

This directory contains an IPython notebook and supporting files working to convert a system (provided in `smarty/data/systems` and `smarty/data/monomers`) into AMBER format after parameterization via SMIRNOFF

## Manifest:
- convert_to_amber.ipynb: IPython notebook working towards converting an example SMIRNOFF system (cyclohexane and ethanol mixture) into AMBER format and validating the energy before and after conversion
- ethanol and cyclohexane files: from `smarty/data/monomers` and `smarty/data/systems`
- `system.prmtop` and `system.crd`: output files from IPython notebook after conversion
