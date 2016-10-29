# Testing examples for SMARTY and ForceField

This directory contains various playgrounds/tests for different aspects of SMARTY, especially the ForceField class structure, as well as for related issues.

## Manifest:
- test_impropers: IPython notebooks and other functionality relating to testing whether impropers are working properly (and consistent with AMBER energetics) for simple cases. Ultimately these were turned into standard tests for smarty
- OE_bondorder_testing.ipynb: IPython notebook testing relating to a possible extension of OpenMM's Topology to store bond orders; this tests whether we can successfully do a round trip of an OEMol where we store just connectivity, elements, atom names, and bond orders, and then re-construct the original molecule (including formal charge) from this information.
