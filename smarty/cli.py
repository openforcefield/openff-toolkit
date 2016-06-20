"""
Command-line driver example for SMARTY.

"""

import sys
import string

from optparse import OptionParser # For parsing of command line arguments

import os
import math
import copy
import re
import numpy
import random

import smarty

def main():
    # Create command-line argument options.
    usage_string = """\
    Sample over atom types, optionally attempting to match atom types in a reference typed set of molecules.

    usage: %prog --basetypes smartsfile --decorators smartsfile [--substitutions smartsfile] --molecules molfile [--reference molfile] --iterations niterations [--temperature temperature]

    example:

    python %prog --basetypes=atomtypes/basetypes.smarts --decorators=atomtypes/decorators.smarts --substitutions=atomtypes/substitutions.smarts \
        --molecules=molecules/zinc-subset-tripos.mol2.gz --reference=molecules/zinc-subset-parm@frosst.mol2.gz --iterations 1000 --temperature=0.1

    """
    version_string = "%prog %__version__"
    parser = OptionParser(usage=usage_string, version=version_string)

    parser.add_option("-b", "--basetypes", metavar='BASETYPES',
                      action="store", type="string", dest='basetypes_filename', default=None,
                      help="Filename defining base atom types as SMARTS atom matches.")

    parser.add_option("-d", "--decorators", metavar='DECORATORS',
                      action="store", type="string", dest='decorators_filename', default=None,
                      help="Filename defining decorator atom types as SMARTS atom matches.")

    parser.add_option("-s", "--substitutions", metavar="SUBSTITUTIONS",
                      action="store", type="string", dest='substitutions_filename', default=None,
                      help="Filename defining substitution definitions for SMARTS atom matches (OPTIONAL).")

    parser.add_option("-r", "--reference", metavar="REFMOL",
                      action="store", type="string", dest='reference_molecules_filename', default=None,
                      help="Reference typed molecules for computing likelihood (must match same molecule and atom ordering in molecules file) (OPTIONAL).")

    parser.add_option("-m", "--molecules", metavar='MOLECULES',
                      action="store", type="string", dest='molecules_filename', default=None,
                      help="Small molecule set (in any OpenEye compatible file format) containing 'dG(exp)' fields with experimental hydration free energies.")

    parser.add_option("-i", "--iterations", metavar='ITERATIONS',
                      action="store", type="int", dest='iterations', default=150,
                      help="MCMC iterations.")

    parser.add_option("-t", "--temperature", metavar='TEMPERATURE',
                      action="store", type="float", dest='temperature', default=0.1,
                      help="Effective temperature for Monte Carlo acceptance, indicating fractional tolerance of mismatched atoms (default: 0.1). If 0 is specified, will behave in a greedy manner.")

    parser.add_option("-l", '--trajectory', metavar="TRAJECTORY_FILE",
            action = "store", dest = "traj_file", default = "trajectory.csv",
            help = "Name for trajectory file output, trajectory saves only changes to the list of 'atomtypes' for each iteration. For now, if the file name already exists, it just won't create a trajectory file")

    verbose = True

    # Parse command-line arguments.
    (options,args) = parser.parse_args()

    # Ensure all required options have been specified.
    if (options.basetypes_filename is None) or (options.decorators_filename is None) or (options.molecules_filename is None):
        parser.print_help()
        parser.error("All input files must be specified.")

    # Load and type all molecules in the specified dataset.
    import smarty.utils
    molecules = smarty.utils.read_molecules(options.molecules_filename, verbose=True)

    # Read reference typed molecules, if specified.
    reference_typed_molecules = None
    if options.reference_molecules_filename is not None:
        reference_typed_molecules = smarty.utils.read_molecules(options.reference_molecules_filename, verbose=True)

    # Construct atom type sampler.
    atomtype_sampler = smarty.AtomTypeSampler(molecules, options.basetypes_filename, options.decorators_filename, replacements_filename=options.substitutions_filename, reference_typed_molecules=reference_typed_molecules, verbose=verbose, temperature=options.temperature)

    # Start sampling atom types.
    atomtype_sampler.run(options.iterations, options.traj_file)
