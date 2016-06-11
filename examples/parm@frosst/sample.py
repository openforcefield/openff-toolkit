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

if __name__=="__main__":

    # Create command-line argument options.
    usage_string = """\
    usage: %prog --basetypes smartsfile --decorators smartsfile --molecules molfile --iterations niterations

    example: %prog --basetypes atomtypes/basetypes.smarts --decorators atomtypes/decorators.smarts --molecules datasets/solvation.sdf --iterations 150

    """
    version_string = "%prog %__version__"
    parser = OptionParser(usage=usage_string, version=version_string)

    parser.add_option("-b", "--basetypes", metavar='BASETYPES',
                      action="store", type="string", dest='basetypes_filename', default='',
                      help="Filename defining base atom types as SMARTS atom matches.")

    parser.add_option("-d", "--decorators", metavar='DECORATORS',
                      action="store", type="string", dest='decorators_filename', default='',
                      help="Filename defining decorator atom types as SMARTS atom matches.")

    parser.add_option("-m", "--molecules", metavar='MOLECULES',
                      action="store", type="string", dest='molecules_filename', default='',
                      help="Small molecule set (in any OpenEye compatible file format) containing 'dG(exp)' fields with experimental hydration free energies.")

    parser.add_option("-i", "--iterations", metavar='ITERATIONS',
                      action="store", type="int", dest='iterations', default=150,
                      help="MCMC iterations.")

    # Parse command-line arguments.
    (options,args) = parser.parse_args()

    # Ensure all required options have been specified.
    if options.basetypes_filename=='' or options.decorators_filename=='' or options.molecules_filename=='':
        parser.print_help()
        parser.error("All input files must be specified.")

    # Load and type all molecules in the specified dataset.
    import smarty.utils
    molecules = smarty.utils.read_molecules(options.molecules_filename, verbose=True)

    # Construct atom type sampler.
    atomtype_sampler = smarty.AtomTypeSampler(options.basetypes_filename, options.decorators_filename, molecules)

    # Start sampling atom types.
    atomtype_sampler.run(options.iterations, verbose=True)
