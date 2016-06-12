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
    atomtype_sampler = smarty.AtomTypeSampler(molecules, options.basetypes_filename, options.decorators_filename, replacements_filename=options.substitutions_filename, reference_typed_molecules=reference_typed_molecules, verbose=verbose)

    # Start sampling atom types.
    atomtype_sampler.run(options.iterations)
