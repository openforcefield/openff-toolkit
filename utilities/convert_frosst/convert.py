#!/bin/env python

import lxml.etree as etree
import numpy as np
from smarty import ForceField

# Define what file we will parse and where we will write to
infile = 'example.frcmod'
inxml = 'template.ffxml'
outxml = 'example.ffxml'

# Function definitions for parsing sections within parameter file
def parse_nonbon_line( line ):
    """Parse an AMBER frcmod nonbon line and return relevant parameters in a dictionary. AMBER uses rmin_half and epsilon in angstroms and kilocalories per mole."""
    tmp = line.split()
    params = {}
    params['smirks'] = tmp[0]
    params['rmin_half'] = tmp[1]
    params['epsilon'] = tmp[2]

    return params

def parse_bond_line( line ):
    """Parse an AMBER frcmod BOND line and return relevant parameters in a dictionary. AMBER uses length and force constant, with the factor of two dropped. Here we multiply by the factor of two before returning. Units are angstroms and kilocalories per mole per square angstrom."""

    tmp = line.split()
    params = {}
    params['smirks'] = tmp[0]
    params['k'] = str(2*float(tmp[1]))
    params['length'] = tmp[2]
    return params

def parse_angl_line( line ):
    """Parse an AMBER frcmod ANGL line and return relevant parameters in a dictionary. AMBER uses angle and force constant, with the factor of two dropped. Here we multiply by the factor of two before returning. Units are degrees and kilocalories per mole."""

    tmp = line.split()
    params = {}
    params['smirks'] = tmp[0]
    params['k'] = str(2*float(tmp[1]))
    params['angle'] = tmp[2]
    return params

def parse_dihe_line( line ):
    """Parse an AMBER frcmod DIHE line and return relevant parameters in a dictionary. Units for k are kilocalories per mole."""

    tmp = line.split()
    params = {}
    params['smirks'] = tmp[0]
    params['idivf1'] = tmp[1]
    params['k1'] = tmp[2]
    params['phase1'] = tmp[3]
    params['periodicity1'] = str(int(np.abs(float(tmp[4]))))
    return params


def parse_impr_line( line ):
    """Parse an AMBER frcmod DIHE line and return relevant parameters in a dictionary. Units for k are kilocalories per mole."""

    tmp = line.split()
    params = {}
    params['smirks'] = tmp[0]
    params['k1'] = tmp[1]
    params['phase1'] = tmp[2]
    params['periodicity1'] = str(int(np.abs(float(tmp[3]))))
    return params

# Obtain sections from target file
file = open(infile, 'r')
text = file.readlines()
file.close()
sections = {}
# Section names from frcmod which we will parse
secnames = ['NONBON', 'BOND', 'ANGL', 'IMPR', 'DIHE']
# Tags that will be used in the FFXML for these (same order)
tag = ['Atom', 'Bond', 'Angle', 'Improper', 'Proper']
# Force names in the FFXML (same order)
force_section = ['NonbondedForce', 'HarmonicBondForce', 'HarmonicAngleForce', 'PeriodicTorsionForce', 'PeriodicTorsionForce']
ct = 0
thissec = None
while ct < len(text):
    line = text[ct]
    tmp = line.split()

    # Skip lines starting with comment or which are blank
    if line[0]=='#' or len(tmp) < 1:
        ct+=1
        continue

    # Check first entry to see if it's a section name, if so initialize storage
    if tmp[0] in secnames:
        print(tmp[0])
        thissec = tmp[0]
        sections[thissec] = []
    # Otherwise store
    else:
        sections[thissec].append(line)

    ct+=1


# Read template forcefield file
ff = ForceField(inxml)
# Use functions to parse sections from target file and add parameters to force field
for (idx, name) in enumerate(secnames):
    for line in sections[name]:
        # Parse line for parameters
        if name=='NONBON':
            params = parse_nonbon_line(line)
        elif name=='BOND':
            params = parse_bond_line(line)
        elif name=='DIHE':
            params = parse_dihe_line(line)
        elif name=='IMPR':
            params = parse_impr_line(line)
        elif name=='ANGL':
            params = parse_angl_line(line)

        smirks = params['smirks']

        # If it's not a torsion, just store in straightforward way
        if not (name=='IMPR' or name=='DIHE'):
            ff.addParameter( params, smirks, force_section[idx], tag[idx] )
        # If it's a torsion, check to see if there are already parameters and
        # if so, add a new term to this torsion
        else:
            # If we have parameters already
            oldparams = ff.getParameter(smirks, force_type=force_section[idx])
            if oldparams:
                # Find what number to use
                idnr = 1
                paramtag = 'k%s' % idnr
                while paramtag in params:
                    idnr+=1
                    paramtag = 'k%s' % idnr
                # Construct new param object with updated numbers
                for paramtag in ('periodicity1', 'phase1', 'idivf1', 'k1'):
                    if paramtag in params:
                        val = params.pop(paramtag)
                        # Update old parameters with this addition
                        # CURRENTLY THIS DOESN'T WORK SINCE SETPARAMETER REQUIRES SAME KEYS AS OLD PARAMETER
                        oldparams[paramtag[:-1]+str(idnr) ] = val
                # Store
                ff.setParameter( oldparams, smirks=smirks, force_type=force_section[idx])
            # Otherwise, just store new parameters
            ff.addParameter( params, smirks, force_section[idx], tag[idx])

  # TO DO: ADD CHECKS TO MAKE SURE NO DUPLICATE PARAMS

# Write SMIRFF XML file
ff.writeFile(outxml)
