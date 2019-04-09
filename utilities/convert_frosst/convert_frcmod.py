#!/bin/env python

import numpy as np
from openforcefield.typing.engines.smirnoff.forcefield import ForceField
from openforcefield.typing.chemistry import environment

# Function definitions for parsing sections within parameter file
def _parse_nonbon_line( line ):
    """Parse an AMBER frcmod nonbon line and return relevant parameters in a dictionary. AMBER uses rmin_half and epsilon in angstroms and kilocalories per mole."""
    tmp = line.split()
    params = {}
    params['smirks'] = tmp[0]
    params['rmin_half'] = tmp[1]
    params['epsilon'] = tmp[2]

    return params

def _parse_bond_line( line ):
    """Parse an AMBER frcmod BOND line and return relevant parameters in a dictionary. AMBER uses length and force constant, with the factor of two dropped. Here we multiply by the factor of two before returning. Units are angstroms and kilocalories per mole per square angstrom."""

    tmp = line.split()
    params = {}
    params['smirks'] = tmp[0]
    params['k'] = str(2*float(tmp[1]))
    params['length'] = tmp[2]
    return params

def _parse_angl_line( line ):
    """Parse an AMBER frcmod ANGL line and return relevant parameters in a dictionary. AMBER uses angle and force constant, with the factor of two dropped. Here we multiply by the factor of two before returning. Units are degrees and kilocalories per mole."""

    tmp = line.split()
    params = {}
    params['smirks'] = tmp[0]
    params['k'] = str(2*float(tmp[1]))
    params['angle'] = tmp[2]
    return params

def _parse_dihe_line( line ):
    """Parse an AMBER frcmod DIHE line and return relevant parameters in a dictionary. Units for k are kilocalories per mole."""

    tmp = line.split()
    params = {}
    params['smirks'] = tmp[0]
    params['idivf1'] = tmp[1]
    params['k1'] = tmp[2]
    params['phase1'] = tmp[3]
    params['periodicity1'] = str(int(np.abs(float(tmp[4]))))
    return params

def _parse_impr_line( line ):
    """Parse an AMBER frcmod DIHE line and return relevant parameters in a dictionary. Units for k are kilocalories per mole."""

    tmp = line.split()
    params = {}
    params['smirks'] = tmp[0]
    params['k1'] = tmp[1]
    params['phase1'] = tmp[2]
    params['periodicity1'] = str(int(np.abs(float(tmp[3]))))
    return params

def add_date_and_author(inxml, date, author):
    """
    Updates the template xml file with the date and authors in the
    input Frcmodish file.
    Parameters
    ----------
    inxml: str, template xml file
    date: str, date from input Frcmod file
    author: str, author list from input Frcmod file
    """
    # read input file
    f = open(inxml,'r')
    input_lines = f.readlines()
    f.close()

    output_lines = list()
    # Save lines, only change those for Date and Author
    for l in input_lines:
        start = l.strip().split('>')[0]
        if start == '<Date':
            output_lines.append("<Date>%s</Date>\n" % date.strip())
        elif start == '<Author':
            output_lines.append("<Author>%s</Author>\n" % author.strip())
        else:
            output_lines.append(l)

    # write fixed lines to ffxml tempate
    f = open(inxml,'w')
    f.writelines(output_lines)
    f.close()


# Main conversion functionality
def convert_frcmod_to_ffxml( infile, inxml, outxml ):
    """Convert a modified AMBER frcmod (with SMIRKS replacing atom types) to SMIRNOFF ffxml format by inserting parameters into a template ffxml file.

    Parameters
    ----------
    infile : str
        File name of input SMIRKS-ified frcmod file containing parameters
    inxml : str
        File name of template SMIRNOFF FFXML file into which to insert these parameters.
    outxml : str
        File name of resulting output SMIRNOFF FFXML

    Notes:
    -------
    Input XML file will normally be the template of a SMIRNOFF XML file without any parameters present (but with requisite force types already specified).
    """

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
    # Why is this a while loop and not a for line in text loop?
    while ct < len(text):
        line = text[ct]
        tmp = line.split()

        # Skip lines starting with comment or which are blank
        if line[0]=='#' or len(tmp) < 1:
            ct+=1
            continue

        # Check first entry to see if it's a section name, if so initialize storage
        if tmp[0] in secnames:
            thissec = tmp[0]
            sections[thissec] = []

        elif tmp[0] in ['DATE','AUTHOR']:
            thissec = tmp[0]
        # Otherwise store
        else:
            if thissec == 'DATE':
                date = line.strip()
            elif thissec == 'AUTHOR':
                author = line.strip()
            else:
                sections[thissec].append(line)

        ct+=1

    # fix date and author in inxml file:
    add_date_and_author(inxml, date, author)
    # Read template forcefield file
    ff = ForceField(inxml)
    # Use functions to parse sections from target file and add parameters to force field
    param_id_by_section={}
    param_prefix_by_sec = {'NONBON':'n' , 'BOND':'b', 'ANGL':'a', 'DIHE':'t', 'IMPR':'i'}
    env_method_by_sec = {'NONBON': environment.AtomChemicalEnvironment,
            'BOND': environment.BondChemicalEnvironment,
            'ANGL': environment.AngleChemicalEnvironment,
            'DIHE': environment.TorsionChemicalEnvironment,
            'IMPR': environment.ImproperChemicalEnvironment}
    for (idx, name) in enumerate(secnames):
        param_id_by_section[name] = 1
        env_method = env_method_by_sec[name]
        for line in sections[name]:
            # Parse line for parameters
            if name=='NONBON':
                params = _parse_nonbon_line(line)
            elif name=='BOND':
                params = _parse_bond_line(line)
            elif name=='DIHE':
                params = _parse_dihe_line(line)
            elif name=='IMPR':
                params = _parse_impr_line(line)
            elif name=='ANGL':
                params = _parse_angl_line(line)

            # Add parameter ID
            params['id'] = param_prefix_by_sec[name]+str( param_id_by_section[name] )

            smirks = params['smirks']
            #Check smirks is valid for chemical enviroment parsing:
            env = env_method(smirks)

            # If it's not a torsion, just store in straightforward way
            if not (name=='IMPR' or name=='DIHE'):
                # Check for duplicates first
                if ff.getParameter( smirks, force_type = force_section[idx] ):
                    raise ValueError("Error: parameter for %s is already present in forcefield." % smirks )
                else:
                    ff.addParameter( params, smirks, force_section[idx], tag[idx] )

                # Increment parameter id
                param_id_by_section[name] +=1
            # If it's a torsion, check to see if there are already parameters and
            # if so, add a new term to this torsion
            else:
                # If we have parameters already
                oldparams = ff.getParameter(smirks, force_type=force_section[idx])
                if oldparams:
                    # Find what number to use
                    idnr = 1
                    paramtag = 'k%s' % idnr
                    # This was "while paramtag in params" so it was overwriting k2 etc.
                    while paramtag in oldparams:
                        idnr+=1
                        paramtag = 'k%s' % idnr
                    # Construct new param object with updated numbers
                    for paramtag in ('periodicity1', 'phase1', 'idivf1', 'k1'):
                        if paramtag in params:
                            val = params.pop(paramtag)
                            oldparams[paramtag[:-1]+str(idnr) ] = val
                    # Store
                    ff.setParameter( oldparams, smirks=smirks, force_type=force_section[idx])
                else:
                    # Otherwise, just store new parameters
                    ff.addParameter( params, smirks, force_section[idx], tag[idx])
                    # Increment parameter ID
                    param_id_by_section[name] += 1


    # Write SMIRNOFF XML file
    ff.writeFile(outxml)

    # Roundtrip to fix formatting (for some reason etree won't format it properly on first write after modification)
    tmp = ForceField(outxml)
    tmp.writeFile(outxml)


if __name__=="__main__":
    from optparse import OptionParser
    usage_string="""\
    Convert specified SMIRKS-ified AMBER frcmod file into SMIRNOFF FFXML format, inserting converted parameters into a template FFXML file and writing to a new output file.

    usage: convert_frcmod.py --frcmod test.frcmod --template template.offxml --xml test.offxml
    """
    parser = OptionParser(usage=usage_string)

    parser.add_option('-f', '--frcmod', type = "string", dest='infile', default = None, action="store", help="Name of input smirks-ified frcmod file.")
    parser.add_option('-t', '--template', type="string", dest='inxml', default = None, action ="store", help="Name of template SMIRNOFF offxml file.")
    parser.add_option('-o', '--xml', type="string", dest='outxml', default =None, action="store", help="Name of output SMIRNOFF offxml file.")
    (options,args) = parser.parse_args()

    if (options.infile is None) or (options.inxml is None) or (options.outxml is None):
        parser.print_help()
        parser.error("Input frcmod and template files and output FFXML file must be specified.")


    convert_frcmod_to_ffxml( options.infile, options.inxml, options.outxml )
