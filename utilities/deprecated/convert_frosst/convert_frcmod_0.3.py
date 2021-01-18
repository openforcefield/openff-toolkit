#!/bin/env python

import numpy as np
from openff.toolkit.typing.engines.smirnoff.forcefield import ForceField

# Function definitions for parsing sections within parameter file
def _parse_nonbon_line( line ):
    """Parse an AMBER frcmod nonbon line and return relevant parameters in a dictionary. AMBER uses rmin_half and epsilon in angstroms and kilocalories per mole."""
    tmp = line.split()
    params = {}
    params['smirks'] = tmp[0]
    params['rmin_half'] = tmp[1] + ' * angstrom'
    params['epsilon'] = tmp[2] + ' * kilocalories_per_mole'

    return params

def _parse_bond_line( line ):
    """Parse an AMBER frcmod BOND line and return relevant parameters in a dictionary. AMBER uses length and force constant, with the factor of two dropped. Here we multiply by the factor of two before returning. Units are angstroms and kilocalories per mole per square angstrom."""

    tmp = line.split()
    params = {}
    params['smirks'] = tmp[0]
    params['k'] = str(2*float(tmp[1])) + "* kilocalories_per_mole/angstrom**2"
    params['length'] = tmp[2] + " * angstrom"
    return params

def _parse_angl_line( line ):
    """Parse an AMBER frcmod ANGL line and return relevant parameters in a dictionary. AMBER uses angle and force constant, with the factor of two dropped. Here we multiply by the factor of two before returning. Units are degrees and kilocalories per mole."""

    tmp = line.split()
    params = {}
    params['smirks'] = tmp[0]
    params['k'] = str(2*float(tmp[1])) + " * kilocalories_per_mole/radian**2"
    params['angle'] = tmp[2] + " * degree"
    return params

def _parse_dihe_line( line ):
    """Parse an AMBER frcmod DIHE line and return relevant parameters in a dictionary. Units for k are kilocalories per mole."""

    tmp = line.split()
    params = {}
    params['smirks'] = tmp[0]
    params['idivf1'] = tmp[1]
    params['k1'] = tmp[2] + " * kilocalories_per_mole"
    params['phase1'] = tmp[3] + " * degree"
    params['periodicity1'] = str(int(np.abs(float(tmp[4]))))
    return params

def _parse_impr_line( line ):
    """Parse an AMBER frcmod DIHE line and return relevant parameters in a dictionary. Units for k are kilocalories per mole."""

    tmp = line.split()
    params = {}
    params['smirks'] = tmp[0]
    params['k1'] = tmp[1] + " * kilocalories_per_mole"
    params['phase1'] = tmp[2] + " * degree"
    params['periodicity1'] = str(int(np.abs(float(tmp[3]))))
    return params


def parse_frcmod(frcmod_file_name):
    """
    Parse a smirnoffish FRCMOD file, returning a dict of lists containing the individual parameters that would go
    into a "SMIRNOFF data" dictionary. When combined with section headers, this dictionary is suitable to be fed
    into the ForceField initializer

    Parameters
    ----------
    frcmod_file_name : str
        path to a smirnoffish FRCMOD file

    Returns
    -------
    parameter_lists : dict of list
        Hierarchical dict of lists, containing parameter entries compliant with the SMIRNOFF spec version 0.2
    author : str or None
        Contents of the Author field, if any
    date: str or None
        Contents of the Date field, if any

    """

    # Obtain sections from target file
    file = open(frcmod_file_name, 'r')
    text = file.readlines()
    file.close()
    sections = {}
    # Section names from frcmod which we will parse
    sec_names = ['NONBON', 'BOND', 'ANGL', 'IMPR', 'DIHE']
    sec_name_2_param_prefix = {'NONBON':'n' , 'BOND':'b', 'ANGL':'a', 'DIHE':'t', 'IMPR':'i'}
    sec_name_2_param_index = dict((sn, 1) for sn in sec_names)
    # Tags that will be used in the FFXML for these (same order)
    force_tags = ['Atom', 'Bond', 'Angle', 'Improper', 'Proper']
    # Force names in the FFXML (same order)
    force_sections = ['vdW', 'Bonds', 'Angles', 'ImproperTorsions', 'ProperTorsions']
    sec_name_2_force_tag = dict((sn, t) for sn, t in zip(sec_names, force_tags))
    sec_name_2_force_section = dict((sn, fs) for sn, fs in zip(sec_names, force_sections))

    date = None
    author = None
    for line in text:
    #while ct < len(text):
        #line = text[ct]
        linesp = line.split()

        # Skip lines starting with comment or which are blank
        if line[0]=='#' or len(linesp) < 1:
            continue

        # Check first entry to see if it's a section name, if so initialize storage
        if linesp[0] in sec_names:
            this_sec = linesp[0]
            sections[this_sec] = []

        elif linesp[0] in ['DATE', 'AUTHOR']:
            this_sec = linesp[0]
        # Otherwise store
        else:
            if this_sec == 'DATE':
                date = line.strip()
            elif this_sec == 'AUTHOR':
                author = line.strip()
            else:
                sections[this_sec].append(line)


    parameter_lists = {}

    # Use functions to parse sections from target file and add parameters to force field
    for (sec_idx, sec_name) in enumerate(sec_names):
        param_list = []
        # sections[sec_name] = []
        for line in sections[sec_name]:
            # Parse line for parameters
            if sec_name=='NONBON':
                param = _parse_nonbon_line(line)
            elif sec_name=='BOND':
                param = _parse_bond_line(line)
            elif sec_name=='DIHE':
                param = _parse_dihe_line(line)
            elif sec_name=='IMPR':
                param = _parse_impr_line(line)
            elif sec_name=='ANGL':
                param = _parse_angl_line(line)

            # Add parameter ID
            param_prefix = sec_name_2_param_prefix[sec_name]
            param_index = sec_name_2_param_index[sec_name]
            param['id'] = param_prefix + str( param_index )

            # If we're dealing with a simple parameter, just append it to the list
            if not (sec_name in ['IMPR', 'DIHE']):
                param_list.append(param)
                # Increment parameter index
                sec_name_2_param_index[sec_name] = sec_name_2_param_index[sec_name] + 1

            # If we're dealing with a potentially multi-term parameter, check if this SMIRKS is already in the list
            else:
                param_matched = False
                for existing_param in param_list:
                    if existing_param['smirks'] != param['smirks']:
                        continue
                    param_matched = True
                    # Find the lowest unoccupied torsion term index
                    term_index = 1
                    while f'k{term_index}' in existing_param:
                        term_index += 1
                    existing_param[f'k{term_index}'] = param['k1']
                    existing_param[f'phase{term_index}'] = param['phase1']
                    existing_param[f'periodicity{term_index}'] = param['periodicity1']
                    existing_param[f'idivf{term_index}'] = param['idivf1']
                    break

                # If the SMIRKS isn't already known, initialize this as a new parameter, since it could be the
                # first term of a multiterm torsion
                if not(param_matched):
                    param_list.append(param)
                    # Increment parameter index ONLY for the first term of a torsion that's found.
                    sec_name_2_param_index[sec_name] = sec_name_2_param_index[sec_name] + 1

        force_section = sec_name_2_force_section[sec_name]
        force_tag = sec_name_2_force_tag[sec_name]
        parameter_lists[force_section] = {force_tag: param_list.copy()}

    return parameter_lists, author, date

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
    from openff.toolkit.typing.engines.smirnoff import XMLParameterIOHandler



    io_handler = XMLParameterIOHandler()
    smirnoff_data = io_handler.parse_file(inxml)
    parameter_lists, author, date = parse_frcmod(infile)
    for section_tag in parameter_lists:
        # Nest the actual parameter lists one level deeper, by their parameter_list_tag (eg Bonds > Bond > List)
        parameter_list_tag = list(parameter_lists[section_tag].keys())[0]
        smirnoff_data['SMIRNOFF'][section_tag][parameter_list_tag] = parameter_lists[section_tag][parameter_list_tag]
    ff=ForceField()
    ff._load_smirnoff_data(smirnoff_data)
    print(ff.to_string())

    # TODO: Add author and date

    ff.to_file(outxml, io_format='XML')
    # # Write SMIRNOFF XML file
    # ff.writeFile(outxml)
    #
    # # Roundtrip to fix formatting (for some reason etree won't format it properly on first write after modification)
    # tmp = ForceField(outxml)
    # tmp.writeFile(outxml)


if __name__=="__main__":
    from optparse import OptionParser
    usage_string="""\
    Convert specified SMIRKS-ified AMBER frcmod file into SMIRNOFF FFXML format, inserting converted parameters into a template FFXML file and writing to a new output file.

    usage: convert_frcmod_0.2.py --frcmod test.frcmod --template template.offxml --xml test.offxml
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
