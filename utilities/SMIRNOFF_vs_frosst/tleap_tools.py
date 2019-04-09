import mdtraj.utils
import os
import shutil
import logging
logger = logging.getLogger(__name__)
from openmoltools.utils import getoutput

# Functionality below adapted from openmoltools.amber.run_tleap

def run_tleap(molecule_name, mol2_filename, frcmod_filename, prmtop_filename=None, inpcrd_filename=None, log_debug_output=False, verbose = True):
    """Run AmberTools tleap to create simulation files for AMBER
    Parameters
    ----------
    molecule_name : str
        The name of the molecule
    mol2_filename : str
        mol2 filename with target atom type names
    frcmod_filename : str
        Amber frcmod file produced by prmchk
    prmtop_filename : str, optional, default=None
        Amber prmtop file produced by tleap, defaults to molecule_name
    inpcrd_filename : str, optional, default=None
        Amber inpcrd file produced by tleap, defaults to molecule_name
    log_debug_output : bool, optional, default=False
        If true, will send output of tleap to logger.
    verbose : bool
        Be verbose
    Returns
    -------
    prmtop_filename : str
        Amber prmtop file produced by tleap
    inpcrd_filename : str
        Amber inpcrd file produced by tleap
    """
    if prmtop_filename is None:
        prmtop_filename = "%s.prmtop" % molecule_name
    if inpcrd_filename is None:
        inpcrd_filename = "%s.inpcrd" % molecule_name

    #Get absolute paths for input/output
    mol2_filename = os.path.abspath( mol2_filename )
    frcmod_filename = os.path.abspath( frcmod_filename )
    prmtop_filename = os.path.abspath( prmtop_filename )
    inpcrd_filename = os.path.abspath( inpcrd_filename )

    #Work in a temporary directory, on hard coded filenames, to avoid any issues AMBER may have with spaces and other special characters in filenames
    with mdtraj.utils.enter_temp_directory():
        shutil.copy( mol2_filename, 'file.mol2' )
        shutil.copy( frcmod_filename, 'file.frcmod' )

        tleap_input = """
    source oldff/leaprc.ff99
    LIG = loadmol2 file.mol2
    loadamberparams file.frcmod
    check LIG
    saveamberparm LIG out.prmtop out.inpcrd
    quit
"""

        file_handle = open('tleap_commands', 'w')
        file_handle.writelines(tleap_input)
        file_handle.close()

        cmd = "tleap -f %s " % file_handle.name
        if log_debug_output: logger.debug(cmd)
        if verbose: print(cmd)

        output = getoutput(cmd)
        if log_debug_output: logger.debug(output)
        if verbose: print(output)

        check_for_errors( output, other_errors = ['Improper number of arguments'], ignore_errors=['Perhaps a format error'])

        #Copy back target files
        shutil.copy( 'out.prmtop', prmtop_filename )
        shutil.copy( 'out.inpcrd', inpcrd_filename )

    return prmtop_filename, inpcrd_filename


def check_for_errors( outputtext, other_errors = None, ignore_errors = None ):
    """Check AMBER package output for the string 'ERROR' (upper or lowercase) and (optionally) specified other strings and raise an exception if it is found (to avoid silent failures which might be noted to log but otherwise ignored).
    Parameters
    ----------
    outputtext : str
        String listing output text from an (AMBER) command which should be checked for errors.
    other_errors : list(str), default None
        If specified, provide strings for other errors which will be chcked for, such as "improper number of arguments", etc.
    ignore_errors: list(str), default None
        If specified, AMBER output lines containing errors but also containing any of the specified strings will be ignored (because, for example, AMBER issues an "ERROR" for non-integer charges in some cases when only a warning is needed).
    Notes
    -----
    If error(s) are found, raise a RuntimeError and attept to print the appropriate errors from the processed text."""
    lines = outputtext.split('\n')
    error_lines = []
    for line in lines:
        if 'ERROR' in line.upper():
            error_lines.append( line )
        if not other_errors is None:
            for err in other_errors:
                if err.upper() in line.upper():
                    error_lines.append( line )

    if not ignore_errors is None and len(error_lines)>0:
        new_error_lines = []
        for ign in ignore_errors:
            ignore = False
            for err in error_lines:
                if ign in err:
                    ignore = True
            if not ignore:
                new_error_lines.append( err )
        error_lines = new_error_lines

    if len(error_lines) > 0:
        print("Unexpected errors encountered running AMBER tool. Offending output:")
        for line in error_lines: print(line)
        raise(RuntimeError("Error encountered running AMBER tool. Exiting."))

    return
