#!/usr/bin/env python

"""
Utility subroutines

"""
#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================
import os
import tempfile
import shutil
import contextlib

#=============================================================================================
# UTILITY SUBROUTINES
#=============================================================================================

def all_subclasses(cls):
    return cls.__subclasses__() + [ g for s in cls.__subclasses__() for g in all_subclasses(s) ]

@contextlib.contextmanager
def temporary_cd(dir_path):
    """Context to temporary change the working directory."""
    prev_dir = os.getcwd()
    os.chdir(os.path.abspath(dir_path))
    try:
        yield
    finally:
        os.chdir(prev_dir)

@contextlib.contextmanager
def temporary_directory():
    """Context for safe creation of temporary directories."""
    tmp_dir = tempfile.mkdtemp()
    try:
        yield tmp_dir
    finally:
        shutil.rmtree(tmp_dir)

def get_data_filename(relative_path):
    """Get the full path to one of the reference files in testsystems.
    In the source distribution, these files are in ``openforcefield/data/``,
    but on installation, they're moved to somewhere in the user's python
    site-packages directory.
    Parameters
    ----------
    name : str
        Name of the file to load (with respect to the repex folder).
    """

    from pkg_resources import resource_filename
    fn = resource_filename('openforcefield', os.path.join('data', relative_path))

    if not os.path.exists(fn):
        raise ValueError("Sorry! %s does not exist. If you just added it, you'll have to re-install" % fn)

    return fn
