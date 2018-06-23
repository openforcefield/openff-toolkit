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
def temporary_directory():
    """Context for safe creation of temporary directories."""
    tmp_dir = tempfile.mkdtemp()
    try:
        yield tmp_dir
    finally:
        shutil.rmtree(tmp_dir)
