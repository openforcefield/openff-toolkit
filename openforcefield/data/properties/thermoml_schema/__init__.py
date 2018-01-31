"""
ThermoML Schema

PyXB generated schemas from the ThermoML Schema spec XML downloaded from NIST

Because different pyxbgen versions are incompatible with one another, this module tries to import all of them
and returns the one that works

"""
from pyxb import PyXBVersionError
done = False
try:
    from .thermoml_schema124 import *
    done = True
except PyXBVersionError:
    pass

if not done:
    try:
        from .thermoml_schema125 import *
        done = True
    except PyXBVersionError:
        pass

if not done:
    try:
        from .thermoml_schema126 import *
        done = True
    except PyXBVersionError as final_error:
        # Final error
        raise
