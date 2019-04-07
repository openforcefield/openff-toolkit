try:
    import openeye
except Exception as e:
    print(e)
    print('Warning: Cannot import openeye toolkit; not all functionality will be available.')

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

#from . import typing
#from . import topology
#from . import utils
