"""Contains a set of 'plugin' classes to enable testing of the plugin system."""
from openforcefield.typing.engines.smirnoff import (ParameterHandler,
                                                    ParameterIOHandler)


class CustomHandler(ParameterHandler):
    _TAGNAME = "CustomHandler"


class CustomIOHandler(ParameterIOHandler):
    _FORMAT = "JSON"
