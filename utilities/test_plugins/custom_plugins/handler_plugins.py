from openff.toolkit.typing.engines.smirnoff import ParameterHandler, ParameterIOHandler


class CustomHandler(ParameterHandler):
    _TAGNAME = "CustomHandler"


class WrongSubclass(list):
    _TAGNAME = "CustomHandler"


class CustomIOHandler(ParameterIOHandler):
    _FORMAT = "JSON"
