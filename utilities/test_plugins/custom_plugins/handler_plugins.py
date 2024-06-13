from openff.toolkit.typing.engines.smirnoff import ParameterHandler, ParameterIOHandler
from openff.toolkit.typing.engines.smirnoff.parameters import ParameterAttribute, ParameterType, _NonbondedHandler
from openff.toolkit import unit

class CustomHandler(ParameterHandler):
    _TAGNAME = "CustomHandler"


class WrongSubclass(list):
    _TAGNAME = "CustomHandler"


class CustomIOHandler(ParameterIOHandler):
    _FORMAT = "JSON"

class BuckinghamHandler(_NonbondedHandler):
    """A custom parameter handler for buckingham interactions."""

    class BuckinghamType(ParameterType):
        """A custom parameter type for buckingham interactions."""

        _ELEMENT_NAME = "Atom"

        # Define unit as a Unit object
        a = ParameterAttribute(default=None,
                               unit=unit.kilojoule_per_mole,
            )

        # Define using a string
        b = ParameterAttribute(default=None,
                               unit="nanometer**-1",
                               )

        # Define as None, since that's supported
        c = ParameterAttribute(
            default=None,
            unit=None,
        )

    _TAGNAME = "TestBuckingham"
    _INFOTYPE = BuckinghamType
