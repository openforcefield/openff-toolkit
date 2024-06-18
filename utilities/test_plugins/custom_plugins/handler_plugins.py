# noqa: INP001
from openff.toolkit.typing.engines.smirnoff import ParameterHandler, ParameterIOHandler
from openff.toolkit.typing.engines.smirnoff.parameters import ParameterAttribute, ParameterType, _NonbondedHandler
from openff.toolkit import unit

class CustomHandler(ParameterHandler):
    _TAGNAME = "CustomHandler"


class WrongSubclass(list):
    _TAGNAME = "CustomHandler"


class CustomIOHandler(ParameterIOHandler):
    _FORMAT = "JSON"

class FOOBuckinghamHandler(_NonbondedHandler):
    """A custom parameter handler for buckingham interactions."""

    class FOOBuckinghamType(ParameterType):
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

        c = ParameterAttribute(
            default=None,
            unit=unit.kilojoule_per_mole * unit.nanometer ** 6,
        )

    _TAGNAME = "FOOBuckingham"
    _INFOTYPE = FOOBuckinghamType
