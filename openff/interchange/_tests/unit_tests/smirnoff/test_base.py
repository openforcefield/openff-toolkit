import pytest
from openff.toolkit.topology import Topology
from openff.toolkit.typing.engines.smirnoff.parameters import (
    AngleHandler,
    ParameterHandler,
)

from openff.interchange._tests import _BaseTest
from openff.interchange.exceptions import InvalidParameterHandlerError
from openff.interchange.smirnoff._base import SMIRNOFFCollection
from openff.interchange.smirnoff._valence import SMIRNOFFAngleCollection


class TestSMIRNOFFCollection(_BaseTest):
    def test_allowed_parameter_handler_types(self):
        class DummyParameterHandler(ParameterHandler):
            pass

        class DummySMIRNOFFCollection(SMIRNOFFCollection):
            type = "Bonds"
            expression = "1+1"

            @classmethod
            def allowed_parameter_handlers(cls):
                return [DummyParameterHandler]

            @classmethod
            def supported_parameters(cls):
                return list()

        dummy_handler = DummySMIRNOFFCollection()
        angle_Handler = AngleHandler(version=0.3)

        assert DummyParameterHandler in dummy_handler.allowed_parameter_handlers()
        assert AngleHandler not in dummy_handler.allowed_parameter_handlers()
        assert (
            DummyParameterHandler
            not in SMIRNOFFAngleCollection.allowed_parameter_handlers()
        )

        dummy_handler = DummyParameterHandler(version=0.3)

        with pytest.raises(InvalidParameterHandlerError):
            SMIRNOFFAngleCollection.create(
                parameter_handler=dummy_handler,
                topology=Topology(),
            )

        with pytest.raises(InvalidParameterHandlerError):
            DummySMIRNOFFCollection.create(
                parameter_handler=angle_Handler,
                topology=Topology(),
            )
