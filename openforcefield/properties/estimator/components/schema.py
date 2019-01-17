# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Calculation Template API.

Authors
-------
* Simon Boothroyd <simon.boothroyd@choderalab.org>

"""
# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

from pydantic import BaseModel

from typing import Dict, Optional

from openforcefield.properties.estimator.components.protocols import ProtocolSchema, \
    available_protocols, ProtocolPath
from openforcefield.properties.estimator.components.groups import ProtocolGroupSchema


# =============================================================================================
# Property Calculation Schema
# =============================================================================================

class CalculationSchema(BaseModel):
    """Defines the set of protocols required to calculate a certain property.
    """
    property_type: str = None
    id: str = None

    protocols: Dict[str, ProtocolSchema] = {}
    groups: Dict[str, ProtocolGroupSchema] = {}

    final_value_source: Optional[ProtocolPath] = None
    final_uncertainty_source: Optional[ProtocolPath] = None

    def validate_interfaces(self):
        """Validates the flow of the data between protocols, ensuring
        that inputs and outputs correctly match up.
        """

        for protocol_name in self.protocols:

            protocol = self.protocols[protocol_name]

            for input_dependence in protocol.input_dependencies:

                if input_dependence.source.is_global:
                    # We handle global inputs separately
                    continue

                # Make sure the other protocol whose output we are interested
                # in actually exists.
                if input_reference.output_protocol_id not in self.protocols:

                    raise Exception('The {} protocol of the {} schema tries to take input from a non-existent '
                                    'protocol: {}'.format(protocol.id, self.id, input_reference.protocol_id))

                other_protocol = self.protocols[input_reference.output_protocol_id]

                other_protocol_object = available_protocols[other_protocol.type]()

                # Make sure the other protocol definitely has the requested output.
                if input_reference.output_property_name not in other_protocol_object.provided_outputs:

                    raise Exception('The {} protocol does not provide an {} output.'.format(
                        other_protocol.id, input_reference.output_property_name))
