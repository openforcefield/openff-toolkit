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

from openforcefield.properties.estimator.components.protocols import ProtocolInputReference, \
    ProtocolSchema, available_protocols
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

    final_value_reference: Optional[ProtocolInputReference] = None
    final_uncertainty_reference: Optional[ProtocolInputReference] = None

    def validate_interfaces(self):
        """Validates the flow of the data between protocols, ensuring
        that inputs and outputs correctly match up.
        """

        for protocol_name in self.protocols:

            protocol = self.protocols[protocol_name]

            for input_reference in protocol.input_references:

                if input_reference.output_protocol_id == 'global':
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

# class CalculationSchema:
#     """Defines the set of protocols required to calculate a certain property.
#     """
#
#     def __init__(self, property_type, schema_id=None):
#
#         self.property_type = property_type
#         self.id = schema_id or '{}Template'.format(str(self.property_type))
#
#         self.protocols = {}
#         self.groups = {}
#
#         self.final_value_reference = None
#         self.final_uncertainty_reference = None
#
#         # A list of protocols which have zero or only global inputs.
#         # These will be the first protocols to be executed.
#         self.starting_protocols = []
#
#         # A key-value pair dictionary of dependants, where
#         # each key represents a protocol and each value is a list
#         # of the protocols which take said protocol as an input.
#         self.dependants_graph = {}
#
#     def build(self):
#         """Builds the schema dependency graph and applies any
#         groupings.
#
#         Notes
#         -----
#
#         This method should only be called after all protocols and groups have
#         been added.
#         """
#         self._validate_inputs()
#
#         self._apply_groups()
#         self._build_dependants_graph()
#
#     @classmethod
#     def from_json(cls):
#         pass
#
#     def to_json(self):
#
#         #         self.property_type = property_type
#         #         self.id = schema_id or '{}Template'.format(str(self.property_type))
#         #
#         #         self.protocols = {}
#         #         self.groups = {}
#         #
#         #         self.final_value_reference = None
#         #         self.final_uncertainty_reference = None
