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

import json

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

    final_value_source: Optional[ProtocolPath] = None
    final_uncertainty_source: Optional[ProtocolPath] = None

    class Config:
        arbitrary_types_allowed = True

        json_encoders = {
            ProtocolPath: lambda v: v.full_path
        }

    def validate_interfaces(self):
        """Validates the flow of the data between protocols, ensuring
        that inputs and outputs correctly match up.
        """

        if self.final_value_source.start_protocol not in self.protocols:
            raise ValueError('The value source {} does not exist.'.format(self.final_value_source))

        if self.final_uncertainty_source.start_protocol not in self.protocols:
            raise ValueError('The uncertainty source {} does not exist.'.format(self.final_uncertainty_source))

        for protocol_name in self.protocols:

            protocol_schema = self.protocols[protocol_name]

            protocol_object = available_protocols[protocol_schema.type](protocol_schema.id)
            protocol_object.schema = protocol_schema

            if protocol_name == self.final_value_source.start_protocol:
                protocol_object.get_value(self.final_value_source)
            if protocol_name == self.final_uncertainty_source.start_protocol:
                protocol_object.get_value(self.final_uncertainty_source)

            for input_path in protocol_object.required_inputs:

                input_value = protocol_object.get_value(input_path)

                if input_value is None:

                    raise Exception('The {} required input of protocol {} in the {} schema was '
                                    'not set.'.format(input_path, protocol_name, self.id))

                if input_value.is_global:
                    # We handle global input validation separately
                    continue

                # Make sure the other protocol whose output we are interested
                # in actually exists.
                if input_value.start_protocol not in self.protocols:

                    raise Exception('The {} protocol of the {} schema tries to take input from a non-existent '
                                    'protocol: {}'.format(protocol_object.id, self.id, input_value.start_protocol))

                other_protocol_schema = self.protocols[input_value.start_protocol]

                other_protocol_object = available_protocols[other_protocol_schema.type](other_protocol_schema.id)
                other_protocol_object.schema = other_protocol_schema

                # Will throw the correct exception if missing.
                other_protocol_object.get_value(input_value)
