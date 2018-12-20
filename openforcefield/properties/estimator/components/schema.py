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
#     @staticmethod
#     def _parse_protocols_node(existing_schemas, calculation_schema, xml_node):
#         """
#
#         Parameters
#         ----------
#         existing_schemas : list(DirectPropertyCalculationTemplate)
#             A list of already loaded calculation schemas.
#         calculation_schema: CalculationSchema
#             The calculation to add the parsed protocols to.
#         xml_node : xml.etree.Element
#             The element containing the xml to read from.
#
#         Returns
#         -------
#
#         """
#
#         protocol_nodes = xml_node.find('protocols')
#
#         if protocol_nodes is None:
#             raise XmlNodeMissingException('protocols')
#
#         # Load in all of the protocols.
#         for protocol_node in protocol_nodes:
#
#             # Handle the case where we import the protocols
#             # from an already defined calculation.
#             if protocol_node.tag == 'import':
#
#                 if protocol_node.text not in existing_schemas:
#                     raise RuntimeError('The {} schema tries to inherit from a non-existent '
#                                        'schema: {}'.format(calculation_schema.id, protocol_node.text))
#
#                 calculation_schema.protocols.update(copy.deepcopy(existing_schemas[protocol_node.text].protocols))
#
#                 continue
#
#             # Try to automatically find which class corresponds to the xml tag.
#             protocol_class = getattr(protocols, protocol_node.tag)
#
#             if protocol_class is None:
#                 raise RuntimeError('The {} calculation contains an undefined protocol of '
#                                    'type: {}'.format(calculation_schema.id, protocol_node.tag))
#
#             # Make sure the found class is actually a protocol.
#             if isinstance(protocol_class, protocols.BaseProtocol):
#                 raise RuntimeError('{} is not a child of Protocol'.format(protocol_node.tag))
#
#             protocol = protocol_class.from_xml(protocol_node)
#
#             if protocol.id in calculation_schema.protocols:
#                 raise Exception('The {} calculation defines two protocols with the same '
#                                 'id: {}'.format(calculation_schema.id, protocol.id))
#
#             calculation_schema.protocols[protocol.id] = protocol
#
#     @staticmethod
#     def _parse_groups_node(calculation_schema, xml_node):
#         """Used to parse the xml <groups> node.
#
#         Parameters
#         ----------
#         calculation_schema: DirectCalculationTemplate
#             The schema to add the parsed groups to.
#         xml_node
#             The xml nodes that contains the groups node to parse.
#         """
#         protocol_group_nodes = xml_node.find('groups')
#
#         # Load in all of the protocol groups.
#         for protocol_group_node in [] if protocol_group_nodes is None else protocol_group_nodes:
#
#             # Try to automatically find which class corresponds to the xml tag.
#             protocol_group_class = getattr(groups, protocol_group_node.tag)
#
#             if protocol_group_class is None:
#                 raise RuntimeError('The {} calculation contains an undefined protocol '
#                                    'group of type: {}'.format(calculation_schema.id, protocol_group_node.tag))
#
#             # Make sure the found class is actually a protocol_group.
#             if isinstance(protocol_group_class, groups.ProtocolGroup):
#                 raise RuntimeError('{} is not a child of ProtocolGroup'.format(protocol_group_node.tag))
#
#             protocol_group = protocol_group_class.from_xml(protocol_group_node, calculation_schema.protocols)
#
#             if protocol_group.id in calculation_schema.groups:
#                 raise Exception('The {} calculation defines two protocol groups with the '
#                                 'same id: {}'.format(calculation_schema.id, protocol_group.id))
#
#             if protocol_group.id in calculation_schema.protocols:
#                 raise Exception('The {} calculation defines a protocol group {} with the same id as an existing '
#                                 'protocol.'.format(calculation_schema.id, protocol_group.id))
#
#             # Make sure protocols haven't been added to more than one group.
#             for existing_group_id in calculation_schema.groups:
#
#                 existing_group = calculation_schema.groups[existing_group_id]
#
#                 for existing_group_protocol_id in existing_group:
#
#                     if existing_group_protocol_id not in protocol_group:
#                         continue
#
#                     raise Exception('The {} protocol has been added to multiple groups in the '
#                                     '{} calculation'.format(calculation_schema.id, existing_group_protocol_id))
#
#             calculation_schema.groups[protocol_group.id] = protocol_group
#
#     @staticmethod
#     def _parse_output_reference(output_type, schema, xml_node):
#
#         # Find the output value of this calculation.
#         output_value_node = xml_node.find('output-{}'.format(output_type))
#
#         if output_value_node is None:
#             raise XmlNodeMissingException('output-{}'.format(output_type))
#
#         output_value_split = output_value_node.text.split(':')
#
#         if len(output_value_split) != 2:
#             raise Exception('The format of the output-{} node '
#                             'should be protocol_id:property_name'.format(output_type))
#
#         protocol_id = output_value_split[0]
#         property_name = output_value_split[1]
#
#         if protocol_id not in schema.protocols:
#             raise Exception('The {} calculation does not contain '
#                             'a protocol with an id {}'.format(schema.id, protocol_id))
#
#         if property_name not in schema.protocols[protocol_id].provided_outputs:
#             raise Exception('The {} protocol does not output a '
#                             'property named {}'.format(protocol_id, property_name))
#
#         return protocols.ProtocolInputReference('', protocol_id, property_name)
#
#     @classmethod
#     def from_xml(cls, xml_node, existing_schemas):
#         """ Imports a set of protocols from an xml definition.
#
#         Parameters
#         ----------
#         xml_node : xml.etree.Element
#             The element containing the xml to read from.
#         existing_schemas : list(DirectPropertyCalculationTemplate)
#             A list of already loaded calculation schemas
#
#         Returns
#         ----------
#         DirectCalculationTemplate
#             The calculated schema created from the xml node.
#         """
#
#         # Determine which property this calculation will yeild.
#         property_node = xml_node.find('property')
#
#         if property_node is None:
#             raise XmlNodeMissingException('property')
#
#         if property_node.text not in all_properties:
#             raise ValueError('{} is not a valid property type'.format(property_node.text))
#
#         return_value = cls(property_node.text)
#
#         # Find the unique id of this calculation.
#         id_node = xml_node.find('id')
#
#         if id_node is None:
#             raise XmlNodeMissingException('id')
#
#         return_value.id = id_node.text
#
#         # Make sure the id is really unique.
#         if return_value.id in existing_schemas:
#             raise Exception('A calculation with id {} has been '
#                             'defined more than once.'.format(return_value.id))
#
#         return_value.id = id_node.text
#
#         # Load in the protocol definitions.
#         cls._parse_protocols_node(existing_schemas, return_value, xml_node)
#         # Load in the group definitions.
#         cls._parse_groups_node(return_value, xml_node)
#
#         return_value.final_value_reference = cls._parse_output_reference('value',
#                                                                          return_value,
#                                                                          xml_node)
#
#         return_value.final_uncertainty_reference = cls._parse_output_reference('uncertainty',
#                                                                                return_value,
#                                                                                xml_node)
#
#         if not graph.is_acyclic(return_value.dependants_graph):
#             raise Exception('The {} contains a cycle. Only acyclic workflows '
#                             'are supported'.format(return_value.id))
#
#         return return_value
#
#     def _validate_inputs(self):
#         """Validates the flow of the data between protocols, ensuring
#         that inputs and outputs correctly match up.
#         """
#
#         for protocol_name in self.protocols:
#
#             protocol = self.protocols[protocol_name]
#
#             for input_reference in protocol.input_references:
#
#                 if input_reference.output_protocol_id == 'global':
#                     # We handle global inputs separately
#                     continue
#
#                 # Make sure the other protocol whose output we are interested
#                 # in actually exists.
#                 if input_reference.output_protocol_id not in self.protocols:
#                     raise Exception('The {} protocol of the {} schema tries to take input from a non-existent '
#                                     'protocol: {}'.format(protocol.id, self.id, input_reference.protocol_id))
#
#                 other_protocol = self.protocols[input_reference.output_protocol_id]
#
#                 # Make sure the other protocol definitely has the requested output.
#                 if input_reference.output_property_name not in other_protocol.provided_outputs:
#                     raise Exception('The {} protocol does not provide an {} output.'.format(
#                         other_protocol.id, input_reference.output_property_name))
#
#     def _build_dependants_graph(self):
#         """Builds a dictionary of key value pairs where each key represents the id of a
#         protocol to be executed in this calculation, and each value a list ids of protocols
#         which must be ran after the protocol identified by the key.
#         """
#
#         for protocol_name in self.protocols:
#             self.dependants_graph[protocol_name] = []
#
#         for dependant_protocol_name in self.protocols:
#
#             dependant_protocol = self.protocols[dependant_protocol_name]
#
#             for input_reference in dependant_protocol.input_references:
#
#                 if input_reference.output_protocol_id == 'global':
#                     # Global inputs are outside the scope of the
#                     # schema dependency graph.
#                     continue
#
#                 if dependant_protocol.id in self.dependants_graph[input_reference.output_protocol_id]:
#                     continue
#
#                 self.dependants_graph[input_reference.output_protocol_id].append(dependant_protocol.id)
#
#         self.starting_protocols = graph.find_root_nodes(self.dependants_graph)
#
#         # Remove all implicit dependencies from the dependants graph.
#         graph.apply_transitive_reduction(self.dependants_graph)
#
#     def _apply_groups(self):
#         """Groups protocols together into a set of user defined groups."""
#
#         if len(self.groups) == 0:
#             # Nothing to do here.
#             return
#
#         for group_id in self.groups:
#
#             group = self.groups[group_id]
#
#             # Remove all grouped protocols from the protocol list.
#             for grouped_protocol_id in group.protocols:
#                 self.protocols.pop(grouped_protocol_id)
#
#             # Point the other protocols to the groups rather than
#             # the removed protocols
#             for grouped_protocol_id in group.protocols:
#
#                 for protocol_id in self.protocols:
#
#                     protocol = self.protocols[protocol_id]
#                     protocol.replace_protocol(grouped_protocol_id, group.id)
#
#                     for input_reference in protocol.input_references:
#
#                         if input_reference.output_protocol_id != group.id:
#                             continue
#
#                         input_reference.grouped_protocol_id = grouped_protocol_id
#
#                 if self.final_value_reference.output_protocol_id == grouped_protocol_id:
#                     self.final_value_reference.output_protocol_id = group.id
#                     self.final_value_reference.grouped_protocol_id = grouped_protocol_id
#
#                 if self.final_uncertainty_reference.output_protocol_id == grouped_protocol_id:
#                     self.final_uncertainty_reference.output_protocol_id = group.id
#                     self.final_uncertainty_reference.grouped_protocol_id = grouped_protocol_id
#
#             # Add the group in their place.
#             self.protocols[group.id] = group
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
