#!/usr/bin/env python

# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Groups API.

Authors
-------
* Simon Boothroyd <simon.boothroyd@choderalab.org>
"""


# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

import copy
import logging
import os
from enum import Enum, unique
from os import path
from typing import List

from openforcefield.properties.estimator.utils import PropertyEstimatorException
from openforcefield.properties.estimator.workflow import protocol_input
from openforcefield.properties.estimator.workflow.decorators import MergeBehaviour
from openforcefield.utils import graph
from openforcefield.utils.serialization import PolymorphicDataType
from .protocols import BaseProtocol, ProtocolSchema, ProtocolPath, available_protocols, register_calculation_protocol


# =============================================================================================
# Groups
# =============================================================================================

class ProtocolGroupSchema(ProtocolSchema):
    """A json serializable representation of a protocol
    definition.
    """
    grouped_protocol_schemas: List[ProtocolSchema] = []


@register_calculation_protocol()
class ProtocolGroup(BaseProtocol):
    """A collection of protocols to be executed in one batch.

    This may be used for example to cluster together multiple protocols
    that will execute in a linear chain so that multiple scheduler
    execution calls are reduced into a single one.

    Additionally, a group may provide enhanced behaviour, for example
    running all protocols within the group self consistently until
    a given condition is met (e.g run a simulation until a given observable
    has converged).
    """

    @property
    def root_protocols(self):
        """List[str]: The ids of the protocols in the group which do not take
                      input from the other grouped protocols."""

        return self._root_protocols

    @property
    def execution_order(self):
        """List[str]: The ids of the protocols in the group, in the order in which
                      they will be internally executed."""

        return self._execution_order

    @property
    def dependants_graph(self):
        """Dict[str, str]: TODO: Docstring."""

        return self._dependants_graph

    @property
    def protocols(self):
        """Dict[str, BaseProtocol]: A dictionary of the protocols in this groups, where the dictionary
                                    key is the protocol id, and the value the protocol itself."""
        return self._protocols

    def __init__(self, protocol_id):
        """Constructs a new ProtocolGroup.
        """
        super().__init__(protocol_id)

        self._dependants_graph = {}

        self._root_protocols = []
        self._execution_order = []

        self._protocols = {}

    def _get_schema(self):

        base_schema = super(ProtocolGroup, self)._get_schema()
        # Convert the base schema to a group one.
        schema = ProtocolGroupSchema.parse_obj(base_schema.dict())

        for protocol_id in self._protocols:
            schema.grouped_protocol_schemas.append(self._protocols[protocol_id].schema)

        return schema

    def _set_schema(self, schema_value):

        super(ProtocolGroup, self)._set_schema(schema_value)

        protocols_to_create = []

        for protocol_schema in schema_value.grouped_protocol_schemas:

            if protocol_schema.id in self._protocols:

                self._protocols[protocol_schema.id].schema = protocol_schema
                continue

            # Recreate the protocol from scratch.
            protocol = available_protocols[protocol_schema.type](protocol_schema.id)
            protocol.schema = protocol_schema

            protocols_to_create.append(protocol)

        if len(protocols_to_create) > 0:
            self.add_protocols(*protocols_to_create)

    def add_protocols(self, *protocols):

        for protocol in protocols:

            if protocol.id in self._protocols:

                raise ValueError('The {} group already contains a protocol '
                                 'with id {}.'.format(self.id, protocol.id))

            self._protocols[protocol.id] = protocol
            self._dependants_graph[protocol.id] = []

        # Pull each of an individual protocols inputs up so that they
        # become a required input of the group.
        for protocol_id in self._protocols:

            protocol = self._protocols[protocol_id]

            for input_path in protocol.required_inputs:

                grouped_path = ProtocolPath.from_string(input_path.full_path)

                if grouped_path.start_protocol != protocol.id:
                    grouped_path.prepend_protocol_id(protocol.id)

                grouped_path.prepend_protocol_id(self.id)

                if grouped_path in self.required_inputs:
                    continue

                input_value = protocol.get_value(input_path)

                if not isinstance(input_value, ProtocolPath):
                    continue

                if input_value.start_protocol not in self._protocols:

                    self.required_inputs.append(grouped_path)
                    continue

                if protocol_id in self._dependants_graph[input_value.start_protocol]:
                    continue

                self._dependants_graph[input_value.start_protocol].append(protocol_id)

        # Figure out the order in which grouped protocols should be executed.
        self._root_protocols = graph.find_root_nodes(self._dependants_graph)
        self._execution_order = graph.topological_sort(self._dependants_graph)

    def set_uuid(self, value):
        """Store the uuid of the calculation this protocol belongs to

        Parameters
        ----------
        value : str
            The uuid of the parent calculation.
        """
        for index in range(len(self._root_protocols)):
            self._root_protocols[index] = graph.append_uuid(self._root_protocols[index], value)

        for index in range(len(self._execution_order)):
            self._execution_order[index] = graph.append_uuid(self._execution_order[index], value)

        new_dependants_graph = {}

        for protocol_id in self._dependants_graph:

            new_protocol_id = graph.append_uuid(protocol_id, value)
            new_dependants_graph[new_protocol_id] = []

            for dependant in self._dependants_graph[protocol_id]:

                new_dependant_id = graph.append_uuid(dependant, value)
                new_dependants_graph[new_protocol_id].append(new_dependant_id)

        self._dependants_graph = new_dependants_graph

        new_protocols = {}

        for protocol_id in self._protocols:

            protocol = self._protocols[protocol_id]
            protocol.set_uuid(value)

            new_protocols[protocol.id] = protocol

        self._protocols = new_protocols
        super(ProtocolGroup, self).set_uuid(value)

    def replace_protocol(self, old_id, new_id):
        """Finds each input which came from a given protocol
         and redirects it to instead take input from a different one.

        Parameters
        ----------
        old_id : str
            The id of the old input protocol.
        new_id : str
            The id of the new input protocol.
        """
        super(ProtocolGroup, self).replace_protocol(old_id, new_id)

        for index in range(len(self._root_protocols)):
            self._root_protocols[index] = self._root_protocols[index].replace(old_id, new_id)

        for index in range(len(self._execution_order)):
            self._execution_order[index] = self._execution_order[index].replace(old_id, new_id)

        new_dependants_graph = {}

        for protocol_id in self._dependants_graph:

            new_protocol_id = protocol_id.replace(old_id, new_id)
            new_dependants_graph[new_protocol_id] = []

            for dependant in self._dependants_graph[protocol_id]:

                new_dependant_id = dependant.replace(old_id, new_id)
                new_dependants_graph[new_protocol_id].append(new_dependant_id)

        self._dependants_graph = new_dependants_graph

        new_protocols = {}

        for protocol_id in self._protocols:

            protocol = self._protocols[protocol_id]
            protocol.replace_protocol(old_id, new_id)

            new_protocols[protocol_id.replace(old_id, new_id)] = protocol

        self._protocols = new_protocols

    def execute(self, directory, available_resources):
        """Executes the protocols within this groups

        Parameters
        ----------
        directory : str
            The root directory in which to run the protocols
        available_resources: BackendResources
            The resources available to execute on.
        Returns
        -------
        bool
            True if all the protocols execute correctly.
        """

        output_dictionary = {}

        for protocol_id_to_execute in self._execution_order:

            protocol_to_execute = self._protocols[protocol_id_to_execute]
            working_directory = path.join(directory, protocol_to_execute.id)

            if not path.isdir(working_directory):
                os.makedirs(working_directory)

            for input_path in protocol_to_execute.required_inputs:

                target_path = protocol_to_execute.get_value(input_path)

                if not isinstance(target_path, ProtocolPath):
                    continue

                if target_path.start_protocol not in self._protocols:
                    continue

                input_value = self._protocols[target_path.start_protocol].get_value(target_path)

                protocol_to_execute.set_value(input_path, input_value)

            return_value = protocol_to_execute.execute(working_directory, available_resources)

            if isinstance(return_value, PropertyEstimatorException):
                return return_value

            for output_path in return_value:

                output_path_prepended = ProtocolPath.from_string(output_path)

                if (output_path_prepended.start_protocol != self.id and
                    output_path_prepended.start_protocol != protocol_id_to_execute):

                    output_path_prepended.prepend_protocol_id(protocol_id_to_execute)

                if output_path_prepended.start_protocol != self.id:
                    output_path_prepended.prepend_protocol_id(self.id)

                    output_path_prepended.prepend_protocol_id(self.id)

                output_dictionary[output_path_prepended.full_path] = return_value[output_path]

        return output_dictionary

    def can_merge(self, other):
        """Determines whether this protocol group can be merged with another.

        Parameters
        ----------
        other : ProtocolGroup
            The protocol group to compare against.

        Returns
        ----------
        bool
            True if the two protocols are safe to merge.
        """

        if not super(ProtocolGroup, self).can_merge(other):
            return False

        # Ensure that the starting points in each group can be
        # merged.
        # TODO: Is this too strict / too lenient / just right?
        for self_root_id in self._root_protocols:

            self_protocol = self._protocols[self_root_id]

            can_merge_with_root = False

            for other_root_id in other.root_protocols:

                other_protocol = other.protocols[other_root_id]

                if not self_protocol.can_merge(other_protocol):
                    continue

                can_merge_with_root = True
                break

            if not can_merge_with_root:
                return False

        return True

    def _try_merge_protocol(self, other_protocol_id, other_group, merged_ids, parent_id=None):
        """Recursively inserts a protocol node into the group.

        Parameters
        ----------
        other_protocol_id : str
            The name of the other protocol to attempt to merge.
        other_group : ProtocolGroup
            The other protocol group which the protocol to merge belongs to.
        parent_id : str, optional
            The id of the new parent of the node to be inserted. If None,
            the protocol will be added as a new parent node.
        merged_ids : Dict[str, str]
            A map between any original protocol ids and their new merged values.
        """

        if other_protocol_id in self._dependants_graph:

            raise RuntimeError('A protocol with id {} has already been merged '
                               'into the group.'.format(other_protocol_id))

        protocols = self._root_protocols if parent_id is None else self._dependants_graph[parent_id]

        protocol_to_merge = other_group.protocols[other_protocol_id]
        existing_protocol = None

        # Start by checking to see if the starting node of the calculation graph is
        # already present in the full graph.
        for protocol_id in protocols:

            protocol = self._protocols[protocol_id]

            if not protocol.can_merge(protocol_to_merge):
                continue

            existing_protocol = protocol
            break

        if existing_protocol is not None:
            # Make a note that the existing node should be used in place
            # of this calculations version.

            other_group.protocols[other_protocol_id] = existing_protocol
            merged_ids[other_protocol_id] = existing_protocol.id

            for protocol_to_update in other_group.protocols:
                other_group.protocols[protocol_to_update].replace_protocol(other_protocol_id, existing_protocol.id)

        else:

            # Add the protocol as a new node in the graph.
            self._protocols[other_protocol_id] = protocol_to_merge
            existing_protocol = self._protocols[other_protocol_id]

            self._dependants_graph[other_protocol_id] = []
            protocols.append(other_protocol_id)

        # Add all of the dependants to the existing node
        for dependant_name in other_group.dependants_graph[other_protocol_id]:
            self._try_merge_protocol(dependant_name, other_group, merged_ids, existing_protocol.id)

    def merge(self, other):
        """Merges another ProtocolGroup with this one. The id
        of this protocol will remain unchanged.

        It is assumed that can_merge has already returned that
        these protocol groups are compatible to be merged together.

        Parameters
        ----------
        other: ProtocolGroup
            The protocol to merge into this one.

        Returns
        -------
        Dict[str, str]
            A map between any original protocol ids and their new merged values.
        """

        merged_ids = super(ProtocolGroup, self).merge(other)

        for root_id in other.root_protocols:
            self._try_merge_protocol(root_id, other, merged_ids)

        self._execution_order = graph.topological_sort(self._dependants_graph)

        return merged_ids

    def get_attribute_type(self, reference_path):
        """Returns the type of one of the protocol input/output attributes.

        Parameters
        ----------
        reference_path: ProtocolPath
            The path pointing to the value whose type to return.

        Returns
        ----------
        type:
            The type of the attribute.
        """

        reference_property, reference_ids = ProtocolPath.to_components(reference_path.full_path)

        if reference_path.start_protocol is None or (reference_path.start_protocol == self.id and
                                                     len(reference_ids) == 1):

            return super(ProtocolGroup, self).get_attribute_type(reference_path)

        # Make a copy of the path so we can alter it safely.
        reference_path_clone = copy.deepcopy(reference_path)

        if reference_path.start_protocol == self.id:
            reference_path_clone.pop_next_in_path()

        target_protocol_id = reference_path_clone.pop_next_in_path()

        if target_protocol_id not in self._protocols:
            raise ValueError('The reference path does not target this protocol'
                             'or any of its children.')

        return self._protocols[target_protocol_id].get_attribute_type(reference_path_clone)

    def get_value(self, reference_path):
        """Returns the value of one of this protocols parameters / inputs.

        Parameters
        ----------
        reference_path: ProtocolPath
            The path pointing to the value to return.

        Returns
        ----------
        object:
            The value of the input
        """

        reference_property, reference_ids = ProtocolPath.to_components(reference_path.full_path)

        if reference_path.start_protocol is None or (reference_path.start_protocol == self.id and
                                                     len(reference_ids) == 1):

            return super(ProtocolGroup, self).get_value(reference_path)

        # Make a copy of the path so we can alter it safely.
        reference_path_clone = copy.deepcopy(reference_path)

        if reference_path.start_protocol == self.id:
            reference_path_clone.pop_next_in_path()

        target_protocol_id = reference_path_clone.pop_next_in_path()

        if target_protocol_id not in self._protocols:

            raise ValueError('The reference path does not target this protocol'
                             'or any of its children.')

        return self._protocols[target_protocol_id].get_value(reference_path_clone)

    def set_value(self, reference_path, value):
        """Sets the value of one of this protocols parameters / inputs.

        Parameters
        ----------
        reference_path: ProtocolPath
            The path pointing to the value to return.
        value: Any
            The value to set.
        """

        reference_property, reference_ids = ProtocolPath.to_components(reference_path.full_path)

        if reference_path.start_protocol is None or (reference_path.start_protocol == self.id and
                                                     len(reference_ids) == 1):

            return super(ProtocolGroup, self).set_value(reference_path, value)

        # Make a copy of the path so we can alter it safely.
        reference_path_clone = copy.deepcopy(reference_path)

        if reference_path.start_protocol == self.id:
            reference_path_clone.pop_next_in_path()

        target_protocol_id = reference_path_clone.pop_next_in_path()

        if target_protocol_id not in self._protocols:
            raise ValueError('The reference path does not target this protocol'
                             'or any of its children.')

        return self._protocols[target_protocol_id].set_value(reference_path_clone, value)


@register_calculation_protocol()
class ConditionalGroup(ProtocolGroup):
    """A collection of protocols which are to execute until
    a given condition is met.
    """
    @unique
    class ConditionType(Enum):
        """The acceptable conditions to place on the group"""
        LessThan = 'lessthan'
        GreaterThan = 'greaterthan'

        @classmethod
        def has_value(cls, value):
            """Checks whether an of the enum items matches a given value.

            Parameters
            ----------
            value: str
                The value to check for.

            Returns
            ---------
            bool
                True if the enum contains the value.
            """
            return any(value == item.value for item in cls)

    class Condition:

        def __init__(self):

            self.type = ConditionalGroup.ConditionType.LessThan

            self.left_hand_value = None
            self.right_hand_value = None

        def __getstate__(self):

            return {
                'type': self.type.value,
                'left_hand_value': PolymorphicDataType.serialize(self.left_hand_value),
                'right_hand_value': PolymorphicDataType.serialize(self.right_hand_value)
            }

        def __setstate__(self, state):

            self.type = ConditionalGroup.ConditionType(state['type'])

            self.left_hand_value = PolymorphicDataType.deserialize(state['left_hand_value']).value
            self.right_hand_value = PolymorphicDataType.deserialize(state['right_hand_value']).value

        def __eq__(self, other):

            return (self.left_hand_value == other.left_hand_value and
                    self.right_hand_value == other.right_hand_value and
                    self.type == other.type)

        def __ne__(self, other):
            return not self.__eq__(other)

    @protocol_input(int, merge_behavior=MergeBehaviour.GreatestValue)
    def max_iterations(self):
        """The maximum number of iterations to run for to try and satisfy the
         groups conditions."""
        pass

    @property
    def conditions(self):
        return self._conditions

    @property
    def dependencies(self):
        """list of ProtocolPath: A list of pointers to the protocols which this
        protocol takes input from.
        """

        return_dependencies = super(ConditionalGroup, self).dependencies

        for condition in self._conditions:

            if (isinstance(condition.left_hand_value, ProtocolPath) and
                condition.left_hand_value.start_protocol is not None and
                condition.left_hand_value.start_protocol != self.id and
                condition.left_hand_value not in return_dependencies):
                
                return_dependencies.append(condition.left_hand_value)

            if (isinstance(condition.right_hand_value, ProtocolPath) and
                condition.right_hand_value.start_protocol is not None and
                condition.right_hand_value.start_protocol != self.id and
                condition.right_hand_value not in return_dependencies):
                
                return_dependencies.append(condition.right_hand_value)

        return return_dependencies

    def __init__(self, protocol_id):
        """Constructs a new ConditionalGroup
        """
        super().__init__(protocol_id)

        self._max_iterations = 10
        self._conditions = []

        self.required_inputs.append(ProtocolPath('conditions'))

    def _evaluate_condition(self, condition_type, left_hand_value, right_hand_value):

        if left_hand_value is None or right_hand_value is None:
            return False

        if condition_type == self.ConditionType.LessThan:
            return left_hand_value < right_hand_value
        elif condition_type == self.ConditionType.GreaterThan:
            return left_hand_value > right_hand_value

        raise NotImplementedError()

    def execute(self, directory, available_resources):
        """Executes the protocols within this groups

        Parameters
        ----------
        directory : str
            The root directory in which to run the protocols
        available_resources: BackendResources
            The resources available to execute on.

        Returns
        -------
        bool
            True if all the protocols execute correctly.
        """

        logging.info('Starting conditional while loop: {}'.format(self.id))

        should_continue = True

        current_iteration = 0

        while should_continue:

            current_iteration += 1
            return_value = super(ConditionalGroup, self).execute(directory, available_resources)

            if isinstance(return_value, PropertyEstimatorException):
                # Exit on exceptions.
                return return_value

            for condition in self._conditions:

                evaluated_left_hand_value = None

                if not isinstance(condition.left_hand_value, ProtocolPath):
                    evaluated_left_hand_value = condition.left_hand_value
                else:
                    evaluated_left_hand_value = self.get_value(condition.left_hand_value)

                evaluated_right_hand_value = None

                if not isinstance(condition.right_hand_value, ProtocolPath):
                    evaluated_right_hand_value = condition.right_hand_value
                else:
                    evaluated_right_hand_value = self.get_value(condition.right_hand_value)

                # Check to see if we have reached our goal.
                if self._evaluate_condition(condition.type, evaluated_left_hand_value, evaluated_right_hand_value):

                    logging.info('Conditional while loop finished after {} iterations: {}'.format(current_iteration,
                                                                                                  self.id))
                    return return_value

            if current_iteration >= self._max_iterations:

                return PropertyEstimatorException(directory=directory,
                                                  message='Conditional while loop failed to '
                                                           'converge: {}'.format(self.id))

            logging.info('Conditional criteria not yet met after {} iterations'.format(current_iteration))

    def can_merge(self, other):
        return super(ConditionalGroup, self).can_merge(other)

    def merge(self, other):
        """Merges another ProtocolGroup with this one. The id
        of this protocol will remain unchanged.

        It is assumed that can_merge has already returned that
        these protocol groups are compatible to be merged together.

        Parameters
        ----------
        other: ConditionalGroup
            The protocol to merge into this one.
        """
        merged_ids = super(ConditionalGroup, self).merge(other)

        for condition in other.conditions:
            self.add_condition(condition)

        return merged_ids

    def add_condition(self, condition):

        if condition in self._conditions:
            return

        self._conditions.append(condition)

    def set_uuid(self, value):
        """Store the uuid of the calculation this protocol belongs to

        Parameters
        ----------
        value : str
            The uuid of the parent calculation.
        """
        super(ConditionalGroup, self).set_uuid(value)

        for condition in self._conditions:

            if isinstance(condition.left_hand_value, ProtocolPath):
                condition.left_hand_value.append_uuid(value)

            if isinstance(condition.right_hand_value, ProtocolPath):
                condition.right_hand_value.append_uuid(value)

    def replace_protocol(self, old_id, new_id):
        """Finds each input which came from a given protocol
         and redirects it to instead take input from a different one.

        Parameters
        ----------
        old_id : str
            The id of the old input protocol.
        new_id : str
            The id of the new input protocol.
        """
        super(ConditionalGroup, self).replace_protocol(old_id, new_id)

        for condition in self._conditions:

            if isinstance(condition.left_hand_value, ProtocolPath):
                condition.left_hand_value.replace_protocol(old_id, new_id)

            if isinstance(condition.right_hand_value, ProtocolPath):
                condition.right_hand_value.replace_protocol(old_id, new_id)

    def get_attribute_type(self, reference_path):

        if reference_path.start_protocol is None or (reference_path.start_protocol == self.id and
                                                     reference_path.last_protocol == self.id):

            if reference_path.property_name == 'conditions' or reference_path.property_name.find('condition_') >= 0:
                return None

        return super(ConditionalGroup, self).get_attribute_type(reference_path)

    def get_value(self, reference_path):
        """Returns the value of one of this protocols parameters / inputs.

        Parameters
        ----------
        reference_path: ProtocolPath
            The path pointing to the value to return.

        Returns
        ----------
        object:
            The value of the input
        """

        if reference_path.start_protocol is None or (reference_path.start_protocol == self.id and
                                                     reference_path.last_protocol == self.id):

            if reference_path.property_name == 'conditions':
                return self._conditions

        return super(ConditionalGroup, self).get_value(reference_path)

    def set_value(self, reference_path, value):
        """Sets the value of one of this protocols parameters / inputs.

        Parameters
        ----------
        reference_path: ProtocolPath
            The path pointing to the value to return.
        value: Any
            The value to set.
        """

        if reference_path.start_protocol is None or (reference_path.start_protocol == self.id and
                                                     reference_path.last_protocol == self.id):

            if reference_path.property_name == 'conditions':

                self._conditions = value
                return

        super(ConditionalGroup, self).set_value(reference_path, value)
