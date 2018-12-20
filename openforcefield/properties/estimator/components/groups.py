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
import os
import logging

from os import path
from enum import Enum, unique

from typing import List

from openforcefield.utils import graph

from .protocols import BaseProtocol, ProtocolInputReference, ProtocolSchema


# =============================================================================================
# Registration Decorators
# =============================================================================================

available_groups = {}


def register_calculation_group():
    """A decorator which registers a class as being a
    protocol group which may be used in calculation schemas.
    """

    def decorator(cls):
        
        if cls.__name__ in available_groups:
            raise ValueError('The {} group is already registered.'.format(cls.__name__))

        available_groups[cls.__name__] = cls
        
        return cls

    return decorator


# =============================================================================================
# Groups
# =============================================================================================

class ProtocolGroupSchema(ProtocolSchema):
    """A json serializable representation of a protocol
    definition.
    """
    grouped_protocol_ids: List[str] = []


@register_calculation_group()
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

    def __init__(self, protocol_ids):
        """Constructs a new ProtocolGroup

        Parameters
        ----------
        protocol_ids : list(str)
            The protocols which will be included in this group.
        """
        super().__init__()

        self._dependants_graph = {}

        self._root_protocols = []
        self._execution_order = []

        self._protocol_ids = protocol_ids
        self._protocols = {}

        #
        # # Groups can take additional global
        # # inputs which the grouped protocols themselves
        # # may not require.
        # self.global_inputs = {}
        #
        # for protocol_id in protocols:
        #     self._protocols[protocol_id] = protocols[protocol_id]
        #     self._dependants_graph[protocol_id] = []
        #
        # # Pull each of an individual protocols inputs up so that they
        # # become a required input of the group.
        # for protocol_id in self._protocols:
        #
        #     protocol = self._protocols[protocol_id]
        #
        #     for input_reference in protocol.input_references:
        #
        #         if input_reference in self.input_references:
        #             continue
        #
        #         if input_reference.output_protocol_id not in self._protocols:
        #             self.input_references.append(copy.deepcopy(input_reference))
        #
        #         if input_reference.output_protocol_id not in self._protocols or \
        #            protocol_id in self._dependants_graph[input_reference.output_protocol_id]:
        #             continue
        #
        #         self._dependants_graph[input_reference.output_protocol_id].append(protocol_id)

        # Do the usual to clean up the graph structure and figure out which order
        # the protocols should execute in.
        graph.apply_transitive_reduction(self._dependants_graph)

        self._root_protocols = graph.find_root_nodes(self._dependants_graph)
        self._execution_order = graph.topological_sort(self._dependants_graph)

    @property
    def schema(self):
        """ProtocolSchema: A serializable schema for this object."""

        base_schema = super(ProtocolGroup, self).schema
        schema = ProtocolGroupSchema.parse_obj(base_schema.dict())

        for protocol_id in self._protocol_ids:
            schema.grouped_protocol_ids.append(protocol_id)

        return schema

    @schema.setter
    def schema(self, schema_value):
        """Sets this protocols properties (i.e id and parameters)
        from a ProtocolSchema"""
        raise NotImplementedError()

    @property
    def protocols(self):
        """dict(str, Protocol): A dictionary of the protocols to execute."""
        return self._protocols

    @property
    def root_protocols(self):
        """list(str): The ids of the root protocols of the dependants graph."""
        return self._root_protocols

    @property
    def dependants_graph(self):
        """dict(str, str): Each key in the dictionary represents a grouped protocol, and each
        string in the value list represents a protocol which depends on the protocol defined by the key."""
        return self._dependants_graph

    def set_uuid(self, value):
        """Store the uuid of the calculation this protocol belongs to

        Parameters
        ----------
        value : str
            The uuid of the parent calculation.
        """
        super(ProtocolGroup, self).set_uuid(value)

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

    def set_global_properties(self, global_properties):
        """Set the value of any inputs which takes values
        from the 'global' (i.e property to calculate) scope

        Parameters
        ----------
        global_properties: dict of str to object
            The list of global properties to draw from.
        """
        for global_property in self.global_inputs:

            if global_property not in global_properties:
                raise Exception('Invalid global property: {}'.format(global_property))

            self.global_inputs[global_property] = global_properties[global_property]

        for protocol_id in self._protocols:
            self._protocols[protocol_id].set_global_properties(global_properties)

    def execute(self, directory):
        """Executes the protocols within this groups

        Parameters
        ----------
        directory : str
            The root directory in which to run the protocols

        Returns
        -------
        bool
            True if all the protocols execute correctly.
        """

        for protocol_id_to_execute in self._execution_order:

            protocol_to_execute = self._protocols[protocol_id_to_execute]
            working_directory = path.join(directory, protocol_to_execute.id)

            if not path.isdir(working_directory):
                os.makedirs(working_directory)

            for input_reference in protocol_to_execute.input_references:

                if input_reference.output_protocol_id in self._protocols:

                    input_protocol = self._protocols[input_reference.output_protocol_id]

                    output_value = input_protocol.get_output_value(input_reference)
                    protocol_to_execute.set_input_value(input_reference, output_value)

            if not protocol_to_execute.execute(working_directory):
                return False

        return True

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

        # Ensure that the two groups at the very least share a root that can be merged.
        # TODO: Is this perhaps too strict....
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

    def _try_merge_protocol(self, other_protocol_id, other_group, parent_id=None):
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
            self._try_merge_protocol(dependant_name, other_group, existing_protocol.id)

    def merge(self, other):
        """Merges another ProtocolGroup with this one. The id
        of this protocol will remain unchanged.

        It is assumed that can_merge has already returned that
        these protocol groups are compatible to be merged together.

        Parameters
        ----------
        other: ProtocolGroup
            The protocol to merge into this one.
        """

        for root_id in other.root_protocols:
            self._try_merge_protocol(root_id, other)

        self._execution_order = graph.topological_sort(self._dependants_graph)

    def set_input_value(self, input_reference, value):
        """Set the value of one of the protocols inputs.

        Parameters
        ----------
        input_reference: ProtocolInputReference
            The input to set.
        value
            The value to set the input to.
        """
        for protocol_id in self._protocols:

            protocol = self._protocols[protocol_id]

            if input_reference not in protocol.input_references:
                continue

            protocol.set_input_value(input_reference, value)

        if input_reference.output_protocol_id == 'global':
            self.global_inputs[input_reference.input_property_name] = value

    def get_input_value(self, input_reference):
        """Gets the value that was set on one of this protocols inputs.

        Parameters
        ----------
        input_reference: ProtocolInputReference
            The input to get.
        Returns
        ----------
        object:
            The value of the input
        """

        for protocol_id in self._protocols:

            protocol = self._protocols[protocol_id]

            if input_reference not in protocol.input_references:
                continue

            return protocol.get_input_value(input_reference)

        raise Exception('The input that was set on a group ({}) was not found. '
                        'This should never happen. {}'.format(self.id, input_reference))

    def get_output_value(self, input_reference):
        """Returns the value of one of this protocols outputs.

        Parameters
        ----------
        input_reference: ProtocolInputReference
            An input reference which points to the output to return.

        Returns
        ----------
        object:
            The value of the output
        """

        if input_reference.grouped_protocol_id is None:

            raise ValueError('Protocols which take input from a group must '
                             'have a grouped_protocol_id set on their input_references.')

        return self._protocols[input_reference.grouped_protocol_id].get_output_value(input_reference)


@register_calculation_group()
class ConditionalGroup(ProtocolGroup):
    """A collection of protocols which are to execute until
    a given condition is met.
    """
    class Condition:
        """Represents a condition placed on a ConditionalGroup"""

        def __init__(self):
            """Constructs a new ConditionalGroup"""
            self.condition_type = None

            self.left_hand_reference = None
            self.right_hand_reference = None

    @unique
    class ConditionType(Enum):
        """The acceptable condtions to place on the group"""
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

    def __init__(self, protocols):
        """Constructs a new ConditionalGroup

        Parameters
        ----------
        protocols : dict(str, Protocol)
            The protocols to include in this group.
        """
        super().__init__(protocols)
        self.conditions = []

        self._max_iterations = 10
        self._stackable = True

        self._condition_type = ConditionalGroup.ConditionType.LessThan

        self._left_hand_value = None
        self._right_hand_value = None

        self._stacked_left_hand_values = []
        self._stacked_right_hand_values = []

    @BaseProtocol.Parameter
    def max_iterations(self):
        pass

    @BaseProtocol.Parameter
    def stackable(self):
        pass

    @BaseProtocol.Parameter
    def condition_type(self):
        pass

    @BaseProtocol.InputPipe
    def left_hand_value(self):
        pass

    @BaseProtocol.InputPipe
    def right_hand_value(self):
        pass

    def _evaluate_condition(self, left_hand_value, right_hand_value):

        if left_hand_value is None or right_hand_value is None:
            return False

        if self.condition_type is self.ConditionType.LessThan:
            return left_hand_value <= right_hand_value
        elif self.condition_type is self.ConditionType.GreaterThan:
            return left_hand_value >= right_hand_value

        raise NotImplementedError()

    def execute(self, directory):
        """Executes the protocols within this groups

        Parameters
        ----------
        directory : str
            The root directory in which to run the protocols

        Returns
        -------
        bool
            True if all the protocols execute correctly.
        """

        logging.info('Starting conditional while loop: {}'.format(self.id))

        should_continue = True

        current_iteration = 0

        while should_continue:

            conditions_met = self._evaluate_condition(self._left_hand_value,
                                                      self._right_hand_value)

            for left_value, right_value in zip(self._stacked_left_hand_values,
                                               self._stacked_right_hand_values):

                conditions_met = conditions_met & self._evaluate_condition(left_value,
                                                                           right_value)

            if conditions_met:

                should_continue = False
                break

            super(ConditionalGroup, self).execute(directory)

            current_iteration += 1

            if self._current_iteration >= self._max_iterations:

                logging.info('Conditional while loop failed to converge: {}'.format(self.id))
                return False

            logging.info('Conditional criteria not yet met')

        logging.info('Conditional while loop finished: {}'.format(self.id))
        return True

    def can_merge(self, other):
        """Determines whether this protocol group can be merged with another.

        Parameters
        ----------
        other : ConditionalGroup
            The protocol group to compare against.

        Returns
        ----------
        bool
            True if the two protocols are safe to merge.
        """
        if not self._stackable and not super(ConditionalGroup, self).can_merge(other):
            return False

        # TODO: Implement stackable group logic.
        return True

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
        #TODO: Implement this for stackables.
        super(ConditionalGroup, self).merge(other)
        self.conditions.extend(other.conditions)

    def set_uuid(self, value):
        """Store the uuid of the calculation this protocol belongs to

        Parameters
        ----------
        value : str
            The uuid of the parent calculation.
        """
        super(ConditionalGroup, self).set_uuid(value)

        for condition in self.conditions:

            condition.left_hand_reference.set_uuid(value)
            condition.right_hand_reference.set_uuid(value)

    def replace_protocol(self, old_id, new_id):
        """Finds each input which came from a given protocol
         and changes it to instead take input from a different one.

        .. todo:
            Implement for conditional groups.

        Parameters
        ----------
        old_id : str
            The id of the old input protocol.
        new_id : str
            The id of the new input protocol.
        """
        super(ConditionalGroup, self).replace_protocol(old_id, new_id)

        for condition in self.conditions:

            condition.left_hand_reference.replace_protocol(old_id, new_id)
            condition.right_hand_reference.replace_protocol(old_id, new_id)