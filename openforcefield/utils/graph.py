# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
A set of utilities for manipulating and validating graph-like structures.

Authors
-------
* Simon Boothroyd <simon.boothroyd@choderalab.org>
"""
# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

import copy


# =============================================================================================
# Utilities
# =============================================================================================

def apply_transitive_reduction(graph):
    """Attempts to remove any implicit dependencies from the
    dependant graph.

    For example, if C depends on B which depends on A (A->B->C), the A->C
    dependency should not be stored as it is already covered by storing that
    A->B and B->C.

    Notes
    -----
    The graph must be directed and acyclic.

    Parameters
    ----------
    graph: dict(str, list(str))
        The graph to reduce. Each key in the dictionary represents a node in the graph, and each
        string in the value list represents a node which depends on the node defined by the key.
    """
    closed_list = []
    closure = {}

    for node_key in graph:
        closure[node_key] = []

    for node_key in graph:
        _visit_protocol(graph, node_key, closed_list, closure)


def _visit_protocol(graph, current_key, closed_list, closure):
    """Visits each node in a graph in turn and tries to remove any
    of its implicit dependencies.

    Parameters
    ----------
    graph: dict(str, list(str))
        The graph to traverse. Each key in the dictionary represents a node in the graph,
        and each string in the value list represents a child of the node defined by the key.
    current_key: str
        The node that is currently being visited.
    closed_list: list(str)
        A list of already visited nodes.
    closure: dict(str, list(str)
        A dictionary storing the transitive closure of each node. This has the same
        structure as the graph parameter.
    """

    if current_key in closed_list:
        # Don't visit protocols more than once.
        return

    closed_list.append(current_key)

    # Build a list of this protocols indirect dependencies.
    indirect_dependencies = []

    for dependent in graph[current_key]:

        _visit_protocol(graph, dependent, closed_list, closure)
        indirect_dependencies.extend(closure[dependent])

    closure[current_key].extend(indirect_dependencies)

    for i in range(len(graph[current_key]) - 1, -1, -1):

        dependent = graph[current_key][i]
        closure[current_key].append(dependent)

        if dependent in indirect_dependencies:
            graph[current_key].pop(i)


def find_root_nodes(graph):
    """Returns the nodes in a graph that do not have inputs.

    Notes
    -----
    The graph must be directed and acyclic.

    Parameters
    ----------
    graph: dict(str, list(str))
        The graph to explore. Each key in the dictionary represents a node in the graph, and each
        string in the value list represents a node which depends on the node defined by the key.

    Returns
    -------
    list(str)
        The root nodes of the graph.
    """

    dependencies = dependants_to_dependencies(graph)
    root_nodes = []

    for node_key in dependencies:

        if len(dependencies[node_key]) > 0:
            continue

        root_nodes.append(node_key)

    return root_nodes


def topological_sort(graph):
    """Performs a topological sort on a graph using Kahn's algorithm.
    The resulting array contains a list of nodes ordered in such a way
    that 'dependant' nodes always come after their dependencies.

    Notes
    -----
    The graph must be directed and acyclic.

    Parameters
    ----------
    graph: dict(str, list(str))
        The graph to explore. Each key in the dictionary represents a node in the graph, and each
        string in the value list represents a node which depends on the node defined by the key.
        
    Returns
    -------
    list(str)
        An ordered list of node keys. The will be empty if a cycle is found.
    """

    sorted_order = []
    open_list = find_root_nodes(graph)

    if len(open_list) == 0:
        return sorted_order

    # Make a copy of the graph as this is a destructive process.
    graph_copy = copy.deepcopy(graph)
    apply_transitive_reduction(graph_copy)

    while len(open_list) > 0:

        current_node_key = open_list.pop()
        sorted_order.append(current_node_key)

        for i in range(len(graph_copy[current_node_key]) - 1, -1, -1):

            comparison_node_key = graph_copy[current_node_key].pop(i)

            has_incoming_edges = False

            for node_key in graph_copy:

                if comparison_node_key not in graph_copy[node_key]:
                    continue

                has_incoming_edges = True
                break

            if not has_incoming_edges:
                open_list.append(comparison_node_key)

    remaining_dependant_count = 0

    # Test to make sure the graph wasn't cyclic.
    for node_key in graph_copy:
        remaining_dependant_count += len(graph_copy[node_key])

    if remaining_dependant_count > 0:
        return[]

    return sorted_order


def dependants_to_dependencies(graph):
    """Inverts a dependant's graph to yield a dependency graph.

    Notes
    -----
    The graph must be directed and acyclic.

    Parameters
    ----------
    graph: dict(str, list(str))
        The graph to invert. Each key in the dictionary represents a node in the graph, and each
        string in the value list represents a node which depends on the node defined by the key.

    Returns
    -------
    dict(str, list(str))
        The inverted graph. Each key in the dictionary represents a node in the graph, and each
        string in the value list represents a node which the node defined by the key depends on.
    """

    dependencies = {}

    for node in graph:

        if node not in dependencies:
            dependencies[node] = []

        for dependant in graph[node]:

            if dependant not in dependencies:
                dependencies[dependant] = []

            if node not in dependencies[dependant]:
                dependencies[dependant].append(node)

    return dependencies


def is_acyclic(graph):
    """Determines if a directed graph is acyclic.

    Parameters
    ----------
    graph: dict(str, list(str))
        The graph to explore. Each key in the dictionary represents a node in the graph, and each
        string in the value list represents a node which depends on the node defined by the key.

    Returns
    -------
    bool
        True if the graph is acyclic
    """

    return len(graph) == 0 or len(topological_sort(graph)) != 0


def append_uuid(base_id, uuid):
    """Appends a uuid to a base id

    Parameters
    ----------
    base_id : str
        The id to append the uuid to.
    uuid : str
        The uuid to append.

    Returns
    -------
    str
        The uuid appended id of the form uuid|base_id
    """
    if base_id.find('|') >= 0:

        base_id_split = base_id.split('|')

        if len(base_id_split) != 2:
            raise ValueError('A uuid appended protocol id should be of the form uuid|base_id')

        base_id = base_id_split[1]

    return '{}|{}'.format(uuid, base_id)


def retrieve_uuid(uuid_appended_id):
    """Retrieves the uuid of a uuid appended id

    Parameters
    ----------
    uuid_appended_id : str
        A uuid appended id of the form uuid|base_id.

    Returns
    -------
    str
        The uuid portion of the appended id.
    """

    if uuid_appended_id.find('|') < 0:
        return None

    id_split = uuid_appended_id.split('|')

    if len(id_split) != 2:
        raise ValueError('A uuid appended protocol id should be of the form uuid|base_id')

    return_uuid = id_split[0]

    return return_uuid
