r"""
This module defines elements by the dictionaries:

 - :py:data:`R_at_nodes_of_element`
 - :py:data:`exponents_of_shape_functions_of_element`
 - :py:data:`R_at_integration_points_of_element`
 - :py:data:`weights_of_integration_points_of_element`
 - :py:data:`corner_nodes_of_element`
 - :py:data:`adjacent_nodes_of_element`

The integration points and weights are computed using Abaqus.
The other dictionaries are defined manually.

.. seealso::
    The module :py:mod:`conforce.expressions`
    describes how to use the nodes, exponents, etc.


Access dictionaries
-------------------

To acces the nodes, exponents, etc. for a specific element,
the element type can be looked up in the dictionaries.
For example,

>>> typ = CPE4R

is a *four*-noded *2* D element with *one* integration point.

>>> ref_nodes = R_at_nodes_of_element[typ]
>>> ref_nodes.shape
(4, 2)
>>> shape_exponents = exponents_of_shape_functions_of_element[typ]
>>> shape_exponents.shape
(4, 2)
>>> int_points = R_at_integration_points_of_element[typ]
>>> int_points.shape
(1, 2)
>>> int_weights = weights_of_integration_points_of_element[typ]
>>> int_weights.shape
(1,)
>>> corner_nodes = corner_nodes_of_element[typ]
>>> corner_nodes.shape
(4,)
>>> adjacent_nodes = adjacent_nodes_of_element[typ]
>>> adjacent_nodes.shape
(4, 4)


Examples
--------

The corner nodes of an element are saved as boolean array.

>>> R_at_nodes = R_at_nodes_of_element[CPE6]
>>> corner_nodes = corner_nodes_of_element[CPE6]
>>> corner_nodes
array([ True,  True,  True, False, False, False])

This boolean array can be used as slice to get the coordinates of the corner nodes.

>>> R_at_nodes[corner_nodes]
array([[0., 0.],
       [1., 0.],
       [0., 1.]])

The i-th node is also accesible by an index.

>>> i = 0
>>> R_at_nodes[i, :]
array([0., 0.])

The adjacent nodes of this i-th node are saved as a boolean array.

>>> adjacent_nodes = adjacent_nodes_of_element[CPE6]
>>> adjacent_nodes[i, :]
array([False, False, False,  True, False,  True])

Consequently, the boolean array can be used as slice to
get the coordinates of nodes adjacent to the i-th node.

>>> R_at_nodes[adjacent_nodes[i, :], :]
array([[0.5, 0. ],
       [0. , 0.5]])

"""
import os
import numpy as np
import json

from conforce.one_element_runner import simulate_one_element


# Declare dictionaries that define each element type

R_at_nodes_of_element = dict()
"""Coordinates of the element nodes in the reference coordinate system"""

exponents_of_shape_functions_of_element = dict()
"""Exponents of powers used by the shape functions defined in the reference space"""

R_at_integration_points_of_element = dict()
"""Coordinates of integration points in the reference coordinates system"""

weights_of_integration_points_of_element = dict()
"""Weights corresponding to the integration points"""

corner_nodes_of_element = dict()
"""Nodes lying on the intersection of three element edges for 3D elements and two element edges for 2D elements."""

adjacent_nodes_of_element = dict()
"""
Adjacency matrices for element nodes.
Two nodes i, j are adjacent, if they are on the 
same outer element edge and there is no other node between them.
"""


# Element types derived from a 20-node quadratic brick:

C3D20 = "C3D20"
"""3D 20-node quadratic brick"""
R_at_nodes_of_element[C3D20] = np.array([
    [-1., -1., -1.], [1.0, -1., -1.], [1.0, 1.0, -1.], [-1., 1.0, -1.],
    [-1., -1., 1.0], [1.0, -1., 1.0], [1.0, 1.0, 1.0], [-1., 1.0, 1.0],
    [0.0, -1., -1.], [1.0, 0.0, -1.], [0.0, 1.0, -1.], [-1., 0.0, -1.],
    [0.0, -1., 1.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0], [-1., 0.0, 1.0],
    [-1., -1., 0.0], [1.0, -1., 0.0], [1.0, 1.0, 0.0], [-1., 1.0, 0.0]
])
exponents_of_shape_functions_of_element[C3D20] = np.array([
    [0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0],
    [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1],
    [2, 0, 0], [0, 2, 0], [2, 1, 0], [1, 2, 0],
    [2, 0, 1], [0, 2, 1], [2, 1, 1], [1, 2, 1],
    [0, 0, 2], [1, 0, 2], [0, 1, 2], [1, 1, 2]
])
adjacent_nodes_of_element[C3D20] = np.array([
    [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1],
    [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
], dtype=bool)
corner_nodes_of_element[C3D20] = np.sum(adjacent_nodes_of_element[C3D20].astype(int), axis=1) == 3

C3D20R = "C3D20R"
"""3D 20-node quadratic brick, reduced integration"""
R_at_nodes_of_element[C3D20R] = R_at_nodes_of_element[C3D20]
exponents_of_shape_functions_of_element[C3D20R] = exponents_of_shape_functions_of_element[C3D20]
adjacent_nodes_of_element[C3D20R] = adjacent_nodes_of_element[C3D20]
corner_nodes_of_element[C3D20R] = corner_nodes_of_element[C3D20]

CPE8 = "CPE8"
"""2D 8-node bilinear quadrilateral"""
R_at_nodes_of_element[CPE8] = R_at_nodes_of_element[C3D20][[0, 1, 2, 3, 8, 9, 10, 11], :2]
exponents_of_shape_functions_of_element[CPE8] = exponents_of_shape_functions_of_element[C3D20][[0, 1, 2, 3, 8, 9, 10, 11], :2]
adjacent_nodes_of_element[CPE8] = adjacent_nodes_of_element[C3D20][[0, 1, 2, 3, 8, 9, 10, 11], :][:, [0, 1, 2, 3, 8, 9, 10, 11]]
corner_nodes_of_element[CPE8] = corner_nodes_of_element[C3D20][[0, 1, 2, 3, 8, 9, 10, 11]]

CPE8R = "CPE8R"
"""2D 8-node bilinear quadrilateral, reduced integration"""
R_at_nodes_of_element[CPE8R] = R_at_nodes_of_element[CPE8]
exponents_of_shape_functions_of_element[CPE8R] = exponents_of_shape_functions_of_element[CPE8]
adjacent_nodes_of_element[CPE8R] = adjacent_nodes_of_element[CPE8]
corner_nodes_of_element[CPE8R] = corner_nodes_of_element[CPE8]


def _add_connections_bypassing_mid_nodes(graph):
    """
    Add connections between nodes i, j,
    if i has a connection to an arbitrary node k
    and this node k has a connection to node j.

    :param graph: numpy array of shape (n, n) containing an adjacency matrix
    :return: numpy array of shape (n, n) containing an adjacency matrix
    """
    new_graph = np.array([
        np.logical_or.reduce(graph[row, :], axis=0)
        for row in graph.astype(bool)
    ], dtype=int) - np.eye(len(graph), dtype=int)
    return new_graph.astype(graph.dtype)


C3D8 = "C3D8"
"""3D 8-node quadratic brick"""
R_at_nodes_of_element[C3D8] = R_at_nodes_of_element[C3D20][:8, :]
exponents_of_shape_functions_of_element[C3D8] = exponents_of_shape_functions_of_element[C3D20][:8, :]
adjacent_nodes_of_element[C3D8] = _add_connections_bypassing_mid_nodes(
    adjacent_nodes_of_element[C3D20]
)[:8, :8]
corner_nodes_of_element[C3D8] = corner_nodes_of_element[C3D20][:8]

C3D8R = "C3D8R"
"""3D 8-node quadratic brick, reduced integration"""
R_at_nodes_of_element[C3D8R] = R_at_nodes_of_element[C3D8]
exponents_of_shape_functions_of_element[C3D8R] = exponents_of_shape_functions_of_element[C3D8]
adjacent_nodes_of_element[C3D8R] = adjacent_nodes_of_element[C3D8]
corner_nodes_of_element[C3D8R] = corner_nodes_of_element[C3D8]

CPE4 = "CPE4"
"""2D 4-node bilinear quadrilateral"""
R_at_nodes_of_element[CPE4] = R_at_nodes_of_element[C3D8][:4, :2]
exponents_of_shape_functions_of_element[CPE4] = exponents_of_shape_functions_of_element[C3D8][:4, :2]
adjacent_nodes_of_element[CPE4] = adjacent_nodes_of_element[C3D8][:4, :4]
corner_nodes_of_element[CPE4] = corner_nodes_of_element[C3D8][:4]

CPE4R = "CPE4R"
"""2D 4-node bilinear quadrilateral, reduced integration"""
R_at_nodes_of_element[CPE4R] = R_at_nodes_of_element[CPE4]
exponents_of_shape_functions_of_element[CPE4R] = exponents_of_shape_functions_of_element[CPE4]
adjacent_nodes_of_element[CPE4R] = adjacent_nodes_of_element[CPE4]
corner_nodes_of_element[CPE4R] = corner_nodes_of_element[CPE4]


# Element types derived from a 10-node quadratic tetrahedron:

C3D10 = "C3D10"
"""3D 10-node quadratic tetrahedron"""
R_at_nodes_of_element[C3D10] = np.array([
    [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
    [0.5, 0.0, 0.0], [0.5, 0.5, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]
])
exponents_of_shape_functions_of_element[C3D10] = np.array([
    [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1],
    [2, 0, 0], [1, 1, 0], [0, 2, 0], [0, 0, 2], [1, 0, 1], [0, 1, 1],
])
adjacent_nodes_of_element[C3D10] = np.array([
    [0, 0, 0, 0, 1, 0, 1, 1, 0, 0],
    [0, 0, 0, 0, 1, 1, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 1, 1, 0, 0, 1],
    [0, 0, 0, 0, 0, 0, 0, 1, 1, 1],
    [1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 1, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 1, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 0, 1, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 1, 0, 0, 0, 0, 0, 0]
], dtype=bool)
corner_nodes_of_element[C3D10] = np.sum(adjacent_nodes_of_element[C3D10].astype(int), axis=1) == 3

C3D4 = "C3D4"
"""3D 4-node linear tetrahedron"""
R_at_nodes_of_element[C3D4] = R_at_nodes_of_element[C3D10][:4, :]
exponents_of_shape_functions_of_element[C3D4] = exponents_of_shape_functions_of_element[C3D10][:4, :]
adjacent_nodes_of_element[C3D4] = _add_connections_bypassing_mid_nodes(
    adjacent_nodes_of_element[C3D10])[:4, :4]
corner_nodes_of_element[C3D4] = corner_nodes_of_element[C3D10][:4]


# Element types derived from a 15-node quadratic triangular prism

C3D15 = "C3D15"
"""3D 15-node quadratic triangular prism"""
R_at_nodes_of_element[C3D15] = np.array([
    [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0],
    [0.5, 0.0, 0.0], [0.5, 0.5, 0.0], [0.0, 0.5, 0.0],
    [0.5, 0.0, 1.0], [0.5, 0.5, 1.0], [0.0, 0.5, 1.0],
    [0.0, 0.0, 0.5], [1.0, 0.0, 0.5], [0.0, 1.0, 0.5],
])
exponents_of_shape_functions_of_element[C3D15] = np.array([
    [0, 0, 0], [1, 0, 0], [0, 1, 0],
    [0, 0, 1], [1, 0, 1], [0, 1, 1],
    [2, 0, 0], [1, 1, 0], [0, 2, 0],
    [2, 0, 1], [1, 1, 1], [0, 2, 1],
    [0, 0, 2], [1, 0, 2], [0, 1, 2],
])
adjacent_nodes_of_element[C3D15] = np.array([
    [0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1],
    [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
], dtype=bool)
corner_nodes_of_element[C3D15] = np.sum(adjacent_nodes_of_element[C3D15].astype(int), axis=1) == 3

CPE6 = "CPE6"
"""2D 6-node quadratic triangle"""
R_at_nodes_of_element[CPE6] = R_at_nodes_of_element[C3D15][[0, 1, 2, 6, 7, 8], :2]
exponents_of_shape_functions_of_element[CPE6] = exponents_of_shape_functions_of_element[C3D15][[0, 1, 2, 6, 7, 8], :2]
adjacent_nodes_of_element[CPE6] = adjacent_nodes_of_element[C3D15][[0, 1, 2, 6, 7, 8], :][:, [0, 1, 2, 6, 7, 8]]
corner_nodes_of_element[CPE6] = corner_nodes_of_element[C3D15][[0, 1, 2, 6, 7, 8]]

C3D6 = "C3D6"
"""3D 6-node quadratic triangular prism"""
R_at_nodes_of_element[C3D6] = R_at_nodes_of_element[C3D15][:6, :]
exponents_of_shape_functions_of_element[C3D6] = exponents_of_shape_functions_of_element[C3D15][:6, :]
adjacent_nodes_of_element[C3D6] = _add_connections_bypassing_mid_nodes(
    adjacent_nodes_of_element[C3D15])[:6, :6]
corner_nodes_of_element[C3D6] = corner_nodes_of_element[C3D15][:6]

CPE3 = "CPE3"
"""2D 3-node linear triangle"""
R_at_nodes_of_element[CPE3] = R_at_nodes_of_element[C3D6][:3, :2]
exponents_of_shape_functions_of_element[CPE3] = exponents_of_shape_functions_of_element[C3D6][:3, :2]
adjacent_nodes_of_element[CPE3] = adjacent_nodes_of_element[C3D6][:3, :3]
corner_nodes_of_element[CPE3] = corner_nodes_of_element[C3D6][:3]


# compute integration points and weights

def _read_integration_data_file():
    """
    Read the _element_integration_data.json file
    that contains integration points and weights
    of the element types and update
     - R_at_integration_points_of_element
     - weights_of_integration_points_of_element

    If the file does not exist, the method
    :py:func`_write_integration_data_file`
    write is executed.
    """
    integration_data_file_path = os.path.abspath(os.path.join(
        __file__, os.pardir, "_element_integration_data.json"))

    if not os.path.exists(integration_data_file_path):
        _write_integration_data_file(integration_data_file_path)

    with open(integration_data_file_path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
        int_points = data["X_at_integration_points_for_element"]
        int_points = {
            str(element_type): np.array(data, dtype=float)
            for element_type, data in int_points.items()
        }

        int_weights = data["weights_of_integration_points_of_element"]
        int_weights = {
            str(element_type): np.array(data, dtype=float)
            for element_type, data in int_weights.items()
        }

    # update integration data in this module
    R_at_integration_points_of_element.update(int_points)
    weights_of_integration_points_of_element.update(int_weights)


def _write_integration_data_file(integration_data_file_path):
    """
    Simulate a one element model for each element type
    defined in :py:attribute:`R_at_nodes_of_element`
    using Abaqus.
    The integration points and weights are extracted out
    of the Abaqus odb-file and are written to the
    given file path in the json format.

    :param integration_data_file_path: str, Write results into this file
    """
    simulation_results = {
        element_type: simulate_one_element(
            X_at_nodes=R_at_nodes,
            U_at_nodes=1e-100 * np.ones_like(R_at_nodes),
            element_type=element_type,
            load_name="dummy_load",
            folder="res/element_definitions"
        )
        for element_type, R_at_nodes in R_at_nodes_of_element.items()
    }

    integration_points = dict()
    integration_weights = dict()

    for element_type, data in simulation_results.items():
        int_point_data = data["integration_points"]

        # integration points
        int_point_for_element_type = np.array(int_point_data["COORD"], dtype=float)
        int_point_for_element_type = np.around(int_point_for_element_type, 14)
        integration_points[element_type] = int_point_for_element_type.tolist()

        # integration weights
        int_weights_for_element_type = np.array(int_point_data["IVOL"], dtype=float).flatten()
        int_weights_for_element_type = np.around(int_weights_for_element_type, 15)
        integration_weights[element_type] = int_weights_for_element_type.tolist()

    with open(integration_data_file_path, "w", encoding="utf-8") as fh:
        json.dump(
            dict(
                README="COMPUTER GENERATED FILE. DO NOT MODIFY! "
                       f"This file is used by {os.path.basename(__file__)} to manager integration points and weights "
                       f"of various element types.",
                X_at_integration_points_for_element=integration_points,
                weights_of_integration_points_of_element=integration_weights
            ),
            fh,
            indent=4
        )


# update R_at_integration_points_of_element and weights_of_integration_points_of_element
_read_integration_data_file()
