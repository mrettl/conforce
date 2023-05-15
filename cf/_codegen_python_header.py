"""
# DO NOT IMPORT OR RUN THIS FILE. THIS IS JUST A TEMPLATE USED BY THE CODE GENERATION!
This module is a C-extension and provides a fast implementation for the computation of:
 - the deformation gradient
 - the first Piola-Kirchhoff stress tensor
 - the configurational stresses
 - the configurational forces


Examples
--------

>>> import numpy
>>> element_type = "CPE4R"
>>> method = "mbf"
>>> X_at_nodes = numpy.array([
...     [0., 0.],
...     [2., 0.],
...     [2., 2.],
...     [0., 2.],
... ])
>>> F_expected = numpy.array([
...     [1.1, 0.0],
...     [0.0, 0.9],
... ])
>>> U_at_nodes = X_at_nodes @ (F_expected - np.eye(2))
>>> S_at_int_points = numpy.array([[
...     [200., 0.],
...     [0., 0.]
... ]])
>>> e_at_int_points = numpy.array([20.])

The deformation gradient for the first (and only) element is evaluated
at the first integration point in the element.

>>> element_id = 0
>>> int_point_id = 0
>>> compute_F(
...     [X_at_nodes],
...     [U_at_nodes],
...     element_type
... )[element_id, int_point_id]
array([[1.1, 0. ],
       [0. , 0.9]])

The first Piola-Kirchhoff stress tensor is computed at the same position
as the deformation gradient.

>>> compute_P(
...     [X_at_nodes],
...     [U_at_nodes],
...     [S_at_int_points],
...     element_type
... )[element_id, int_point_id]
array([[180.,   0.],
       [  0.,   0.]])

The configuration stresses are computed using the motion based formulation "mbf" method
(see :py:func:`cf.expressions.eval_CS_mbf`).

>>> compute_CS(
...     [e_at_int_points],
...     [X_at_nodes],
...     [U_at_nodes],
...     [S_at_int_points],
...     element_type,
...     method
... )[element_id, int_point_id]
array([[-178.,   -0.],
       [  -0.,   20.]])

The configuration forces are evaluated at the nodes of each element.

>>> compute_CF(
...     [e_at_int_points],
...     [X_at_nodes],
...     [U_at_nodes],
...     [S_at_int_points],
...     element_type,
...     method
... )[element_id]
array([[ 178.,  -20.],
       [-178.,  -20.],
       [-178.,   20.],
       [ 178.,   20.]])
"""
import ctypes
import os
from collections import namedtuple

import numpy as np


# load c library
if os.name == 'nt':
    # windows
    lib = ctypes.cdll.LoadLibrary(os.path.abspath(os.path.join(
        __file__, os.path.pardir, 'REPLACE_THIS_BY_LIBRARY_FILE_NAME.dll'
    )))
else:
    # linux
    lib = ctypes.cdll.LoadLibrary(os.path.abspath(os.path.join(
        __file__, os.path.pardir, 'REPLACE_THIS_BY_LIBRARY_FILE_NAME.so'
    )))


# function lookup dictionaries
ElementInfo = namedtuple("ElementInfo", [
    "number_of_dimensions",
    "number_of_nodes",
    "number_of_integration_points"
])
map_type_to_info = dict()
map_type_to_F_function = dict()
map_type_to_P_function = dict()
map_type_and_method_to_CS_function = dict()
map_type_and_method_to_CF_function = dict()


def compute_F(
        X_at_nodes,
        U_at_nodes,
        element_type
):
    """
    Computes the deformation gradients for num_elem elements of the given element type.
    This element type has n nodes and ips integration points.
    2d elements have d=2 dimensions. 3d elements have d=3 dimensions.

    .. seealso::
        Module :py:mod:`cf.element_definitions` for a list of supported element types

    :param X_at_nodes: Array of shape (num_elem, n, d) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, d) containing the displacements at n nodes of num_elem elements.
    :param element_type: str, name of an element type
    :return: F_at_int_points: Array of shape (num_elem, ips, d, d) containing the deformation gradients
        evaluated on ips integration points for num_elem element.
    """
    fun = map_type_to_F_function[str(element_type)]
    return fun(
        X_at_nodes,
        U_at_nodes
    )


def compute_P(
        X_at_nodes,
        U_at_nodes,
        S_at_int_points,
        element_type
):
    """
    Computes the first Piola-Kirchhoff stress tensors for num_elem elements of the given element type.
    This element type has n nodes and ips integration points.
    2d elements have d=2 dimensions. 3d elements have d=3 dimensions.

    .. seealso::
        Module :py:mod:`cf.element_definitions` for a list of supported element types

    :param X_at_nodes: Array of shape (num_elem, n, d) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, d) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, d, d)
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :param element_type: str, name of an element type
    :return: P_at_int_points: Array of shape (num_elem, ips, d, d) containing the 1. Piola-Kirchhoff stress tensors
        evaluated on ips integration points for num_elem element.
    """
    fun = map_type_to_P_function[str(element_type)]
    return fun(
        X_at_nodes,
        U_at_nodes,
        S_at_int_points
    )


def compute_CS(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points,
        element_type,
        method
):
    """
    Computes the configurational stresses for num_elem elements of the given element type.
    This element type has n nodes and ips integration points.
    2d elements have d=2 dimensions. 3d elements have d=3 dimensions.

    .. seealso::
        - Module :py:mod:`cf.element_definitions` for a list of supported element types.
        - Supported methods are:

            - "mbf" :py:func:`cf.expressions.eval_CS_mbf`
            - "dbf" :py:func:`cf.expressions.eval_CS_dbf`

    :param e_at_int_points: Array of shape (num_elem,) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, d) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, d) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, d, d)
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :param element_type: str, name of an element type
    :param method: str, name of the method for the computation of the configurational stress.
    :return: CS_at_int_points: Array of shape (num_elem, ips, d, d) containing the configurational stresses
        evaluated on ips integration points for num_elem element.
    """
    fun = map_type_and_method_to_CS_function[(str(element_type), str(method))]
    return fun(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points
    )


def compute_CF(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points,
        element_type,
        method
):
    """
    Computes the configurational forces for num_elem elements of the given element type.
    This element type has n nodes and ips integration points.
    2d elements have d=2 dimensions. 3d elements have d=3 dimensions.

    .. seealso::
        - Module :py:mod:`cf.element_definitions` for a list of supported element types.
        - Supported methods are:

            - "mbf" :py:func:`cf.expressions.eval_CS_mbf`
            - "dbf" :py:func:`cf.expressions.eval_CS_dbf`

    :param e_at_int_points: Array of shape (num_elem,) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, d) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, d) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, d, d)
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :param element_type: str, name of an element type
    :param method: str, name of the method for the computation of the configurational stress.
    :return: CF_at_nodes: Array of shape (num_elem, n, d) containing the configurational forces
        evaluated on n nodes for num_elem element.
    """
    fun = map_type_and_method_to_CF_function[(str(element_type), str(method))]
    return fun(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points
    )

