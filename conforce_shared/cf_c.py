"""
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
(see :py:func:`conforce.expressions.eval_CS_mbf`).

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
        __file__, os.path.pardir, 'cf_c.dll'
    )))
else:
    # linux
    lib = ctypes.cdll.LoadLibrary(os.path.abspath(os.path.join(
        __file__, os.path.pardir, 'cf_c.so'
    )))


# function lookup dictionaries
ElementInfo = namedtuple("ElementInfo", [
    "number_of_dimensions",
    "number_of_nodes",
    "number_of_integration_points"
])

map_type_to_info = dict()
"""mapping of element type name to :py:class:`ElementInfo`"""

map_type_to_F_function = dict()
"""mapping of element type name to the function that computes the deformation gradient"""

map_type_to_P_function = dict()
"""mapping of element type name to the function that computes the first Piola-Kirchhoff stresses"""

map_type_and_method_to_CS_function = dict()
"""mapping of element type name and method to the function that computes the configurational stresses"""

map_type_and_method_to_CF_function = dict()
"""mapping of element type name and method  to the function that computes the configurational forces"""


def compute_F(
        X_at_nodes,
        U_at_nodes,
        element_type
):
    """
    Compute the deformation gradients for num_elem elements of the given element type.
    This element type has n nodes and ips integration points.
    2d elements have d=2 dimensions. 3d elements have d=3 dimensions.

    .. seealso::
        Module :py:mod:`conforce.element_definitions` for a list of supported element types

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
    Compute the first Piola-Kirchhoff stress tensors for num_elem elements of the given element type.
    This element type has n nodes and ips integration points.
    2d elements have d=2 dimensions. 3d elements have d=3 dimensions.

    .. seealso::
        Module :py:mod:`conforce.element_definitions` for a list of supported element types

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
    Compute the configurational stresses for num_elem elements of the given element type.
    This element type has n nodes and ips integration points.
    2d elements have d=2 dimensions. 3d elements have d=3 dimensions.

    .. seealso::
        - Module :py:mod:`conforce.element_definitions` for a list of supported element types.
        - Supported methods are:

            - "mbf" :py:func:`conforce.expressions.eval_CS_mbf`
            - "dbf" :py:func:`conforce.expressions.eval_CS_dbf`

    :param e_at_int_points: Array of shape (num_elem, ips)
        containing the internal energy densities at ips integration points for num_elem elements.
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
    Compute the configurational forces for num_elem elements of the given element type.
    This element type has n nodes and ips integration points.
    2d elements have d=2 dimensions. 3d elements have d=3 dimensions.

    .. seealso::
        - Module :py:mod:`conforce.element_definitions` for a list of supported element types.
        - Supported methods are:

            - "mbf" :py:func:`conforce.expressions.eval_CS_mbf`
            - "dbf" :py:func:`conforce.expressions.eval_CS_dbf`

    :param e_at_int_points: Array of shape (num_elem, ips)
        containing the internal energy densities at ips integration points for num_elem elements.
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


map_type_to_info['C3D20'] = ElementInfo(
    number_of_dimensions=3,
    number_of_nodes=20,
    number_of_integration_points=27
)


def compute_F_for_C3D20(
        X_at_nodes,
        U_at_nodes):
    """
    Compute the deformation gradients for num_elem elements of typ C3D20.
    Each element has n=20 nodes and ips=27 integration points.

    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :return: F_at_int_points: Array of shape (num_elem, ips, 3, 3) containing the deformation gradients
        evaluated on ips integration points for num_elem element.
    """

    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    F_at_int_points = np.zeros((num_elem, 27, 3, 3), dtype=np.float64)

    assert X_at_nodes.shape == (num_elem, 20, 3)
    assert U_at_nodes.shape == (num_elem, 20, 3)

    lib.compute_F_for_C3D20(
        ctypes.c_size_t(num_elem),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        F_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return F_at_int_points


map_type_to_F_function['C3D20'] = compute_F_for_C3D20


def compute_P_for_C3D20(
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the first Piola-Kirchhoff stress tensors for num_elem elements of typ C3D20.
    Each element has n=20 nodes and ips=27 integration points.

    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3)
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: P_at_int_points: Array of shape (num_elem, ips, 3, 3) containing the 1. Piola-Kirchhoff stress tensors
        evaluated on ips integration points for num_elem element.
    """

    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    P_at_int_points = np.zeros((num_elem, 27, 3, 3), dtype=np.float64)

    assert X_at_nodes.shape == (num_elem, 20, 3)
    assert U_at_nodes.shape == (num_elem, 20, 3)
    assert S_at_int_points.shape == (num_elem, 27, 3, 3)

    lib.compute_P_for_C3D20(
        ctypes.c_size_t(num_elem),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        P_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return P_at_int_points


map_type_to_P_function['C3D20'] = compute_P_for_C3D20


def compute_CS_for_C3D20_using_dbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational stresses for num_elem elements of typ C3D20.
    Each element has n=20 nodes and ips=27 integration points.

    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CS_at_int_points: Array of shape (num_elem, ips, 3, 3) containing the configurational stresses
        evaluated on ips integration points for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CS_at_int_points = np.zeros((num_elem, 27, 3, 3), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 27)
    assert X_at_nodes.shape == (num_elem, 20, 3)
    assert U_at_nodes.shape == (num_elem, 20, 3)
    assert S_at_int_points.shape == (num_elem, 27, 3, 3)

    lib.compute_CS_for_C3D20_using_dbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CS_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return CS_at_int_points


map_type_and_method_to_CS_function[('C3D20', 'dbf')] = compute_CS_for_C3D20_using_dbf


def compute_CF_for_C3D20_using_dbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational forces for num_elem elements of typ C3D20.
    Each element has n=20 nodes and ips=27 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CF_at_nodes: Array of shape (num_elem, n, 3) containing the configurational forces
        evaluated on n nodes for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CF_at_nodes = np.zeros((num_elem, 20, 3), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 27)
    assert X_at_nodes.shape == (num_elem, 20, 3)
    assert U_at_nodes.shape == (num_elem, 20, 3)
    assert S_at_int_points.shape == (num_elem, 27, 3, 3)

    lib.compute_CF_for_C3D20_using_dbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CF_at_nodes.ctypes.data_as(ctypes.c_void_p)
    )

    return CF_at_nodes


map_type_and_method_to_CF_function[('C3D20', 'dbf')] = compute_CF_for_C3D20_using_dbf

map_type_to_info['C3D20'] = ElementInfo(
    number_of_dimensions=3,
    number_of_nodes=20,
    number_of_integration_points=27
)


def compute_CS_for_C3D20_using_mbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational stresses for num_elem elements of typ C3D20.
    Each element has n=20 nodes and ips=27 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CS_at_int_points: Array of shape (num_elem, ips, 3, 3) containing the configurational stresses
        evaluated on ips integration points for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CS_at_int_points = np.zeros((num_elem, 27, 3, 3), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 27)
    assert X_at_nodes.shape == (num_elem, 20, 3)
    assert U_at_nodes.shape == (num_elem, 20, 3)
    assert S_at_int_points.shape == (num_elem, 27, 3, 3)

    lib.compute_CS_for_C3D20_using_mbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CS_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return CS_at_int_points


map_type_and_method_to_CS_function[('C3D20', 'mbf')] = compute_CS_for_C3D20_using_mbf


def compute_CF_for_C3D20_using_mbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational forces for num_elem elements of typ C3D20.
    Each element has n=20 nodes and ips=27 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CF_at_nodes: Array of shape (num_elem, n, 3) containing the configurational forces
        evaluated on n nodes for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CF_at_nodes = np.zeros((num_elem, 20, 3), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 27)
    assert X_at_nodes.shape == (num_elem, 20, 3)
    assert U_at_nodes.shape == (num_elem, 20, 3)
    assert S_at_int_points.shape == (num_elem, 27, 3, 3)

    lib.compute_CF_for_C3D20_using_mbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CF_at_nodes.ctypes.data_as(ctypes.c_void_p)
    )

    return CF_at_nodes


map_type_and_method_to_CF_function[('C3D20', 'mbf')] = compute_CF_for_C3D20_using_mbf

map_type_to_info['C3D20R'] = ElementInfo(
    number_of_dimensions=3,
    number_of_nodes=20,
    number_of_integration_points=8
)


def compute_F_for_C3D20R(
        X_at_nodes,
        U_at_nodes):
    """
    Compute the deformation gradients for num_elem elements of typ C3D20R.
    Each element has n=20 nodes and ips=8 integration points.

    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :return: F_at_int_points: Array of shape (num_elem, ips, 3, 3) containing the deformation gradients 
        evaluated on ips integration points for num_elem element.
    """

    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    F_at_int_points = np.zeros((num_elem, 8, 3, 3), dtype=np.float64)

    assert X_at_nodes.shape == (num_elem, 20, 3)
    assert U_at_nodes.shape == (num_elem, 20, 3)

    lib.compute_F_for_C3D20R(
        ctypes.c_size_t(num_elem),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        F_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return F_at_int_points


map_type_to_F_function['C3D20R'] = compute_F_for_C3D20R


def compute_P_for_C3D20R(
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the first Piola-Kirchhoff stress tensors for num_elem elements of typ C3D20R.
    Each element has n=20 nodes and ips=8 integration points.
    
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: P_at_int_points: Array of shape (num_elem, ips, 3, 3) containing the 1. Piola-Kirchhoff stress tensors
        evaluated on ips integration points for num_elem element.
    """

    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    P_at_int_points = np.zeros((num_elem, 8, 3, 3), dtype=np.float64)

    assert X_at_nodes.shape == (num_elem, 20, 3)
    assert U_at_nodes.shape == (num_elem, 20, 3)
    assert S_at_int_points.shape == (num_elem, 8, 3, 3)

    lib.compute_P_for_C3D20R(
        ctypes.c_size_t(num_elem),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        P_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return P_at_int_points


map_type_to_P_function['C3D20R'] = compute_P_for_C3D20R


def compute_CS_for_C3D20R_using_dbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational stresses for num_elem elements of typ C3D20R.
    Each element has n=20 nodes and ips=8 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CS_at_int_points: Array of shape (num_elem, ips, 3, 3) containing the configurational stresses
        evaluated on ips integration points for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CS_at_int_points = np.zeros((num_elem, 8, 3, 3), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 8)
    assert X_at_nodes.shape == (num_elem, 20, 3)
    assert U_at_nodes.shape == (num_elem, 20, 3)
    assert S_at_int_points.shape == (num_elem, 8, 3, 3)

    lib.compute_CS_for_C3D20R_using_dbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CS_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return CS_at_int_points


map_type_and_method_to_CS_function[('C3D20R', 'dbf')] = compute_CS_for_C3D20R_using_dbf


def compute_CF_for_C3D20R_using_dbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational forces for num_elem elements of typ C3D20R.
    Each element has n=20 nodes and ips=8 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CF_at_nodes: Array of shape (num_elem, n, 3) containing the configurational forces
        evaluated on n nodes for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CF_at_nodes = np.zeros((num_elem, 20, 3), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 8)
    assert X_at_nodes.shape == (num_elem, 20, 3)
    assert U_at_nodes.shape == (num_elem, 20, 3)
    assert S_at_int_points.shape == (num_elem, 8, 3, 3)

    lib.compute_CF_for_C3D20R_using_dbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CF_at_nodes.ctypes.data_as(ctypes.c_void_p)
    )

    return CF_at_nodes


map_type_and_method_to_CF_function[('C3D20R', 'dbf')] = compute_CF_for_C3D20R_using_dbf

map_type_to_info['C3D20R'] = ElementInfo(
    number_of_dimensions=3,
    number_of_nodes=20,
    number_of_integration_points=8
)


def compute_CS_for_C3D20R_using_mbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational stresses for num_elem elements of typ C3D20R.
    Each element has n=20 nodes and ips=8 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CS_at_int_points: Array of shape (num_elem, ips, 3, 3) containing the configurational stresses
        evaluated on ips integration points for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CS_at_int_points = np.zeros((num_elem, 8, 3, 3), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 8)
    assert X_at_nodes.shape == (num_elem, 20, 3)
    assert U_at_nodes.shape == (num_elem, 20, 3)
    assert S_at_int_points.shape == (num_elem, 8, 3, 3)

    lib.compute_CS_for_C3D20R_using_mbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CS_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return CS_at_int_points


map_type_and_method_to_CS_function[('C3D20R', 'mbf')] = compute_CS_for_C3D20R_using_mbf


def compute_CF_for_C3D20R_using_mbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational forces for num_elem elements of typ C3D20R.
    Each element has n=20 nodes and ips=8 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CF_at_nodes: Array of shape (num_elem, n, 3) containing the configurational forces
        evaluated on n nodes for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CF_at_nodes = np.zeros((num_elem, 20, 3), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 8)
    assert X_at_nodes.shape == (num_elem, 20, 3)
    assert U_at_nodes.shape == (num_elem, 20, 3)
    assert S_at_int_points.shape == (num_elem, 8, 3, 3)

    lib.compute_CF_for_C3D20R_using_mbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CF_at_nodes.ctypes.data_as(ctypes.c_void_p)
    )

    return CF_at_nodes


map_type_and_method_to_CF_function[('C3D20R', 'mbf')] = compute_CF_for_C3D20R_using_mbf

map_type_to_info['CPE8'] = ElementInfo(
    number_of_dimensions=2,
    number_of_nodes=8,
    number_of_integration_points=9
)


def compute_F_for_CPE8(
        X_at_nodes,
        U_at_nodes):
    """
    Compute the deformation gradients for num_elem elements of typ CPE8.
    Each element has n=8 nodes and ips=9 integration points.

    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :return: F_at_int_points: Array of shape (num_elem, ips, 2, 2) containing the deformation gradients 
        evaluated on ips integration points for num_elem element.
    """

    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    F_at_int_points = np.zeros((num_elem, 9, 2, 2), dtype=np.float64)

    assert X_at_nodes.shape == (num_elem, 8, 2)
    assert U_at_nodes.shape == (num_elem, 8, 2)

    lib.compute_F_for_CPE8(
        ctypes.c_size_t(num_elem),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        F_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return F_at_int_points


map_type_to_F_function['CPE8'] = compute_F_for_CPE8


def compute_P_for_CPE8(
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the first Piola-Kirchhoff stress tensors for num_elem elements of typ CPE8.
    Each element has n=8 nodes and ips=9 integration points.
    
    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 2, 2) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: P_at_int_points: Array of shape (num_elem, ips, 2, 2) containing the 1. Piola-Kirchhoff stress tensors
        evaluated on ips integration points for num_elem element.
    """

    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    P_at_int_points = np.zeros((num_elem, 9, 2, 2), dtype=np.float64)

    assert X_at_nodes.shape == (num_elem, 8, 2)
    assert U_at_nodes.shape == (num_elem, 8, 2)
    assert S_at_int_points.shape == (num_elem, 9, 2, 2)

    lib.compute_P_for_CPE8(
        ctypes.c_size_t(num_elem),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        P_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return P_at_int_points


map_type_to_P_function['CPE8'] = compute_P_for_CPE8


def compute_CS_for_CPE8_using_dbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational stresses for num_elem elements of typ CPE8.
    Each element has n=8 nodes and ips=9 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 2, 2) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CS_at_int_points: Array of shape (num_elem, ips, 2, 2) containing the configurational stresses
        evaluated on ips integration points for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CS_at_int_points = np.zeros((num_elem, 9, 2, 2), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 9)
    assert X_at_nodes.shape == (num_elem, 8, 2)
    assert U_at_nodes.shape == (num_elem, 8, 2)
    assert S_at_int_points.shape == (num_elem, 9, 2, 2)

    lib.compute_CS_for_CPE8_using_dbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CS_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return CS_at_int_points


map_type_and_method_to_CS_function[('CPE8', 'dbf')] = compute_CS_for_CPE8_using_dbf


def compute_CF_for_CPE8_using_dbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational forces for num_elem elements of typ CPE8.
    Each element has n=8 nodes and ips=9 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 2, 2) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CF_at_nodes: Array of shape (num_elem, n, 2) containing the configurational forces
        evaluated on n nodes for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CF_at_nodes = np.zeros((num_elem, 8, 2), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 9)
    assert X_at_nodes.shape == (num_elem, 8, 2)
    assert U_at_nodes.shape == (num_elem, 8, 2)
    assert S_at_int_points.shape == (num_elem, 9, 2, 2)

    lib.compute_CF_for_CPE8_using_dbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CF_at_nodes.ctypes.data_as(ctypes.c_void_p)
    )

    return CF_at_nodes


map_type_and_method_to_CF_function[('CPE8', 'dbf')] = compute_CF_for_CPE8_using_dbf

map_type_to_info['CPE8'] = ElementInfo(
    number_of_dimensions=2,
    number_of_nodes=8,
    number_of_integration_points=9
)


def compute_CS_for_CPE8_using_mbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational stresses for num_elem elements of typ CPE8.
    Each element has n=8 nodes and ips=9 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 2, 2) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CS_at_int_points: Array of shape (num_elem, ips, 2, 2) containing the configurational stresses
        evaluated on ips integration points for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CS_at_int_points = np.zeros((num_elem, 9, 2, 2), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 9)
    assert X_at_nodes.shape == (num_elem, 8, 2)
    assert U_at_nodes.shape == (num_elem, 8, 2)
    assert S_at_int_points.shape == (num_elem, 9, 2, 2)

    lib.compute_CS_for_CPE8_using_mbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CS_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return CS_at_int_points


map_type_and_method_to_CS_function[('CPE8', 'mbf')] = compute_CS_for_CPE8_using_mbf


def compute_CF_for_CPE8_using_mbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational forces for num_elem elements of typ CPE8.
    Each element has n=8 nodes and ips=9 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 2, 2) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CF_at_nodes: Array of shape (num_elem, n, 2) containing the configurational forces
        evaluated on n nodes for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CF_at_nodes = np.zeros((num_elem, 8, 2), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 9)
    assert X_at_nodes.shape == (num_elem, 8, 2)
    assert U_at_nodes.shape == (num_elem, 8, 2)
    assert S_at_int_points.shape == (num_elem, 9, 2, 2)

    lib.compute_CF_for_CPE8_using_mbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CF_at_nodes.ctypes.data_as(ctypes.c_void_p)
    )

    return CF_at_nodes


map_type_and_method_to_CF_function[('CPE8', 'mbf')] = compute_CF_for_CPE8_using_mbf

map_type_to_info['CPE8R'] = ElementInfo(
    number_of_dimensions=2,
    number_of_nodes=8,
    number_of_integration_points=4
)


def compute_F_for_CPE8R(
        X_at_nodes,
        U_at_nodes):
    """
    Compute the deformation gradients for num_elem elements of typ CPE8R.
    Each element has n=8 nodes and ips=4 integration points.

    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :return: F_at_int_points: Array of shape (num_elem, ips, 2, 2) containing the deformation gradients 
        evaluated on ips integration points for num_elem element.
    """

    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    F_at_int_points = np.zeros((num_elem, 4, 2, 2), dtype=np.float64)

    assert X_at_nodes.shape == (num_elem, 8, 2)
    assert U_at_nodes.shape == (num_elem, 8, 2)

    lib.compute_F_for_CPE8R(
        ctypes.c_size_t(num_elem),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        F_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return F_at_int_points


map_type_to_F_function['CPE8R'] = compute_F_for_CPE8R


def compute_P_for_CPE8R(
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the first Piola-Kirchhoff stress tensors for num_elem elements of typ CPE8R.
    Each element has n=8 nodes and ips=4 integration points.
    
    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 2, 2) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: P_at_int_points: Array of shape (num_elem, ips, 2, 2) containing the 1. Piola-Kirchhoff stress tensors
        evaluated on ips integration points for num_elem element.
    """

    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    P_at_int_points = np.zeros((num_elem, 4, 2, 2), dtype=np.float64)

    assert X_at_nodes.shape == (num_elem, 8, 2)
    assert U_at_nodes.shape == (num_elem, 8, 2)
    assert S_at_int_points.shape == (num_elem, 4, 2, 2)

    lib.compute_P_for_CPE8R(
        ctypes.c_size_t(num_elem),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        P_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return P_at_int_points


map_type_to_P_function['CPE8R'] = compute_P_for_CPE8R


def compute_CS_for_CPE8R_using_dbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational stresses for num_elem elements of typ CPE8R.
    Each element has n=8 nodes and ips=4 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 2, 2) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CS_at_int_points: Array of shape (num_elem, ips, 2, 2) containing the configurational stresses
        evaluated on ips integration points for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CS_at_int_points = np.zeros((num_elem, 4, 2, 2), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 4)
    assert X_at_nodes.shape == (num_elem, 8, 2)
    assert U_at_nodes.shape == (num_elem, 8, 2)
    assert S_at_int_points.shape == (num_elem, 4, 2, 2)

    lib.compute_CS_for_CPE8R_using_dbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CS_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return CS_at_int_points


map_type_and_method_to_CS_function[('CPE8R', 'dbf')] = compute_CS_for_CPE8R_using_dbf


def compute_CF_for_CPE8R_using_dbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational forces for num_elem elements of typ CPE8R.
    Each element has n=8 nodes and ips=4 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 2, 2) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CF_at_nodes: Array of shape (num_elem, n, 2) containing the configurational forces
        evaluated on n nodes for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CF_at_nodes = np.zeros((num_elem, 8, 2), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 4)
    assert X_at_nodes.shape == (num_elem, 8, 2)
    assert U_at_nodes.shape == (num_elem, 8, 2)
    assert S_at_int_points.shape == (num_elem, 4, 2, 2)

    lib.compute_CF_for_CPE8R_using_dbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CF_at_nodes.ctypes.data_as(ctypes.c_void_p)
    )

    return CF_at_nodes


map_type_and_method_to_CF_function[('CPE8R', 'dbf')] = compute_CF_for_CPE8R_using_dbf

map_type_to_info['CPE8R'] = ElementInfo(
    number_of_dimensions=2,
    number_of_nodes=8,
    number_of_integration_points=4
)


def compute_CS_for_CPE8R_using_mbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational stresses for num_elem elements of typ CPE8R.
    Each element has n=8 nodes and ips=4 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 2, 2) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CS_at_int_points: Array of shape (num_elem, ips, 2, 2) containing the configurational stresses
        evaluated on ips integration points for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CS_at_int_points = np.zeros((num_elem, 4, 2, 2), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 4)
    assert X_at_nodes.shape == (num_elem, 8, 2)
    assert U_at_nodes.shape == (num_elem, 8, 2)
    assert S_at_int_points.shape == (num_elem, 4, 2, 2)

    lib.compute_CS_for_CPE8R_using_mbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CS_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return CS_at_int_points


map_type_and_method_to_CS_function[('CPE8R', 'mbf')] = compute_CS_for_CPE8R_using_mbf


def compute_CF_for_CPE8R_using_mbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational forces for num_elem elements of typ CPE8R.
    Each element has n=8 nodes and ips=4 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 2, 2) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CF_at_nodes: Array of shape (num_elem, n, 2) containing the configurational forces
        evaluated on n nodes for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CF_at_nodes = np.zeros((num_elem, 8, 2), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 4)
    assert X_at_nodes.shape == (num_elem, 8, 2)
    assert U_at_nodes.shape == (num_elem, 8, 2)
    assert S_at_int_points.shape == (num_elem, 4, 2, 2)

    lib.compute_CF_for_CPE8R_using_mbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CF_at_nodes.ctypes.data_as(ctypes.c_void_p)
    )

    return CF_at_nodes


map_type_and_method_to_CF_function[('CPE8R', 'mbf')] = compute_CF_for_CPE8R_using_mbf

map_type_to_info['C3D8'] = ElementInfo(
    number_of_dimensions=3,
    number_of_nodes=8,
    number_of_integration_points=8
)


def compute_F_for_C3D8(
        X_at_nodes,
        U_at_nodes):
    """
    Compute the deformation gradients for num_elem elements of typ C3D8.
    Each element has n=8 nodes and ips=8 integration points.

    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :return: F_at_int_points: Array of shape (num_elem, ips, 3, 3) containing the deformation gradients 
        evaluated on ips integration points for num_elem element.
    """

    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    F_at_int_points = np.zeros((num_elem, 8, 3, 3), dtype=np.float64)

    assert X_at_nodes.shape == (num_elem, 8, 3)
    assert U_at_nodes.shape == (num_elem, 8, 3)

    lib.compute_F_for_C3D8(
        ctypes.c_size_t(num_elem),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        F_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return F_at_int_points


map_type_to_F_function['C3D8'] = compute_F_for_C3D8


def compute_P_for_C3D8(
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the first Piola-Kirchhoff stress tensors for num_elem elements of typ C3D8.
    Each element has n=8 nodes and ips=8 integration points.
    
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: P_at_int_points: Array of shape (num_elem, ips, 3, 3) containing the 1. Piola-Kirchhoff stress tensors
        evaluated on ips integration points for num_elem element.
    """

    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    P_at_int_points = np.zeros((num_elem, 8, 3, 3), dtype=np.float64)

    assert X_at_nodes.shape == (num_elem, 8, 3)
    assert U_at_nodes.shape == (num_elem, 8, 3)
    assert S_at_int_points.shape == (num_elem, 8, 3, 3)

    lib.compute_P_for_C3D8(
        ctypes.c_size_t(num_elem),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        P_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return P_at_int_points


map_type_to_P_function['C3D8'] = compute_P_for_C3D8


def compute_CS_for_C3D8_using_dbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational stresses for num_elem elements of typ C3D8.
    Each element has n=8 nodes and ips=8 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CS_at_int_points: Array of shape (num_elem, ips, 3, 3) containing the configurational stresses
        evaluated on ips integration points for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CS_at_int_points = np.zeros((num_elem, 8, 3, 3), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 8)
    assert X_at_nodes.shape == (num_elem, 8, 3)
    assert U_at_nodes.shape == (num_elem, 8, 3)
    assert S_at_int_points.shape == (num_elem, 8, 3, 3)

    lib.compute_CS_for_C3D8_using_dbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CS_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return CS_at_int_points


map_type_and_method_to_CS_function[('C3D8', 'dbf')] = compute_CS_for_C3D8_using_dbf


def compute_CF_for_C3D8_using_dbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational forces for num_elem elements of typ C3D8.
    Each element has n=8 nodes and ips=8 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CF_at_nodes: Array of shape (num_elem, n, 3) containing the configurational forces
        evaluated on n nodes for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CF_at_nodes = np.zeros((num_elem, 8, 3), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 8)
    assert X_at_nodes.shape == (num_elem, 8, 3)
    assert U_at_nodes.shape == (num_elem, 8, 3)
    assert S_at_int_points.shape == (num_elem, 8, 3, 3)

    lib.compute_CF_for_C3D8_using_dbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CF_at_nodes.ctypes.data_as(ctypes.c_void_p)
    )

    return CF_at_nodes


map_type_and_method_to_CF_function[('C3D8', 'dbf')] = compute_CF_for_C3D8_using_dbf

map_type_to_info['C3D8'] = ElementInfo(
    number_of_dimensions=3,
    number_of_nodes=8,
    number_of_integration_points=8
)


def compute_CS_for_C3D8_using_mbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational stresses for num_elem elements of typ C3D8.
    Each element has n=8 nodes and ips=8 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CS_at_int_points: Array of shape (num_elem, ips, 3, 3) containing the configurational stresses
        evaluated on ips integration points for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CS_at_int_points = np.zeros((num_elem, 8, 3, 3), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 8)
    assert X_at_nodes.shape == (num_elem, 8, 3)
    assert U_at_nodes.shape == (num_elem, 8, 3)
    assert S_at_int_points.shape == (num_elem, 8, 3, 3)

    lib.compute_CS_for_C3D8_using_mbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CS_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return CS_at_int_points


map_type_and_method_to_CS_function[('C3D8', 'mbf')] = compute_CS_for_C3D8_using_mbf


def compute_CF_for_C3D8_using_mbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational forces for num_elem elements of typ C3D8.
    Each element has n=8 nodes and ips=8 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CF_at_nodes: Array of shape (num_elem, n, 3) containing the configurational forces
        evaluated on n nodes for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CF_at_nodes = np.zeros((num_elem, 8, 3), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 8)
    assert X_at_nodes.shape == (num_elem, 8, 3)
    assert U_at_nodes.shape == (num_elem, 8, 3)
    assert S_at_int_points.shape == (num_elem, 8, 3, 3)

    lib.compute_CF_for_C3D8_using_mbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CF_at_nodes.ctypes.data_as(ctypes.c_void_p)
    )

    return CF_at_nodes


map_type_and_method_to_CF_function[('C3D8', 'mbf')] = compute_CF_for_C3D8_using_mbf

map_type_to_info['C3D8R'] = ElementInfo(
    number_of_dimensions=3,
    number_of_nodes=8,
    number_of_integration_points=1
)


def compute_F_for_C3D8R(
        X_at_nodes,
        U_at_nodes):
    """
    Compute the deformation gradients for num_elem elements of typ C3D8R.
    Each element has n=8 nodes and ips=1 integration points.

    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :return: F_at_int_points: Array of shape (num_elem, ips, 3, 3) containing the deformation gradients 
        evaluated on ips integration points for num_elem element.
    """

    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    F_at_int_points = np.zeros((num_elem, 1, 3, 3), dtype=np.float64)

    assert X_at_nodes.shape == (num_elem, 8, 3)
    assert U_at_nodes.shape == (num_elem, 8, 3)

    lib.compute_F_for_C3D8R(
        ctypes.c_size_t(num_elem),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        F_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return F_at_int_points


map_type_to_F_function['C3D8R'] = compute_F_for_C3D8R


def compute_P_for_C3D8R(
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the first Piola-Kirchhoff stress tensors for num_elem elements of typ C3D8R.
    Each element has n=8 nodes and ips=1 integration points.
    
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: P_at_int_points: Array of shape (num_elem, ips, 3, 3) containing the 1. Piola-Kirchhoff stress tensors
        evaluated on ips integration points for num_elem element.
    """

    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    P_at_int_points = np.zeros((num_elem, 1, 3, 3), dtype=np.float64)

    assert X_at_nodes.shape == (num_elem, 8, 3)
    assert U_at_nodes.shape == (num_elem, 8, 3)
    assert S_at_int_points.shape == (num_elem, 1, 3, 3)

    lib.compute_P_for_C3D8R(
        ctypes.c_size_t(num_elem),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        P_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return P_at_int_points


map_type_to_P_function['C3D8R'] = compute_P_for_C3D8R


def compute_CS_for_C3D8R_using_dbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational stresses for num_elem elements of typ C3D8R.
    Each element has n=8 nodes and ips=1 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CS_at_int_points: Array of shape (num_elem, ips, 3, 3) containing the configurational stresses
        evaluated on ips integration points for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CS_at_int_points = np.zeros((num_elem, 1, 3, 3), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 1)
    assert X_at_nodes.shape == (num_elem, 8, 3)
    assert U_at_nodes.shape == (num_elem, 8, 3)
    assert S_at_int_points.shape == (num_elem, 1, 3, 3)

    lib.compute_CS_for_C3D8R_using_dbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CS_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return CS_at_int_points


map_type_and_method_to_CS_function[('C3D8R', 'dbf')] = compute_CS_for_C3D8R_using_dbf


def compute_CF_for_C3D8R_using_dbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational forces for num_elem elements of typ C3D8R.
    Each element has n=8 nodes and ips=1 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CF_at_nodes: Array of shape (num_elem, n, 3) containing the configurational forces
        evaluated on n nodes for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CF_at_nodes = np.zeros((num_elem, 8, 3), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 1)
    assert X_at_nodes.shape == (num_elem, 8, 3)
    assert U_at_nodes.shape == (num_elem, 8, 3)
    assert S_at_int_points.shape == (num_elem, 1, 3, 3)

    lib.compute_CF_for_C3D8R_using_dbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CF_at_nodes.ctypes.data_as(ctypes.c_void_p)
    )

    return CF_at_nodes


map_type_and_method_to_CF_function[('C3D8R', 'dbf')] = compute_CF_for_C3D8R_using_dbf

map_type_to_info['C3D8R'] = ElementInfo(
    number_of_dimensions=3,
    number_of_nodes=8,
    number_of_integration_points=1
)


def compute_CS_for_C3D8R_using_mbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational stresses for num_elem elements of typ C3D8R.
    Each element has n=8 nodes and ips=1 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CS_at_int_points: Array of shape (num_elem, ips, 3, 3) containing the configurational stresses
        evaluated on ips integration points for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CS_at_int_points = np.zeros((num_elem, 1, 3, 3), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 1)
    assert X_at_nodes.shape == (num_elem, 8, 3)
    assert U_at_nodes.shape == (num_elem, 8, 3)
    assert S_at_int_points.shape == (num_elem, 1, 3, 3)

    lib.compute_CS_for_C3D8R_using_mbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CS_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return CS_at_int_points


map_type_and_method_to_CS_function[('C3D8R', 'mbf')] = compute_CS_for_C3D8R_using_mbf


def compute_CF_for_C3D8R_using_mbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational forces for num_elem elements of typ C3D8R.
    Each element has n=8 nodes and ips=1 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CF_at_nodes: Array of shape (num_elem, n, 3) containing the configurational forces
        evaluated on n nodes for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CF_at_nodes = np.zeros((num_elem, 8, 3), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 1)
    assert X_at_nodes.shape == (num_elem, 8, 3)
    assert U_at_nodes.shape == (num_elem, 8, 3)
    assert S_at_int_points.shape == (num_elem, 1, 3, 3)

    lib.compute_CF_for_C3D8R_using_mbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CF_at_nodes.ctypes.data_as(ctypes.c_void_p)
    )

    return CF_at_nodes


map_type_and_method_to_CF_function[('C3D8R', 'mbf')] = compute_CF_for_C3D8R_using_mbf

map_type_to_info['CPE4'] = ElementInfo(
    number_of_dimensions=2,
    number_of_nodes=4,
    number_of_integration_points=4
)


def compute_F_for_CPE4(
        X_at_nodes,
        U_at_nodes):
    """
    Compute the deformation gradients for num_elem elements of typ CPE4.
    Each element has n=4 nodes and ips=4 integration points.

    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :return: F_at_int_points: Array of shape (num_elem, ips, 2, 2) containing the deformation gradients 
        evaluated on ips integration points for num_elem element.
    """

    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    F_at_int_points = np.zeros((num_elem, 4, 2, 2), dtype=np.float64)

    assert X_at_nodes.shape == (num_elem, 4, 2)
    assert U_at_nodes.shape == (num_elem, 4, 2)

    lib.compute_F_for_CPE4(
        ctypes.c_size_t(num_elem),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        F_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return F_at_int_points


map_type_to_F_function['CPE4'] = compute_F_for_CPE4


def compute_P_for_CPE4(
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the first Piola-Kirchhoff stress tensors for num_elem elements of typ CPE4.
    Each element has n=4 nodes and ips=4 integration points.
    
    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 2, 2) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: P_at_int_points: Array of shape (num_elem, ips, 2, 2) containing the 1. Piola-Kirchhoff stress tensors
        evaluated on ips integration points for num_elem element.
    """

    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    P_at_int_points = np.zeros((num_elem, 4, 2, 2), dtype=np.float64)

    assert X_at_nodes.shape == (num_elem, 4, 2)
    assert U_at_nodes.shape == (num_elem, 4, 2)
    assert S_at_int_points.shape == (num_elem, 4, 2, 2)

    lib.compute_P_for_CPE4(
        ctypes.c_size_t(num_elem),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        P_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return P_at_int_points


map_type_to_P_function['CPE4'] = compute_P_for_CPE4


def compute_CS_for_CPE4_using_dbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational stresses for num_elem elements of typ CPE4.
    Each element has n=4 nodes and ips=4 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 2, 2) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CS_at_int_points: Array of shape (num_elem, ips, 2, 2) containing the configurational stresses
        evaluated on ips integration points for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CS_at_int_points = np.zeros((num_elem, 4, 2, 2), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 4)
    assert X_at_nodes.shape == (num_elem, 4, 2)
    assert U_at_nodes.shape == (num_elem, 4, 2)
    assert S_at_int_points.shape == (num_elem, 4, 2, 2)

    lib.compute_CS_for_CPE4_using_dbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CS_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return CS_at_int_points


map_type_and_method_to_CS_function[('CPE4', 'dbf')] = compute_CS_for_CPE4_using_dbf


def compute_CF_for_CPE4_using_dbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational forces for num_elem elements of typ CPE4.
    Each element has n=4 nodes and ips=4 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 2, 2) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CF_at_nodes: Array of shape (num_elem, n, 2) containing the configurational forces
        evaluated on n nodes for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CF_at_nodes = np.zeros((num_elem, 4, 2), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 4)
    assert X_at_nodes.shape == (num_elem, 4, 2)
    assert U_at_nodes.shape == (num_elem, 4, 2)
    assert S_at_int_points.shape == (num_elem, 4, 2, 2)

    lib.compute_CF_for_CPE4_using_dbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CF_at_nodes.ctypes.data_as(ctypes.c_void_p)
    )

    return CF_at_nodes


map_type_and_method_to_CF_function[('CPE4', 'dbf')] = compute_CF_for_CPE4_using_dbf

map_type_to_info['CPE4'] = ElementInfo(
    number_of_dimensions=2,
    number_of_nodes=4,
    number_of_integration_points=4
)


def compute_CS_for_CPE4_using_mbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational stresses for num_elem elements of typ CPE4.
    Each element has n=4 nodes and ips=4 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 2, 2) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CS_at_int_points: Array of shape (num_elem, ips, 2, 2) containing the configurational stresses
        evaluated on ips integration points for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CS_at_int_points = np.zeros((num_elem, 4, 2, 2), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 4)
    assert X_at_nodes.shape == (num_elem, 4, 2)
    assert U_at_nodes.shape == (num_elem, 4, 2)
    assert S_at_int_points.shape == (num_elem, 4, 2, 2)

    lib.compute_CS_for_CPE4_using_mbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CS_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return CS_at_int_points


map_type_and_method_to_CS_function[('CPE4', 'mbf')] = compute_CS_for_CPE4_using_mbf


def compute_CF_for_CPE4_using_mbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational forces for num_elem elements of typ CPE4.
    Each element has n=4 nodes and ips=4 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 2, 2) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CF_at_nodes: Array of shape (num_elem, n, 2) containing the configurational forces
        evaluated on n nodes for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CF_at_nodes = np.zeros((num_elem, 4, 2), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 4)
    assert X_at_nodes.shape == (num_elem, 4, 2)
    assert U_at_nodes.shape == (num_elem, 4, 2)
    assert S_at_int_points.shape == (num_elem, 4, 2, 2)

    lib.compute_CF_for_CPE4_using_mbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CF_at_nodes.ctypes.data_as(ctypes.c_void_p)
    )

    return CF_at_nodes


map_type_and_method_to_CF_function[('CPE4', 'mbf')] = compute_CF_for_CPE4_using_mbf

map_type_to_info['CPE4R'] = ElementInfo(
    number_of_dimensions=2,
    number_of_nodes=4,
    number_of_integration_points=1
)


def compute_F_for_CPE4R(
        X_at_nodes,
        U_at_nodes):
    """
    Compute the deformation gradients for num_elem elements of typ CPE4R.
    Each element has n=4 nodes and ips=1 integration points.

    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :return: F_at_int_points: Array of shape (num_elem, ips, 2, 2) containing the deformation gradients 
        evaluated on ips integration points for num_elem element.
    """

    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    F_at_int_points = np.zeros((num_elem, 1, 2, 2), dtype=np.float64)

    assert X_at_nodes.shape == (num_elem, 4, 2)
    assert U_at_nodes.shape == (num_elem, 4, 2)

    lib.compute_F_for_CPE4R(
        ctypes.c_size_t(num_elem),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        F_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return F_at_int_points


map_type_to_F_function['CPE4R'] = compute_F_for_CPE4R


def compute_P_for_CPE4R(
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the first Piola-Kirchhoff stress tensors for num_elem elements of typ CPE4R.
    Each element has n=4 nodes and ips=1 integration points.
    
    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 2, 2) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: P_at_int_points: Array of shape (num_elem, ips, 2, 2) containing the 1. Piola-Kirchhoff stress tensors
        evaluated on ips integration points for num_elem element.
    """

    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    P_at_int_points = np.zeros((num_elem, 1, 2, 2), dtype=np.float64)

    assert X_at_nodes.shape == (num_elem, 4, 2)
    assert U_at_nodes.shape == (num_elem, 4, 2)
    assert S_at_int_points.shape == (num_elem, 1, 2, 2)

    lib.compute_P_for_CPE4R(
        ctypes.c_size_t(num_elem),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        P_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return P_at_int_points


map_type_to_P_function['CPE4R'] = compute_P_for_CPE4R


def compute_CS_for_CPE4R_using_dbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational stresses for num_elem elements of typ CPE4R.
    Each element has n=4 nodes and ips=1 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 2, 2) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CS_at_int_points: Array of shape (num_elem, ips, 2, 2) containing the configurational stresses
        evaluated on ips integration points for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CS_at_int_points = np.zeros((num_elem, 1, 2, 2), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 1)
    assert X_at_nodes.shape == (num_elem, 4, 2)
    assert U_at_nodes.shape == (num_elem, 4, 2)
    assert S_at_int_points.shape == (num_elem, 1, 2, 2)

    lib.compute_CS_for_CPE4R_using_dbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CS_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return CS_at_int_points


map_type_and_method_to_CS_function[('CPE4R', 'dbf')] = compute_CS_for_CPE4R_using_dbf


def compute_CF_for_CPE4R_using_dbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational forces for num_elem elements of typ CPE4R.
    Each element has n=4 nodes and ips=1 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 2, 2) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CF_at_nodes: Array of shape (num_elem, n, 2) containing the configurational forces
        evaluated on n nodes for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CF_at_nodes = np.zeros((num_elem, 4, 2), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 1)
    assert X_at_nodes.shape == (num_elem, 4, 2)
    assert U_at_nodes.shape == (num_elem, 4, 2)
    assert S_at_int_points.shape == (num_elem, 1, 2, 2)

    lib.compute_CF_for_CPE4R_using_dbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CF_at_nodes.ctypes.data_as(ctypes.c_void_p)
    )

    return CF_at_nodes


map_type_and_method_to_CF_function[('CPE4R', 'dbf')] = compute_CF_for_CPE4R_using_dbf

map_type_to_info['CPE4R'] = ElementInfo(
    number_of_dimensions=2,
    number_of_nodes=4,
    number_of_integration_points=1
)


def compute_CS_for_CPE4R_using_mbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational stresses for num_elem elements of typ CPE4R.
    Each element has n=4 nodes and ips=1 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 2, 2) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CS_at_int_points: Array of shape (num_elem, ips, 2, 2) containing the configurational stresses
        evaluated on ips integration points for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CS_at_int_points = np.zeros((num_elem, 1, 2, 2), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 1)
    assert X_at_nodes.shape == (num_elem, 4, 2)
    assert U_at_nodes.shape == (num_elem, 4, 2)
    assert S_at_int_points.shape == (num_elem, 1, 2, 2)

    lib.compute_CS_for_CPE4R_using_mbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CS_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return CS_at_int_points


map_type_and_method_to_CS_function[('CPE4R', 'mbf')] = compute_CS_for_CPE4R_using_mbf


def compute_CF_for_CPE4R_using_mbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational forces for num_elem elements of typ CPE4R.
    Each element has n=4 nodes and ips=1 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 2, 2) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CF_at_nodes: Array of shape (num_elem, n, 2) containing the configurational forces
        evaluated on n nodes for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CF_at_nodes = np.zeros((num_elem, 4, 2), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 1)
    assert X_at_nodes.shape == (num_elem, 4, 2)
    assert U_at_nodes.shape == (num_elem, 4, 2)
    assert S_at_int_points.shape == (num_elem, 1, 2, 2)

    lib.compute_CF_for_CPE4R_using_mbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CF_at_nodes.ctypes.data_as(ctypes.c_void_p)
    )

    return CF_at_nodes


map_type_and_method_to_CF_function[('CPE4R', 'mbf')] = compute_CF_for_CPE4R_using_mbf

map_type_to_info['C3D10'] = ElementInfo(
    number_of_dimensions=3,
    number_of_nodes=10,
    number_of_integration_points=4
)


def compute_F_for_C3D10(
        X_at_nodes,
        U_at_nodes):
    """
    Compute the deformation gradients for num_elem elements of typ C3D10.
    Each element has n=10 nodes and ips=4 integration points.

    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :return: F_at_int_points: Array of shape (num_elem, ips, 3, 3) containing the deformation gradients 
        evaluated on ips integration points for num_elem element.
    """

    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    F_at_int_points = np.zeros((num_elem, 4, 3, 3), dtype=np.float64)

    assert X_at_nodes.shape == (num_elem, 10, 3)
    assert U_at_nodes.shape == (num_elem, 10, 3)

    lib.compute_F_for_C3D10(
        ctypes.c_size_t(num_elem),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        F_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return F_at_int_points


map_type_to_F_function['C3D10'] = compute_F_for_C3D10


def compute_P_for_C3D10(
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the first Piola-Kirchhoff stress tensors for num_elem elements of typ C3D10.
    Each element has n=10 nodes and ips=4 integration points.
    
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: P_at_int_points: Array of shape (num_elem, ips, 3, 3) containing the 1. Piola-Kirchhoff stress tensors
        evaluated on ips integration points for num_elem element.
    """

    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    P_at_int_points = np.zeros((num_elem, 4, 3, 3), dtype=np.float64)

    assert X_at_nodes.shape == (num_elem, 10, 3)
    assert U_at_nodes.shape == (num_elem, 10, 3)
    assert S_at_int_points.shape == (num_elem, 4, 3, 3)

    lib.compute_P_for_C3D10(
        ctypes.c_size_t(num_elem),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        P_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return P_at_int_points


map_type_to_P_function['C3D10'] = compute_P_for_C3D10


def compute_CS_for_C3D10_using_dbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational stresses for num_elem elements of typ C3D10.
    Each element has n=10 nodes and ips=4 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CS_at_int_points: Array of shape (num_elem, ips, 3, 3) containing the configurational stresses
        evaluated on ips integration points for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CS_at_int_points = np.zeros((num_elem, 4, 3, 3), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 4)
    assert X_at_nodes.shape == (num_elem, 10, 3)
    assert U_at_nodes.shape == (num_elem, 10, 3)
    assert S_at_int_points.shape == (num_elem, 4, 3, 3)

    lib.compute_CS_for_C3D10_using_dbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CS_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return CS_at_int_points


map_type_and_method_to_CS_function[('C3D10', 'dbf')] = compute_CS_for_C3D10_using_dbf


def compute_CF_for_C3D10_using_dbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational forces for num_elem elements of typ C3D10.
    Each element has n=10 nodes and ips=4 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CF_at_nodes: Array of shape (num_elem, n, 3) containing the configurational forces
        evaluated on n nodes for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CF_at_nodes = np.zeros((num_elem, 10, 3), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 4)
    assert X_at_nodes.shape == (num_elem, 10, 3)
    assert U_at_nodes.shape == (num_elem, 10, 3)
    assert S_at_int_points.shape == (num_elem, 4, 3, 3)

    lib.compute_CF_for_C3D10_using_dbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CF_at_nodes.ctypes.data_as(ctypes.c_void_p)
    )

    return CF_at_nodes


map_type_and_method_to_CF_function[('C3D10', 'dbf')] = compute_CF_for_C3D10_using_dbf

map_type_to_info['C3D10'] = ElementInfo(
    number_of_dimensions=3,
    number_of_nodes=10,
    number_of_integration_points=4
)


def compute_CS_for_C3D10_using_mbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational stresses for num_elem elements of typ C3D10.
    Each element has n=10 nodes and ips=4 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CS_at_int_points: Array of shape (num_elem, ips, 3, 3) containing the configurational stresses
        evaluated on ips integration points for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CS_at_int_points = np.zeros((num_elem, 4, 3, 3), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 4)
    assert X_at_nodes.shape == (num_elem, 10, 3)
    assert U_at_nodes.shape == (num_elem, 10, 3)
    assert S_at_int_points.shape == (num_elem, 4, 3, 3)

    lib.compute_CS_for_C3D10_using_mbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CS_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return CS_at_int_points


map_type_and_method_to_CS_function[('C3D10', 'mbf')] = compute_CS_for_C3D10_using_mbf


def compute_CF_for_C3D10_using_mbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational forces for num_elem elements of typ C3D10.
    Each element has n=10 nodes and ips=4 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CF_at_nodes: Array of shape (num_elem, n, 3) containing the configurational forces
        evaluated on n nodes for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CF_at_nodes = np.zeros((num_elem, 10, 3), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 4)
    assert X_at_nodes.shape == (num_elem, 10, 3)
    assert U_at_nodes.shape == (num_elem, 10, 3)
    assert S_at_int_points.shape == (num_elem, 4, 3, 3)

    lib.compute_CF_for_C3D10_using_mbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CF_at_nodes.ctypes.data_as(ctypes.c_void_p)
    )

    return CF_at_nodes


map_type_and_method_to_CF_function[('C3D10', 'mbf')] = compute_CF_for_C3D10_using_mbf

map_type_to_info['C3D4'] = ElementInfo(
    number_of_dimensions=3,
    number_of_nodes=4,
    number_of_integration_points=1
)


def compute_F_for_C3D4(
        X_at_nodes,
        U_at_nodes):
    """
    Compute the deformation gradients for num_elem elements of typ C3D4.
    Each element has n=4 nodes and ips=1 integration points.

    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :return: F_at_int_points: Array of shape (num_elem, ips, 3, 3) containing the deformation gradients 
        evaluated on ips integration points for num_elem element.
    """

    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    F_at_int_points = np.zeros((num_elem, 1, 3, 3), dtype=np.float64)

    assert X_at_nodes.shape == (num_elem, 4, 3)
    assert U_at_nodes.shape == (num_elem, 4, 3)

    lib.compute_F_for_C3D4(
        ctypes.c_size_t(num_elem),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        F_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return F_at_int_points


map_type_to_F_function['C3D4'] = compute_F_for_C3D4


def compute_P_for_C3D4(
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the first Piola-Kirchhoff stress tensors for num_elem elements of typ C3D4.
    Each element has n=4 nodes and ips=1 integration points.
    
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: P_at_int_points: Array of shape (num_elem, ips, 3, 3) containing the 1. Piola-Kirchhoff stress tensors
        evaluated on ips integration points for num_elem element.
    """

    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    P_at_int_points = np.zeros((num_elem, 1, 3, 3), dtype=np.float64)

    assert X_at_nodes.shape == (num_elem, 4, 3)
    assert U_at_nodes.shape == (num_elem, 4, 3)
    assert S_at_int_points.shape == (num_elem, 1, 3, 3)

    lib.compute_P_for_C3D4(
        ctypes.c_size_t(num_elem),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        P_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return P_at_int_points


map_type_to_P_function['C3D4'] = compute_P_for_C3D4


def compute_CS_for_C3D4_using_dbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational stresses for num_elem elements of typ C3D4.
    Each element has n=4 nodes and ips=1 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CS_at_int_points: Array of shape (num_elem, ips, 3, 3) containing the configurational stresses
        evaluated on ips integration points for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CS_at_int_points = np.zeros((num_elem, 1, 3, 3), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 1)
    assert X_at_nodes.shape == (num_elem, 4, 3)
    assert U_at_nodes.shape == (num_elem, 4, 3)
    assert S_at_int_points.shape == (num_elem, 1, 3, 3)

    lib.compute_CS_for_C3D4_using_dbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CS_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return CS_at_int_points


map_type_and_method_to_CS_function[('C3D4', 'dbf')] = compute_CS_for_C3D4_using_dbf


def compute_CF_for_C3D4_using_dbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational forces for num_elem elements of typ C3D4.
    Each element has n=4 nodes and ips=1 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CF_at_nodes: Array of shape (num_elem, n, 3) containing the configurational forces
        evaluated on n nodes for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CF_at_nodes = np.zeros((num_elem, 4, 3), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 1)
    assert X_at_nodes.shape == (num_elem, 4, 3)
    assert U_at_nodes.shape == (num_elem, 4, 3)
    assert S_at_int_points.shape == (num_elem, 1, 3, 3)

    lib.compute_CF_for_C3D4_using_dbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CF_at_nodes.ctypes.data_as(ctypes.c_void_p)
    )

    return CF_at_nodes


map_type_and_method_to_CF_function[('C3D4', 'dbf')] = compute_CF_for_C3D4_using_dbf

map_type_to_info['C3D4'] = ElementInfo(
    number_of_dimensions=3,
    number_of_nodes=4,
    number_of_integration_points=1
)


def compute_CS_for_C3D4_using_mbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational stresses for num_elem elements of typ C3D4.
    Each element has n=4 nodes and ips=1 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CS_at_int_points: Array of shape (num_elem, ips, 3, 3) containing the configurational stresses
        evaluated on ips integration points for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CS_at_int_points = np.zeros((num_elem, 1, 3, 3), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 1)
    assert X_at_nodes.shape == (num_elem, 4, 3)
    assert U_at_nodes.shape == (num_elem, 4, 3)
    assert S_at_int_points.shape == (num_elem, 1, 3, 3)

    lib.compute_CS_for_C3D4_using_mbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CS_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return CS_at_int_points


map_type_and_method_to_CS_function[('C3D4', 'mbf')] = compute_CS_for_C3D4_using_mbf


def compute_CF_for_C3D4_using_mbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational forces for num_elem elements of typ C3D4.
    Each element has n=4 nodes and ips=1 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CF_at_nodes: Array of shape (num_elem, n, 3) containing the configurational forces
        evaluated on n nodes for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CF_at_nodes = np.zeros((num_elem, 4, 3), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 1)
    assert X_at_nodes.shape == (num_elem, 4, 3)
    assert U_at_nodes.shape == (num_elem, 4, 3)
    assert S_at_int_points.shape == (num_elem, 1, 3, 3)

    lib.compute_CF_for_C3D4_using_mbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CF_at_nodes.ctypes.data_as(ctypes.c_void_p)
    )

    return CF_at_nodes


map_type_and_method_to_CF_function[('C3D4', 'mbf')] = compute_CF_for_C3D4_using_mbf

map_type_to_info['C3D15'] = ElementInfo(
    number_of_dimensions=3,
    number_of_nodes=15,
    number_of_integration_points=9
)


def compute_F_for_C3D15(
        X_at_nodes,
        U_at_nodes):
    """
    Compute the deformation gradients for num_elem elements of typ C3D15.
    Each element has n=15 nodes and ips=9 integration points.

    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :return: F_at_int_points: Array of shape (num_elem, ips, 3, 3) containing the deformation gradients 
        evaluated on ips integration points for num_elem element.
    """

    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    F_at_int_points = np.zeros((num_elem, 9, 3, 3), dtype=np.float64)

    assert X_at_nodes.shape == (num_elem, 15, 3)
    assert U_at_nodes.shape == (num_elem, 15, 3)

    lib.compute_F_for_C3D15(
        ctypes.c_size_t(num_elem),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        F_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return F_at_int_points


map_type_to_F_function['C3D15'] = compute_F_for_C3D15


def compute_P_for_C3D15(
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the first Piola-Kirchhoff stress tensors for num_elem elements of typ C3D15.
    Each element has n=15 nodes and ips=9 integration points.
    
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: P_at_int_points: Array of shape (num_elem, ips, 3, 3) containing the 1. Piola-Kirchhoff stress tensors
        evaluated on ips integration points for num_elem element.
    """

    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    P_at_int_points = np.zeros((num_elem, 9, 3, 3), dtype=np.float64)

    assert X_at_nodes.shape == (num_elem, 15, 3)
    assert U_at_nodes.shape == (num_elem, 15, 3)
    assert S_at_int_points.shape == (num_elem, 9, 3, 3)

    lib.compute_P_for_C3D15(
        ctypes.c_size_t(num_elem),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        P_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return P_at_int_points


map_type_to_P_function['C3D15'] = compute_P_for_C3D15


def compute_CS_for_C3D15_using_dbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational stresses for num_elem elements of typ C3D15.
    Each element has n=15 nodes and ips=9 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CS_at_int_points: Array of shape (num_elem, ips, 3, 3) containing the configurational stresses
        evaluated on ips integration points for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CS_at_int_points = np.zeros((num_elem, 9, 3, 3), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 9)
    assert X_at_nodes.shape == (num_elem, 15, 3)
    assert U_at_nodes.shape == (num_elem, 15, 3)
    assert S_at_int_points.shape == (num_elem, 9, 3, 3)

    lib.compute_CS_for_C3D15_using_dbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CS_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return CS_at_int_points


map_type_and_method_to_CS_function[('C3D15', 'dbf')] = compute_CS_for_C3D15_using_dbf


def compute_CF_for_C3D15_using_dbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational forces for num_elem elements of typ C3D15.
    Each element has n=15 nodes and ips=9 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CF_at_nodes: Array of shape (num_elem, n, 3) containing the configurational forces
        evaluated on n nodes for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CF_at_nodes = np.zeros((num_elem, 15, 3), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 9)
    assert X_at_nodes.shape == (num_elem, 15, 3)
    assert U_at_nodes.shape == (num_elem, 15, 3)
    assert S_at_int_points.shape == (num_elem, 9, 3, 3)

    lib.compute_CF_for_C3D15_using_dbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CF_at_nodes.ctypes.data_as(ctypes.c_void_p)
    )

    return CF_at_nodes


map_type_and_method_to_CF_function[('C3D15', 'dbf')] = compute_CF_for_C3D15_using_dbf

map_type_to_info['C3D15'] = ElementInfo(
    number_of_dimensions=3,
    number_of_nodes=15,
    number_of_integration_points=9
)


def compute_CS_for_C3D15_using_mbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational stresses for num_elem elements of typ C3D15.
    Each element has n=15 nodes and ips=9 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CS_at_int_points: Array of shape (num_elem, ips, 3, 3) containing the configurational stresses
        evaluated on ips integration points for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CS_at_int_points = np.zeros((num_elem, 9, 3, 3), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 9)
    assert X_at_nodes.shape == (num_elem, 15, 3)
    assert U_at_nodes.shape == (num_elem, 15, 3)
    assert S_at_int_points.shape == (num_elem, 9, 3, 3)

    lib.compute_CS_for_C3D15_using_mbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CS_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return CS_at_int_points


map_type_and_method_to_CS_function[('C3D15', 'mbf')] = compute_CS_for_C3D15_using_mbf


def compute_CF_for_C3D15_using_mbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational forces for num_elem elements of typ C3D15.
    Each element has n=15 nodes and ips=9 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CF_at_nodes: Array of shape (num_elem, n, 3) containing the configurational forces
        evaluated on n nodes for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CF_at_nodes = np.zeros((num_elem, 15, 3), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 9)
    assert X_at_nodes.shape == (num_elem, 15, 3)
    assert U_at_nodes.shape == (num_elem, 15, 3)
    assert S_at_int_points.shape == (num_elem, 9, 3, 3)

    lib.compute_CF_for_C3D15_using_mbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CF_at_nodes.ctypes.data_as(ctypes.c_void_p)
    )

    return CF_at_nodes


map_type_and_method_to_CF_function[('C3D15', 'mbf')] = compute_CF_for_C3D15_using_mbf

map_type_to_info['CPE6'] = ElementInfo(
    number_of_dimensions=2,
    number_of_nodes=6,
    number_of_integration_points=3
)


def compute_F_for_CPE6(
        X_at_nodes,
        U_at_nodes):
    """
    Compute the deformation gradients for num_elem elements of typ CPE6.
    Each element has n=6 nodes and ips=3 integration points.

    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :return: F_at_int_points: Array of shape (num_elem, ips, 2, 2) containing the deformation gradients 
        evaluated on ips integration points for num_elem element.
    """

    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    F_at_int_points = np.zeros((num_elem, 3, 2, 2), dtype=np.float64)

    assert X_at_nodes.shape == (num_elem, 6, 2)
    assert U_at_nodes.shape == (num_elem, 6, 2)

    lib.compute_F_for_CPE6(
        ctypes.c_size_t(num_elem),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        F_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return F_at_int_points


map_type_to_F_function['CPE6'] = compute_F_for_CPE6


def compute_P_for_CPE6(
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the first Piola-Kirchhoff stress tensors for num_elem elements of typ CPE6.
    Each element has n=6 nodes and ips=3 integration points.
    
    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 2, 2) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: P_at_int_points: Array of shape (num_elem, ips, 2, 2) containing the 1. Piola-Kirchhoff stress tensors
        evaluated on ips integration points for num_elem element.
    """

    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    P_at_int_points = np.zeros((num_elem, 3, 2, 2), dtype=np.float64)

    assert X_at_nodes.shape == (num_elem, 6, 2)
    assert U_at_nodes.shape == (num_elem, 6, 2)
    assert S_at_int_points.shape == (num_elem, 3, 2, 2)

    lib.compute_P_for_CPE6(
        ctypes.c_size_t(num_elem),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        P_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return P_at_int_points


map_type_to_P_function['CPE6'] = compute_P_for_CPE6


def compute_CS_for_CPE6_using_dbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational stresses for num_elem elements of typ CPE6.
    Each element has n=6 nodes and ips=3 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 2, 2) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CS_at_int_points: Array of shape (num_elem, ips, 2, 2) containing the configurational stresses
        evaluated on ips integration points for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CS_at_int_points = np.zeros((num_elem, 3, 2, 2), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 3)
    assert X_at_nodes.shape == (num_elem, 6, 2)
    assert U_at_nodes.shape == (num_elem, 6, 2)
    assert S_at_int_points.shape == (num_elem, 3, 2, 2)

    lib.compute_CS_for_CPE6_using_dbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CS_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return CS_at_int_points


map_type_and_method_to_CS_function[('CPE6', 'dbf')] = compute_CS_for_CPE6_using_dbf


def compute_CF_for_CPE6_using_dbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational forces for num_elem elements of typ CPE6.
    Each element has n=6 nodes and ips=3 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 2, 2) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CF_at_nodes: Array of shape (num_elem, n, 2) containing the configurational forces
        evaluated on n nodes for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CF_at_nodes = np.zeros((num_elem, 6, 2), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 3)
    assert X_at_nodes.shape == (num_elem, 6, 2)
    assert U_at_nodes.shape == (num_elem, 6, 2)
    assert S_at_int_points.shape == (num_elem, 3, 2, 2)

    lib.compute_CF_for_CPE6_using_dbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CF_at_nodes.ctypes.data_as(ctypes.c_void_p)
    )

    return CF_at_nodes


map_type_and_method_to_CF_function[('CPE6', 'dbf')] = compute_CF_for_CPE6_using_dbf

map_type_to_info['CPE6'] = ElementInfo(
    number_of_dimensions=2,
    number_of_nodes=6,
    number_of_integration_points=3
)


def compute_CS_for_CPE6_using_mbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational stresses for num_elem elements of typ CPE6.
    Each element has n=6 nodes and ips=3 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 2, 2) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CS_at_int_points: Array of shape (num_elem, ips, 2, 2) containing the configurational stresses
        evaluated on ips integration points for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CS_at_int_points = np.zeros((num_elem, 3, 2, 2), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 3)
    assert X_at_nodes.shape == (num_elem, 6, 2)
    assert U_at_nodes.shape == (num_elem, 6, 2)
    assert S_at_int_points.shape == (num_elem, 3, 2, 2)

    lib.compute_CS_for_CPE6_using_mbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CS_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return CS_at_int_points


map_type_and_method_to_CS_function[('CPE6', 'mbf')] = compute_CS_for_CPE6_using_mbf


def compute_CF_for_CPE6_using_mbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational forces for num_elem elements of typ CPE6.
    Each element has n=6 nodes and ips=3 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 2, 2) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CF_at_nodes: Array of shape (num_elem, n, 2) containing the configurational forces
        evaluated on n nodes for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CF_at_nodes = np.zeros((num_elem, 6, 2), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 3)
    assert X_at_nodes.shape == (num_elem, 6, 2)
    assert U_at_nodes.shape == (num_elem, 6, 2)
    assert S_at_int_points.shape == (num_elem, 3, 2, 2)

    lib.compute_CF_for_CPE6_using_mbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CF_at_nodes.ctypes.data_as(ctypes.c_void_p)
    )

    return CF_at_nodes


map_type_and_method_to_CF_function[('CPE6', 'mbf')] = compute_CF_for_CPE6_using_mbf

map_type_to_info['C3D6'] = ElementInfo(
    number_of_dimensions=3,
    number_of_nodes=6,
    number_of_integration_points=2
)


def compute_F_for_C3D6(
        X_at_nodes,
        U_at_nodes):
    """
    Compute the deformation gradients for num_elem elements of typ C3D6.
    Each element has n=6 nodes and ips=2 integration points.

    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :return: F_at_int_points: Array of shape (num_elem, ips, 3, 3) containing the deformation gradients 
        evaluated on ips integration points for num_elem element.
    """

    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    F_at_int_points = np.zeros((num_elem, 2, 3, 3), dtype=np.float64)

    assert X_at_nodes.shape == (num_elem, 6, 3)
    assert U_at_nodes.shape == (num_elem, 6, 3)

    lib.compute_F_for_C3D6(
        ctypes.c_size_t(num_elem),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        F_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return F_at_int_points


map_type_to_F_function['C3D6'] = compute_F_for_C3D6


def compute_P_for_C3D6(
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the first Piola-Kirchhoff stress tensors for num_elem elements of typ C3D6.
    Each element has n=6 nodes and ips=2 integration points.
    
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: P_at_int_points: Array of shape (num_elem, ips, 3, 3) containing the 1. Piola-Kirchhoff stress tensors
        evaluated on ips integration points for num_elem element.
    """

    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    P_at_int_points = np.zeros((num_elem, 2, 3, 3), dtype=np.float64)

    assert X_at_nodes.shape == (num_elem, 6, 3)
    assert U_at_nodes.shape == (num_elem, 6, 3)
    assert S_at_int_points.shape == (num_elem, 2, 3, 3)

    lib.compute_P_for_C3D6(
        ctypes.c_size_t(num_elem),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        P_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return P_at_int_points


map_type_to_P_function['C3D6'] = compute_P_for_C3D6


def compute_CS_for_C3D6_using_dbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational stresses for num_elem elements of typ C3D6.
    Each element has n=6 nodes and ips=2 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CS_at_int_points: Array of shape (num_elem, ips, 3, 3) containing the configurational stresses
        evaluated on ips integration points for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CS_at_int_points = np.zeros((num_elem, 2, 3, 3), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 2)
    assert X_at_nodes.shape == (num_elem, 6, 3)
    assert U_at_nodes.shape == (num_elem, 6, 3)
    assert S_at_int_points.shape == (num_elem, 2, 3, 3)

    lib.compute_CS_for_C3D6_using_dbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CS_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return CS_at_int_points


map_type_and_method_to_CS_function[('C3D6', 'dbf')] = compute_CS_for_C3D6_using_dbf


def compute_CF_for_C3D6_using_dbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational forces for num_elem elements of typ C3D6.
    Each element has n=6 nodes and ips=2 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CF_at_nodes: Array of shape (num_elem, n, 3) containing the configurational forces
        evaluated on n nodes for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CF_at_nodes = np.zeros((num_elem, 6, 3), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 2)
    assert X_at_nodes.shape == (num_elem, 6, 3)
    assert U_at_nodes.shape == (num_elem, 6, 3)
    assert S_at_int_points.shape == (num_elem, 2, 3, 3)

    lib.compute_CF_for_C3D6_using_dbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CF_at_nodes.ctypes.data_as(ctypes.c_void_p)
    )

    return CF_at_nodes


map_type_and_method_to_CF_function[('C3D6', 'dbf')] = compute_CF_for_C3D6_using_dbf

map_type_to_info['C3D6'] = ElementInfo(
    number_of_dimensions=3,
    number_of_nodes=6,
    number_of_integration_points=2
)


def compute_CS_for_C3D6_using_mbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational stresses for num_elem elements of typ C3D6.
    Each element has n=6 nodes and ips=2 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CS_at_int_points: Array of shape (num_elem, ips, 3, 3) containing the configurational stresses
        evaluated on ips integration points for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CS_at_int_points = np.zeros((num_elem, 2, 3, 3), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 2)
    assert X_at_nodes.shape == (num_elem, 6, 3)
    assert U_at_nodes.shape == (num_elem, 6, 3)
    assert S_at_int_points.shape == (num_elem, 2, 3, 3)

    lib.compute_CS_for_C3D6_using_mbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CS_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return CS_at_int_points


map_type_and_method_to_CS_function[('C3D6', 'mbf')] = compute_CS_for_C3D6_using_mbf


def compute_CF_for_C3D6_using_mbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational forces for num_elem elements of typ C3D6.
    Each element has n=6 nodes and ips=2 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 3) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 3) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 3, 3) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CF_at_nodes: Array of shape (num_elem, n, 3) containing the configurational forces
        evaluated on n nodes for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CF_at_nodes = np.zeros((num_elem, 6, 3), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 2)
    assert X_at_nodes.shape == (num_elem, 6, 3)
    assert U_at_nodes.shape == (num_elem, 6, 3)
    assert S_at_int_points.shape == (num_elem, 2, 3, 3)

    lib.compute_CF_for_C3D6_using_mbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CF_at_nodes.ctypes.data_as(ctypes.c_void_p)
    )

    return CF_at_nodes


map_type_and_method_to_CF_function[('C3D6', 'mbf')] = compute_CF_for_C3D6_using_mbf

map_type_to_info['CPE3'] = ElementInfo(
    number_of_dimensions=2,
    number_of_nodes=3,
    number_of_integration_points=1
)


def compute_F_for_CPE3(
        X_at_nodes,
        U_at_nodes):
    """
    Compute the deformation gradients for num_elem elements of typ CPE3.
    Each element has n=3 nodes and ips=1 integration points.

    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :return: F_at_int_points: Array of shape (num_elem, ips, 2, 2) containing the deformation gradients 
        evaluated on ips integration points for num_elem element.
    """

    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    F_at_int_points = np.zeros((num_elem, 1, 2, 2), dtype=np.float64)

    assert X_at_nodes.shape == (num_elem, 3, 2)
    assert U_at_nodes.shape == (num_elem, 3, 2)

    lib.compute_F_for_CPE3(
        ctypes.c_size_t(num_elem),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        F_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return F_at_int_points


map_type_to_F_function['CPE3'] = compute_F_for_CPE3


def compute_P_for_CPE3(
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the first Piola-Kirchhoff stress tensors for num_elem elements of typ CPE3.
    Each element has n=3 nodes and ips=1 integration points.
    
    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 2, 2) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: P_at_int_points: Array of shape (num_elem, ips, 2, 2) containing the 1. Piola-Kirchhoff stress tensors
        evaluated on ips integration points for num_elem element.
    """

    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    P_at_int_points = np.zeros((num_elem, 1, 2, 2), dtype=np.float64)

    assert X_at_nodes.shape == (num_elem, 3, 2)
    assert U_at_nodes.shape == (num_elem, 3, 2)
    assert S_at_int_points.shape == (num_elem, 1, 2, 2)

    lib.compute_P_for_CPE3(
        ctypes.c_size_t(num_elem),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        P_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return P_at_int_points


map_type_to_P_function['CPE3'] = compute_P_for_CPE3


def compute_CS_for_CPE3_using_dbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational stresses for num_elem elements of typ CPE3.
    Each element has n=3 nodes and ips=1 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 2, 2) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CS_at_int_points: Array of shape (num_elem, ips, 2, 2) containing the configurational stresses
        evaluated on ips integration points for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CS_at_int_points = np.zeros((num_elem, 1, 2, 2), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 1)
    assert X_at_nodes.shape == (num_elem, 3, 2)
    assert U_at_nodes.shape == (num_elem, 3, 2)
    assert S_at_int_points.shape == (num_elem, 1, 2, 2)

    lib.compute_CS_for_CPE3_using_dbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CS_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return CS_at_int_points


map_type_and_method_to_CS_function[('CPE3', 'dbf')] = compute_CS_for_CPE3_using_dbf


def compute_CF_for_CPE3_using_dbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational forces for num_elem elements of typ CPE3.
    Each element has n=3 nodes and ips=1 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 2, 2) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CF_at_nodes: Array of shape (num_elem, n, 2) containing the configurational forces
        evaluated on n nodes for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CF_at_nodes = np.zeros((num_elem, 3, 2), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 1)
    assert X_at_nodes.shape == (num_elem, 3, 2)
    assert U_at_nodes.shape == (num_elem, 3, 2)
    assert S_at_int_points.shape == (num_elem, 1, 2, 2)

    lib.compute_CF_for_CPE3_using_dbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CF_at_nodes.ctypes.data_as(ctypes.c_void_p)
    )

    return CF_at_nodes


map_type_and_method_to_CF_function[('CPE3', 'dbf')] = compute_CF_for_CPE3_using_dbf

map_type_to_info['CPE3'] = ElementInfo(
    number_of_dimensions=2,
    number_of_nodes=3,
    number_of_integration_points=1
)


def compute_CS_for_CPE3_using_mbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational stresses for num_elem elements of typ CPE3.
    Each element has n=3 nodes and ips=1 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 2, 2) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CS_at_int_points: Array of shape (num_elem, ips, 2, 2) containing the configurational stresses
        evaluated on ips integration points for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CS_at_int_points = np.zeros((num_elem, 1, 2, 2), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 1)
    assert X_at_nodes.shape == (num_elem, 3, 2)
    assert U_at_nodes.shape == (num_elem, 3, 2)
    assert S_at_int_points.shape == (num_elem, 1, 2, 2)

    lib.compute_CS_for_CPE3_using_mbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CS_at_int_points.ctypes.data_as(ctypes.c_void_p)
    )

    return CS_at_int_points


map_type_and_method_to_CS_function[('CPE3', 'mbf')] = compute_CS_for_CPE3_using_mbf


def compute_CF_for_CPE3_using_mbf(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points):
    """
    Compute the configurational forces for num_elem elements of typ CPE3.
    Each element has n=3 nodes and ips=1 integration points.
    
    :param e_at_int_points: Array of shape (num_elem, ips) containing the internal energy densities of num_elem elements.
    :param X_at_nodes: Array of shape (num_elem, n, 2) containing the coordinates at n nodes of num_elem elements.
    :param U_at_nodes: Array of shape (num_elem, n, 2) containing the displacements at n nodes of num_elem elements.
    :param S_at_int_points: Array of shape (num_elem, ips, 2, 2) 
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: CF_at_nodes: Array of shape (num_elem, n, 2) containing the configurational forces
        evaluated on n nodes for num_elem element.
    """

    e_at_int_points = np.ascontiguousarray(e_at_int_points, dtype=np.float64)
    X_at_nodes = np.ascontiguousarray(X_at_nodes, dtype=np.float64)
    U_at_nodes = np.ascontiguousarray(U_at_nodes, dtype=np.float64)
    S_at_int_points = np.ascontiguousarray(S_at_int_points, dtype=np.float64)

    num_elem = X_at_nodes.shape[0]

    CF_at_nodes = np.zeros((num_elem, 3, 2), dtype=np.float64)

    assert e_at_int_points.shape == (num_elem, 1)
    assert X_at_nodes.shape == (num_elem, 3, 2)
    assert U_at_nodes.shape == (num_elem, 3, 2)
    assert S_at_int_points.shape == (num_elem, 1, 2, 2)

    lib.compute_CF_for_CPE3_using_mbf(
        ctypes.c_size_t(num_elem),
        e_at_int_points.ctypes.data_as(ctypes.c_void_p),
        X_at_nodes.ctypes.data_as(ctypes.c_void_p),
        U_at_nodes.ctypes.data_as(ctypes.c_void_p),
        S_at_int_points.ctypes.data_as(ctypes.c_void_p),
        CF_at_nodes.ctypes.data_as(ctypes.c_void_p)
    )

    return CF_at_nodes


map_type_and_method_to_CF_function[('CPE3', 'mbf')] = compute_CF_for_CPE3_using_mbf
