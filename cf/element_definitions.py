r"""
The elements are defined in a reference coordinate system.
For the application, they have to be transformed to the real coordinate system.
This module defines elements by four dictionaries:

 - :py:data:`integration_points`: Integration points in the reference coordinate system
 - :py:data:`integration_weights`: Integration weights corresponding to the integration points
 - :py:data:`reference_nodes`: Coordinates of the element nodes in the reference coordinate system
 - :py:data:`shape_function_exponents`: Exponents of polynomial powers used inside the shape functions


ACCESS TO DICTIONARIES
----------------------

To acces the nodes, exponents, etc. for a specific element,
the element typ can be looked up in the dictionaries.
For example,

>>> typ = "CPE4R"

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


Shape functions in reference space
----------------------------------

>>> import sympy as sy
>>> from cf.expressions import eval_H, eval_dH_dR, eval_J, eval_dH_dX

The shape functions are defined in the reference coordinate system.
For a 2D element this coordinate system has 2 variables:

>>> r, s = sy.symbols("r s", real=True)

A shape function consists out of coefficients

>>> C = sy.MatrixSymbol("C", len(shape_exponents), 1).as_explicit()

and powers of the reference space variables

>>> P = sy.zeros(len(shape_exponents), 1)
>>> for i in range(len(shape_exponents)):
...     P[i, 0] = r**shape_exponents[i, 0] * s**shape_exponents[i, 1]
>>> P
Matrix([
[  1],
[  r],
[  s],
[r*s]])

The shape function is a matrix multiplication of the coefficients and the powers

>>> h0 = (C.T * P)[0]
>>> h0
r*s*C[3, 0] + r*C[1, 0] + s*C[2, 0] + C[0, 0]

The coefficients are computed out of the coordinates of the element nodes
and the powers of the shape functions.

>>> H = eval_H(sy.Matrix([r, s]), ref_nodes, shape_exponents)
>>> H[0]
0.25*r*s - 0.25*r - 0.25*s + 0.25

Each shape function corresponds to one node.
At this node the shape function is one, whereas at other nodes
it is zero.
E.g. The i-th shape function is one at the i-th node.

>>> i = 1
>>> node_i = ref_nodes[i]
>>> float(H.xreplace(dict(zip([r, s], node_i)))[i])
1.0

but is zero at any j-th node (j != i).

>>> j = 0
>>> node_j = ref_nodes[j]
>>> float(H.xreplace(dict(zip([r, s], node_j)))[i])
0.0


Interpolation of nodal values
-----------------------------

Shape functions are used to interpolate quantities inside the element.
Values for a quantity are provided at the element nodes.
E.g. The coordinates are provided at the element nodes.

>>> quantity = sy.Matrix(ref_nodes)

The interpolation is done by a matrix multiplication with the shape functions.

>>> interpolation = quantity.T * H
>>> interpolation
Matrix([
[1.0*r],
[1.0*s]])

In the simple example the interpolated quantities are :math:`(1*r, 1*s)`.


Derivative of interpolation in reference space
----------------------------------------------

To derive the interpolation of a quantity,
it is sufficient to derive the shape function
with respect to the reference space coordinates.

>>> dH_dR = eval_dH_dR(H, sy.Matrix([r, s]))

The derivative of the quantities with respect to
the reference space coordinates is again a
matrix multiplication.

>>> dQ_dR = quantity.T * dH_dR

E.g. the i-th quantity derived by the j-th reference space coordinate is

>>> i, j = (0, 0)
>>> dQ_dR[i, j]
1.00000000000000


Integration in the reference space
----------------------------------

Since, the shape functions are multidimensional polynomials, they can be
integrated using Gaussian quadrature.
E.g. the function

>>> f = r*s + 1

is integrated over the element volume in the reference space, by

>>> integral = 0
>>> for int_point, int_weight in zip(int_points, int_weights):
...     integral += int_weight * float(f.xreplace(dict(zip([r, s], int_point))))
>>> integral
4.0

However, reduced elements like CPE4R have less integration points
than it would be necessary to integrate the shape function exactly.


Real space
----------

For an application the element is transformed to the real space.
E.g. the element nodes have now real space coordinates

>>> real_nodes = np.array([
...     [1.0, 2.0],
...     [3.0, 3.0],
...     [4.0, 5.0],
...     [1.0, 4.0]
... ])

The coordiantes in the real space are interpolated
from the real space coordinates of the element nodes.

>>> real_nodes.T @ H
Matrix([
[0.25*r*s + 1.25*r + 0.25*s + 2.25],
[              0.5*r + 1.0*s + 3.5]])

For the space transformation the jacobian is computed.

>>> J = eval_J(real_nodes, dH_dR)

The jacobian derives the real space coordinates
with respect to the reference space coordinates.
The i-th real space coordinate derived by the
j-th reference space coordinate is

>>> i, j = (1, 0)
>>> J[i, j]
0.500000000000000


Derivative of interpolation in real space
-----------------------------------------

Imagine, displacements are given in the real space

>>> strain_tensor = np.array([
...     [0.3, 0.0],
...     [0.0, -0.1]
... ])
>>> real_U = real_nodes @ strain_tensor
>>> real_U
array([[ 0.3, -0.2],
       [ 0.9, -0.3],
       [ 1.2, -0.5],
       [ 0.3, -0.4]])

They can be interpolated using the shape functions

>>> U = real_U.T * H

However, at the moment you can only derive them with
respect to the reference space.

To derive an interpolation like `U` with respect to the
real space, the jacobian is used to transfer `dH_dR` to

>>> dH_dX = eval_dH_dX(dH_dR, J)

Now the displacements can be derived with respect to the
real space.

>>> dU_dX = real_U.T * dH_dX

The derivative is evaluated at the center of element and returns the
strain tensor.

>>> np.array(dU_dX.subs({r: 0., s: 0.}), dtype=float).round(7)
array([[ 0.3,  0. ],
       [-0. , -0.1]])



Integration in the real space
-----------------------------

The jacobian is also used for the integration
in the real space.
E.g. the strain energy density is computed out of strain and stress tensors.

>>> stress_tensor = 2 * 400 * strain_tensor + 600 * np.trace(strain_tensor) * np.eye(2)
>>> stress_tensor
array([[360.,   0.],
       [  0.,  40.]])

>>> se_density = 0.5 * np.tensordot(strain_tensor, stress_tensor)
>>> se_density
52.0


The strain energy is the integral of the strain energy density
over the element in the real space coordinates.
In addition to the integral in the reference space,
the determinat of the jacobian is added.

>>> integral = 0.
>>> for int_point, int_weight in zip(int_points, int_weights):
...     integral += (
...         int_weight
...         * se_density
...         * float(J.det().xreplace(dict(zip([r, s], int_point))))
...     )
>>> integral
234.0

"""
import numpy as np


R_at_integration_points_of_element = {}
"""coordinates of integration points in the reference coordinates system"""

weights_of_integration_points_of_element = {}
"""weights corresponding to the integration points"""

R_at_nodes_of_element = {}
"""coordinates of the element nodes in the reference coordinate system"""

exponents_of_shape_functions_of_element = {}
"""exponents of powers used by the shape functions defined in the reference space"""

# Calculation of gaussian weights
EdW3 = 1. / np.sqrt(3)
Ed3 = 1. / 3.
W3d5 = np.sqrt(3. / 5.)
tri41 = 1. / 4. + 3 * np.sqrt(5) / 20
tri42 = 1. / 4. - np.sqrt(5) / 20
tri31 = 1. / 6.
tri32 = 4. / 6.
B5d9 = 5. / 9.
B8d9 = 8. / 9.


###################################################################################
# 2D - elements ###################################################################
###################################################################################
# CPE4
R_at_integration_points_of_element['CPE4'] = np.array((
    (-EdW3, -EdW3), (EdW3, -EdW3),
    (-EdW3, EdW3), (EdW3, EdW3)
))
weights_of_integration_points_of_element['CPE4'] = np.array((1., 1., 1., 1.))
R_at_nodes_of_element['CPE4'] = np.array((
    (-1., -1.), (1., -1.),
    (1., 1.), (-1., 1.)
))
exponents_of_shape_functions_of_element['CPE4'] = np.array((
    (0, 0), (1, 0),
    (0, 1), (1, 1)
))

# CPE4R
R_at_integration_points_of_element['CPE4R'] = np.array(((0., 0.),))
weights_of_integration_points_of_element['CPE4R'] = np.array((4.,))
R_at_nodes_of_element['CPE4R'] = R_at_nodes_of_element['CPE4']
exponents_of_shape_functions_of_element['CPE4R'] = exponents_of_shape_functions_of_element['CPE4']

# CPE8
R_at_integration_points_of_element['CPE8'] = np.array((
    (-W3d5, -W3d5), (0., -W3d5), (W3d5, -W3d5),
    (-W3d5, 0.), (0., 0.), (W3d5, 0.),
    (-W3d5, W3d5), (0., W3d5), (W3d5, W3d5)
))
weights_of_integration_points_of_element['CPE8'] = np.array((
    B5d9 * B5d9, B8d9 * B5d9, B5d9 * B5d9,
    B8d9 * B5d9, B8d9 * B8d9, B8d9 * B5d9,
    B5d9 * B5d9, B8d9 * B5d9, B5d9 * B5d9
))
R_at_nodes_of_element['CPE8'] = np.array((
    (-1., -1.), (1., -1.), (1., 1.), (-1., 1.),
    (0., -1.), (1., 0.), (0., 1.), (-1., 0.)
))
exponents_of_shape_functions_of_element['CPE8'] = np.array((
    (0, 0), (1, 0), (0, 1), (1, 1),
    (2, 0), (0, 2), (2, 1), (1, 2)
))

# CPE8R
R_at_integration_points_of_element['CPE8R'] = np.array((
    (-EdW3, -EdW3), (EdW3, -EdW3),
    (-EdW3, EdW3), (EdW3, EdW3)
))
weights_of_integration_points_of_element['CPE8R'] = np.array((1., 1., 1., 1.))
R_at_nodes_of_element['CPE8R'] = R_at_nodes_of_element['CPE8']
exponents_of_shape_functions_of_element['CPE8R'] = exponents_of_shape_functions_of_element['CPE8']

# CPE3
R_at_integration_points_of_element['CPE3'] = np.array(((Ed3, Ed3),))
weights_of_integration_points_of_element['CPE3'] = np.array([0.5])
R_at_nodes_of_element['CPE3'] = np.array((
    (0., 0.),
    (1., 0.),
    (0., 1.)
))
exponents_of_shape_functions_of_element['CPE3'] = np.array((
    (0, 0),
    (1, 0),
    (0, 1)
))

# CPE6
R_at_integration_points_of_element['CPE6'] = np.array((
    (tri31, tri31),
    (tri32, tri31),
    (tri31, tri32)
))
weights_of_integration_points_of_element['CPE6'] = np.array((1. / 6, 1. / 6, 1. / 6))
R_at_nodes_of_element['CPE6'] = np.array((
    (0., 0.), (1., 0.), (0., 1.),
    (0.5, 0.), (0.5, 0.5), (0., 0.5)
))
exponents_of_shape_functions_of_element['CPE6'] = np.array((
    (0, 0),
    (1, 0), (0, 1), (1, 1),
    (0, 2), (2, 0)
))

###################################################################################
# 3D - elements ###################################################################
###################################################################################

# C3D8
R_at_integration_points_of_element['C3D8'] = np.array((
    (-EdW3, -EdW3, -EdW3),
    (EdW3, -EdW3, -EdW3),
    (-EdW3, EdW3, -EdW3),
    (EdW3, EdW3, -EdW3),
    (-EdW3, -EdW3, EdW3),
    (EdW3, -EdW3, EdW3),
    (-EdW3, EdW3, EdW3),
    (EdW3, EdW3, EdW3)
))
weights_of_integration_points_of_element['C3D8'] = np.array((1., 1., 1., 1., 1., 1., 1., 1.))
R_at_nodes_of_element['C3D8'] = np.array((
    (-1., -1., -1.),
    (1., -1., -1.),
    (1., 1., -1.),
    (-1., 1., -1.),
    (-1., -1., 1.),
    (1., -1., 1.),
    (1., 1., 1.),
    (-1., 1., 1.)
))
exponents_of_shape_functions_of_element['C3D8'] = np.array((
    (0, 0, 0),
    (1, 0, 0), (0, 1, 0), (0, 0, 1),
    (1, 1, 0), (1, 0, 1), (0, 1, 1),
    (1, 1, 1)
))

# C3D8R
R_at_integration_points_of_element['C3D8R'] = np.array(((0., 0., 0.),))
weights_of_integration_points_of_element['C3D8R'] = np.zeros(1) + 8

R_at_nodes_of_element['C3D8R'] = R_at_nodes_of_element['C3D8']
exponents_of_shape_functions_of_element['C3D8R'] = exponents_of_shape_functions_of_element['C3D8']

# C3D20
R_at_integration_points_of_element['C3D20'] = np.array((
    (-W3d5, -W3d5, -W3d5),
    (0., -W3d5, -W3d5),
    (W3d5, -W3d5, -W3d5),
    (-W3d5, 0., -W3d5),
    (0., 0., -W3d5),
    (W3d5, 0., -W3d5),
    (-W3d5, W3d5, -W3d5),
    (0., W3d5, -W3d5),
    (W3d5, W3d5, -W3d5),
    (-W3d5, -W3d5, 0.),
    (0., -W3d5, 0.),
    (W3d5, -W3d5, 0.),
    (-W3d5, 0., 0.),
    (0., 0., 0.),
    (W3d5, 0., 0.),
    (-W3d5, W3d5, 0.),
    (0., W3d5, 0.),
    (W3d5, W3d5, 0.),
    (-W3d5, -W3d5, W3d5),
    (0., -W3d5, W3d5),
    (W3d5, -W3d5, W3d5),
    (-W3d5, 0., W3d5),
    (0., 0., W3d5),
    (W3d5, 0., W3d5),
    (-W3d5, W3d5, W3d5),
    (0., W3d5, W3d5),
    (W3d5, W3d5, W3d5)
))
weights_of_integration_points_of_element['C3D20'] = np.array((
    B5d9 * B5d9 * B5d9,
    B8d9 * B5d9 * B5d9,
    B5d9 * B5d9 * B5d9,
    B5d9 * B8d9 * B5d9,
    B8d9 * B8d9 * B5d9,
    B5d9 * B8d9 * B5d9,
    B5d9 * B5d9 * B5d9,
    B8d9 * B5d9 * B5d9,
    B5d9 * B5d9 * B5d9,
    B5d9 * B5d9 * B8d9,
    B8d9 * B5d9 * B8d9,
    B5d9 * B5d9 * B8d9,
    B5d9 * B8d9 * B8d9,
    B8d9 * B8d9 * B8d9,
    B5d9 * B8d9 * B8d9,
    B5d9 * B5d9 * B8d9,
    B8d9 * B5d9 * B8d9,
    B5d9 * B5d9 * B8d9,
    B5d9 * B5d9 * B5d9,
    B8d9 * B5d9 * B5d9,
    B5d9 * B5d9 * B5d9,
    B5d9 * B8d9 * B5d9,
    B8d9 * B8d9 * B5d9,
    B5d9 * B8d9 * B5d9,
    B5d9 * B5d9 * B5d9,
    B8d9 * B5d9 * B5d9,
    B5d9 * B5d9 * B5d9
))
R_at_nodes_of_element['C3D20'] = np.array((
    (-1., -1., -1.),
    (1., -1., -1.),
    (1., 1., -1.),
    (-1., 1., -1.),
    (-1., -1., 1.),
    (1., -1., 1.),
    (1., 1., 1.),
    (-1., 1., 1.),
    (0., -1., -1.),
    (1., 0., -1.),
    (0., 1., -1.),
    (-1., 0., -1.),
    (0., -1., 1.),
    (1., 0., 1.),
    (0., 1., 1.),
    (-1., 0., 1.),
    (-1., -1., 0.),
    (1., -1., 0.),
    (1., 1., 0.),
    (-1., 1., 0.)
))
exponents_of_shape_functions_of_element['C3D20'] = np.array((
    (0, 0, 0),
    (1, 0, 0), (0, 1, 0), (0, 0, 1),
    (1, 1, 0), (1, 0, 1), (0, 1, 1),
    (1, 1, 1),
    (2, 0, 0), (0, 2, 0), (0, 0, 2),
    (2, 1, 0), (2, 0, 1), (1, 2, 0), (0, 2, 1), (1, 0, 2), (0, 1, 2),
    (2, 1, 1), (1, 2, 1), (1, 1, 2)
))

# C3D20R
R_at_integration_points_of_element['C3D20R'] = np.array((
    (-EdW3, -EdW3, -EdW3),
    (EdW3, -EdW3, -EdW3),
    (-EdW3, EdW3, -EdW3),
    (EdW3, EdW3, -EdW3),
    (-EdW3, -EdW3, EdW3),
    (EdW3, -EdW3, EdW3),
    (-EdW3, EdW3, EdW3),
    (EdW3, EdW3, EdW3)
))
weights_of_integration_points_of_element['C3D20R'] = np.array((1., 1., 1., 1., 1., 1., 1., 1.))
R_at_nodes_of_element['C3D20R'] = R_at_nodes_of_element['C3D20']
exponents_of_shape_functions_of_element['C3D20R'] = exponents_of_shape_functions_of_element['C3D20']

# C3D4
R_at_integration_points_of_element['C3D4'] = np.array(((0.25, 0.25, 0.25),))
weights_of_integration_points_of_element['C3D4'] = np.zeros(1) + 1. / 6.
R_at_nodes_of_element['C3D4'] = np.array((
    (0., 0., 0.),
    (1., 0., 0.),
    (0., 1., 0.),
    (0., 0., 1.)
))
exponents_of_shape_functions_of_element['C3D4'] = np.array(((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)))

# C3D10
R_at_integration_points_of_element['C3D10'] = np.array((
    (tri42, tri42, tri42),
    (tri41, tri42, tri42),
    (tri42, tri41, tri42),
    (tri42, tri42, tri41)
))
weights_of_integration_points_of_element['C3D10'] = np.array((1. / 24, 1. / 24, 1. / 24, 1. / 24))
R_at_nodes_of_element['C3D10'] = np.array((
    (0., 0., 0.),
    (1., 0., 0.),
    (0., 1., 0.),
    (0., 0., 1.),
    (0.5, 0., 0.),
    (0.5, 0.5, 0.),
    (0., 0.5, 0.),
    (0., 0., 0.5),
    (0.5, 0., 0.5),
    (0., 0.5, 0.5)
))
exponents_of_shape_functions_of_element['C3D10'] = np.array((
    (0, 0, 0),
    (1, 0, 0), (0, 1, 0), (0, 0, 1),
    (1, 1, 0), (1, 0, 1), (0, 1, 1),
    (2, 0, 0), (0, 2, 0), (0, 0, 2)
))

# C3D10R
# not validated against abaqus
R_at_integration_points_of_element['C3D10R'] = np.array(((0.25, 0.25, 0.25),))
weights_of_integration_points_of_element['C3D10R'] = np.zeros(1) + 1. / 6.
R_at_nodes_of_element['C3D10R'] = R_at_nodes_of_element['C3D10']
exponents_of_shape_functions_of_element['C3D10R'] = exponents_of_shape_functions_of_element['C3D10']
