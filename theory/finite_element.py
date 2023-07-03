r"""
Finite elements are commonly used to approximate a field (e.g. displacement field).
The field is approximated by shape or test functions.
Isoparametric finite elements use the same shape functions for:

- the displacement field
- the interpolation of the undeformed position of the element

>>> import numpy as np
>>> import sympy as sy
>>> import conforce_3.element_definitions as el_def
>>> import conforce_3.expressions as expr

Define an element type
----------------------

Element types are defined in the reference coordinate system and
are prototypse for elements in the real space.
In this example a 2d plain strain element with four nodes and reduced integration is considere.

>>> typ = el_def.CPE4R

An element type requires at least the reference positions of the nodes,

>>> R_at_nodes = el_def.R_at_nodes_of_element[typ]

the exponents of the polynomial powers of the shape functions,

>>> shape_exponents = el_def.exponents_of_shape_functions_of_element[typ]

the reference position of the integration points,

>>> R_at_int_points = el_def.R_at_integration_points_of_element[typ]

and the integration weights.

>>> int_weights = el_def.weights_of_integration_points_of_element[typ]

This element type has `n` nodes in `d` dimensions and `ips` integration points.

>>> n, d = R_at_nodes.shape
>>> ips = len(R_at_int_points)
>>> n, d, ips
(4, 2, 1)


Reference space `R`
-------------------

For a 2D element the reference space has 2 variables:

>>> R = expr.eval_R(d)
>>> r1, r2 = R

Shape functions in reference space `H`
--------------------------------------

The shape functions are defined in the reference coordinate system.
A shape function consists out of coefficients

>>> COEF = sy.MatrixSymbol("C", len(shape_exponents), 1).as_explicit()

and powers of the reference space variables

>>> POWER = sy.zeros(len(shape_exponents), 1)
>>> for i in range(len(shape_exponents)):
...     POWER[i, 0] = r1**shape_exponents[i, 0] * r2**shape_exponents[i, 1]
>>> POWER
Matrix([
[    1],
[   r0],
[   r1],
[r0*r1]])

The shape function is a matrix multiplication of the coefficients and the powers

>>> h0 = (COEF.T * POWER)[0]
>>> h0
r0*r1*C[3, 0] + r0*C[1, 0] + r1*C[2, 0] + C[0, 0]

The coefficients are computed out of the coordinates of the element nodes
and the powers of the shape functions.

>>> H = expr.eval_H(R, R_at_nodes, shape_exponents)
>>> H[0]
0.25*r0*r1 - 0.25*r0 - 0.25*r1 + 0.25

Each shape function corresponds to one node.
At this node the shape function is one, whereas at other nodes
it is zero.
E.g. the i-th shape function is one at the i-th node.

>>> i = 1
>>> node_i = R_at_nodes[i]
>>> float(H.xreplace(dict(zip([r1, r2], node_i)))[i])
1.0

but is zero at any j-th node (j != i).

>>> j = 0
>>> node_j = R_at_nodes[j]
>>> float(H.xreplace(dict(zip([r1, r2], node_j)))[i])
0.0


Interpolation of nodal values using `H`
---------------------------------------

Shape functions are used to interpolate quantities inside an element.
Values for a quantity are provided at the element nodes.
E.g. the coordinates are provided at the element nodes.

>>> quantity_at_nodes = sy.Matrix(R_at_nodes)

The interpolation is done by a matrix multiplication with the shape functions.

>>> interpolation = quantity_at_nodes.T * H
>>> interpolation
Matrix([
[1.0*r0],
[1.0*r1]])

In this simple example the interpolated quantities are :math:`(1*r, 1*s)`.


Derivative of interpolation in reference space using `dH_dR`
------------------------------------------------------------

To derive the interpolation of a quantity,
it is sufficient to derive the shape function
with respect to the reference space coordinates.

>>> dH_dR = expr.eval_dH_dR(H, sy.Matrix([r1, r2]))

The derivative of the quantities with respect to
the reference space coordinates is again a
matrix multiplication.

>>> dQ_dR = quantity_at_nodes.T * dH_dR

E.g. the i-th quantity derived by the j-th reference space coordinate is

>>> i, j = (0, 0)
>>> dQ_dR[i, j] # doctest: +ELLIPSIS
1.0...


Integration in the reference space using integration points
-----------------------------------------------------------

Since, the shape functions are multidimensional polynomials, they can be
integrated using Gaussian quadrature.
E.g. the function

>>> f = r1*r2 + 1

is integrated over the element volume in the reference space, by

>>> integral = 0
>>> for R_at_int_point, int_weight in zip(R_at_int_points, int_weights):
...     integral += int_weight * float(f.xreplace(dict(zip([r1, r2], R_at_int_point))))
>>> integral
4.0

However, reduced elements like CPE4R have less integration points
than it would be necessary to integrate the shape function exactly.


Real space `X`
--------------

For an application the element is transformed to the real space.
E.g. the element nodes have now real space coordinates

>>> X_at_nodes = np.array([
...     [1.0, 2.0],
...     [3.0, 3.0],
...     [4.0, 5.0],
...     [1.0, 4.0]
... ])

The coordiantes in the real space are interpolated
from the real space coordinates of the element nodes.

>>> X_at_nodes.T @ H
Matrix([
[0.25*r0*r1 + 1.25*r0 + 0.25*r1 + 2.25],
[                0.5*r0 + 1.0*r1 + 3.5]])


Jacobian of real space `X` with respect to reference space `R`
--------------------------------------------------------------

For the space transformation the jacobian is computed.

>>> dX_dR = expr.eval_dX_dR(X_at_nodes, dH_dR)

The jacobian derives the real space coordinates
with respect to the reference space coordinates.
The i-th real space coordinate derived by the
j-th reference space coordinate is

>>> i, j = (1, 0)
>>> dX_dR[i, j] # doctest: +ELLIPSIS
0.5...

Derivatives in real space `dH_dX` and `dU_dX`
---------------------------------------------

Imagine, displacements are given in the real space

>>> STRAIN = np.array([
...     [0.3, 0.0],
...     [0.0, -0.1]
... ])
>>> U_at_nodes = X_at_nodes @ STRAIN
>>> U_at_nodes
array([[ 0.3, -0.2],
       [ 0.9, -0.3],
       [ 1.2, -0.5],
       [ 0.3, -0.4]])

They can be interpolated using the shape functions

>>> U = U_at_nodes.T * H

However, at the moment you can only derive them with
respect to the reference space.

To derive an interpolation like `U` with respect to the
real space, the jacobian `dX_dR` is used to transfer `dH_dR` to

>>> dH_dX = expr.eval_dH_dX(dH_dR, dX_dR)

Now the displacements can be derived with respect to the
real space.

>>> dU_dX = U_at_nodes.T * dH_dX

This is exactly the same as

>>> dU_dX = expr.eval_dU_dX(U_at_nodes, dH_dX)


The derivative is evaluated at the center of element

>>> dU_dX_at_00 = np.array(dU_dX.subs({r1: 0., r2: 0.}), dtype=float).round(7)
>>> dU_dX_at_00
array([[ 0.3,  0. ],
       [-0. , -0.1]])

and the strain tensor is computed.

>>> 0.5 * (dU_dX_at_00 + dU_dX_at_00.T)
array([[ 0.3,  0. ],
       [ 0. , -0.1]])

.. note::
    In the small deformation theory the strain tensor is computed as

    `STRAIN = 0.5 * (dU_dX + dU_dX.T)`

    Note, that if `dU_dX` is symmetric, `dU_dX` equals the strain tensor.


Integration in the real space
-----------------------------

The jacobian `dX_dR` is also used for the integration in the real space.
E.g. the strain energy density is computed out of strain and stress tensors.
The Cauchy stress tensor is computed using Lamé parameters :math:`\lambda=600` and :math:`\mu=400`.

>>> S = 2 * 400 * STRAIN + 600 * np.trace(STRAIN) * np.eye(2)
>>> S
array([[360.,   0.],
       [  0.,  40.]])

The strain energy density is

>>> e = 0.5 * np.tensordot(STRAIN, S)
>>> e
52.0

The strain energy is the integral of the strain energy density
over the element in the real space coordinates.
In addition to the integral in the reference space,
the determinat of the jacobian `dX_dR` is added.

>>> integral = 0.
>>> for R_at_int_point, int_weight in zip(R_at_int_points, int_weights):
...     integral += (
...         int_weight
...         * e
...         * float(dX_dR.det().xreplace(dict(zip([r1, r2], R_at_int_point))))
...     )
>>> integral
234.0

"""
import doctest

if __name__ == '__main__':
    doctest.testmod()
