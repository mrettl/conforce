r"""
This module contains symbolic expressions of various physical quantities.

Abbreviations
-------------

    - `d`: number of space dimensions (2D, 3D)
    - `n`: number of nodes
    - `ips`: number of integration points

    - `R`: (d, 1)-Matrix of reference space coordinates
    - `X`: (d, 1)-Matrix of real space coordinates

    - `H`: (n, 1)-Matrix of shape functions
    - `dH_dR`: (n, d)-Matrix: :math:`dh\_dr_{ik} = \partial h_{i} / \partial r_{k}`
    - `dH_dX`: (n, d)-Matrix: :math:`dh\_dx_{ik} = \partial h_{i} / \partial x_{k}`
    - `J`: (d, d)-Jacobian matrix: :math:`j_{ik}= \partial x_{i} / \partial r_{k}`
    - `U`: (d,)-Matrix of displacements in the real space
    - `U_at_nodes`: (n, d)-Matrix of displacements in the real space at the nodes
    - `dU_dX`: (d, d)-Matrix: :math:`du\_dx_{ik} = \partial u_{i} / \partial x_{k}`
    - `S`: (d, d)-Cauchy stress tensor
    - `F`: (d, d)-Deformation gradient
    - `P`: (d, d)-First Piola-Kirchhoff stress tensor
    - `e`: internal energy density
    - `CS`: (d, d)-Matrix of configurational stresses. Implemented formulations:

        - `mbf`: motion based formulation (original formulation)
        - `dbf`: deformation based formulation

    - `CF`: (d,)-Matrix of a configurational force in the real space
    - `CF_at_nodes`: (n, d)-Matrix of configurational forces in the real space at the nodes


Define an element
-----------------

>>> import cf.element_definitions as el_def
>>> typ = el_def.CPE4R
>>> R_at_nodes = el_def.R_at_nodes_of_element[typ]
>>> shape_exponents = el_def.exponents_of_shape_functions_of_element[typ]
>>> R_at_int_points = el_def.R_at_integration_points_of_element[typ]
>>> int_weights = el_def.weights_of_integration_points_of_element[typ]

>>> n, d = R_at_nodes.shape
>>> ips = len(R_at_int_points)
>>> n, d, ips
(4, 2, 1)


Reference space `R`
-------------------

For a 2D element the reference space has 2 variables:

>>> R = eval_R(2)
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

>>> H = eval_H(R, R_at_nodes, shape_exponents)
>>> H[0]
0.25*r0*r1 - 0.25*r0 - 0.25*r1 + 0.25

Each shape function corresponds to one node.
At this node the shape function is one, whereas at other nodes
it is zero.
E.g. The i-th shape function is one at the i-th node.

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

Shape functions are used to interpolate quantities inside the element.
Values for a quantity are provided at the element nodes.
E.g. The coordinates are provided at the element nodes.

>>> quantity = sy.Matrix(R_at_nodes)

The interpolation is done by a matrix multiplication with the shape functions.

>>> interpolation = quantity.T * H
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

>>> dH_dR = eval_dH_dR(H, sy.Matrix([r1, r2]))

The derivative of the quantities with respect to
the reference space coordinates is again a
matrix multiplication.

>>> dQ_dR = quantity.T * dH_dR

E.g. the i-th quantity derived by the j-th reference space coordinate is

>>> i, j = (0, 0)
>>> dQ_dR[i, j]
1.00000000000000


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


Jacobian matrix `J`
-------------------

For the space transformation the jacobian is computed.

>>> J = eval_J(X_at_nodes, dH_dR)

The jacobian derives the real space coordinates
with respect to the reference space coordinates.
The i-th real space coordinate derived by the
j-th reference space coordinate is

>>> i, j = (1, 0)
>>> J[i, j]
0.500000000000000

.. note::

    In some conventions the determinant of the deformation gradient :math:`det(F)`
    is called `J` or the Jacobian determinant.
    However, this is a completly other Jacobian matrix as defined here.
    This Jacobian referes to the trasformation from the reference space into the real space,
    whereas the determinant of the deformation gradient considers the transformation from
    the undeformed real space into the deformed real space.

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
real space, the jacobian is used to transfer `dH_dR` to

>>> dH_dX = eval_dH_dX(dH_dR, J)

Now the displacements can be derived with respect to the
real space.

>>> dU_dX = U_at_nodes.T * dH_dX

This is exactly the same as

>>> dU_dX = eval_dU_dX(U_at_nodes, dH_dX)


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

The jacobian is also used for the integration in the real space.
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
the determinat of the jacobian is added.

>>> integral = 0.
>>> for R_at_int_point, int_weight in zip(R_at_int_points, int_weights):
...     integral += (
...         int_weight
...         * e
...         * float(J.det().xreplace(dict(zip([r1, r2], R_at_int_point))))
...     )
>>> integral
234.0


Deformation gradient `F`
------------------------

The deformation gradient

.. math::
    F = \frac{\partial (X + U)}{\partial X} = I + \frac{\partial U}{\partial X}

is the derivative of the deformed coordinates with respect to the undeformed coordinates.

.. note::

    In the 2D case, a plane strain state is assumed.
    This implies :math:`F_{zz}=1`.

Undeformed
^^^^^^^^^^

The deformation gradient of the undeformed state in 2D is

>>> F_undeformed = eval_F(2, sy.ZeroMatrix(2, 2))
>>> F_undeformed
Matrix([
[1, 0],
[0, 1]])

The determinant of the deformation gradient states the volume change and is

>>> sy.det(F_undeformed)
1

for the undeformed state. So the volume does not change.

Deformed
^^^^^^^^

However, for a deformation of :math:`\varepsilon_{xx}=0.2; \gamma_{xy}=0.1`,
the deformation gradient is

>>> F_deformed = eval_F(2, sy.Matrix([[0.2, 0.1], [0, 0]]))
>>> F_deformed
Matrix([
[1.2, 0.1],
[  0,   1]])

and its determinat

>>> float(sy.det(F_deformed))
1.2

states that the volume increases by a factor of 1.2.
Note, that only the tensile strain :math:`\varepsilon_{xx}=0.2` results in a volume change.
The simple shear :math:`\gamma_{xy}=0.1` does not influence the volume.

Pure shear
^^^^^^^^^^

However, a pure shear deformation of :math:`\gamma_{xy}=0.1` without a tensile strain
would decrease the volume by a factor of

>>> float(sy.det(eval_F(2, sy.Matrix([[0.0, 0.1], [0.1, 0]]))))
0.99

The deformation gradient is a linear transformation from the
undeformed to the deformation state.
Considering a simple shear deformation,

>>> F = eval_F(2, sy.Matrix([[0.0, 0.1], [0.0, 0]]))

the points of the undeformed unit square can be transformed to
the deformed state by a vector matrix multiplication

>>> X_unit_square = np.array([
...     [0, 0],
...     [0, 1],
...     [1, 1],
...     [1, 0],
... ])
>>> X_deformed_unit_square = X_unit_square @ F.T
>>> X_deformed_unit_square
Matrix([
[  0, 0],
[0.1, 1],
[1.1, 1],
[  1, 0]])

For the back transformation from the deformed to the undeformed state,
the inverse of the deformation gradient is used.

>>> np.array(X_deformed_unit_square @ F.inv().T, dtype=float)
array([[0., 0.],
       [0., 1.],
       [1., 1.],
       [1., 0.]])


First Piola-Kirchhoff stresses `P`
----------------------------------

The Cauchy stress `S` considers the deformed state,
whereas the first Piola-Kirchhoff stress `P` is defined on the undeformed state.

Tensile Test
^^^^^^^^^^^^

Assuming a tensile load with the reaction force

>>> RF = sy.Matrix([
...     [10_000.],
...     [0.]])

acting on the upper surface of a specimen with
an undeformed area and undeformed length of

>>> width_undeformed = 10.
>>> length_undeformed = 100.

This load causes a strain of :math:`\varepsilon_{x} = 0.1`
and :math:`\varepsilon_{y} = -0.03`.
The deformation gradient is

>>> F = eval_F(2, sy.Matrix([
...     [0.1, 0.0],
...     [0.0, -0.03]]))

The traction (force per area with depth=1)

>>> Tx_undeformed = RF / (1 * width_undeformed)
>>> Tx_undeformed
Matrix([
[1000.0],
[     0]])

acts on the undeformed upper surface with normal vector

>>> Nx_undeformed = sy.Matrix([
...     [1.],
...     [0.]])

Furthermore, no force acts on the side of the specimen
with normal vector

>>> Ny_undeformed = sy.Matrix([
...     [0.],
...     [1.]])

and hence the corresponding traction is

>>> Ty_undeformed = sy.Matrix([
...     [0.],
...     [0.]
... ]) / (1*length_undeformed)

The first Piola-Kirchhoff stress is

>>> T_undeformed = sy.Matrix.hstack(Tx_undeformed, Ty_undeformed)
>>> N_undeformed = sy.Matrix.hstack(Nx_undeformed, Ny_undeformed)
>>> P = T_undeformed * N_undeformed.inv()
>>> P
Matrix([
[1000.0, 0],
[     0, 0]])

Compare this to the computation of the true Cauchy stress
that corresponds to the deformed state with 

>>> (length_deformed,), (width_deformed,) = np.around(np.array(
...     F @ sy.Matrix([
...         [length_undeformed],
...         [width_undeformed]
... ]), dtype=float), 6)
>>> length_deformed
110.0
>>> width_deformed
9.7

For pure tension, the normal vectors of the surfaces
stay the same.

>>> Nx_deformed = Nx_undeformed
>>> Ny_deformed = Ny_undeformed
>>> N_deformed = N_undeformed

The tractions of the Cauchy stress are computed with respect to
the deformed width and length.

>>> Tx_deformed = RF / (1 * width_deformed)
>>> TY_deformed = sy.Matrix([
...     [0.],
...     [0.]
... ]) / (1*length_deformed)
>>> T_deformed = sy.Matrix.hstack(Tx_deformed, TY_deformed)
>>> S = T_deformed * N_deformed.inv()
>>> np.around(np.array(S, dtype=float), 6)
array([[1030.927835,    0.      ],
       [   0.      ,    0.      ]])

The first Piola-Kirchhoff stress can be also computed
directly with the Cauchy stress and the deformation gradient.

>>> eval_P(F, S)
Matrix([
[1000.0, 0],
[     0, 0]])


Simple shear and tension
^^^^^^^^^^^^^^^^^^^^^^^^

Considering a combined load case (:math:`\gamma_{yx} = 0.1; \varepsilon_{xx}=0.2`)
with a deformation gradient

>>> F = eval_F(2, sy.Matrix([
...     [0.2, 0.1],
...     [0.0, 0.0]]))

and a Cauchy stress of

>>> S = sy.Matrix([
...     [280., 40.],
...     [40.0, 0.0]])

the corresponding first Piola-Kirchhoff tensor is

>>> eval_P(F, S)
Matrix([
[276.0, 48.0],
[ 40.0,    0]])

not symmetric anymore.


Configurational stresses `CS`
-----------------------------

.. todo: CS


Configurational forces `CF`
---------------------------

.. todo: CF

Computation
-----------

.. todo:
    Computation

"""
from typing import Union, List, Dict, Any

import numpy as np
import sympy as sy

from cf.symbolic_util import create_symbolic_matrix, inverse, create_replacement_rules, apply_replacement_rules

R_3d = create_symbolic_matrix("r{row}", 3, 1)
"""symbolic reference space coordinates for the 3d space"""


def eval_R(d_: int):
    return sy.Matrix(R_3d[:d_])


def eval_H(R: sy.MatrixBase, R_at_nodes_: np.ndarray, exponents_: np.ndarray):
    R_at_nodes_ = np.array(R_at_nodes_, dtype=float)
    exponents_ = np.array(exponents_, dtype=int)

    num_points_ = R_at_nodes_.shape[0]
    num_powers_ = exponents_.shape[0]
    assert num_points_ == num_powers_

    # symbolic powers of shape functions
    POWERS_ = sy.Matrix([
        sy.Mul(*[
            sy.Pow(ri, power)
            for ri, power in zip(R.T.tolist()[0], powers)
        ]).nsimplify()
        for powers in exponents_
    ])

    # A_[i, j] = i-th shape power evaluated at j-th point
    A_ = sy.Matrix([
        POWERS_.T.subs({
            ri: xi
            for ri, xi
            in zip(R.T.tolist()[0], point)
        })
        for point in R_at_nodes_
    ]).T

    # coefficients solve the system COEFFICIENTS_ * A_ = I,
    # such that every shape function is one at exactly one point and zero at the other points
    COEFFICIENTS_ = A_.inv()

    # shape functions are the combination of coefficients and shapes powers
    H_ = COEFFICIENTS_ * POWERS_

    return H_


def eval_dH_dR(H: sy.MatrixBase, R: sy.MatrixBase):
    n_ = H.shape[0]
    d_ = R.shape[0]

    dH_dR_ = sy.zeros(n_, d_)
    for i in range(n_):
        for j in range(d_):
            dH_dR_[i, j] = H[i].diff(R[j])

    return dH_dR_


def eval_J(X_at_nodes: Union[sy.MatrixBase, np.ndarray], dH_dR: sy.MatrixBase):
    return X_at_nodes.T * dH_dR


def eval_dH_dX(dH_dR: sy.MatrixBase, J: sy.MatrixBase):
    return dH_dR * inverse(J)


def eval_dU_dX(U_at_nodes: Union[sy.MatrixBase, np.ndarray], dH_dX: sy.MatrixBase):
    return sy.Matrix(U_at_nodes.T) * dH_dX


def eval_F(d_: int, dU_dX: sy.MatrixBase):
    return sy.eye(d_) + dU_dX


def eval_P(F: sy.MatrixBase, S: sy.MatrixBase):
    return S * inverse(F).T * F.det()


def eval_CS_mbf(d_: int, e: sy.Expr, F: sy.MatrixBase, P: sy.MatrixBase):
    return e * sy.eye(d_) - F.T * P


def eval_CS_dbf(d_: int, e: sy.Expr, dU_dX: sy.MatrixBase, P: sy.MatrixBase):
    return e * sy.eye(d_) - dU_dX.T * P


def eval_CF_at_nodes(
        dH_dX: sy.MatrixBase,
        CS: sy.MatrixBase,
        J: sy.MatrixBase,
        int_weights: np.ndarray,
        int_points_replacements: List[Dict[sy.Expr, Any]]):

    CF_contributions = list()
    for w, int_point_replacement in zip(int_weights, int_points_replacements):
        expr_ = dH_dX * (CS.T * J.det() * w)

        expr_ = expr_.xreplace(int_point_replacement)
        CF_contributions.append(expr_)

    CF_at_nodes_ = sy.MatAdd(*CF_contributions).doit()
    return CF_at_nodes_


class Computation(object):
    def __init__(
        self,
        R_at_nodes_, exponents_,
        R_at_int_points_, int_weights_,
        X_at_nodes_, U_at_nodes_,
        S_at_int_points_, e_at_int_points_,
        is_dbf: bool = True
    ):
        """
        Compute `F`, `P`, `CS`, and `CF`.

        :param R_at_nodes_: array of shape (n, d) of reference space nodal coordinates
        :param exponents_: array of shape (n, d) of integer exponents
        :param R_at_int_points_: array of shape (ips, d) of reference space coordinates
        :param int_weights_: array of shape (ips,)
        :param X_at_nodes_: array of shape (n, d) of real space nodal coordinates
        :param U_at_nodes_: array of shape (n, d) of real space nodal displacements
        :param S_at_int_points_: array of shape (ips, d, d) of real space Cauchy stresses
        :param e_at_int_points_: array of shape (ips,) of energy densities
        :param is_dbf: use "dbf" if True else "mbf"
        """

        # check shapes
        n_, d_ = R_at_nodes_.shape
        ips_ = len(int_weights_)

        assert d_ in {2, 3}
        assert R_at_nodes_.shape == (n_, d_)
        assert exponents_.shape == (n_, d_)
        assert R_at_int_points_.shape == (ips_, d_)
        assert int_weights_.shape == (ips_,)
        assert X_at_nodes_.shape == (n_, d_)
        assert U_at_nodes_.shape == (n_, d_)
        assert S_at_int_points_.shape == (ips_, d_, d_)
        assert e_at_int_points_.shape == (ips_,)

        #
        symbols_to_expressions = dict()
        n_, d_ = R_at_nodes_.shape
        dim_n = list(range(n_))

        R = eval_R(d_)
        dim_r = [r.name for r in R]

        X = create_symbolic_matrix("x{row}", d_, 1, *R)
        dim_x = [x.name for x in X]

        U = create_symbolic_matrix("U{row}", dim_x, 1, *R)
        S = create_symbolic_matrix("S{row}{col}", dim_x, dim_x, *R, is_symmetric=True)

        e = sy.Function("e", real=True)(*R)
        symbols_to_expressions[e] = e_at_int_points_

        #
        replacements_by_nodes = create_replacement_rules(R, *R_at_nodes_)
        replacements_by_int_points = create_replacement_rules(R, *R_at_int_points_)

        X_at_nodes = sy.Matrix(apply_replacement_rules(X, *replacements_by_nodes)[:, :, 0]).as_immutable()
        symbols_to_expressions[X_at_nodes] = X_at_nodes_

        U_at_nodes = sy.Matrix(apply_replacement_rules(U, *replacements_by_nodes)[:, :, 0]).as_immutable()
        symbols_to_expressions[U_at_nodes] = U_at_nodes_

        S_at_int_points = apply_replacement_rules(S, *replacements_by_int_points)
        symbols_to_expressions[S_at_int_points] = S_at_int_points_

        e_at_int_points = apply_replacement_rules(e, *replacements_by_int_points)
        symbols_to_expressions[e_at_int_points] = e_at_int_points_

        #
        H = create_symbolic_matrix("H{row}", dim_n, 1, *R)
        H_ = eval_H(R, R_at_nodes_, exponents_).doit()
        symbols_to_expressions[H] = H_

        dH_dR = create_symbolic_matrix("dH{row}_d{col}", dim_n, dim_r, *R)
        dH_dR_ = eval_dH_dR(H_, R).doit()
        symbols_to_expressions[dH_dR] = dH_dR_

        J = create_symbolic_matrix("J{row}{col}", dim_x, dim_r, *R)
        J_ = eval_J(X_at_nodes, dH_dR_).doit()
        symbols_to_expressions[J] = J_

        dH_dX = create_symbolic_matrix("dH{row}_d{col}", dim_n, dim_x, *R)
        dH_dX_ = eval_dH_dX(dH_dR_, J).doit()
        symbols_to_expressions[dH_dX] = dH_dX_

        dU_dX = create_symbolic_matrix("dU{row}_d{col}", dim_x, dim_x, *R)
        dU_dX_ = eval_dU_dX(U_at_nodes, dH_dX).doit()
        symbols_to_expressions[dU_dX] = dU_dX_

        #
        F = create_symbolic_matrix("F{row}{col}", dim_x, dim_x, *R)
        F_ = eval_F(d_, dU_dX).doit()
        symbols_to_expressions[F] = F_

        P = create_symbolic_matrix("P{row}{col}", dim_x, dim_x, *R)
        P_ = eval_P(F, S).doit()
        symbols_to_expressions[P] = P_

        CS = create_symbolic_matrix("CS{row}{col}", dim_x, dim_x, *R)
        if is_dbf:
            CS_ = eval_CS_dbf(d_, e, dU_dX, P).doit()
        else:
            CS_ = eval_CS_mbf(d_, e, F, P).doit()

        symbols_to_expressions[CS] = CS_

        #
        CF_at_nodes = sy.MatrixSymbol("CF_at_nodes", n_, d_)
        CF_at_nodes_ = eval_CF_at_nodes(
            dH_dX,
            CS,
            J,
            int_weights_,
            replacements_by_int_points
        ).doit()
        symbols_to_expressions[CF_at_nodes] = CF_at_nodes_

        #
        self.symbols_to_expressions = symbols_to_expressions
        self.replacements_by_int_points = replacements_by_int_points
        self.R = R
        self.F = F
        self.P = P
        self.CS = CS
        self.CF_at_nodes = CF_at_nodes

    def map_symbolic_to_expression(self, symbolic: sy.Expr, expr: sy.Expr):
        """
        Add the mapping `symbolic -> expr` to :py:attr:`symbols_to_expressions`.

        :param symbolic: symbol that is mapped to the expression
        :param expr: An expression
        """
        self.symbols_to_expressions[symbolic] = expr
