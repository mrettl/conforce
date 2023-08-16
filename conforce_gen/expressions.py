r"""
This module contains symbolic expressions of various physical quantities.
"""
from typing import Union, List, Dict, Any

import numpy as np
import sympy as sy

from conforce_gen.symbolic_util import create_symbolic_matrix, inverse, create_replacement_rules, apply_replacement_rules

R_3d = create_symbolic_matrix("r{row}", 3, 1)
"""symbolic reference space coordinates for the 3d space"""


def eval_R(d: int):
    """
    Evaluate symbolic reference space coordinates.

    :param d: number of dimensions (2 or 3)
    :return: symbolic matrix of shape (d, 1)
    """
    return sy.Matrix(R_3d[:d])


def eval_H(R: sy.MatrixBase, R_at_nodes: np.ndarray, exponents: np.ndarray):
    """
    Evaluate all shape functions of an element.

    :param R: matrix of shape (d, 1); symbolic reference space coordinates
    :param R_at_nodes: array of shape (n, d); reference space coordinates at nodes
    :param exponents: matrix of shape (n, d); exponents of powers of the shape functions
    :return: symbolic matrix of shape (n, 1)
    """
    R_at_nodes = np.array(R_at_nodes, dtype=float)
    exponents = np.array(exponents, dtype=int)

    num_points_ = R_at_nodes.shape[0]
    num_powers_ = exponents.shape[0]
    assert num_points_ == num_powers_

    # symbolic powers of shape functions
    POWERS_ = sy.Matrix([
        sy.Mul(*[
            sy.Pow(ri, power)
            for ri, power in zip(R.T.tolist()[0], powers)
        ]).nsimplify()
        for powers in exponents
    ])

    # A_[i, j] = i-th shape power evaluated at j-th point
    A_ = sy.Matrix([
        POWERS_.T.subs({
            ri: xi
            for ri, xi
            in zip(R.T.tolist()[0], point)
        })
        for point in R_at_nodes
    ]).T

    # coefficients solve the system COEFFICIENTS_ * A_ = I,
    # such that every shape function is one at exactly one point and zero at the other points
    COEFFICIENTS_ = A_.inv()

    # shape functions are the combination of coefficients and shapes powers
    H_ = COEFFICIENTS_ * POWERS_

    return H_


def eval_dH_dR(H: sy.MatrixBase, R: sy.MatrixBase):
    """
    Evaluate derivatives of shape functions with respect to reference space coordinates.

    :param H: matrix of shape (n, 1); shape functions
    :param R: matrix of shape (d, 1); symbolic reference space coordinates
    :return: symbolic matrix of shape (n, d)
    """
    n_ = H.shape[0]
    d_ = R.shape[0]

    dH_dR_ = sy.zeros(n_, d_)
    for i in range(n_):
        for j in range(d_):
            dH_dR_[i, j] = H[i].diff(R[j])

    return dH_dR_


def eval_dX_dR(X_at_nodes: Union[sy.MatrixBase, np.ndarray], dH_dR: sy.MatrixBase):
    """
    Evaluate derivatives of undeformed real space coordinates with respect to reference space coordinates.

    :param X_at_nodes: matrix of shape (n, d); undeformed real space coordinates at nodes
    :param dH_dR: matrix of shape (n, d); derivatives of shape functions with respect to reference space coordinates
    :return: symbolic matrix with shape (d, d)
    """
    return X_at_nodes.T * dH_dR


def eval_dH_dX(dH_dR: sy.MatrixBase, dX_dR: sy.MatrixBase):
    """
    Evaluate derivatives of shape functions with respect to undeformed real space coordinates.

    :param dH_dR: matrix of shape (n, d); derivatives of shape functions with respect to reference space coordinates
    :param dX_dR: matrix of shape (d, d);
        derivatives of undeformed real space coordinates with respect to reference space coordinates

    :return: symbolic matrix of shape (n, d)
    """
    return dH_dR * inverse(dX_dR)


def eval_dU_dX(U_at_nodes: Union[sy.MatrixBase, np.ndarray], dH_dX: sy.MatrixBase):
    """
    Evaluate derivatives of displacements with respect to undeformed real space coordinates.

    :param U_at_nodes: matrix of shape (n, d); displacements at nodes
    :param dH_dX: matrix of shape (n, d);
        derivatives of shape functions with respect to undeformed real space coordinates

    :return: symbolic matrix of shape (d, d)
    """
    return sy.Matrix(U_at_nodes.T) * dH_dX


def eval_F(d: int, dU_dX: sy.MatrixBase):
    """
    Evaluate the deformation gradient.

    :param d: number of dimensions
    :param dU_dX: matrix of shape (d, d);
        derivatives of displacements with respect to undeformed real space coordinates

    :return: symbolic matrix of shape (d, d)
    """
    return sy.eye(d) + dU_dX


def eval_P(F: sy.MatrixBase, S: sy.MatrixBase):
    """
    Evaluate the First Piola-Kirchhoff stress tensor.

    :param F: matrix of shape (d, d); deformation gradient
    :param S: matrix of shape (d, d); symmetric Cauchy stress tensor
    :return: symbolic matrix of shape (d, d)
    """
    return S * inverse(F).T * F.det()


def eval_CS_mbf(d: int, e: sy.Expr, F: sy.MatrixBase, P: sy.MatrixBase):
    """
    Evaluate motion based configurational stress tensor.

    :param d: number of dimensions
    :param e: energy density
    :param F: matrix of shape (d, d); deformation gradient
    :param P: matrix of shape (d, d); First Piola-Kirchhoff stress tensor
    :return: symbolic matrix of shape (d, d)
    """
    return e * sy.eye(d) - F.T * P


def eval_CS_dbf(d: int, e: sy.Expr, dU_dX: sy.MatrixBase, P: sy.MatrixBase):
    """
    Evaluate displacement based configurational stress tensor.

    :param d: number of dimensions
    :param e: energy density
    :param dU_dX: matrix of shape (d, d);
        derivatives of displacements with respect to undeformed real space coordinates

    :param P: matrix of shape (d, d); First Piola-Kirchhoff stress tensor
    :return: symbolic matrix of shape (d, d)
    """
    return e * sy.eye(d) - dU_dX.T * P


def eval_CF_at_nodes(
        dH_dX: sy.MatrixBase,
        CS: sy.MatrixBase,
        dX_dR: sy.MatrixBase,
        int_weights: np.ndarray,
        int_points_replacements: List[Dict[sy.Expr, Any]]):
    """
    Evaluate the configurational forces at the element nodes.

    :param dH_dX: derivatives of shape functions with respect to undeformed real space coordinates
    :param CS: matrix of shape (d, d); configurational stress tensor
    :param dX_dR: matrix of shape (d, d);
        derivatives of undeformed real space coordinates with respect to reference space coordinates

    :param int_weights: matrix of shape (ips, );
        integration weights corresponding to the integration points

    :param int_points_replacements: list of length (ips, );
        Each list entry corresponds to an integration point
        and contains a dictionary that maps the reference space coordinate symbols
        to the position of an integration point. (E.g. `[{r: 0., s: 0.}]`)

    :return: symbolic matrix of shape (n, d)
    """

    CF_contributions = list()
    for w, int_point_replacement in zip(int_weights, int_points_replacements):
        expr_ = dH_dX * (CS.T * dX_dR.det() * w)

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

        dX_dR = create_symbolic_matrix("d{row}_d{col}", dim_x, dim_r, *R)
        dX_dR_ = eval_dX_dR(X_at_nodes, dH_dR_).doit()
        symbols_to_expressions[dX_dR] = dX_dR_

        dH_dX = create_symbolic_matrix("dH{row}_d{col}", dim_n, dim_x, *R)
        dH_dX_ = eval_dH_dX(dH_dR_, dX_dR).doit()
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
            dX_dR,
            int_weights_,
            replacements_by_int_points
        ).doit()
        symbols_to_expressions[CF_at_nodes] = CF_at_nodes_

        #
        self.symbols_to_expressions = symbols_to_expressions
        """dictionary mapping symbols to (symbolic) expressions"""

        self.replacements_by_int_points = replacements_by_int_points
        """list of dictionary in the form [{r0: 0, r1: 0}, {r0: 1, r1: 0}, ...]"""

        self.R = R
        """symbolic reference coordinates"""

        self.F = F
        """deformation gradient"""

        self.P = P
        """first Piola-Kirchhoff stress tensor"""

        self.CS = CS
        """configurational stresses"""

        self.CF_at_nodes = CF_at_nodes
        """configurational forces at element nodes"""

    def map_symbolic_to_expression(self, symbolic: sy.Expr, expr: sy.Expr):
        """
        Add the mapping `symbolic -> expr` to :py:attr:`symbols_to_expressions`.

        :param symbolic: symbol that is mapped to the expression
        :param expr: An expression
        """
        self.symbols_to_expressions[symbolic] = expr
