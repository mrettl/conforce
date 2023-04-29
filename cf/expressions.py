from itertools import product, chain
from typing import Union, List, Dict, Any

import numpy as np
import sympy as sy

from cf.math_util import create_symbolic_matrix, inverse, create_replacement_rules, apply_replacement_rules

R_3d = sy.Matrix(sy.symbols("r s t", real=True))


def eval_R(d_: int):
    return sy.Matrix(R_3d[:d_])


def eval_H(R: sy.MatrixBase, R_at_nodes_: np.ndarray, exponents_: np.ndarray):
    R_at_nodes_ = np.array(R_at_nodes_, dtype=float)
    exponents_ = np.array(exponents_, dtype=int)

    num_points_ = R_at_nodes_.shape[0]
    num_powers_ = exponents_.shape[0]
    assert num_points_ == num_powers_

    POWERS_ = sy.Matrix([
        sy.Mul(*[
            sy.Pow(rst, power)
            for rst, power in zip(R.T.tolist()[0], powers)
        ]).nsimplify()
        for powers in exponents_
    ])
    """symbolic powers of shape functions"""

    A_ = sy.Matrix([
        POWERS_.T.subs({
            rst: xyz
            for rst, xyz
            in zip(R.T.tolist()[0], point)
        })
        for point in R_at_nodes_
    ]).T
    """A_[i, j] = i-th shape power evaluated at j-th point"""

    COEFFICIENTS_ = A_.inv()
    """
    coefficients solve the system COEFFICIENTS_ * A_ = I,
    such that every shape function is one at exactly one point and zero at the other points
    """

    H_ = COEFFICIENTS_ * POWERS_
    """shape functions are the combination of coefficients and shapes powers"""

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


def eval_CF(
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

        X = create_symbolic_matrix("{row}", ["x", "y", "z"][:d_], 1, *R)
        dim_x = [x.name for x in X]

        U = create_symbolic_matrix("U{row}", dim_x, 1, *R)
        S = create_symbolic_matrix("S{row}{col}", dim_x, dim_x, *R, is_symmetric=True)

        e = sy.Function("e", real=True)(*R)
        symbols_to_expressions[e] = e_at_int_points_

        #
        nodes_replacements = create_replacement_rules(R, R_at_nodes_)
        int_points_replacements = create_replacement_rules(R, R_at_int_points_)

        X_at_nodes = sy.Matrix(apply_replacement_rules(X, nodes_replacements)[:, :, 0]).as_immutable()
        symbols_to_expressions[X_at_nodes] = X_at_nodes_

        U_at_nodes = sy.Matrix(apply_replacement_rules(U, nodes_replacements)[:, :, 0]).as_immutable()
        symbols_to_expressions[U_at_nodes] = U_at_nodes_

        S_at_int_points = apply_replacement_rules(S, int_points_replacements)
        symbols_to_expressions[S_at_int_points] = S_at_int_points_

        e_at_int_points = apply_replacement_rules(e, int_points_replacements)
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
        CF_at_nodes_ = eval_CF(
            dH_dX,
            CS,
            J,
            int_weights_,
            int_points_replacements
        ).doit()
        symbols_to_expressions[CF_at_nodes] = CF_at_nodes_

        #
        self.symbols_to_expressions = symbols_to_expressions
        self.int_points_replacements = int_points_replacements
        self.R = R
        self.F = F
        self.P = P
        self.CS = CS
        self.CF_at_nodes = CF_at_nodes

    def expand_matrices_in_symbols_to_expressions(self):
        self.symbols_to_expressions.update(chain(
            *[
                {
                    symbols[idx]: expressions[idx]
                    for idx in product(*[range(int(dim)) for dim in symbols.shape])
                }.items()
                for symbols, expressions in self.symbols_to_expressions.items()
                if hasattr(symbols, "shape")
            ]
        ))

    def map_symbolic_to_expression(self, symbolic: sy.Expr, expr: sy.Expr):
        self.symbols_to_expressions[symbolic] = expr

    def at_int_point(self, symbolic: Union[sy.MatrixBase, sy.MatrixExpr]):
        return apply_replacement_rules(
            symbolic,
            self.int_points_replacements
        )
