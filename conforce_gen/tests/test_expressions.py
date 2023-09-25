import unittest

import numpy as np

from conforce_gen.expressions import *
from conforce_gen.element_definitions import (
    R_at_nodes_of_element,
    exponents_of_shape_functions_of_element,
    weights_of_integration_points_of_element,
    R_at_integration_points_of_element
)
from conforce_gen.symbolic_util import create_replacement_rules, apply_replacement_rules


class TestExpressions(unittest.TestCase):
    def test_eval_H(self):
        """
        test if each shape function is one at exactly one node and zero at other nodes
        """
        for element_type in R_at_nodes_of_element.keys():
            with self.subTest(element_type):
                nodes = R_at_nodes_of_element[element_type]
                exponents = exponents_of_shape_functions_of_element[element_type]

                n, d = nodes.shape
                R = eval_R(d)
                H = eval_H(R, nodes, exponents)
                H_at_nodes = apply_replacement_rules(
                    H, *create_replacement_rules(R, *nodes))

                for i, H_at_node_i in enumerate(H_at_nodes):
                    np.testing.assert_array_almost_equal(
                        np.eye(n)[i],
                        np.array(H_at_node_i, dtype=float).flat
                    )

    def test_eval_F(self):
        """
        Deforms an element as defined by a deformation gradient.
        Computes the deformation gradients of the deformed element.
        Check if the computed deformation gradient matches the defined deformation gradient.
        """
        for element_type in R_at_nodes_of_element.keys():
            with self.subTest(element_type):
                R_at_nodes = R_at_nodes_of_element[element_type]
                exponents = exponents_of_shape_functions_of_element[element_type]

                n, d = R_at_nodes.shape
                R = eval_R(d)
                R_at_origin = {r: v for r, v in zip(R, [0.]*len(R))}
                R_at_1 = {r: v for r, v in zip(R, [1.]*len(R))}

                X_at_nodes = R_at_nodes @ np.array([
                    [2.0, 0.3, 0.2],
                    [0.1, 3.0, 0.1],
                    [0.2, 0.3, 0.9]
                ])[:d, :d] + np.array([123., 456., -789.])[:d]

                F_expected = np.array([
                    [2., 0.1, 0.01],
                    [0., 0.7, 0.02],
                    [-0.1, 0.03, 0.7]
                ])[:d, :d]

                U_at_nodes = X_at_nodes @ F_expected.T - X_at_nodes
                H = eval_H(R, R_at_nodes, exponents)
                dH_dR = eval_dH_dR(H, R)
                dX_dR = eval_dX_dR(X_at_nodes, dH_dR)
                dH_dX = eval_dH_dX(dH_dR, dX_dR)
                dU_dX = eval_dU_dX(U_at_nodes, dH_dX)
                F = eval_F(d, dU_dX)

                for R_values in [R_at_origin, R_at_1]:
                    np.testing.assert_array_almost_equal(
                        F_expected,
                        np.array(F.subs(R_values), dtype=float)
                    )

    def test_eval_P(self):
        """
        validate first piola kirchhoff tensor.
        Bergstroem [1]_ provides a similar example

        .. [1] Bergstroem, mechanics of solid polymers, page 169
        """
        for typ in R_at_nodes_of_element.keys():
            with self.subTest(typ):
                R_at_nodes = R_at_nodes_of_element[typ]
                n, d = R_at_nodes.shape

                F_expected = np.array([
                    [2., 0., 0.],
                    [0., 0.7, 0.],
                    [0., 0., 1.]
                ])[:d, :d]
                sxx = sy.symbols("sxx", real=True)
                S = sy.Matrix([
                    [sxx, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0]
                ])[:d, :d]
                P_expected = sy.Matrix([
                    [0.7 * sxx, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0]
                ])[:d, :d]

                P = eval_P(sy.Matrix(F_expected), S[:d, :d])

                for sxx_ in [0, 1, 2]:
                    np.testing.assert_array_almost_equal(
                        np.array(P_expected.replace(sxx, sxx_), dtype=float),
                        np.array(P.replace(sxx, sxx_), dtype=float)
                    )


if __name__ == '__main__':
    unittest.main()
