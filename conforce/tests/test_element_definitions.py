import unittest

import sympy as sy

from conforce.element_definitions import *
from conforce.one_element_runner import simulate_one_element
from conforce.expressions import eval_R, eval_H, eval_dH_dR, eval_dX_dR, R_3d
from conforce.symbolic_util import create_replacement_rules, apply_replacement_rules


class TestElementDefinitions(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        cls.abaqus_data = {
            element_type: simulate_one_element(
                X_at_nodes=R_at_nodes,
                U_at_nodes=np.ones_like(R_at_nodes),
                element_type=element_type,
                load_name="translation",
                folder="res/tests/test_element_definitions"
            )
            for element_type, R_at_nodes in R_at_nodes_of_element.items()
        }

    def test_check_shapes(self):
        for element_type, _ in self.abaqus_data.items():
            with self.subTest(element_type):
                R_at_nodes = R_at_nodes_of_element[element_type]
                R_at_int_points = R_at_integration_points_of_element[element_type]
                weights_of_int_points = weights_of_integration_points_of_element[element_type]
                exponents = exponents_of_shape_functions_of_element[element_type]
                corner_nodes = corner_nodes_of_element[element_type]
                adjacency_matrix = adjacent_nodes_of_element[element_type]

                n, d = R_at_nodes.shape
                ips = R_at_int_points.shape[0]

                # check shapes
                self.assertEqual((n, d), R_at_nodes.shape)
                self.assertEqual((ips, d), R_at_int_points.shape)
                self.assertEqual((ips,), weights_of_int_points.shape)
                self.assertEqual((n, d), exponents.shape)
                self.assertEqual((n,), corner_nodes.shape)
                self.assertEqual((n, n), adjacency_matrix.shape)

    def test_sum_of_shape_functions(self):
        for element_type in R_at_nodes_of_element.keys():
            with self.subTest(element_type):
                R_at_nodes = R_at_nodes_of_element[element_type]
                exponents = exponents_of_shape_functions_of_element[element_type]

                # computation
                _, d = R_at_nodes.shape
                R = eval_R(d)
                H = eval_H(R, R_at_nodes, exponents)

                # test
                self.assertAlmostEqual(
                    1.,
                    sum(H)
                )

    def test_validation_against_abaqus(self):
        for element_type, data in self.abaqus_data.items():
            with self.subTest(element_type):
                # current
                R_at_nodes_current = R_at_nodes_of_element[element_type]
                R_at_int_points_current = R_at_integration_points_of_element[element_type]
                weights_of_int_points_current = weights_of_integration_points_of_element[element_type]
                exponents = exponents_of_shape_functions_of_element[element_type]

                # abaqus data
                X_at_nodes_abaqus = np.array(data["nodes"]["COORD"])
                X_at_int_points_abaqus = np.array(data["integration_points"]["COORD"])
                vol_abaqus = data["element"]["EVOL"]

                # computation
                _, d = X_at_nodes_abaqus.shape
                R = sy.Matrix(R_3d[:d])
                H = eval_H(R, R_at_nodes_current, exponents)
                X = X_at_nodes_abaqus.T @ H
                dH_dR = eval_dH_dR(H, R)
                dX_dR = eval_dX_dR(X_at_nodes_abaqus, dH_dR)
                detJ = dX_dR.det()

                int_repl_rules = create_replacement_rules(R, *R_at_int_points_current)

                # check volume
                vol_current = np.sum([
                    weight_at_int_point * detJ.xreplace(int_repl_rule)
                    for weight_at_int_point, int_repl_rule
                    in zip(weights_of_int_points_current, int_repl_rules)
                ])
                self.assertAlmostEqual(
                    vol_abaqus,
                    vol_current,
                    places=6
                )

                # check integration point position
                X_at_int_points_current = np.array(apply_replacement_rules(
                    X, *int_repl_rules
                ), dtype=float)[:, :, 0]
                np.testing.assert_array_almost_equal(
                    X_at_int_points_current,
                    X_at_int_points_abaqus,
                    decimal=6
                )


if __name__ == '__main__':
    unittest.main()
