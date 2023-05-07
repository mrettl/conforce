import unittest

import sympy as sy

from cf.element_definitions import *
from cf import one_element_runner
from cf.expressions import eval_H, eval_dH_dR, eval_J, R_3d
from cf.symbolic_util import create_replacement_rules, apply_replacement_rules


class TestElementDefinitions(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        cls.abaqus_data = one_element_abaqus_runner.simulate_all_element_types()

    def test_validation_against_abaqus(self):
        for element_type, data in self.abaqus_data.items():
            with self.subTest(element_type):
                # current
                R_at_nodes_current = R_at_nodes_of_element[element_type]
                R_at_int_points_current = R_at_integration_points_of_element[element_type]
                weights_of_int_points_current = weights_of_integration_points_of_element[element_type]
                exponents = exponents_of_shape_functions_of_element[element_type]

                # abaqus data
                model = data["origin"]
                X_at_nodes_abaqus = np.array(model["nodes"]["COORD"])
                X_at_int_points_abaqus = np.array(model["integration_points"]["COORD"])
                vol_abaqus = model["element"]["EVOL"]

                # computation
                n, d = X_at_nodes_abaqus.shape
                R = sy.Matrix(R_3d[:d])
                H = eval_H(R, R_at_nodes_current, exponents)
                X = X_at_nodes_abaqus.T @ H
                dH_dR = eval_dH_dR(H, R)
                J = eval_J(X_at_nodes_abaqus, dH_dR)
                detJ = J.det()

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

                # TODO: compare IVOL with integration weights


if __name__ == '__main__':
    unittest.main()
