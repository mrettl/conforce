import unittest

import numpy as np

from conforce_gen import element_definitions as el_def
from conforce import cf_c
from conforce_gen.one_element_runner import simulate_one_element
from conforce.tensor_util import tensor_from_abaqus_notation


class TestEnergyGradient(unittest.TestCase):
    def test_energy_gradient(self):
        STRAIN = 1e-3*np.array([
            [1., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.]
        ])
        delta = 1e-3

        for element_type in el_def.R_at_nodes_of_element.keys():
            with self.subTest(element_type):
                X_at_nodes = el_def.R_at_nodes_of_element[element_type]
                n, d = X_at_nodes.shape

                origin = simulate_one_element(
                    X_at_nodes=X_at_nodes,
                    U_at_nodes=X_at_nodes @ STRAIN[:d, :d],
                    element_type=element_type,
                    load_name="origin",
                    folder=f"res/tests/test_energy_gradient/{element_type}"
                )

                increments = list()
                for i in range(n):
                    increment = list()
                    increments.append(increment)
                    for j in range(d):
                        X_at_nodes_inc = np.array(X_at_nodes)
                        X_at_nodes_inc[i, j] += delta

                        increment.append(simulate_one_element(
                            X_at_nodes=X_at_nodes_inc,
                            U_at_nodes=X_at_nodes_inc @ STRAIN[:d, :d],
                            element_type=element_type,
                            load_name=f"increment_{i}_{j}",
                            folder=f"res/tests/test_energy_gradient/{element_type}"
                        ))

                dALLSE_dX = np.array([
                    [
                        (
                                increment_axis["model"]["ALLSE"]
                                - origin["model"]["ALLSE"]
                        ) / delta
                        for increment_axis in increment
                    ]
                    for increment in increments
                ])

                CF_at_nodes = cf_c.compute_CF(
                    e_at_int_points=np.array(origin["integration_points"]["SENER"]).reshape((1, -1)),
                    X_at_nodes=[X_at_nodes],
                    U_at_nodes=[X_at_nodes @ STRAIN[:d, :d]],
                    S_at_int_points=[tensor_from_abaqus_notation(origin["integration_points"]["S"])[:, :d, :d]],
                    element_type=element_type,
                    method="dbf"
                )[0]

                np.testing.assert_array_almost_equal(
                    dALLSE_dX,
                    CF_at_nodes,
                    decimal=2
                )


if __name__ == '__main__':
    unittest.main()
