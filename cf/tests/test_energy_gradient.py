import unittest

import numpy as np

from cf.tests import one_element_abaqus_runner
from cf import element_definitions
from cf import cf_c


class TestEnergyGradient(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        cls.abaqus_data = one_element_abaqus_runner.simulate_all_element_types()


    def test_something(self):
        element_type = "C3D8"
        data = self.abaqus_data[element_type]
        origin = data["origin"]
        gradient = data["gradient"]

        n, d = np.shape(origin["nodes"]["COORD"])

        expected_CF = np.zeros((n, d), dtype=float)
        for derivative in gradient:
            dALLSE = np.array(derivative["dALLSE"])
            dX = np.array(derivative["dX"])
            i, j = np.argwhere(dX != 0)[0]
            dALLSE_dXij = dALLSE / dX[i, j]

            expected_CF[i, j] = dALLSE_dXij

        S = np.array(origin["integration_points"]["S"])
        S_at_integration_points = np.zeros((len(S), 3, 3), dtype=float)
        S_at_integration_points[:, 0, 0] = S[:, 0]
        S_at_integration_points[:, 1, 1] = S[:, 1]
        S_at_integration_points[:, 2, 2] = S[:, 2]
        S_at_integration_points[:, 0, 1] = S[:, 3]
        S_at_integration_points[:, 0, 2] = S[:, 4]
        S_at_integration_points[:, 1, 2] = S[:, 5]

        actual_CF = cf_c.compute_static_dbf_CF_for_multiple_C3D8(
            e=[origin["element"]["ESEDEN"]],
            X_at_nodes=[origin["nodes"]["COORD"]],
            U_at_nodes=[origin["nodes"]["U"]],
            S_at_int_points=[S_at_integration_points[:, :3, :3]]
        ).reshape((n, d))

        self.assertEqual(True, False)  # add assertion here


if __name__ == '__main__':
    unittest.main()
