import unittest

import numpy as np

import cf.symbolic_util
from cf import element_definitions as el_def, one_element_runner
from cf import cf_c
from cf import expressions as exp
import Conf_Forces_py as cf_m


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
        S_at_integration_points[:, 0, 1] = S_at_integration_points[:, 1, 0] = S[:, 3]
        S_at_integration_points[:, 0, 2] = S_at_integration_points[:, 2, 0] = S[:, 4]
        S_at_integration_points[:, 1, 2] = S_at_integration_points[:, 2, 1] = S[:, 5]

        a = cf_c.compute_static_dbf_CF_for_multiple_C3D8(
            e_at_int_points=[np.array(origin["integration_points"]["SENER"]).flatten()],
            X_at_nodes=[origin["nodes"]["COORD"]],
            U_at_nodes=[origin["nodes"]["U"]],
            S_at_int_points=[np.zeros_like(S_at_integration_points[:, :3, :3])]
        ).reshape((n, d))

        actual_CF = cf_c.compute_static_mbf_CF_for_multiple_C3D8(
            e_at_int_points=[np.array(origin["integration_points"]["SENER"]).flatten()],
            X_at_nodes=[origin["nodes"]["COORD"]],
            U_at_nodes=[origin["nodes"]["U"]],
            S_at_int_points=[S_at_integration_points[:, :3, :3]]
        ).reshape((n, d))

        markus_CF = cf_m.calc_Conf_Force_C3D8_static(
            Coords=np.array([origin["nodes"]["COORD"]]),
            Element_U=np.array([origin["nodes"]["U"]]),
            S_vec=np.array([S]),
            PENER=np.array([[0.]*8]),
            SENER=np.array([[origin["element"]["ESEDEN"]]*8]),
            method="dbf"
        )
        #
        R = cf.math_util.create_symbolic_matrix("{row}", ["r", "s", "t"], 1)
        S = cf.math.create_symbolic_matrix("S{row}{col}", 3, 3, *R, is_symmetric=True)
        H_ = exp.eval_H(
            R,
            el_def.R_at_nodes_of_element[element_type],
            el_def.exponents_of_shape_functions_of_element[element_type]).doit()

        dH_dR_ = exp.eval_dH_dR(H_, R).doit()
        J_ = exp.eval_J(np.array(origin["nodes"]["COORD"]), dH_dR_).doit()
        dH_dX_ = exp.eval_dH_dX(dH_dR_, J_).doit()
        dU_dX_ = exp.eval_dU_dX(np.array(origin["nodes"]["U"]), dH_dX_).doit()
        F_ = exp.eval_F(3, dU_dX_).doit()
        P_ = exp.eval_P(F_, S)
        P_at_integration_points = [
            exp.eval_P(F_, S_at_integration_point)
            for S_at_integration_point in S_at_integration_points
        ]

        e = origin["element"]["ESEDEN"] * np.eye(3, 3)
        CS = e - F_.T * P_at_integration_points[0]

        self.assertEqual(True, False)  # add assertion here


if __name__ == '__main__':
    unittest.main()
