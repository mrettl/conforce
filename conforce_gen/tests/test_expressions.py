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
from conforce import cf_c


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

    def test_issue(self):
        el_type = "CPE4"
        X_at_nodes = np.array([
            [0.671303629875183, 23.3493576049805],  # 37
            [0.669628739356995, 23.3450508117676],  # 4266
            [0.675243616104126, 23.3419837951660],  # 5516
            [0.678554773330688, 23.3468837738037],  # 4447
        ])
        U_at_nodes = np.array([
            [-0.15737548, 0.0021416706],  # 37
            [-0.15733641, 0.0021185826],  # 4266
            [-0.15731329, 0.0021713173],  # 5516
            [-0.15735458, 0.0022020650],  # 4447
        ])
        S_at_integration_points = np.array([
            [[3.10283732414246, 17.6723880767822], [17.6723880767822, 267.967407226562]],
            [[4.19220495223999, 77.5216522216797], [77.5216522216797, 266.878051757812]],
            [[75.1460189819336, 21.1991348266602], [21.1991348266602, 195.924224853516]],
            [[82.9048919677734, 73.6524200439453], [73.6524200439453, 188.16535949707]]
        ])
        e_at_integration_points = np.array([
            0.159019857645035,
            0.19316029548645,
            0.0722203999757767,
            0.098083920776844
        ])

        # expected values: computed by an alternative implementation
        F_0_expected = np.array([[0.999775991585619, -0.0088190636027492], [0.00905579548622659, 1.00144575849556]])
        P_0_expected = np.array([[3.26317719, 17.64033065], [20.06115969, 267.74734274]])
        CS_dbf_0_expected = np.array([[-0.0219189225734775, -2.42071359533884], [-0.000225424848754423, -0.0725069397837795]])

        F_1_expected = np.array([[0.999508925913042, -0.0086973048636878], [0.00968958831308582, 1.00115680404222]])
        P_1_expected = np.array([[4.8712839547833, 77.4429626069841], [79.9324493598944, 265.995841966838]])
        CS_dbf_1_expected = np.array([[-0.578960070347097, -2.53935996949774], [-0.0500991388920703, 0.559000285624483]])

        F_2_expected = np.array([[0.999900593692058, -0.00848802828964601], [0.00876009323249827, 1.00066015469683]])
        P_2_expected = np.array([[75.3755658354376, 20.5387413666026], [22.8761398982684, 195.719042352157]])
        CS_dbf_2_expected = np.array([[-0.120683911623952, -1.71247537793113], [0.624688143960002, 0.117348972660648]])

        F_3_expected = np.array([[0.999679862066307, -0.00835093860164093], [0.00928392738801309, 1.00033481726182]])
        P_3_expected = np.array([[83.5477167943428, 72.8591581132411], [75.2484375096666, 187.421336908049]])
        CS_dbf_3_expected = np.array([[-0.573770315705013, -1.71668110248966], [0.672507377453637, 0.643774377917609]])

        CF_dbf_37_expected = np.array([-0.00999102603740184, 0.000419372408536019])
        CF_dbf_4266_expected = np.array([0.00643790591343163, -0.00110200119515269])
        CF_dbf_5516_expected = np.array([0.00797998576036596, -0.00163124659990473])
        CF_dbf_4447_expected = np.array([-0.00442686563639576, 0.0023138753865214])

        # test deformation gradients
        ((F_0_actual, F_1_actual, F_2_actual, F_3_actual),) = cf_c.compute_F_for_CPE4(
            X_at_nodes=(X_at_nodes,),
            U_at_nodes=(U_at_nodes,)
        )
        np.testing.assert_array_almost_equal(F_0_actual, F_0_expected, decimal=5)
        np.testing.assert_array_almost_equal(F_1_actual, F_1_expected, decimal=5)
        np.testing.assert_array_almost_equal(F_2_actual, F_2_expected, decimal=5)
        np.testing.assert_array_almost_equal(F_3_actual, F_3_expected, decimal=5)

        # test First Piola-Kirchhoff stress tensor
        ((P_0_actual, P_1_actual, P_2_actual, P_3_actual),) = cf_c.compute_P_for_CPE4(
            X_at_nodes=(X_at_nodes,),
            U_at_nodes=(U_at_nodes,),
            S_at_int_points=(S_at_integration_points,)
        )
        np.testing.assert_array_almost_equal(P_0_actual, P_0_expected, decimal=3)
        np.testing.assert_array_almost_equal(P_1_actual, P_1_expected, decimal=3)
        np.testing.assert_array_almost_equal(P_2_actual, P_2_expected, decimal=3)
        np.testing.assert_array_almost_equal(P_3_actual, P_3_expected, decimal=3)

        # test configurational stress
        ((CS_dbf_0_actual, CS_dbf_1_actual, CS_dbf_2_actual, CS_dbf_3_actual),) = cf_c.compute_CS_for_CPE4_using_dbf(
            e_at_int_points=(e_at_integration_points,),
            X_at_nodes=(X_at_nodes,),
            U_at_nodes=(U_at_nodes,),
            S_at_int_points=(S_at_integration_points,)
        )
        np.testing.assert_array_almost_equal(CS_dbf_0_actual, CS_dbf_0_expected, decimal=3)
        np.testing.assert_array_almost_equal(CS_dbf_1_actual, CS_dbf_1_expected, decimal=3)
        np.testing.assert_array_almost_equal(CS_dbf_2_actual, CS_dbf_2_expected, decimal=3)
        np.testing.assert_array_almost_equal(CS_dbf_3_actual, CS_dbf_3_expected, decimal=3)

        # test configurational force
        ((CF_dbf_37_actual, CF_dbf_4266_actual, CF_dbf_5516_actual, CF_dbf_4447_actual),) = cf_c.compute_CF_for_CPE4_using_dbf(
            e_at_int_points=(e_at_integration_points,),
            X_at_nodes=(X_at_nodes,),
            U_at_nodes=(U_at_nodes,),
            S_at_int_points=(S_at_integration_points,)
        )
        np.testing.assert_array_almost_equal(CF_dbf_37_actual, CF_dbf_37_expected, decimal=3)
        np.testing.assert_array_almost_equal(CF_dbf_4266_actual, CF_dbf_4266_expected, decimal=3)
        np.testing.assert_array_almost_equal(CF_dbf_5516_actual, CF_dbf_5516_expected, decimal=3)
        np.testing.assert_array_almost_equal(CF_dbf_4447_actual, CF_dbf_4447_expected, decimal=3)



if __name__ == '__main__':
    unittest.main()
