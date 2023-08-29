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

    def test_surface_integral(self):
        element_type = "CPE3"
        nodes = R_at_nodes_of_element[element_type]
        exponents = exponents_of_shape_functions_of_element[element_type]

        n, d = nodes.shape
        R = eval_R(d)
        H = eval_H(R, nodes, exponents)

        X_at_nodes = np.array([
            [np.cos(angle_rad), np.sin(angle_rad)]
            for angle_rad in np.deg2rad([90, 210, 330])
        ])

        import sympy as sy

        t = sy.symbols("t", real=True)

        r0 = 0*(1-t) + 1*t
        r1 = 0
        H_s1 = H.xreplace({R[0]: r0, R[1]: r1})
        R_s1 = nodes.T @ H_s1

        content = sy.integrate(1, (t, 0, 1))
        first_momentum = sy.integrate(R_s1, (t, 0, 1))

        t_ip = first_momentum / content
        t_weight = content

        X_s1 = X_at_nodes.T @ H_s1
        dx, dy = X_s1.replace(t, 1) - X_s1.replace(t, 0)
        lt_real = (dx**2 + dy**2)**0.5

        tangential_vector = X_s1.diff(t).replace(t, 0.5)
        tangential_vector = sy.Matrix([
            tangential_vector[0],
            tangential_vector[1],
            0
        ])
        tangential_vector_2 = sy.Matrix([
            0, 0, 1
        ])

        normal_vector = tangential_vector.cross(tangential_vector_2)
        normal_vector = normal_vector / (normal_vector[0]**2 + normal_vector[1]**2)**0.5


        integral = normal_vector * H_s1[0].replace(t, 0.5) * lt_real

        print("ok")
        self.assertTrue(False)

    def test_surface_integral2(self):
        import sympy as sy

        element_type = "C3D4"
        nodes = R_at_nodes_of_element[element_type]
        exponents = exponents_of_shape_functions_of_element[element_type]

        n, d = nodes.shape
        R = eval_R(d)
        H = eval_H(R, nodes, exponents)

        X_at_nodes = 2 * nodes

        # reference face coordinate system
        v, w = sy.symbols("v w", real=True)

        p0 = nodes[0, :]
        p1 = nodes[1, :]
        p2 = nodes[2, :]

        p = sy.Matrix(
            (1 - v) * p0 + v * p1
            + (1 - w) * p1 + w * p2
        )

        # orthogonality
        vec_0 = p0

        vec_v = p1 - p0
        vec_v = vec_v / (vec_v @ vec_v)**0.5

        vec_w = p2 - p0
        vec_w - (vec_v @ vec_w) / (vec_v @ vec_v) * vec_v
        vec_w = vec_w / (vec_w @ vec_w)**0.5

        # project points to v, w base
        pvw0, pvw1, pvw2 = [
            np.array([
                (p - vec_0) @ vec_v,
                (p - vec_0) @ vec_w,
            ])
            for p in [p0, p1, p2]
        ]

        # area
        vw01_mean_w = (pvw1[1] + pvw0[1]) / 2
        vw01_dv = pvw1[0] - pvw0[0]

        vw12_mean_w = (pvw2[1] + pvw1[1]) / 2
        vw12_dv = pvw2[0] - pvw1[0]

        vw20_mean_w = (pvw0[1] + pvw2[1]) / 2
        vw20_dv = pvw0[0] - pvw2[0]

        area = -(
            vw01_mean_w * vw01_dv
            + vw12_mean_w * vw12_dv
            + vw20_mean_w * vw20_dv
        )

        # line variable
        s = sy.symbols("s", real=True)

        def int_test(pvwi, pvwj):
            lvw_ji = sy.Matrix(pvwj * (1 + s) / 2 + pvwi * (1 - s) / 2)
            ipa_ji = lvw_ji.replace(s, -3 ** -0.5)
            ipb_ji = lvw_ji.replace(s, 3 ** -0.5)
            lvwji_d = lvw_ji.replace(s, 1) - lvw_ji.replace(s, -1)

            int_vdw_ij = (ipa_ji[0] + ipb_ji[0]) * lvwji_d[1] / 2

            int_v2dw_ij = (ipa_ji[0]**2 + ipb_ji[0]**2) * lvwji_d[1] / 2
            int_w2dv_ij = (ipa_ji[1] ** 2 + ipb_ji[1] ** 2) * lvwji_d[0] / 2

            return np.array([
                int_vdw_ij,
                int_v2dw_ij,
                int_w2dv_ij
            ])

        int_vdw, int_v2dw, int_w2dv = (
            int_test(pvw0, pvw1)
            + int_test(pvw1, pvw2)
            + int_test(pvw2, pvw0)
        )

        area = int_vdw
        vs = 1/(2*area)*int_v2dw
        ws = -1/(2*area)*int_w2dv

        # curve parametrization
        R_of_vw = sy.Matrix(
            vec_0
            + vec_v * v
            + vec_w * w
        )

        X = X_at_nodes.T @ H  # type: sy.MatrixExpr
        X_of_vw = X.xreplace({
            R[0]: R_of_vw[0],
            R[1]: R_of_vw[1],
            R[2]: R_of_vw[2]
        })  # type: sy.MatrixExpr

        dX_dv = X_of_vw.diff(v)
        dX_dw = X_of_vw.diff(w)
        dX_dvw = sy.Matrix(np.column_stack([dX_dv, dX_dw]))
        det = sy.sqrt(sy.det(dX_dvw.T @ dX_dvw))

        # normal vector
        N = dX_dv.cross(dX_dw)
        N = N / N.norm()

        #
        CS = sy.Matrix([
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]
        ])
        (CS @ N * H[0] * det).xreplace({
            R[0]: R_of_vw[0],
            R[1]: R_of_vw[1],
            R[2]: R_of_vw[2]
        }).xreplace({
            v: vs,
            w: ws
        }) * area

        print("ok")
        self.assertTrue(False)





if __name__ == '__main__':
    unittest.main()
