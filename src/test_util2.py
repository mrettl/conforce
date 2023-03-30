import unittest

import ele_def
from util2 import *


class MyTestCase(unittest.TestCase):
    def test_reference_elements(self):
        for typ in ele_def.poly_power.keys():
            with self.subTest(typ):
                nodes = ele_def.ref_nodes[typ]
                powers = ele_def.poly_power[typ]
                int_points = ele_def.int_points[typ]
                int_weights = ele_def.int_weights[typ]
                int_scheme = IntegrationScheme(int_points, int_weights)

                ref_element = Element.create_ref_space_element(nodes, powers, int_scheme)

                np.testing.assert_array_almost_equal(nodes, ref_element.nodes)
                self.assertIsNot(nodes, ref_element.nodes)  # copied

                self.assertEqual(ref_element.coordinates, ref_element.reference_coordinates)
                self.assertEqual(ref_element.jacobian, sy.eye(ref_element.count_dimensions()))
                self.assertEqual(ref_element.inv_jacobian, sy.eye(ref_element.count_dimensions()))
                self.assertEqual(ref_element.det, 1.)

                # check if shape is 1 at corresponding node and 0 at other nodes
                shapes = sy.lambdify(
                    ref_element.reference_coordinates,
                    ref_element.reference_shapes
                )
                for i in range(len(nodes)):
                    for j in range(len(nodes)):
                        value = shapes(*nodes[i])[j]

                        if i == j:
                            self.assertAlmostEqual(1., value, 6)
                        else:
                            self.assertAlmostEqual(0., value, 6)

                # test methods
                self.assertEqual(ref_element.count_nodes(), len(nodes))
                random_nodes = np.random.random((5, ref_element.count_dimensions()))

                def function(_nodes):
                    if ref_element.count_dimensions() == 2:
                        return 1 + _nodes[:, 0] * 2. + _nodes[:, 1] * 1.5
                    else:
                        return 1 + _nodes[:, 0] * 2. + _nodes[:, 1] * 1.5 + _nodes[:, 2] * -0.5

                def diff_function_d_ref_coordinates(_nodes, axis):
                    return np.array([[2, 1.5, -0.5][axis]] * len(_nodes))

                interpolation = sy.lambdify(
                    ref_element.reference_coordinates,
                    ref_element(function(nodes))
                )

                interpolated_values = np.array([interpolation(*n) for n in random_nodes]).flatten()
                expected_values = function(random_nodes)
                np.testing.assert_array_almost_equal(interpolated_values, expected_values)

                for axis in range(0, ref_element.count_dimensions()):
                    diff = sy.lambdify(
                        ref_element.reference_coordinates,
                        ref_element.diff(axis, function(nodes))
                    )

                    diff_values = np.array([diff(*n) for n in random_nodes]).flatten()
                    expected_diff_values = diff_function_d_ref_coordinates(random_nodes, axis)
                    np.testing.assert_array_almost_equal(diff_values, expected_diff_values)

    def test_real_element(self):
        nodes = np.array([
            [0, 0],
            [1, 0],
            [1, 1],
            [0, 1]
        ])
        shape_powers = np.array([
            [0, 0],
            [1, 0],
            [0, 1],
            [1, 1],
        ])
        real_nodes = -1 + nodes * [2, 3] + 1. * nodes * nodes[:, ::-1]

        real_element = Element\
            .create_ref_space_element(nodes, shape_powers, None)\
            .create_real_space_element(real_nodes)

        np.testing.assert_array_almost_equal(real_element.nodes, real_nodes)
        r, s = real_element.reference_coordinates

        self.assertEqual(
            real_element.coordinates,
            sy.Matrix([
                [-1. + 2.*r + 1.*r*s],
                [-1. + 3.*s + 1.*r*s]
            ]))
        self.assertEqual(
            real_element.jacobian,
            sy.Matrix([
                [2. + 1.*s, 1.*r],
                [1.*s, 3. + 1.*r]
            ]))
        self.assertAlmostEqual(
            real_element.det.subs({r: 0., s: 0.}),
            2.*3.
        )

        dxy_dx = sy.simplify(real_element.diff(0))
        self.assertEqual(
            dxy_dx,
            sy.Matrix([[1.], [0.]])
        )

        dxy_dy = sy.simplify(real_element.diff(1))
        self.assertEqual(
            dxy_dy,
            sy.Matrix([[0.], [1.]])
        )

    def test_deformation(self):
        typ = "C3D8"
        nodes = ele_def.ref_nodes[typ]
        powers = ele_def.poly_power[typ]
        int_points = ele_def.int_points[typ]
        int_weight = ele_def.int_weights[typ]
        int_scheme = IntegrationScheme(int_points, int_weight)

        real_nodes = np.array([
            (0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
            (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)
        ])

        with self.subTest("deformation gradient"):
            expected_deformation_gradient = np.array([
                [1., 1., 0.],
                [0., 1., 0.],
                [0., 0., 1.]
            ])
            displacement = real_nodes @ expected_deformation_gradient.T - real_nodes

            ref_element = Element.create_ref_space_element(nodes, powers, int_scheme)
            real_element = ref_element.create_real_space_element(real_nodes)
            deformation = Deformation(
                element=real_element,
                internal_energy_density=1.,
                displacement=displacement,
                stress_tensor=np.eye(3)
            )
            np.testing.assert_array_almost_equal(
                deformation.deformation_gradient,
                expected_deformation_gradient
            )

        with self.subTest("bergstroem, mechanics of solid polymers, page 169"):
            expected_deformation_gradient = np.array([
                [2., 0., 0.],
                [0., 0.7, 0.],
                [0., 0., 0.7]
            ])
            displacement = real_nodes @ expected_deformation_gradient.T - real_nodes

            stress = sy.symbols("\\sigma", real=True)
            stress_tensor = sy.Matrix([
                [stress, 0, 0],
                [0, 0, 0],
                [0, 0, 0]
            ])

            ref_element = Element.create_ref_space_element(nodes, powers, int_scheme)
            real_element = ref_element.create_real_space_element(real_nodes)
            deformation = Deformation(
                element=real_element,
                internal_energy_density=1.,
                displacement=displacement,
                stress_tensor=stress_tensor
            )

            np.testing.assert_array_almost_equal(
                deformation.deformation_gradient,
                expected_deformation_gradient
            )

            piola_stress_tensor = sy.Matrix(deformation.piola_stress_tensor)
            self.assertAlmostEqual(
                0.,
                (0.49*stress - piola_stress_tensor[0, 0]) / stress
            )
            piola_stress_tensor[0, 0] = 0
            self.assertEqual(sy.zeros(3), piola_stress_tensor)

    def test_integration_scheme(self):
        typ = "CPE4"
        nodes = ele_def.ref_nodes[typ]
        powers = ele_def.poly_power[typ]
        int_points = ele_def.int_points[typ]
        int_weights = ele_def.int_weights[typ]
        int_scheme = IntegrationScheme(
            int_points,
            int_weights)

        real_nodes = np.array([
            [-1., -1.], [1., -1.],
            [2., 3.], [-1., 2.]
        ])

        ref_element = Element.create_ref_space_element(nodes, powers, None)
        real_element = ref_element.create_real_space_element(real_nodes)

        self.assertEqual(
            int_scheme.integrate(real_element.reference_coordinates, real_element.det),
            8.5
        )
        self.assertEqual(
            sy.Matrix(int_scheme.integrate(real_element.reference_coordinates, real_element.det, sy.Matrix([t, 2.]))),
            sy.Matrix([8.5*t, 2*8.5])
        )

    def test_configurational_force(self):
        assert False  # TODO


if __name__ == '__main__':
    unittest.main()
