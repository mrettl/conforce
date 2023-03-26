import unittest

import numpy as np
import sympy as sy

import util
import ele_def


r, s, t = sy.symbols("r s t")  # bild space
x, y, z = sy.symbols("x y z")  # real space


class MyTestCase(unittest.TestCase):
    def test_shapes(self):
        for typ in ele_def.poly_power.keys():
            with self.subTest(typ):
                points = ele_def.ref_nodes[typ]
                powers = ele_def.poly_power[typ]

                shapes = util.Shapes.from_points(points, powers)
                sym_shapes = shapes.to_symbolic(np.eye(len(points)), r, s, t)

                for i in range(len(points)):
                    for j in range(len(points)):
                        value = shapes.eval_shape(i, points[j])

                        if i == j:
                            self.assertAlmostEqual(1., value, 6)
                        else:
                            self.assertAlmostEqual(0., value, 6)

                    np.testing.assert_array_almost_equal(
                        sym_shapes(*points[j]).flatten(),
                        np.eye(len(points))[j]
                    )

    def test_gen_2d_ten_from_vec(self):
        vector = [1, 2, 3, 4, 5, 6]
        tensor = util.tensor_from_vector_notation(vector)

        self.assertEqual((3, 3), tensor.shape)

        vector_2 = util.vector_from_tensor_notation(tensor)

        self.assertEqual((6, 1), vector_2.shape)
        self.assertEqual(vector_2, sy.Matrix(vector))

    def test_to_symbolic(self):
        typ = "CPE4"
        points = ele_def.ref_nodes[typ]
        powers = ele_def.poly_power[typ]
        i_points = ele_def.int_points[typ]
        i_weights = ele_def.int_weights[typ]
        shapes = util.Shapes.from_points(points, powers)

        a = util.get_jac(
            shapes.coef,
            points,
            r, s, t
        )

        b = shapes.to_symbolic(r, s, t)

        print("ok")


if __name__ == '__main__':
    unittest.main()
