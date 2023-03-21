import unittest

import numpy as np

import util
import ele_def


class MyTestCase(unittest.TestCase):
    def test_shapes_from_points(self):
        for typ in ele_def.poly_power.keys():
            with self.subTest(typ):
                points = ele_def.bild_points[typ]
                powers = ele_def.poly_power[typ]

                shapes = util.Shapes.from_points(points, powers)

                for i in range(len(points)):
                    for j in range(len(points)):
                        value = shapes(i, points[j])

                        if i == j:
                            self.assertAlmostEqual(1., value, 6)
                        else:
                            self.assertAlmostEqual(0., value, 6)

    def test_generate_shape_func_coeff(self):
        util.generate_shape_func_coeff(
            np.array([
                [-1, -1, 0],
                [1, -1, 0],
                [1, 1, 0],
                [-1, 1, 0]
            ]),
            3*np.array([
                [-1, -1, 0],
                [1, -1, 0],
                [1, 1, 0],
                [-1, 1, 0]
            ]),
            np.array([
                [0, 0, 0],
                [1, 0, 0],
                [0, 1, 0],
                [1, 1, 0]
            ])
        )
        print("ok")


if __name__ == '__main__':
    unittest.main()
