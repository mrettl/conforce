import unittest

import numpy as np
import sympy as sy

import auxiliary
import ele_def


r, s, t = sy.symbols("r s t")  # bild space
x, y, z = sy.symbols("x y z")  # real space


class MyTestCase(unittest.TestCase):
    def test_gen_Configurational_Forces_Static(self):
        typ = "CPE4"
        points = ele_def.bild_points[typ]
        powers = ele_def.poly_power[typ]
        int_points = ele_def.int_points[typ]
        int_weights = ele_def.int_weights[typ]

        a = auxiliary.gen_Configurational_Forces_Static(
            powers,
            points,
            int_points,
            int_weights,
            typ
        )
        print("ok")


if __name__ == '__main__':
    unittest.main()
