import unittest

import numpy as np

import Conf_Forces_py as cf


class MyTestCase(unittest.TestCase):
    def test_resulting_nodal_forces(self):
        element_connectivity = [
            [1, 2, 3],
            [2, 4, 5, 3]
        ]
        element_nodal_forces = [
            [[0, -1, 0], [1, 0, 0], [-1, -1, 0]],
            [[-1, 0, 0], [0, 0, 0], [0, 1, 0], [-1, -1, 0]]
        ]

        node_labels, node_forces = cf.resulting_nodal_forces(
            element_connectivity,
            element_nodal_forces)

        actual_result = {
            label: tuple(force)
            for label, force
            in zip(node_labels, node_forces)
        }
        expected_result = {
            1: (0., -1., 0.),
            2: (0., 0., 0.),
            3: (-2., -2., 0.),
            4: (0., 0., 0.),
            5: (0., 1., 0.),
        }

        self.assertDictEqual(expected_result, actual_result)

        self.assertEqual(True, False)  # add assertion here


if __name__ == '__main__':
    unittest.main()
