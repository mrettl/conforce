import unittest

import numpy as np

from cf.one_element_runner import simulate_one_element
from cf import element_definitions as el_def
from cf_shared.element_type_mapping import map_abaqus_element_type_to_supported_element_type as map_el_type
from cf_shared.cf_c import map_type_to_info


class TestElementTypeMapping(unittest.TestCase):
    def test_compatibility(self):
        for el_type, supported_el_type in map_el_type.items():
            with self.subTest(el_type):
                X_at_nodes = el_def.R_at_nodes_of_element[supported_el_type]
                X_at_nodes = X_at_nodes - np.min(X_at_nodes, axis=0)

                result = simulate_one_element(
                    X_at_nodes=X_at_nodes,
                    U_at_nodes=np.ones_like(X_at_nodes),
                    element_type=el_type,
                    load_name=f"test_compatibility",
                    folder="res/tests/test_element_type_mapping"
                )
                n, d = np.shape(result["nodes"]["RF"])
                ips, _ = np.shape(result["integration_points"]["SENER"])

                result_expected = simulate_one_element(
                    X_at_nodes=X_at_nodes,
                    U_at_nodes=np.ones_like(X_at_nodes),
                    element_type=supported_el_type,
                    load_name=f"test_compatibility",
                    folder="res/tests/test_element_type_mapping"
                )

                expected = map_type_to_info[supported_el_type]

                self.assertEqual(n, expected.number_of_nodes)
                self.assertEqual(d, expected.number_of_dimensions)
                self.assertEqual(ips, expected.number_of_integration_points)

                np.testing.assert_array_almost_equal(
                    result_expected["integration_points"]["COORD"],
                    result["integration_points"]["COORD"],
                    decimal=6
                )
                np.testing.assert_array_almost_equal(
                    result_expected["integration_points"]["IVOL"],
                    result["integration_points"]["IVOL"],
                    decimal=6
                )


if __name__ == '__main__':
    unittest.main()
