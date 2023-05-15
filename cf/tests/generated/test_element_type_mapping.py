import unittest

import numpy as np

from cf_shared.element_type_mapping import map_abaqus_element_type_to_supported_element_type
from cf import element_definitions as el_def
from cf.one_element_runner import simulate_one_element


class TestElementMapping(unittest.TestCase):
    def test_compatibility(self):
        for abq_el_type, supported_el_type in map_abaqus_element_type_to_supported_element_type.items():
            with self.subTest(abq_el_type):
                data = simulate_one_element(
                    X_at_nodes=el_def.R_at_nodes_of_element[supported_el_type],
                    U_at_nodes=np.ones_like(el_def.R_at_nodes_of_element[supported_el_type]),
                    element_type=supported_el_type,
                    load_name="translation",
                    folder="res/tests/test_element_type_mapping"
                )

                np.testing.assert_array_almost_equal(
                    el_def.R_at_nodes_of_element[supported_el_type],
                    np.array(data["nodes"]["COORD"])
                )
                np.testing.assert_array_almost_equal(
                    el_def.R_at_integration_points_of_element[supported_el_type],
                    np.array(data["integration_points"]["COORD"])
                )
                np.testing.assert_array_almost_equal(
                    el_def.weights_of_integration_points_of_element[supported_el_type],
                    np.array(data["integration_points"]["IVOL"]).flatten()
                )


if __name__ == '__main__':
    unittest.main()
