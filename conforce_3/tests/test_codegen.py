import unittest
import shutil

from conforce_3 import element_definitions as el_def
from conforce_3.codegen import *

folder = "conforce_3/tests/generated"


class TestCodegen(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        if os.path.exists(folder):
            shutil.rmtree(folder)

    def test_compilation(self):
        name = "test_codegen"
        element_type = "CPE4R"
        os.makedirs(folder, exist_ok=True)

        for ext in [".dll", ".so", ".py"]:
            file_path = os.path.join(folder, name + ext)
            if os.path.exists(file_path):
                os.remove(file_path)

        with CPyCodeCompiler(
                name=name,
                compile_at_exit=True,
                write_header_at_enter=True,
                folder=folder
        ) as compiler:
            write_code_for_element_type(
                element_type=element_type,
                is_dbf=False,
                write_F=True,
                write_P=True,
                write_CS=True,
                write_CF=True,
                compiler=compiler
            )

        from conforce_3.tests.generated import test_codegen

        X_at_nodes = el_def.R_at_nodes_of_element[element_type]
        U_at_nodes = 2 * X_at_nodes

        F_at_int_point = test_codegen.compute_F(
            [X_at_nodes],
            [U_at_nodes],
            element_type
        )[0, 0]
        F_at_int_point_expected = 2 * np.eye(2) + np.eye(2)

        np.testing.assert_array_almost_equal(
            F_at_int_point,
            F_at_int_point_expected
        )

        self.assertEqual(1, test_codegen.map_type_to_info["CPE4R"].number_of_integration_points)
        self.assertEqual(2, test_codegen.map_type_to_info["CPE4R"].number_of_dimensions)
        self.assertEqual(4, test_codegen.map_type_to_info["CPE4R"].number_of_nodes)


if __name__ == '__main__':
    unittest.main()
