import importlib
import unittest

from cf.codegen import *


class TestCodegen(unittest.TestCase):
    def test_something(self):
        name = "test_cf"
        folder = "simulations/tests/test_codegen"
        os.makedirs(folder, exist_ok=True)

        for ext in [".dll", ".so", ".py"]:
            file_path = os.path.join(folder, name + ext)
            if os.path.exists(file_path):
                os.remove(file_path)

        with CPyCodeCompiler(
                name=name,
                compile_at_exit=False,
                write_header_at_enter=True,
                folder=folder
        ) as compiler:
            write_code_for_element_type(
                element_type="CPE4",
                is_dbf=False,
                compiler=compiler
            )

        # importlib.import_module(name, f"{folder}/{name}".replace("/", "."))


if __name__ == '__main__':
    unittest.main()
