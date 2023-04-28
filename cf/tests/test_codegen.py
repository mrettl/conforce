import unittest
import shutil

from cf.codegen import *

folder = "cf/tests/generated"


class TestCodegen(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        shutil.rmtree(folder)

    def test_compilation(self):
        name = "test_codegen"
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
                element_type="CPE4",
                is_dbf=False,
                compiler=compiler
            )

        from cf.tests.generated import test_codegen
        # TODO: call functions
        print("ok")


if __name__ == '__main__':
    unittest.main()
