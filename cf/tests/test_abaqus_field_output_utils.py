import unittest

from cf.tests.abaqus_tensile_specimen_script import build_2d, build_3d
from cf.abaqus_field_output_utils import add_field_outputs


class TestAbaqusFieldOutputUtils(unittest.TestCase):
    def test_something(self):

        # TODO:
        self.assertEqual(True, False)  # add assertion here


if __name__ == '__main__':
    unittest.main()
