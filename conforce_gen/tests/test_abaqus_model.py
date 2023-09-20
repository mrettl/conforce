import os
import subprocess
import json
import unittest


class TestModel(unittest.TestCase):
    def test_abaqus_model(self):
        home = os.path.abspath(".")
        try:
            os.chdir("conforce_gen/tests/abaqus_model")
            if os.path.exists("results.json"):
                os.remove("results.json")

            self.assertEqual(
                subprocess.call([
                    "abaqus",
                    "cae",
                    r"noGUI=abaqus_script.py"
                ], shell=True),
                0
            )

            with open("results.json", "r", encoding="UTF-8") as fh:
                results = json.load(fh)

        finally:
            os.chdir(home)

        expected_results = {
            'CF_INSTANCE': [1.2095028068870306e-05, -4.2431202018633485e-06],
            'CF_ALL': [5.787210102425888e-05, 1.483038067817688e-05],
            'CF_NSET': [2.874352503567934e-06, 1.3584576663561165e-06],
            'CF_ELSET': [-6.5695494413375854e-06, 2.4996697902679443e-06]
        }

        self.assertDictEqual(results, expected_results)


    if __name__ == '__main__':
        unittest.main()
