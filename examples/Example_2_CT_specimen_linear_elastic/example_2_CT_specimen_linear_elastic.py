r"""
Example 2 - CT-specimen with linear elastic material behaviour
==============================================================

Kolednik [1]_

>>> w = 50  # mm
>>> a = 27  # mm
>>> b = 1  # mm
>>> b_n = 1  # mm
>>> nu = 0.3
>>> youngs_modulus = 200_000 # MPa

.. list-table:: Parameters
    :header-rows: 1

    * - Parameter
      - value
    * - :math:`w`
      - 50 mm
    * - :math:`a`
      - 27 mm
    * - :math:`H`
      - 60 mm
    * - :math:`\nu` (Poisson's ratio)
      - 0.3
    * - :math:`E` (Young's modulus)
      - 200 GPa
    * - :math:`u`
      - 0.5 mm

.. image:: example_2_images/00_scheme.png
    :width: 400
    :alt: scheme


Simulation
----------

>>> import os
>>> import numpy as np
>>> HOME_DIR = os.path.abspath(".")
>>> os.chdir(__file__ + "/..")


:py:mod:`example_2_abaqus_script`

>>> import subprocess
>>> if not os.path.exists("results.json"):
...     _ = subprocess.call("abaqus cae noGUI=example_2_abaqus_script.py", shell=True)

>>> import json
>>> with open("results.json", "r", encoding="UTF-8") as fh:
...     results = json.load(fh)


>>> load = results["reaction_force"]
>>> np.around(load, 3)  # N
56319.531

Theory
------

# TODO: wo finde ich die Norm?

>>> geometry_factor = (
...     (2 + a/w)
...     * (
...         0.866
...         + a/w * 4.64
...         - (a/w)**2 * 13.32
...         + (a/w)**3 * 14.72
...         - (a/w)**4 * 5.6
...     )
...     / (1 - a/w)**(3./2)
... )
>>> geometry_factor
10.821389667870234

# TODO: warum muss ich b=bn=25 mm einsetzen, wenn load in 2D N/mm !?!?!?!?

>>> k = (
...     load
...      * geometry_factor
...     / (b * b_n * w) ** 0.5
... )
>>> k
86190.03660880736

>>> j_integral_theory = k**2 * (1 - nu**2) / youngs_modulus
>>> j_integral_theory

Abaqus J-Integral
-----------------

>>> j_integral_abaqus = results["J"][57][-1]
>>> np.around(j_integral_abaqus, 3)
54.635

Configurational forces
----------------------

>>> (cfx_in_a, _) = results["CF"]["A"]
>>> np.around(cfx_in_a, 3)
-54.484

>>> (cfx_in_b, _) = results["CF"]["B"]
>>> np.around(cfx_in_b, 3)
-54.633

>>> (cfx_in_c, _) = results["CF"]["C"]
>>> np.around(cfx_in_c, 3)
-32.306

References
----------

.. [1] Kolednik, Otmar, Ronald Schöngrundner, and Franz Dieter Fischer.
    "A new view on J-integrals in elastic–plastic materials."
    International Journal of Fracture 187.1 (2014): 77-107.

Change to home directory

>>> os.chdir(HOME_DIR)
"""
import doctest

if __name__ == '__main__':
    doctest.testmod(verbose=True)
