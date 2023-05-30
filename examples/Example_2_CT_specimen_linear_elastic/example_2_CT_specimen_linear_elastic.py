r"""
Example 2 - CT-specimen with linear elastic material behaviour
==============================================================

Problem description
-------------------

Kolednik [1]_ suggested the following compact tension (CT) test.
Contrary to Kolednik, in this example a linear elastic material behaviour is used
instead of an elastic plastic behaviour.
This allows a comparisson with theoertical prediction for CT specimens by Anderson [2]_.

.. image:: example_2_images/00_scheme.png
    :width: 400
    :alt: scheme

Geometeric dimensions:
^^^^^^^^^^^^^^^^^^^^^^

>>> w = 50  # mm
>>> a = 27  # mm
>>> b = 25  # mm
>>> h = 60  # mm

Material properties:
^^^^^^^^^^^^^^^^^^^^

>>> nu = 0.3
>>> youngs_modulus = 200_000 # MPa

Boundary condition:
^^^^^^^^^^^^^^^^^^^

>>> u = 0.5  # mm

Simulation
----------

The model is scripted and there is no need to manually open the `conforce`-plugin.
First, change to the directory where the \*.inp file is located.

>>> import os
>>> import numpy as np
>>> HOME_DIR = os.path.abspath(".")
>>> os.chdir(__file__ + "/..")

Next, call the Abaqus script (:py:mod:`example_2_abaqus_script`).
The script simulates the \*.inp file and writes a `results.json` file.
To save time, the script is not exectued if the `results.json` file already exists.

>>> import subprocess
>>> if not os.path.exists("results.json"):
...     _ = subprocess.call("abaqus cae noGUI=example_2_abaqus_script.py", shell=True)

After the simulation, the result file is loaded to fetch the reaction force `load`.
This force is needed to reach the defined displacement `u`.

>>> import json
>>> with open("results.json", "r", encoding="UTF-8") as fh:
...     results = json.load(fh)
>>> load = results["reaction_force"]
>>> np.around(load, 3)  # N
56319.531

Theory
------

Anderson [2]_ defines an analytical solution for this problem.

>>> geometry_factor = (
...     (2 + a/w)
...     * (
...         0.866
...         + 4.64 * (a/w)
...         - 13.32 * (a/w)**2
...         + 14.72 * (a/w)**3
...         - 5.6 * (a/w)**4
...     )
...     / (1 - a/w)**(3./2)
... )
>>> geometry_factor
10.821389667870234

The fracture thoughness is:

>>> k = (
...     load  # N
...     / (b * w**0.5)  # 1/mm**(3/2)
...      * geometry_factor
... )  # N / mm**(3/2) = MPa mm**0.5
>>> np.around(k, 3)  # MPa mm**0.5
3447.601

For a plane strain state, the J-Integral is computed by:

>>> j_integral_theory = k**2 * (1 - nu**2) / youngs_modulus
>>> np.around(j_integral_theory, 3)  # mJ/mm**2
54.081

Abaqus J-Integral
-----------------

Abaqus comptues the J-Integral using the virtual crack extension method by Parks [3]_.
According to Abaqus J-Integral over the region `B` is:

>>> j_integral_abaqus = results["J"][57][-1]
>>> np.around(j_integral_abaqus, 3)  # mJ/mm**2
54.635

This is in good aggreement with the prediction made by Anderson.

Configurational forces
----------------------

Finally. the configurational forces computed by `conforce` are summed up
for the regions `A`, `B`, and `C`.

.. note::
    The configurational forces have a negative sign.

Only the x-component of the configurational force is considered.
The configurational force for region `A` shows a good aggreement with Anderson.

>>> (cfx_in_a, _) = results["CF"]["A"]
>>> np.around(cfx_in_a, 3)
-54.484

The configurational force for region `B` shows an excellent aggreement to the Abaqus results
as they are computed both in the same domain `B`.

>>> (cfx_in_b, _) = results["CF"]["B"]
>>> np.around(cfx_in_b, 3)
-54.633

The configurational force for region `C` is much smaller than predicted by Abaqus and Anderson.
The configurational force is the deviation of the energy if the whole region is shifted rightwards.
For small regions, just the domain near the crack tip is shifted and according to the crack similarity assumption
this cooresponds to the deviation of the energy if the crack grows.
However, the region `C` contains almost the whole specimen.
Consequently, the assumption that the stress field is dominated only by the crack tip is wrong,
as the boundary condition have also a considerable influence.
Hence, the crack similarity assumptioin is invalid and the configurational force of region `C` is underestimated.

>>> (cfx_in_c, _) = results["CF"]["C"]
>>> np.around(cfx_in_c, 3)
-32.306

Like conventional forces, the configurational forces fullfill the equilibrium of forces.
The configuraional forces summed up for the whole model are zero.

>>> (cfx_whole_model, _) = results["CF"][" ALL NODES"]
>>> np.around(cfx_whole_model, 3)
-0.0

References
----------

.. [1] Kolednik, Otmar, Ronald Schöngrundner, and Franz Dieter Fischer.
    "A new view on J-integrals in elastic-plastic materials."
    International Journal of Fracture 187.1 (2014): 77-107.

.. [2] Anderson, T. L., and T. Anderson.
    "Fracture mechanics: fundamentals and applications. 2005."
    CRC press, Taylor and Francis Group, ISBN 10 (2005): 0-8493.

.. [3] Parks, D. M.
    "The virtual crack extension method for nonlinear material behavior."
    Computer methods in applied mechanics and engineering 12.3 (1977): 353-364.

Change to home directory

>>> os.chdir(HOME_DIR)
"""
import doctest

if __name__ == '__main__':
    doctest.testmod(verbose=True)
