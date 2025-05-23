r"""
Problem description
-------------------

:ref:`Figure 1 <example_6_model_image>` depicts the model considered in this example.
We are interested in the gradient of the strain energy with respect to the position of the cylindrical hole.
Configurational forces are used to estimate this energy gradient.
The configurational forces are computed using ConForce.
To validate this approach, the energy gradient is estimated numerically by a difference quotient,
for which two additional FEM simulations with slightly displaced holes are run.

.. _example_6_model_image:

.. figure:: example_6_images/model.png
    :width: 400
    :alt: model

    Figure 1: 3D model with a cylindrical hole.
    The left side is fixed. A displacement :math:`u_{y}` is applied on the right side.
    The regions with radii :math:`r_{1}` and :math:`r_{2}` are used to compute
    the configurational forces. The dimensions are listed in section :ref:`Geometric dimensions <example_6_dimensions>`.

The mesh is shown in :ref:`Figure 2 <example_6_mesh_image>`.
We use 20-node bricks with reduced integration.
The image shows a red-, green-, and cream-colored region.
The cream-colored region has the radius :math:`r_{1}` and the green region has the radius :math:`r_{2}`.
The configurational forces are computed for all nodes inside two domains:

1. The smaller domain :math:`\mathcal{A}` contains all nodes of the cream-colored region.
2. The larger domain :math:`\mathcal{B}` contains all nodes of the green- and cream-colored region.

.. _example_6_mesh_image:

.. figure:: example_6_images/mesh.png
    :width: 400
    :alt: mesh

    Figure 2: The FEM mesh consists of 20-node brick elements with reduced integration (C3D20R).

Simulation
----------

>>> import os
>>> import json
>>> import subprocess
>>> import numpy as np
>>> import scipy.optimize as optimize
>>> HOME_DIR = os.path.abspath(".")
>>> os.chdir(os.path.abspath(__file__ + "/.."))

The model is automatically build and simulated by an Abaqus script.
We provide the name of the model and the position of the hole by cx and cy.

>>> model_name = "bending_model_0"
>>> cx = 50  # mm
>>> cy = 40  # mm
>>> subprocess.call(f"abaqus cae noGui=example_6_abaqus_script.py -- {model_name} {cx:.3f} {cy:.3f}", shell=True)
0

The Abaqus script writes results to a json file.
The command below loads the results.

>>> with open(f"{model_name}_results.json", "r", encoding="UTF-8") as fh:
...     results = json.load(fh)


Material
^^^^^^^^

The compressible Neo-Hookean material model for hyperelastic materials is used.
The first parameter

>>> results["C10"]  # MPa
86.0

defines the deviatoric behavior.
The second parameter

>>> results["D1"]  # 1/MPa
0.0012

controls the behavior under a hydrostatic load.


Geometric dimensions
^^^^^^^^^^^^^^^^^^^^

.. _example_6_dimensions:

The dimensions of the model shown in :ref:`Figure 1 <example_6_model_image>` are:

>>> results["lx"]  # mm
100.0

>>> results["ly"]  # mm
50.0

>>> results["lz"]  # mm
5.0

>>> results["cx"]  # mm
50.0

>>> results["cy"]  # mm
40.0

>>> results["r"]  # mm
5.0

>>> results["r1"]  # mm
6.0

>>> results["r2"]  # mm
8.0


Applied displacement
^^^^^^^^^^^^^^^^^^^^

>>> results["uy"]  # mm
10.0


Configurational forces
^^^^^^^^^^^^^^^^^^^^^^

The configurational forces in x- and y-direction are computed for two domains.
Domain :math:`\mathcal{A}` has configurational forces of

>>> CF_A = results["CF_A"]
>>> # in mJ/mm**2
>>> CF_A[:-1]  # doctest: +ELLIPSIS
[4.90..., -16.89...]

Domain :math:`\mathcal{B}` is larger and configurational forces of

>>> CF_B = results["CF_B"]
>>> # in mJ/mm**2
>>> CF_B[:-1]  # doctest: +ELLIPSIS
[4.89..., -16.90...]

matches the results of :math:`\mathcal{A}` with less than 0.3% deviation.
We consider only the volume integral :eq:`eq_23` and neglect the surface influence on the configurational forces.


Numerical energy gradient
^^^^^^^^^^^^^^^^^^^^^^^^^

The energy gradient with respect to the position of the hole can also be found by a numerical derivation.
For the numerical derivation, we move the hole once by a small dx

>>> model_name_dx = "bending_model_dx"
>>> dx = 0.1  # mm
>>> cx_dx = cx + dx
>>> subprocess.call(f"abaqus cae noGui=example_6_abaqus_script.py -- {model_name_dx} {cx_dx:.3f} {cy:.3f}", shell=True)
0

and once by a small dy:

>>> model_name_dy = "bending_model_dy"
>>> dy = 0.1  # mm
>>> cy_dy = cy + dy
>>> subprocess.call(f"abaqus cae noGui=example_6_abaqus_script.py -- {model_name_dy} {cx:.3f} {cy_dy:.3f}", shell=True)
0

The simulation computes the strain energy (ALLSE) and writes it into json files.

>>> with open(f"{model_name_dx}_results.json", "r", encoding="UTF-8") as fh:
...     results_dx = json.load(fh)

>>> with open(f"{model_name_dy}_results.json", "r", encoding="UTF-8") as fh:
...     results_dy = json.load(fh)

Next, the energies (ALLSE) are read from the json files.

>>> ALLSE_0 = results["ALLSE"]
>>> ALLSE_dx = results_dx["ALLSE"]
>>> ALLSE_dy = results_dy["ALLSE"]

With the energies, we can use a forward difference quotient to numerically estimate the energy gradient.

>>> Gx = (ALLSE_dx - ALLSE_0) / dx
>>> Gy = (ALLSE_dy - ALLSE_0) / dy
>>> # in mJ/mm**2
>>> [Gx, Gy]  # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
[4.9..., -17...]


Conclusion
----------

The energy gradient with respect to x computed by the configurational forces
is almost identical to the result of the numerical derivation.
The energy gradient with respect to y is also in a good agreement,
although the configurational forces method is slightly smaller than
the result obtained by the numerical derivation.

To minimize the strain energy, the hole should be moved in the opposite direction of the gradient.
Consequently, moving the hole up and to the left would decrease the strain energy.
This seems plausible because the bending stiffness (and the energy) will be less
if the hole is at the top compared to the hole being in the middle.
Furthermore, the bending moment (and the energy density) is greater on the left side.
A hole on the left side reduces the total energy more than a hole on the right side.

"""
import doctest


if __name__ == '__main__':
    doctest.testmod(verbose=True)
