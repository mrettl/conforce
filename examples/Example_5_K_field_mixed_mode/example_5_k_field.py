r"""
Problem description
-------------------

:ref:`Figure 1 <example_5_scheme_image>` depicts the model considered in this example,
which is basically Example 4 with an additional mode II loading.
The circular model contains a crack with its tip in the model center.
The green dashed boundary is displaced according to a defined displacement field function.
The blue arrows depict the configurational forces at the nodes.
The configurational forces are summed up for several contours to estimate the energy release rate.
For example, the 3rd contour contains all dark-brown nodes
and the 15th contour contains all orange and dark-brown nodes.

.. _example_5_scheme_image:

.. figure:: example_5_images/K_field_model.png
    :width: 800
    :alt: scheme

    Figure 1: Model of a crack model whose boundaries (green dashed line) are deformed according to
    the theoretical displacement field around a crack (KI and KII).
    The configurational forces are evaluated in the red region.
    The depicted deformation is scaled by a factor of 100.


Simulation
----------

>>> import os
>>> import json
>>> import subprocess
>>> import numpy as np
>>> import scipy.optimize as optimize
>>> HOME_DIR = os.path.abspath(".")
>>> os.chdir(__file__ + "/..")

The model is automatically build and simulated by the Abaqus script.

>>> if not os.path.exists("results.json"):
...     _ = subprocess.call("abaqus cae noGUI=example_5_abaqus_script.py", shell=True)

The Abaqus script writes results to a json file.
The command below loads the results.

>>> with open("results.json", "r", encoding="UTF-8") as fh:
...     results = json.load(fh)

Model parameters
----------------

Linear rectangular (CPE4) and triangular (CPE3) plane strain elements are used.
Furthermore, the flag NLGEOM is turned on.

Geometric parameters
^^^^^^^^^^^^^^^^^^^^

The radius of the model is:

>>> R_mm = results["regions"][-1]["R_mm"]
>>> R_mm
50.0

The red framed region depicted in :ref:`Figure 1 <example_5_scheme_image>`
is a square with side lengths of:

>>> length_region_A_mm = results["regions"][0]["R_mm"]
>>> length_region_A_mm  # doctest: +ELLIPSIS
1.4...

The thickness of the model is:

>>> t_mm = 1

Material properties
^^^^^^^^^^^^^^^^^^^

The Young's modulus is:

>>> E_MPa = results["E_MPa"]
>>> E_MPa
210000

The Poisson's ratio is:

>>> nu = results["nu"]
>>> nu
0.3

Applied displacement
^^^^^^^^^^^^^^^^^^^^

Anderson [1]_ describes displacement fields :math:`u(r, \varphi)`
that correspond to certain stress intensity factors.
In this model, a mode-I stress intensity factor of

>>> KI_MPa_m_05 = results["KI_MPa_m_05"]
>>> KI_MPa_m_05
20

and a mode-II stress intensity factor of

>>> KII_MPa_m_05 = results["KII_MPa_m_05"]
>>> KII_MPa_m_05
10

are chosen.
The displacements of the outer boundary (green dashed line) of the model
are prescribed by the displacement field :math:`u(r=R, \varphi)`
according to Anderson.
The displacement field function can be found in the Abaqus script.

Results
-------

In this example, the energy release rate is predicted by two methods:

1. configurational forces computed by ConForce
2. vectorial J-Integral computed by Abaqus [2]_

The vectorial J-Integral, proposed by Budiansky and Rice [3]_, is obtained by
computing the scalar J-Integral of Rice [4]_ with two q-vectors:
The first q-vector points in the direction of the crack and results in the first component, denoted as J1.
The second  q-vector points normal to the crack face and results in the second component, denoted as J2.
The vectorial J-Integral is related to the configurational force vector, expressed as :math:`\vec{J}=-\vec{CF}`.

Crack tip
^^^^^^^^^

At the crack tip, we find a good agreement between the values computed
by the configurational forces method

>>> CF_x_tip, CF_y_tip = results["contours"][0]["CF_mbf_mJ_mm2"]
>>> CF_x_tip, CF_y_tip  # doctest: +ELLIPSIS
(-1.945..., 1.2289...)

and those obtained by the J-Integral.
The flipped sign is intended.

>>> J1_tip, J2_tip = results["J1_mJ_mm2"][0], results["J2_mJ_mm2"][0]
>>> J1_tip, J2_tip  # doctest: +ELLIPSIS
(1.887..., -1.234...)


Influence of contour size
^^^^^^^^^^^^^^^^^^^^^^^^^

According to Schmitz and Ricoeur [5]_, the J2 component might vary with the size of the contour region.
For this reason, we compute the J-Integral and configurational forces for two regions.

The smaller region (dark brown nodes in :ref:`Figure 1 <example_5_scheme_image>`) has a side length of

>>> r_mm_c3 = results["contours"][3]["R_mm"]
>>> 2 * r_mm_c3  # doctest: +ELLIPSIS
0.307...

and the larger region (orange and dark brown nodes in :ref:`Figure 1 <example_5_scheme_image>`) has a side length of

>>> r_mm_c15 = results["contours"][15]["R_mm"]
>>> 2 * r_mm_c15  # doctest: +ELLIPSIS
1.527...

If the y-component of the configurational force vector and of the J-Integral show a deviation between the smaller and larger region,
Schmitz and Ricoeur would provide a correction method.
However, in our case the configurational forces for the smaller region

>>> CF_x_c3, CF_y_c3 = results["contours"][3]["CF_mbf_mJ_mm2"]
>>> CF_x_c3, CF_y_c3  # doctest: +ELLIPSIS
(-2.149..., 1.686...)

are in good agreement with the configurational forces of the larger region.

>>> CF_x_c15, CF_y_c15 = results["contours"][15]["CF_mbf_mJ_mm2"]
>>> CF_x_c15, CF_y_c15  # doctest: +ELLIPSIS
(-2.155..., 1.697...)

The same holds true for the J-Integral computed by Abaqus. The J2 value for the smaller region

>>> J1_c3, J2_c3 = results["J1_mJ_mm2"][3], results["J2_mJ_mm2"][3]
>>> J1_c3, J2_c3  # doctest: +ELLIPSIS
(2.147..., -1.689...)

are close to the J2 value for the larger region.

>>> J1_c15, J2_c15 = results["J1_mJ_mm2"][15], results["J2_mJ_mm2"][15]
>>> J1_c15, J2_c15  # doctest: +ELLIPSIS
(2.154..., -1.699...)

Consequently, Schmitz and Ricoeur's correction method for J2 is not required in this example.


Crack growth direction from ConForce
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We consider the configurational forces at the 15th contour.
According to the configurational force method, the crack propagates with an angle (in degrees) of

>>> angle_CF_rad = np.arctan2(-CF_y_c15, -CF_x_c15)
>>> np.rad2deg(angle_CF_rad)  # doctest: +ELLIPSIS
-38.2...


Maximum Tangential Stress criterion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Maximum Tangential Stress (MTS) criterion, proposed by Erdogan and Sih [6]_,
is an alternative method to predict the direction of the crack growth under mixed mode I/II loading.
For the given values of KI and KII, the MTS criterion predicts a crack growth angle (in degrees) of

>>> angle_mts_rad = optimize.root_scalar(
...     lambda angle: KI_MPa_m_05*np.sin(angle) + KII_MPa_m_05*(3*np.cos(angle) - 1),
...     bracket=(-np.pi/2, np.pi/2)
... ).root
>>> np.rad2deg(angle_mts_rad)  # doctest: +ELLIPSIS
-40.2...


Conclusion
----------

We showed that the configurational forces computed by ConForce and the J-Integral of Abaqus
lead to the same estimation of the energy release rate.

Furthermore, we computed the crack growth angle using two methods:
The Maximum Tangential Stress criterion and the configurational forces approach.
Both methods lead to similar angles with about 2° deviation.


References
----------

.. [1] T. L. Anderson,
    Fracture mechanics: fundamentals and applications.
    Boca Raton: CRC Press, 1991.

.. [2] D. M. Parks,
    The virtual crack extension method for nonlinear material behavior,
    Computer Methods in Applied Mechanics and Engineering, vol. 12, no. 3, pp. 353–364, Dec. 1977, doi: `10.1016/0045-7825(77)90023-8 <https://doi.org/10.1016/0045-7825(77)90023-8>`_.

.. [3] B. Budiansky, J.R. Rice,
    Conservation Laws and Energy-Release Rates,
    Journal of Applied Mechanics 40 (1973) 201–203. https://doi.org/10.1115/1.3422926.

.. [4] J.R. Rice,
    A Path Independent Integral and the Approximate Analysis of Strain Concentration by Notches and Cracks,
    Journal of Applied Mechanics 35 (1968) 379–386. https://doi.org/10.1115/1.3601206.

.. [5] K. Schmitz, A. Ricoeur,
    Theoretical and computational aspects of configurational forces in three-dimensional crack problems,
    International Journal of Solids and Structures 282 (2023) 112456. https://doi.org/10.1016/j.ijsolstr.2023.112456.


.. [6] F. Erdogan, G.C. Sih,
    On the Crack Extension in Plates Under Plane Loading and Transverse Shear,
    Journal of Basic Engineering 85 (1963) 519–525. https://doi.org/10.1115/1.3656897.


"""
import doctest


if __name__ == '__main__':
    doctest.testmod(verbose=True)
