r"""
Problem description
-------------------

:ref:`Figure 1 <example_4_scheme_image>` depicts the model considered in this example.
The circular model contains a crack with its tip in the model center.
The green dashed boundary is deformed according to a defined displacement field function.
The region :math:`\mathcal{A}` is used to compute the resulting configurational force.
The blue arrows depict the configurational force at the nodes.

.. _example_4_scheme_image:

.. figure:: example_4_images/K_field_model.png
    :width: 800
    :alt: scheme

    Figure 1: FEM mesh of a crack model whose boundaries (green dashed line) are deformed according to
    the theoretical displacement field around a crack.
    The configurational forces are evaluated in the red region.
    The depicted deformation is scaled by a factor of 100.


Simulation
----------

>>> import os
>>> import json
>>> import subprocess
>>> import numpy as np
>>> HOME_DIR = os.path.abspath(".")
>>> os.chdir(__file__ + "/..")

The model is automatically build and simulated by the Abaqus script.

>>> subprocess.call("abaqus cae noGUI=example_4_abaqus_script.py", shell=True)
0

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

>>> R_mm = max(results["R_mm"])
>>> R_mm
50.0

The region :math:`\mathcal{A}` depicted in :ref:`Figure 1 <example_4_scheme_image>`
is a square with side lengths of:

>>> length_region_A_mm = results["inner_region_length_mm"]
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

is choosen.
The displacements of the outer boundary (green dashed line) of the model
are prescribed by the displacement field :math:`u(r=R, \varphi)`
according to Anderson.
The displacement field function can be found in the Abaqus script.

Results
-------

Theoretical energy release rate
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

According to Anderson, the applied energy release rate

>>> G_applied_mJ_mm2 = (
...     1_000  #  mm / m
...     * KI_MPa_m_05**2
...     * (1 - nu**2)
...     / E_MPa
... )
>>> G_applied_mJ_mm2  # doctest: +ELLIPSIS
1.733...

can be computed for a plane strain state from the stress intensity factor,
the Poisson's ratio, and the Young's modulus.


Energy release rate from node closure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Another approach to obtain the energy release rate
is to compare the strain energies for two crack lengths.
For this, we simulate the model twice:

    - In the first simulation, the tip lies in the model center.
      This results in a strain energy of :math:`\Pi_{0}`.
    - In the second simulation, the first two nodes that were open in the previous
      simulation are closed. This reduces the crack length by :math:`da`
      and leads to a strain energy of :math:`\Pi_{-1}`.

The energy release rate

.. math::

    G = \frac{\Pi_{-1} - \Pi_{0}}{t \cdot da}

is computed from the strain energies and the area of the closed crack face.
This results in a energy release rate of:

>>> G_mJ_mm2 = results["G_mJ_mm2"]
>>> G_mJ_mm2  # doctest: +ELLIPSIS
1.735...


Energy release rate from the Abaqus J-Integral
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Abaqus computes the J-Integral using the virtual crack extension method by Parks [2]_.
The J-Integral is evaluated for several contours for a crack growing in x-direction.
For the contour that encloses the region :math:`\mathcal{A}`,
the J_integral is:

>>> J_mJ_mm2 = results["J_mJ_mm2"][-1]
>>> J_mJ_mm2  # doctest: +ELLIPSIS
1.738...

Energy release rate from conforce
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Conforce computes the nodal configurational forces.
The energy release rate

.. math::

    G_{CF} = -\vec{v} \cdot \vec{CF}

is computed from the crack growth direction and the resulting configurational force vector.
The resulting configurational force vector sums up all configurational forces in the region :math:`\mathcal{A}`.

The crack growth direction points in x-direction:

>>> v = np.array([1, 0])

The configurational forces are computed by two formulation.
The motion based formulation results in an energy release rate of:

>>> G_mbf_mJ_mm2 = -np.dot(v, results["CF_mbf_mJ_mm2"][0])
>>> G_mbf_mJ_mm2  # doctest: +ELLIPSIS
1.738...

The displacement based formulation results in an energy release rate of:

>>> G_dbf_mJ_mm2 = -np.dot(v, results["CF_dbf_mJ_mm2"][0])
>>> G_dbf_mJ_mm2  # doctest: +ELLIPSIS
1.738...


Conclusion
----------

All approaches compute similar energy release rates of about 1.7 mJ/mm²
with a deviation of less than 0.3%.

References
----------

.. [1] T. L. Anderson,
    Fracture mechanics: fundamentals and applications.
    Boca Raton: CRC Press, 1991.

.. [2] D. M. Parks,
    “The virtual crack extension method for nonlinear material behavior,”
    Computer Methods in Applied Mechanics and Engineering, vol. 12, no. 3, pp. 353–364, Dec. 1977, doi: `10.1016/0045-7825(77)90023-8 <https://doi.org/10.1016/0045-7825(77)90023-8>`_.



"""