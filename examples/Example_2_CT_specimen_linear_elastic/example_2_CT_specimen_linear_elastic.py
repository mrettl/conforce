r"""
Example 2 - CT-specimen with linear elastic material behaviour
==============================================================

Problem description
-------------------

Kolednik [1]_ suggests the following compact tension (CT) test.
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
>>> E = 200_000 # MPa

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

>>> J_theory = k**2 * (1 - nu**2) / E
>>> np.around(J_theory, 3)  # mJ/mm**2
54.081

Abaqus J-Integral
-----------------

Abaqus computes the J-Integral using the virtual crack extension method by Parks [3]_.
Abaqus computes the J-Integral over several contours.
The regions correspond to the following contour indices
- region `A` corresponds to the contour index 21
- region `B` corresponds to the contour index 57
- region `FAR_FFIELD`

According to Abaqus the J-Integral for region `A` is

>>> J_in_a_abaqus = results["J"]["J at J_NEAR_CRACK_TIP_Contour_21"]
>>> np.around(J_in_a_abaqus, 3)  # mJ/mm**2
54.503

for region `B`

>>> J_in_b_abaqus = results["J"]["J at J_NEAR_CRACK_TIP_Contour_57"]
>>> np.around(J_in_b_abaqus, 3)  # mJ/mm**2
54.635

and for region `FAR_FIELD`

>>> J_in_far_abaqus = results["J"]["J at J_FAR_FAR_FIELD_Contour_1"]
>>> np.around(J_in_far_abaqus, 3)  # mJ/mm**2
54.765

This is in good aggreement with the prediction made by Anderson.


Configurational forces
----------------------

Conforce can predict J-Integrals by summing up configurational forces inside a region.
The resulting configurational force is projected onto the crack extension direction.
In this case the crack extension direction is simply the x-component.

The configurational forces for the regions `A`, `B`, and `FAR_FIELD` show a good aggreement with Anderson and Abaqus.

.. note::
    The configurational forces have a negative sign.

>>> (cfx_in_a, _) = results["CF"]["A"]
>>> np.around(cfx_in_a, 3)  # mJ/mm**2
-54.484

>>> (cfx_in_b, _) = results["CF"]["B"]
>>> np.around(cfx_in_b, 3)  # mJ/mm**2
-54.633

>>> (cfx_in_far, _) = results["CF"]["FAR_FIELD"]
>>> np.around(cfx_in_far, 3)  # mJ/mm**2
-54.794

Like conventional forces, the configurational forces fullfill the equilibrium of forces.
The configuraional forces summed up for the whole model are zero.

>>> (cfx_whole_model, _) = results["CF"]["ALL_NODES"]
>>> np.around(cfx_whole_model, 3)
-0.0

Comparison
----------

The figure compares the J-Integral of Abaqus and Anderson with the negative configurational forces of conforce.
Abaqus and conforce show a good aggreement for all contours in the regions `A` and `B`.

>>> fig = compare_J_and_negative_cfx(results, J_theory)
>>> save_fig(fig, "example_2_images/01_comparison_over_contours.png")

.. image:: example_2_images/01_comparison_over_contours.png
    :width: 400
    :alt: comparison


References
----------

.. [1] Kolednik, Otmar, Ronald Schoengrundner, and Franz Dieter Fischer.
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
import os
import doctest


def save_fig(fig, path):
    if not os.path.exists(path):
        fig.savefig(path)


def compare_J_and_negative_cfx(results, J_theory):
    import matplotlib.pyplot as plt

    contour_ids_to_set_names = {
        int(set_name.split("_")[-1]): set_name
        for set_name in results["J"].keys()
        if "J_NEAR" in set_name.upper()
    }
    cf_regions = {
        int(set_name.split("_")[-1]): set_name
        for set_name in results["CF"].keys()
        if "J_NEAR" in set_name.upper()
    }
    contour_ids = sorted(contour_ids_to_set_names.keys())

    J_map = {k: v for k, v in results["J"].items()}
    cfx_map = {k: v[0] for k, v in results["CF"].items()}

    fig, ax = plt.subplots()  # type: plt.Figure, plt.Axes
    fig.set_size_inches(3.15, 3.15)
    ax.axvline(21, ls=":", c="k")
    ax.text(21/2, 50, "region A", ha="center")
    ax.text(21 + (57-21)/2, 50, "region B", ha="center")
    ax.axhline(J_theory, ls="-", c="k", label="J-Integral (Anderson)")
    ax.plot(contour_ids, [J_map[contour_ids_to_set_names[idx]] for idx in contour_ids], label="J-Integral (Abaqus)")
    ax.plot(contour_ids, [-cfx_map[cf_regions[idx]] for idx in contour_ids], label="-CF (conforce)")
    ax.set_xlabel("contour index")
    ax.set_ylabel("J-Integral and -CF [mJ/mm²]")
    ax.set_xlim(left=0, right=contour_ids[-1])
    ax.legend()
    fig.tight_layout()
    return fig


if __name__ == '__main__':
    doctest.testmod(verbose=True)
