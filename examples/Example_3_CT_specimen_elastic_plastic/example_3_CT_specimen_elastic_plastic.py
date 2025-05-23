r"""

Example 3 shows how `ConForce` can evaluate configurational forces in an elastic-plastic
model with a crack.
Small strain plasticity is assumed.
Import packages and change to the directory where the \*.inp file is located.

>>> import os
>>> import json
>>> import matplotlib.pyplot as plt
>>> HOME_DIR = os.path.abspath(".")
>>> os.chdir(os.path.abspath(__file__ + "/.."))

Problem description
-------------------

Kolednik [1]_ suggested the following compact tension (CT) test.

.. image:: example_3_images/00_scheme.png
    :width: 400
    :alt: scheme


Geometric dimensions:
^^^^^^^^^^^^^^^^^^^^^

>>> w = 50  # mm
>>> a = 27  # mm
>>> b = 25  # mm
>>> h = 60  # mm

Material properties:
^^^^^^^^^^^^^^^^^^^^

Schoengrundner [2]_ provides material data for an elastic-plastic behaviour with isotropic hardening.
Region `D` is modeled without plasticty to prevent large deformations
in the vicinity of constraint nodes.

>>> nu = 0.3
>>> E = 200_000 # MPa
>>> yield_stress_over_plastic_strain = np.array([
...     [    270.,         0.],
...     [ 273.052, 0.00109256],
...     [ 300.462, 0.00837633],
...     [ 329.014,  0.0198685],
...     [ 367.845,  0.0392919],
...     [ 394.113,  0.0580678],
...     [ 415.812,   0.082347],
...     [ 440.557,   0.116176],
...     [ 456.166,   0.145311],
...     [ 472.155,   0.181244],
...     [ 578.368,   0.399433],
...     [   1460.,         3.],
... ])
>>> fig = plot_stress_strain(E, yield_stress_over_plastic_strain)
>>> save_fig(fig, "example_3_images/01_stress_strain.png")

.. image:: example_3_images/01_stress_strain.png
    :width: 400
    :alt: stress strain curve

Applied displacements:
^^^^^^^^^^^^^^^^^^^^^^

>>> u_max = 0.5  # mm

Literature of crack-driving force
---------------------------------

Kolednik [1]_ investigates this problem.
The non-linear elastic J-integrals (j_nl_el\_\*)
define the energy density as sum of elastic and plastic energy densities (e=SENER+PENER).
The incremental plasticity J-integrals (j_inc_pl\_\*)
define the energy density only as the elastic energy densities (e=SENER) and ignore the plastic contribution.

>>> literature_data = {
...     "u": [0.10, 0.25, 0.50],  # mm
...     "j_nl_el_contour_tip": [1.513, 8.45, 31.80],  # mJ/mm²
...     "j_inc_pl_contour_tip": [1.736, 10.15, 38.58],  # mJ/mm²
...     "j_nl_el_contour_1": [2.132, 11.92, 44.72],  # mJ/mm²
...     "j_inc_pl_contour_1": [2.271, 13.90, 51.18],  # mJ/mm²
...     "j_nl_el_contour_3": [2.239, 12.88, 46.36],  # mJ/mm²
...     "j_inc_pl_contour_3": [2.239, 13.05, 47.42],  # mJ/mm²
...     "j_nl_el_contour_7": [2.246, 13.67, 47.03],  # mJ/mm²
...     "j_inc_pl_contour_7": [2.246, 13.20, 47.35],  # mJ/mm²
...     "j_nl_el_contour_25": [2.247, 13.56, 48.13],  # mJ/mm²
...     "j_inc_pl_contour_25": [2.247, 13.56, 47.72],  # mJ/mm²
...     "j_nl_el_contour_far": [2.247, 13.56, 48.21],  # mJ/mm²
...     "j_inc_pl_contour_far": [2.247, 13.10, 30.40],  # mJ/mm²
... }


Simulation
----------

The input file is located in the example folder.
The model evaluation is scripted and there is no need to open the `ConForce`-plugin manually.
First, call the Abaqus script (:py:mod:`example_3_abaqus_script`).
The script simulates the \*.inp file and writes a `results.json` file.

>>> import subprocess
>>> subprocess.call("abaqus cae noGUI=example_3_abaqus_script.py", shell=True)
0

After the simulation, the result file is loaded to fetch the reaction force `load`.
This force is needed to reach the defined displacement `u`.

>>> with open("results.json", "r", encoding="UTF-8") as fh:
...     results = json.load(fh)
>>> load = np.array(results["reaction_force"])[:, 1]
>>> u = np.array(results["u"])[:, 1]

With the obtained results the force-displacement curve is plotted.

>>> fig = plot_force_displacement(u, load)
>>> save_fig(fig, "example_3_images/02_force_displacement.png")

.. image:: example_3_images/02_force_displacement.png
    :width: 400
    :alt: force displacement

Abaqus J-Integral
-----------------

Abaqus computes the J-integral using the virtual crack extension method by Parks [3]_.
This J-integral is compared to the literature values of Kolednik [1]_ that are marked by "x".

>>> fig = compare_abaqus_j_with_literature(results, literature_data)
>>> save_fig(fig, "example_3_images/03_j_integral.png")

.. image:: example_3_images/03_j_integral.png
    :width: 400
    :alt: comparison of J-Integral with literature

The Abaqus J-integrals agree well with the literature values.

Configurational forces
----------------------

Furthermore, the x-components of the configurational forces are summed up in regions sorrounded by the contours.
This corresponds to the J-integral with a negative sign.
There is good agreement between the `ConForce` and the literature values.

>>> fig = compare_conforce_cfx_with_literature(results, literature_data)
>>> save_fig(fig, "example_3_images/04_negative_cfx.png")

.. image:: example_3_images/04_negative_cfx.png
    :width: 400
    :alt: comparison of configurational forces with literature

References
----------

.. [1] O. Kolednik, R. Schöngrundner, and F. D. Fischer,
    “A new view on J-integrals in elastic–plastic materials,”
    Int J Fract, vol. 187, no. 1, pp. 77–107, May 2014, doi: `10.1007/s10704-013-9920-6 <https://doi.org/10.1007/s10704-013-9920-6>`_.

.. [2] Schoengrundner.
    "Numerische Studien zur Ermittlung der risstreibenden Kraft in elastisch-plastischen Materialien bei unterschiedlichen Belastungsbedingungen"
    Ph.D. thesis, University of Leoben (2010).

.. [3] D. M. Parks,
    “The virtual crack extension method for nonlinear material behavior,”
    Computer Methods in Applied Mechanics and Engineering, vol. 12, no. 3, pp. 353–364, Dec. 1977, doi: `10.1016/0045-7825(77)90023-8 <https://doi.org/10.1016/0045-7825(77)90023-8>`_.


After the evaluation, the directory is changed to the home directory.

>>> os.chdir(HOME_DIR)
"""
import doctest
import os.path

import numpy as np


def save_fig(fig, path):
    if not os.path.exists(path):
        fig.savefig(path)


def plot_stress_strain(youngs_modulus, yield_stress_over_plastic_strain):
    import matplotlib.pyplot as plt

    yield_stress_cal = yield_stress_over_plastic_strain[:, 0]
    plastic_strain_cal = yield_stress_over_plastic_strain[:, 1]

    # compute stress and strains
    plastic_strains = np.linspace(0, 0.5, 200)
    stresses = np.interp(plastic_strains, plastic_strain_cal, yield_stress_cal)
    elastic_strains = stresses / youngs_modulus
    total_strains = elastic_strains + plastic_strains

    # compute stress strain curve
    strains = np.append([0.], total_strains)
    stresses = np.append([0.], stresses)

    # plot
    fig, (ax1, ax2) = plt.subplots(nrows=2)  # type: plt.Figure, (plt.Axes, plt.Axes)
    fig.set_size_inches(6., 6.)
    zoom_strain = 0.05
    zoom_stress = np.interp(0.05, strains, stresses)

    ax1.plot(strains, stresses, "k-")
    ax1.plot([zoom_strain, zoom_strain, 0], [0, zoom_stress, zoom_stress], "k:")
    ax1.set_xlim(left=0, right=0.5)
    ax1.set_ylim(bottom=0)
    ax1.set_xlabel(r"$\varepsilon$ [-]")
    ax1.set_ylabel(r"$\sigma$ [MPa]")

    ax2.plot(strains, stresses, "k-")
    ax2.set_xlim(left=0, right=zoom_strain)
    ax2.set_ylim(bottom=0, top=zoom_stress)
    ax2.set_xlabel(r"$\varepsilon$ [-]")
    ax2.set_ylabel(r"$\sigma$ [MPa]")

    fig.tight_layout()
    return fig


def plot_force_displacement(u, load):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()  # type: plt.Figure, plt.Axes
    fig.set_size_inches(6., 6.)
    ax.plot(u, load, "k-")
    ax.set_xlabel("displacement [mm]")
    ax.set_ylabel("load [N]")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    return fig


def compare_abaqus_j_with_literature(results, literature):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()  # type: plt.Figure, plt.Axes
    fig.set_size_inches(6., 6.)

    j = results["J"]
    t, u = np.array(results["u"]).T
    for contour in [1, 3, 7, 25]:
        (line, ) = ax.plot(
            u,
            np.interp(t, *np.array(j[f"J at J_NEAR_CRACK_TIP_Contour_{contour+1:02}"]).T),
            label=f"contour {contour}; Abaqus"
        )
        ax.plot(
            literature["u"],
            literature[f"j_nl_el_contour_{contour}"],
            label=f"contour {contour}; Literature", ls="", marker="x",
            color=line.get_color()
        )

    (line, ) = ax.plot(
        u,
        np.interp(t, *np.array(j["J at J_FAR_FAR_FIELD_Contour_1"]).T),
        label=f"far field; Abaqus"
    )
    ax.plot(
        literature["u"],
        literature[f"j_nl_el_contour_far"],
        label=f"far field; Literature", ls="", marker="x",
        color=line.get_color()
    )

    ax.legend()
    ax.set_xlabel("displacement [mm]")
    ax.set_ylabel("J-integral [mJ/mm²]")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    return fig


def compare_conforce_cfx_with_literature(results, literature):
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(nrows=2)  # type: plt.Figure, (plt.Axes, plt.Axes)
    fig.set_size_inches(6, 6*2)

    t, u = np.array(results["u"]).T
    for ax, cf, lit_fix, energy_expression in (
            (ax1, results["CONF_FORCE_EL_at_frame"], "inc_pl", "SENER"),
            (ax2, results["CONF_FORCE_EL_PL_at_frame"], "nl_el", "SENER+PENER")
    ):
        ax.set_title(f"ConForce with e={energy_expression} and literature with j_{lit_fix} ")

        # crack tip
        (line,) = ax.plot(
            u,
            -np.array(cf["J_NEAR_J_CRACK_TIP_Contour_01"])[:, 0],
            label=f"crack tip; ConForce"
        )
        ax.plot(
            literature["u"],
            literature[f"j_{lit_fix}_contour_tip"],
            label=f"crack tip; Literature", ls="", marker="x",
            color=line.get_color()
        )

        # contours
        for contour in [1, 3, 7, 25]:
            (line, ) = ax.plot(
                u,
                -np.array(cf[f"J_NEAR_J_CRACK_TIP_Contour_{contour+1:02}"])[:, 0],
                label=f"contour {contour}; ConForce"
            )
            ax.plot(
                literature["u"],
                literature[f"j_{lit_fix}_contour_{contour}"],
                label=f"contour {contour}; Literature", ls="", marker="x",
                color=line.get_color()
            )

        # far field
        (line, ) = ax.plot(
            u,
            -np.array(cf["J_FAR_J_FAR_FIELD_Contour_1"])[:, 0],
            label=f"far field; ConForce"
        )
        ax.plot(
            literature["u"],
            literature[f"j_{lit_fix}_contour_far"],
            label=f"far field; Literature", ls="", marker="x",
            color=line.get_color()
        )

        ax.legend()
        ax.set_xlabel("displacement [mm]")
        ax.set_ylabel(r"$-\mathrm{CF}_{x}$ [mJ/mm²]")
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)

    fig.tight_layout()
    return fig


if __name__ == '__main__':
    doctest.testmod(verbose=True)
