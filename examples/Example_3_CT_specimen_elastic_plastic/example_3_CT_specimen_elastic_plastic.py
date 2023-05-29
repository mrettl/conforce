r"""
Example 3 - CT-specimen with linear elastic material behaviour
==============================================================

Problem description
-------------------

Kolednik [1]_ suggested the following compact tension (CT) test.

.. image:: example_3_images/00_scheme.png
    :width: 400
    :alt: scheme


Geometeric dimensions:
^^^^^^^^^^^^^^^^^^^^^^

>>> import os
>>> import numpy as np
>>> import matplotlib.pyplot as plt
>>> w = 50  # mm
>>> a = 27  # mm
>>> b = 25  # mm
>>> h = 60  # mm

Material properties:
^^^^^^^^^^^^^^^^^^^^

>>> nu = 0.3
>>> youngs_modulus = 200_000 # MPa
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
>>> fig = plot_stress_strain(youngs_modulus, yield_stress_over_plastic_strain)
>>> _ = fig.savefig(__file__ + "/../example_3_images/01_stress_strain.svg")

.. image:: example_3_images/01_stress_strain.svg
    :width: 400
    :alt: stress strain curve

Boundary condition:
^^^^^^^^^^^^^^^^^^^

>>> u = 0.5  # mm

Simulation
----------

The model is scripted and there is no need to manually open the `conforce`-plugin.
First, change to the directory where the \*.inp file is located.

>>> HOME_DIR = os.path.abspath(".")
>>> os.chdir(__file__ + "/..")

Next, call the Abaqus script (:py:mod:`example_3_abaqus_script`).
The script simulates the \*.inp file and writes a `results.json` file.
To save time, the script is not exectued if the `results.json` file already exists.

>>> import subprocess
>>> if not os.path.exists("results.json"):
...     _ = subprocess.call("abaqus cae noGUI=example_3_abaqus_script.py", shell=True)

After the simulation, the result file is loaded to fetch the reaction force `load`.
This force is needed to reach the defined displacement `u`.

>>> import matplotlib.pyplot as plt import json
>>> with open("results.json", "r", encoding="UTF-8") as fh:
...     results = json.load(fh)
>>> load = np.array(results["reaction_force"])[:, 1]
>>> u = np.array(results["u"])[:, 1]

>>> fig, ax = plt.subplots()  # type: plt.Figure, plt.Axes
>>> _ = ax.plot(u, load, "k-")
>>> _ = ax.set_xlabel("u [mm]"), ax.set_ylabel("load [N]")
>>> _ = ax.set_xlim(left=0), ax.set_ylim(bottom=0)
>>> _ = fig.savefig("example_3_images/02_force_displacement.svg")

.. image:: example_3_images/02_force_displacement.svg
    :width: 400
    :alt: force displacement


Abaqus J-Integral
-----------------

.. todo: describe results

Configurational forces
----------------------

.. todo: describe results

References
----------

.. [1] Kolednik, Otmar, Ronald SchÃ¶ngrundner, and Franz Dieter Fischer.
    "A new view on J-integrals in elasticâ€“plastic materials."
    International Journal of Fracture 187.1 (2014): 77-107.

.. [2] Parks, D. M.
    "The virtual crack extension method for nonlinear material behavior."
    Computer methods in applied mechanics and engineering 12.3 (1977): 353-364.

Change to home directory

>>> os.chdir(HOME_DIR)
"""
import doctest


def plot_stress_strain(youngs_modulus, yield_stress_over_plastic_strain):
    import numpy as np
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


if __name__ == '__main__':
    doctest.testmod(verbose=True)
