r"""
Example 1 - Two-phase bar
=========================

Considering a bar that consists of a stiffer and a softer material
with yoiungs_modulus :math:`E_{1}`, and :math:`E_{2}`.
The bar is fixed on the left side and a displacement :math:`u` is applied to the right side.
The configurational forces on the interface should be computed using the Abaqus plugin.

.. list-table:: Parameters
    :header-rows: 1

    * - Parameter
      - value
    * - :math:`l`
      - 10 mm
    * - :math:`h`
      - 10 mm
    * - :math:`u`
      - 0.1 mm
    * - :math:`\nu` (Poisson's ratio)
      - 0
    * - :math:`E_{1}`
      - 210 GPa
    * - :math:`E_{2}`
      - 105 GPa

.. image:: example_1_images/00_scheme.png
    :width: 400
    :alt: scheme

Theory
------

Kollnig et. al. [1]_ provide a theoretical solution for the configurational forces on the interface.
The solution can be computed symbolically using sympy.
First, define symbolic and real quantities.

>>> import sympy as sy
>>> (E1, E2, A, l, l1, e, e1, e2) = sy.symbols("E1 E2 A l l1 e e1 e2", real=True)
>>> E1_val = 210_000  # MPa
>>> E2_val = E1_val / 2  # MPa
>>> l1_val = 10  # mm
>>> l_val = 2 * l1_val
>>> A_val = 10  # mm^2
>>> u_val = 0.1  # mm

Next, the energy of the beam is formulated as

>>> ALLSE = 0.5 * E1 * e1**2 * A * l1 + 0.5 * E2 * e2**2 * A * (l - l1)

The strain :math:`e = \frac{u}/{l} = e_{1} + e_{2}` is the sum of the strain
`e1` in the left material and the strain `e2` in the right material.
This kinematic relationship between `e1`, `e2` and `e` is used to eliminated `e2`.

>>> e2_expr = sy.solve(sy.Eq(l*e, l1*e1 + (l - l1) * e2), e2)[0]
>>> ALLSE = ALLSE.replace(e2, e2_expr)
>>> e2_expr
(e*l - e1*l1)/(l - l1)

Furthermore, minimizing the strain energy `ALLSE` leads to the quasi-static solution of `e1`.
For this reason, the derivative of `ALLSE` with respect to the not constrained strain `e1` is set to zero.
Consequently, `e1` is eliminated.

>>> e1_expr = sy.solve(sy.Eq(0, ALLSE.diff(e1)), e1)[0]
>>> ALLSE = ALLSE.replace(e1, e1_expr)
>>> e1_expr
E2*e*l/(E1*l - E1*l1 + E2*l1)

Finally, the total derivative of the strain energy with respect to the position of the interface is computed.
The symbolic values are replaced by the numeric values.
This results in the expected configurational force in x-direction at the interface.

>>> dAllSE_dl1 = ALLSE.diff(l1)
>>> CFx_at_interface_theory = dAllSE_dl1.xreplace({
...     A: A_val,
...     E1: E1_val,
...     E2: E2_val,
...     e: u_val / l_val,
...     l: l_val,
...     l1: l1_val
... })
>>> CFx_at_interface_theory  # N
11.6666666666667

Apply Plug-in
-------------

In this section the configurational forces are computed using the Abaqus Plug-in.
First, open the folder that contains the \*.inp file.
This folder is the working directory.

.. image:: example_1_images/01_folder.png
    :alt: working directory

Open a shell in the working directory.

.. image:: example_1_images/02_start_cmd_in_folder.png
    :alt: start cmd in working directory

Type `abaqus cae` into the shell to start Abaqus in the working directory.

.. image:: example_1_images/03_cmd_in_folder.png
    :alt: start abaqus cae in working directory

Navigate to `File -> Import -> Model...`

.. image:: example_1_images/04_open_model.png
    :alt: import model from inp

Set the File Filter, such that \*.inp files are displayed and
select `Two_phase_bar.inp`. Click OK

.. image:: example_1_images/05_open_model.png
    :alt: select inp

Abaqus shows the model.
The Load module displays the boundary conditions and loads.

.. image:: example_1_images/06_model_load.png
    :alt: show imported model

Create a new Job.

.. image:: example_1_images/07_create_job.png
    :alt: new job

Use the default settings and click OK.

.. image:: example_1_images/08_job.png
    :alt: job settings

Submit the job and wait for completion.

.. image:: example_1_images/09_submit.png
    :alt: submit button

Open the \*.odb. Assure that the file is **not** opened Read-only.

.. image:: example_1_images/10_open_odb.png
    :alt: open odb

Navigate to `Plug-ins -> Conf. Force`

.. image:: example_1_images/11_plugin.png
    :alt: plugin toolbar

The Plug-In GUI opens.

- The odb drop-down lists all opened odb files.
- Configurational stresses and forces can be computed using various methods.
  Select your desired method.
- Configurational stresses and forces need a energy density.
  Set the energy density to `SENER` to only consider the elastic strain energy density.
  Or use `SENER+PENER` to build the sum of elastic and plastic energy density.
- Select quantities in the field output section that are computed and saved into the odb.
- The name of the configurational stress and force can be modified
  to compute them with different methods.

Click Apply to start the computation.

.. image:: example_1_images/12_plugin.png
    :alt: plugin gui

After all requested quantities are written into the odb,
a log summary is printed.

.. image:: example_1_images/13_log_summary.png
    :alt: log summary

The requested field outputs can be selected and plotted.

.. image:: example_1_images/14_CF.png
    :alt: result visualization


Verify Results
--------------

For the comparison with the theoretical values, the x-component of the
configurational forces on the interface are summed up.
This can be done by probing values in abaqus.
Navigate to `Tools -> Query -> Prove values`.

.. image:: example_1_images/15_probe_values.png
    :alt: probe values

Select `Node` as Probe position and activate the x-component of the
configurational forces in the viewport.
Select the two nodes on the interface in the viewport.

.. image:: example_1_images/16_CF_at_interface.png
    :alt: CF at interface nodes

The plugin computes configurational forces of

>>> CFx_at_interface_plugin = 5.8591 + 5.8591  # N
>>> CFx_at_interface_plugin
11.7182

The theoretical and computed values differ by a factor of

>>> CFx_at_interface_plugin / CFx_at_interface_theory
1.00441714285714

References
----------

.. [1] Kolling, S., and R. Mueller.
    "On configurational forces in short-time dynamics and their computation with an explicit solver."
    Computational Mechanics 35 (2005): 392-399.
"""
import doctest

if __name__ == '__main__':
    doctest.testmod(verbose=True)
