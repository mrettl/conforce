r"""
The first Piola-Kirchhoff stress `P` is defined on the undeformed state.
This is contrary to the Cauchy stress `S` that considers the deformed state,

>>> import numpy as np
>>> import sympy as sy
>>> import conforce_3.expressions as expr

Tensile Test
------------

This example demonstrates the influence of the lateral strain
of a tensile specimen on the (true) Cauchy stress
and the Piola-Kirchhoff stress `P`.

Assuming a tensile load with the reaction force

>>> RF = sy.Matrix([
...     [10_000.],
...     [0.]])

acting on the upper surface of a specimen with
an undeformed area and undeformed length of

>>> width_undeformed = 10.
>>> length_undeformed = 100.

This causes a strain of :math:`\varepsilon_{\textrm{x}} = 0.1`
and a lateral strain of :math:`\varepsilon_{\textrm{y}} = -0.03`.
The deformation gradient is

>>> F = expr.eval_F(2, sy.Matrix([
...     [0.1, 0.0],
...     [0.0, -0.03]]))

The traction (force per area with depth=1)

>>> Tx_undeformed = RF / (1 * width_undeformed)
>>> Tx_undeformed
Matrix([
[1000.0],
[     0]])

acts on the undeformed upper surface with normal vector

>>> Nx_undeformed = sy.Matrix([
...     [1.],
...     [0.]])

Furthermore, no force acts on the side of the specimen
with normal vector

>>> Ny_undeformed = sy.Matrix([
...     [0.],
...     [1.]])

and hence the corresponding tractionon the specimen sides is

>>> Ty_undeformed = sy.Matrix([
...     [0.],
...     [0.]
... ]) / (1*length_undeformed)

The traction vector of the undeformed state is defined as

.. math::

    T_{\textrm{undeformed}} = P \cdot N_{\textrm{undeformed}}

Consequently, the first Piola-Kirchhoff stress can be computed
by multiplying the base of the undeformed traction vectors
times the inverse base of the undeformed normal vectors.

>>> T_undeformed = sy.Matrix.hstack(Tx_undeformed, Ty_undeformed)
>>> N_undeformed = sy.Matrix.hstack(Nx_undeformed, Ny_undeformed)
>>> P = T_undeformed * N_undeformed.inv()
>>> P
Matrix([
[1000.0, 0],
[     0, 0]])

Compare this to the computation of the true Cauchy stress
that corresponds to the deformed state with

>>> (length_deformed,), (width_deformed,) = np.around(np.array(
...     F @ sy.Matrix([
...         [length_undeformed],
...         [width_undeformed]
... ]), dtype=float), 6)
>>> length_deformed
110.0
>>> width_deformed
9.7

For pure tension, the normal vectors of the surfaces
stay the same.

>>> Nx_deformed = Nx_undeformed
>>> Ny_deformed = Ny_undeformed
>>> N_deformed = N_undeformed

The tractions of the Cauchy stress are computed with respect to
the deformed width and length.

>>> Tx_deformed = RF / (1 * width_deformed)
>>> TY_deformed = sy.Matrix([
...     [0.],
...     [0.]
... ]) / (1*length_deformed)
>>> T_deformed = sy.Matrix.hstack(Tx_deformed, TY_deformed)

The Cauchy stress is

>>> S = T_deformed * N_deformed.inv()
>>> np.around(np.array(S, dtype=float), 6)
array([[1030.927835,    0.      ],
       [   0.      ,    0.      ]])

higher than the first Piola-Kirchhoff stress.

The first Piola-Kirchhoff stress can be also computed
directly with the Cauchy stress and the deformation gradient.

>>> expr.eval_P(F, S)
Matrix([
[1000.0, 0],
[     0, 0]])


Simple shear and tension
------------------------

Considering a combined load case (:math:`\gamma_{\textrm{yx}} = 0.1; \varepsilon_{\textrm{xx}}=0.2`)
with a deformation gradient

>>> F = expr.eval_F(2, sy.Matrix([
...     [0.2, 0.1],
...     [0.0, 0.0]]))

and a Cauchy stress of

>>> S = sy.Matrix([
...     [280., 40.],
...     [40.0, 0.0]])

the corresponding first Piola-Kirchhoff tensor is

>>> expr.eval_P(F, S)
Matrix([
[276.0, 48.0],
[ 40.0,    0]])

not symmetric anymore.


Rotation
--------

This example demonstrates how a rigid body rotation affects the
first Piola-Kirchhoff stress.
A rigid body rotation of 45° degree leads to a deformation gradient of

>>> angle = np.deg2rad(45)
>>> dU_dX = sy.Matrix([
...     [sy.cos(angle) - 1, sy.sin(angle)],
...     [-sy.sin(angle), sy.cos(angle) - 1],
... ])
>>> F = expr.eval_F(2, dU_dX)
>>> F  # doctest: +ELLIPSIS
Matrix([
[ 0.7071..., 0.7071...],
[-0.7071..., 0.7071...]])

A (true) Cauchy-stress appears in the specimen.
The Cauchy-stress corresponds to the deformed and hence rotated state.
Consequently, the Cauchy-stress

>>> S = sy.Matrix([
...     [1., 0.],
...     [0., 0.],
... ])

points normal to the rotated state.

However, the first Piola-Kirchhoff stress tensor is defined in the
undeformed state.
This results in

>>> P = expr.eval_P(F, S)
>>> np.around(np.array(P, dtype=float), 3)
array([[0.707, 0.707],
       [0.   , 0.   ]])

Nevertheless, the stress magnitude does not change since there
is no change of the surface area.

>>> np.linalg.norm(np.array(P[0, :], dtype=float))
1.0

"""
import doctest

if __name__ == '__main__':
    doctest.testmod()
