r"""

The deformation gradient is the derivative of the deformed coordinates with respect to the undeformed coordinates.

.. math::
    F = \frac{\partial (X + U)}{\partial X} = I + \frac{\partial U}{\partial X}

.. note::

    In the 2D case, a plane strain state is assumed.
    This implies :math:`F_{zz}=1`.

>>> import numpy as np
>>> import sympy as sy
>>> import conforce_gen.expressions as expr

Undeformed
----------

The deformation gradient of the undeformed state in 2D is the identity matrix.

>>> F_undeformed = expr.eval_F(2, sy.ZeroMatrix(2, 2))
>>> F_undeformed  # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
Matrix([
[1, 0],
[0, 1]])

In order to integrate over a volume, the volume change is required.
The change in volume is computed using the determinant of the deformation gradient,
which is

>>> sy.det(F_undeformed)
1

for the undeformed state. So the volume does not change.

Deformed
--------

However, for a deformation of :math:`\varepsilon_{xx}=0.2; \gamma_{xy}=0.1`,
the deformation gradient is

>>> F_deformed = expr.eval_F(2, sy.Matrix([[0.2, 0.1], [0, 0]]))
>>> F_deformed  # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
Matrix([
[1.2, 0.1],
[  0,   1]])

and its determinat

>>> float(sy.det(F_deformed))
1.2

states that the volume increases by a factor of 1.2.
Note, that the tensile strain :math:`\varepsilon_{xx}=0.2` results in a volume change.
The simple shear :math:`\gamma_{xy}=0.1` does not contribute to the volume change.

Pure shear
----------

However, a pure shear deformation of :math:`\gamma_{xy}=0.1` without a tensile strain
would decrease the volume by a factor of

>>> float(sy.det(expr.eval_F(2, sy.Matrix([[0.0, 0.1], [0.1, 0]]))))
0.99

Transformation
--------------

The deformation gradient is a linear transformation from the
undeformed to the deformation state.
Considering a simple shear deformation,

>>> F = expr.eval_F(2, sy.Matrix([[0.0, 0.1], [0.0, 0]]))
>>> F  # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
Matrix([
[1.0, 0.1],
[  0,   1]])

the points of the undeformed unit square can be transformed to
the deformed state by a vector matrix multiplication.
The "@"-symbol stands for a matrix multiplication.

>>> X_unit_square = np.array([
...     [0, 0],
...     [0, 1],
...     [1, 1],
...     [1, 0],
... ])
>>> X_deformed_unit_square = X_unit_square @ F.T
>>> X_deformed_unit_square
Matrix([
[  0, 0],
[0.1, 1],
[1.1, 1],
[1.0, 0]])

For the back transformation from the deformed to the undeformed state,
the inverse of the deformation gradient is used.

>>> np.array(X_deformed_unit_square @ F.inv().T, dtype=float)
array([[0., 0.],
       [0., 1.],
       [1., 1.],
       [1., 0.]])

"""
import doctest

if __name__ == '__main__':
    doctest.testmod()
