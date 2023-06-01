r"""
Configurational forces are the derivative of the energy with respect to a change in the geometry.

>>> import numpy as np
>>> import sympy as sy
>>> import conforce.expressions as expr
>>> import conforce.element_definitions as el_def


Helmholtz free energy density `e`
---------------------------------

The Helmholtz free energy density `e` is defined as

.. math::
    e = \psi - temperature \cdot entropy

where :math:`\psi` is the internal energy density.
For a hyper-elastic material the Helmholtz free energy density is the strain energy density.
Consequently, the partial derivative of `e` with respect to the deformation gradient `F` is
the first Piola-Kirchhoff stress tensor.

.. math::
    :label: eq_partial_e_partial_F

    p_{ij} = \frac{\partial e}{\partial f_{ij}}

Furthermore, the Helmholtz free energy density for a hyper-elastic material is only
a function of

.. math::
    :label: eq_e_func

    e = e(X, F)

the two independent quantities `X` and `F`.
`X` is a position of the material in the undeformed state and
`F` is the deformation gradient.


Static equilibrium in the reference space
-----------------------------------------

According to the static equilibrium the sum of all forces is zero.
In the deformed state, the static equilibrium is known as

.. math::
    :label: eq_stat_force_balance_deformed

    \frac{\partial \sigma_{ij}}{\partial x_{j}}
    = -b_{deformed, j}

where :math:`\sigma` is the true Cauchy stress tensor and :math:`B_{deformed}`
is a body force acting on the deformed shape.

The balance of forces can also be written in terms of the undeformed state.
Then the first Piola-Kirchhoff stress `P` is used and the equilibrium is

.. math::
    :label: eq_stat_force_balance_undeformed

    \frac{\partial p_{ij}}{\partial x_{j}}
    = -b_{undeformed, i}

with the body force :math:`B_{undeformed}`.

Although, :math:`B_{deformed}` and  :math:`B_{undeformed}` are defined on either
the deformed or undeformed state, they both act on the **displacement** field.
An alternative idea is to define a body force that acts on the
**geometry** instead on the **displacements**.

Configurational body force `G`
------------------------------

This is exactly, what the configurational body force `G` does.
Like for the other two balance laws
:eq:`eq_stat_force_balance_deformed` and :eq:`eq_stat_force_balance_deformed`,
the divergence of some stress tensor (we call it configurational stress `CS`)
corresponds to the negative configurational body force `G`.

.. math::
    :label: eq_cf_body_forces_balance

    \frac{d {cs}_{ij}}{d x_{j}} = -g_{i}

The static configurational body force is defined as

.. math::
    :label: eq_cf_body_forces

    g_{i}
    = - \frac{\partial e}{\partial x_{i}}
    - f_{ji} \cdot b_{undeformed, j}

If no body force :math:`B_{undeformed}` like a pressure or gravity acts on the **displacement** field,
the configurational body forces are the partial derivative of the Helmholtz energy with
respect to a **geometry** change.

However, the term :math:`\frac{\partial e}{\partial x_{i}}` can not be computed in a straight forward way.
Instead, equation :eq:`eq_cf_body_forces_balance` is used to compute `G` with the help
of the configurational stress tensor `CS`.

Configurational stress tensor `CS` (mbf)
----------------------------------------

The configurational stress tensor is also known as energy-momentum tensor
and was first introduced by Eshelby.
This section explains how the configurational stress tensor `CS`
and subsequently the configurational body force `G` are computed
using the **motion based formulation**.

First, derive the Helmholtz energy density `e` with respect to a geometry change `X`.
According to the chain rule, the total derivative with repsect to `X` equals
the partial derivatives with respect to the independent quantities `X` and `F`
times the total derivatives of those quantities with respect to `X`.

.. math::
    :label: eq_de_dx

    \frac{d e}{d x_{i}}
    = \frac{\partial e}{\partial x_{i}}
    + \frac{\partial e}{\partial f_{jk}} \cdot \frac{d f_{jk}}{d x_{i}}

The term :math:`\frac{\partial e}{\partial f_{jk}}` can be replaced
by the first Piola-Kirchhoff stress tensor using equation :eq:`eq_partial_e_partial_F`.

.. math::
    :label: eq_de_dx_2

    \frac{d e}{d x_{i}}
    = \frac{\partial e}{\partial x_{i}}
    + p_{jk} \cdot \frac{d f_{jk}}{d x_{i}}


.. note::

    **Auxiliary calculation**

    .. math::

        \frac{d (p_{jk} \cdot f_{ji}) }{d x_{i}}
        = f_{ji} \cdot \frac{d p_{jk} }{d x_{i}}
        + p_{ji} \cdot \frac{d f_{jk} }{d x_{i}}

    .. math::
        :label: eq_p_partial_f_partial_x

        p_{jk} \cdot \frac{d f_{jk} }{d x_{i}}
        = \frac{d (p_{jk} \cdot f_{ji}) }{d x_{i}}
        - f_{ji} \cdot \frac{d p_{jk} }{d x_{k}}

.. todo:: fix indices


Next, inserting :eq:`eq_p_partial_f_partial_x` into :eq:`eq_de_dx_2` yields

.. math::
    :label: eq_de_dx_3

    \frac{d e}{d x_{i}}
    = \frac{\partial e}{\partial x_{i}}
    + \frac{d (p_{jk} \cdot f_{ji}) }{d x_{i}}
    - f_{ji} \cdot \frac{d p_{jk} }{d x_{i}}


The term :math:`\frac{d p_{jk} }{d x_{i}}` occurs in the force balance
and can be replaced by :math:`B_{undeformed}` using
equation :eq:`eq_stat_force_balance_undeformed`.

.. math::
    :label: eq_de_dx_4

    \frac{d e}{d x_{i}}
    = \frac{\partial e}{\partial x_{i}}
    + \frac{d (p_{jk} \cdot f_{jk}) }{d x_{i}}
    + b_{undeformed, i} \cdot f_{jk}

Collecting the total derivatives yields

.. math::
    :label: eq_de_dx_5

    \left(
        \frac{d e}{d x_{i}}
        - \frac{d (p_{jk} \cdot f_{jk}) }{d x_{i}}
    \right)
    + \left(
        - \frac{\partial e}{\partial x_{i}}
        - b_{undeformed, i} \cdot f_{jk}
    \right)
    = 0

    \frac{ d (e - p_{jk} \cdot f_{jk})}{d x_{i}}
    + \left(
        - \frac{\partial e}{\partial x_{i}}
        - b_{undeformed, i} \cdot f_{jk}
    \right)
    = 0

Note, that the right bracket is exactly the configurational body force from equation :eq:`eq_cf_body_forces`.
Inserting the equation yields

.. math::
    :label: eq_de_dx_6

    \frac{ d (e - p_{jk} \cdot f_{jk})}{d x_{i}}
    = -g

This equation is the (static) configurational force balance known from equation :eq:`eq_cf_body_forces_balance`.
Consequently, the configurational stress has to be

.. math::
    :label: eq_cs_explicit

    cs_{jk} = e - p_{jk} \cdot f_{jk}

Finally, the configurational body force `G` is computed

.. math::
    :label: eq_de_dx_7

    \frac{d {cs}_{jk}}{d x_{i}} = -g_{i}


Configurational stress tensor `CS` (dbf)
----------------------------------------

An alternative formulation is the **displacement based formulation**
of the configurational stress tensor

.. math::
    :label: eq_cs_dbf_explicit

    cs_{dbf, jk} = e - \frac{u_{k}}{x_{j}} \cdot f_{jk}

The dbf and mbf quantities can be transformed to each other.

.. math::

    cs_{mbf, jk} = cs_{dbf, jk} - \textrm{second-Piola-Kirchhoff tensor}

.. math::

    g_{mbf, i} = g_{dbf, i}

Inside a material the configurational body forces stay the same.
The configurational stresses differ by the second Piola-Kirchhoff stress tensor.

However, at the boundaries of a material the configurational body forces might deviate.


Configurational forces `CF`
---------------------------

The configurational forces `CF` are defined as the integral over the
configurational body forces weighted by a shape function.
We use the element shape function to compute configurational forces for each node.
This leads to the integral

.. math::

    \textrm{cf_at_node}_{j}
    = \int h \cdot g_{j} \; dV
    = \int \frac{d h}{d x_{i}} \cdot cs_{ij} \; dV



"""
import doctest

if __name__ == '__main__':
    doctest.testmod()
