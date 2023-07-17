r"""
Configurational forces are the derivative of the Helmholtz energy with respect to a change in the geometry.

Helmholtz free energy density :math:`e`
---------------------------------------

The Helmholtz free energy density :math:`e` is defined as

.. math::
    e = \psi - \textrm{temperature} \cdot \textrm{entropy}

where :math:`\psi` is the internal energy density.
For a hyper-elastic material, the Helmholtz free energy density is the strain energy density.
The partial derivative of :math:`e` with respect to the deformation gradient :math:`F=(f_{ij})`
at a fixed position :math:`X=(x_{i})` is the first Piola-Kirchhoff stress tensor. [1]_

.. math::
    :label: eq_partial_e_partial_F

    p_{ij} = \frac{\partial e}{\partial f_{ij}}

Furthermore, the Helmholtz free energy density for a hyper-elastic material is
a function of

.. math::
    :label: eq_e_func

    e = e(X, F)

the two independent quantities :math:`X` and :math:`F`.
:math:`X` is a location in the undeformed state and
:math:`F` is the deformation gradient. [2]_


Static force equilibrium
------------------------

According to the static equilibrium, the sum of all forces is zero.
In the deformed state, static equilibrium is given by

.. math::
    :label: eq_stat_force_balance_deformed

    \frac{d \sigma_{ij}}{d x_{\textrm{deformed}, j}}
    = -b_{\textrm{deformed}, i}

where :math:`\sigma` is the true Cauchy stress tensor and
:math:`B_{\textrm{deformed}} = (b_{\textrm{deformed}, i})`
is a body force acting on the deformed shape.

The balance of forces can also be written in terms of the undeformed state.
Then the first Piola-Kirchhoff stress :math:`P=(p_{ij})` is used and the equilibrium is

.. math::
    :label: eq_stat_force_balance_undeformed

    \frac{d p_{ij}}{d x_{j}}
    = -b_{\textrm{undeformed}, i}

with the body force :math:`B_{\textrm{undeformed}}=(b_{\textrm{undeformed}, i})`.

Although, :math:`B_{\textrm{deformed}}` and  :math:`B_{\textrm{undeformed}}` are defined on either
the deformed or undeformed state, they both act on the **displacement** field.
An alternative idea is to define a body force that acts on the
**geometry** instead on the **displacements**.

Configurational body force `G`
------------------------------

This is exactly what the configurational body force :math:`G=(g_{i})` does.
As with the other two balance laws
:eq:`eq_stat_force_balance_deformed` and :eq:`eq_stat_force_balance_undeformed`,
the divergence of a stress tensor (we call it the configurational or Eshelby stress :math:`CS=({cs}_{ij})` [3]_)
corresponds to the negative configurational body force :math:`G`.

.. math::
    :label: eq_cf_body_forces_balance

    \frac{d {cs}_{ij}}{d x_{j}} = -g_{i}

The static configurational body force is defined as

.. math::
    :label: eq_cf_body_forces

    g_{i}
    = -\frac{\partial e}{\partial x_{i}}
    -f_{ji} \cdot b_{\textrm{undeformed}, j}

.. note:

    Gurtin [4]_ defines configurational body forces in an alternative way.
    He splits the configurational force into an internal contribution :math:`-\frac{\partial e}{\partial x_{i}}`
    and an external contribution :math:`-f_{ji} \cdot b_{\textrm{undeformed}, j}`.

If not body force :math:`B_{\textrm{undeformed}}` like a pressure or gravity is present,
the configurational body force is :math:`-\frac{\partial e}{\partial x_{i}}`,
the partial derivative of the Helmholtz energy density :math:`e`
with respect to the position :math:`X`
under a constant deformation gradient :math:`F`.

This is the derivative with respect to a position,
but not with respect to a change in **geometry**.
However, with a little trick, :math:`dx` can be interpreted as :math:`-dl`
where :math:`dl` stands for geometry change of e.g. the position of the crack tip or the position of an inclusion.

The image below demonstrates in which cases the derivatives can be replaced :math:`dx = -dl`.
The eye observes a region near a geometric feature (ellipse).
The left side shows the original state and the right side shows that something has moved.
In (a) it is not possible to tell whether the eye has moved by :math:`dx`
or the ellipse moved by :math:`-dl`.
In this case :math:`dx = -dl` holds true and the configurational force acting on the ellipse
can be interpreted as the negative derivative of the energy density with respect to a change in geometry.

.. image:: theory_images/similarity.png
    :width: 600
    :alt: similarity

Unlike (a), in (b) there is a boundary condition in the region near the ellipse.
The boundary condition pins a piece of the material to this fixed position.
In this case we can see whether the ellipse or the eye is moving.
If the eye moves, the boundary condition moves with the ellipse,
if the ellipse moves, the boundary conditions remain in the same place.
In this case :math:`dx` **cannot** be interpreted as a change in geometry :math:`-dl`.

Even though the configurational body force :math:`G` corresponds to a geometry change,
the derivative :math:`\frac{\partial e}{\partial x_{i}}` cannot be computed in a straightforward way.
Instead, equation :eq:`eq_cf_body_forces_balance` is used to compute :math:`G` using the configurational stress tensor :math:`CS`.

Configurational stress tensor `CS`
----------------------------------

Motion based formulation (mbf)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The configurational stress tensor is also known as the energy-momentum tensor
and was first introduced by Eshelby [3]_.
This section explains how the configurational stress tensor :math:`CS`
and subsequently the configurational body force :math:`G` are computed
using the **motion based formulation**.

First, derive the Helmholtz energy density :math:`e` with respect to a :math:`X`.
According to the chain rule, the total derivative with respect to :math:`X` is equal to
the partial derivatives with respect to the independent quantities :math:`X` and :math:`F`
times the total derivatives of these quantities with respect to :math:`X`.

.. math::
    :label: eq_de_dx

    \frac{d e}{d x_{i}}
    = \frac{\partial e}{\partial x_{i}}
    + \frac{\partial e}{\partial f_{jk}} \cdot \frac{d f_{jk}}{d x_{i}}

The term :math:`\frac{\partial e}{\partial f_{jk}}` can be replaced
by the first Piola-Kirchhoff stress tensor :math:`P` using equation :eq:`eq_partial_e_partial_F`.

.. math::
    :label: eq_de_dx_2

    \frac{d e}{d x_{i}}
    = \frac{\partial e}{\partial x_{i}}
    + p_{jk} \cdot \frac{d f_{jk}}{d x_{i}}


.. note::

    **Identity from Gurtin [page 16, equation (1-26)]** [4]_

    .. math::

        P : \textrm{grad}(F)
        = \textrm{div}(F^{T} P)
        - F^{T} \textrm{div}(P)

    The identity in the summation notation is

    .. math::
        :label: eq_p_df_dx

        p_{jk} \frac{d f_{jk}}{d x_{i}}
        = \frac{d (f^{T}_{kj} \cdot p_{ij}) }{d x_{i}}
        - f^{T}_{kj} \frac{d p_{ji}}{d x_{i}}


Insert the above identity :eq:`eq_p_df_dx` into :eq:`eq_de_dx_2`.

.. math::
    :label: eq_de_dx_3

    \frac{d e}{d x_{i}}
    = \frac{\partial e}{\partial x_{i}}
    + \frac{d (f^{T}_{kj} \cdot p_{ij}) }{d x_{i}}
    - f^{T}_{kj} \frac{d p_{ji}}{d x_{i}}


The term :math:`\frac{d p_{ji}}{d x_{i}}` occurs in the force balance
and can be replaced by :math:`B_{\textrm{undeformed}}` using
equation :eq:`eq_stat_force_balance_undeformed`.

.. math::
    :label: eq_de_dx_4

    \frac{d e}{d x_{i}}
    = \frac{\partial e}{\partial x_{i}}
    + \frac{d (f^{T}_{kj} \cdot p_{ij}) }{d x_{i}}
    + b_{\textrm{undeformed}, j} \cdot f^{T}_{kj}

Collecting the total derivatives yields

.. math::
    :label: eq_de_dx_5

    \left(
        \frac{d e}{d x_{i}}
        - \frac{d (f^{T}_{kj} \cdot p_{ij}) }{d x_{i}}
    \right)
    + \left(
        - \frac{\partial e}{\partial x_{i}}
        - b_{\textrm{undeformed}, j} \cdot f^{T}_{kj}
    \right)
    = 0

Pull out the total derivative in the left bracket.
Note the kronecker delta :math:`\delta_{ik}` nearby the Helmholtz energy density.

.. math::
    :label: eq_de_dx_5_b

    \frac{ d (e \delta_{ik} - f^{T}_{kj} \cdot p_{ij})}{d x_{i}}
    + \left(
        - \frac{\partial e}{\partial x_{i}}
        - b_{\textrm{undeformed}, i} \cdot f^{T}_{kj}
    \right)
    = 0

The right bracket is exactly the definition of the configurational body force from equation :eq:`eq_cf_body_forces`.
Inserting this equation yields

.. math::
    :label: eq_de_dx_6

    \frac{ d (e \cdot \delta_{ik} - f^{T}_{kj} \cdot p_{ij})}{d x_{i}}
    = -g_{i}

This equation is the (static) configurational force balance known from equation :eq:`eq_cf_body_forces_balance`.
Consequently, the configurational stress has to be

.. math::
    :label: eq_cs_explicit

    cs_{ik} = e \cdot \delta_{ik} - f^{T}_{kj} \cdot p_{ij}

Finally, the configurational body force :math:`G` is computed

.. math::
    :label: eq_de_dx_7

    \frac{d {cs}_{ik}}{d x_{k}} = -g_{k}


Displacement based formulation (dbf)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An alternative formulation is the **displacement based formulation**
of the configurational stress tensor [4]_

.. math::
    :label: eq_cs_dbf_explicit

    cs_{dbf, ik} = e \cdot \delta_{ik} - \frac{u_{k}}{x_{j}} p_{ji}

According to Gurtin [4]_, the dbf and mbf quantities can be transformed to each other using the
first Piola-Kirchhoff stress :math:`P` and the body force :math:`G_{\textrm{undeformed}}`.

.. math::

    cs_{mbf, ik} = cs_{dbf, ik} - P_{ik}

.. math::

    g_{mbf, i} = g_{dbf, i} - b_{\textrm{undeformed}, i}

Within a material, the configurational body forces remain the same.
The configurational stresses differ by the second Piola-Kirchhoff stress tensor.
In the absence of body forces, the configurational body forces are identical.


Plasticity
^^^^^^^^^^

The implemented formulation does **not** support plasticity in general,
since the Helmholtz energy density not only depends on :math:`X` and :math:`F`,
but also on plastic hardening parameters which are not considered here.
Furthermore, it would be necessary to split the deformation gradient :math:`F`
into an elastic and a plastic part.

However, under the assumption of small strain plasticity, the formulation is valid
and has already been used by Kolednik [5]_.
Kolednik proposes two modifications for the configurational stress:

- incremental plasticity, which considers only the elastic strain energy density (SENER)

.. math::

    cs_{\textrm{ep}, ik}
    = e_{\textrm{elastic}} \cdot \delta_{ik}
    - f^{T}_{kj} \cdot p_{ij}

- deformation plasticity that considers both elastic and plastic strain energy densities (SENER+PENER)

.. math::

    cs_{\textrm{nlel}, ik}
    = (e_{\textrm{elastic}} + e_{\textrm{plastic}}) \cdot \delta_{ik}
    - f^{T}_{kj} \cdot p_{ij}


Configurational forces `CF`
---------------------------

The configurational force :math:`CF` corresponds to movement of a volume.
In our implementation, :math:`CF` is associated to element nodes and therefore their corresponding volume.
The (nodal) configurational force is therefore the volume integral

.. math::
    :label: q_cf_0

    \textrm{cf_at_node}_{ij}
    = \int_{\textrm{body}} g_{j} \cdot h_{i} \; dV

over the configurational body force weighted by the shape function of the :math:`i`-th node.
Since :math:`G` is unknown, the configurational force balance from equation :eq:`eq_cf_body_forces_balance`
is used to evaluate the integral.
The configurational force balance is true at each position.
Consequently, the (weighted) integral over the balance is

.. math::
    :label: q_cf_2

    \int_{\textrm{body}} g_{j} \cdot h_{i} \; dV
    = \int_{\textrm{body}} \frac{d {cs}_{jk}}{d x_{k}} \cdot h_{i} \; dV

true for each test/shape function :math:`H=(h_{i})`.
As mentioned above, we use the element shape functions :math:`H` as test/shape function.
The left side is exactly the integral we want to evaluate, but contains the unknown configurational body force.
The right side contains the configurational stress that can be computed.
A partial integration transfers the derivative of the configurational stress :math:`CS`
to the shape function.

.. math::
    :label: q_cf_3

    \int_{\textrm{body}} g_{j} \cdot h_{i} \; dV
    = \int_{\textrm{body}} {cs}_{jk} \cdot \frac{d h_{i}}{d x_{k}} \; dV
    - \int_{\textrm{surface}} {cs}_{jk} \cdot n_{k} \cdot h_{i} \; dA

However, we assume that the surface integral is zero [2]_ [5]_.
This is a valid approach if the shape functions are zero at the surface.
This is the case for all nodes that are not directly on the surface.
We compute the configurational force purely from the volume integral

.. math::
    :label: q_cf_4

    \int_{\textrm{body}} g_{j} \cdot h_{i} \; dV
    = \int_{\textrm{body}} {cs}_{jk} \cdot \frac{d h_{i}}{d x_{k}} \; dV

If you also want to evaluate configurational forces on the surface,
we provide methods to compute the configurational stresses.
Just use them to evaluate the surface integral and subtract it
from our configurational forces.


References
----------

.. [1] J. S. Bergstrom,
    “Mechanics of Solid Polymers: Theory and Computational Modeling”.
    Elsevier, 2015.


.. [2] R. Mueller and G. A. Maugin,
    “On material forces and finite element discretizations,”
    Computational Mechanics, vol. 29, no. 1, pp. 52–60, Jul. 2002, doi: `10.1007/s00466-002-0322-2 <https://doi.org/10.1007/s00466-002-0322-2>`_.

.. [3] J. D. Eshelby,
    “The force on an elastic singularity,”
    Philosophical Transactions of the Royal Society of London. Series A, Mathematical and Physical Sciences, vol. 244, no. 877, pp. 87–112, 1951.

.. [4] M. E. Gurtin,
    Configurational forces as basic concepts of continuum physics.
    in Applied mathematical sciences, no. 137. New York: Springer, 2000.

.. [5] O. Kolednik, R. Schöngrundner, and F. D. Fischer,
    “A new view on J-integrals in elastic–plastic materials,”
    Int J Fract, vol. 187, no. 1, pp. 77–107, May 2014, doi: `10.1007/s10704-013-9920-6 <https://doi.org/10.1007/s10704-013-9920-6>`_.

"""
import doctest

if __name__ == '__main__':
    doctest.testmod()
