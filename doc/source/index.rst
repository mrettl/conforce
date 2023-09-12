Welcome to the conforce documentation!
======================================

.. image:: conforce_logo.png
    :width: 400
    :alt: conforce logo

Conforce computes computes nodal configurational (or material) forces from the displacement, stress, and energy density fields
for finite elements. It is implemented as a Python 3 package and as an Abaqus plug-in.
The implementation is based on the paper "On material forces and finite element discretizations" by Mueller [1]_.
Configurational forces correspond to the change in energy density when some geometry measure is changed.
Common application of configurational forces is in fracture mechanics,
where they can be used to obtain the energy release rate of cracks.

How to read the documentation
=============================

The :doc:`README` section
describes how to install and use the Python 3 package and the Abaqus plug-in.
Conforce supports 3D elements and 2D plane strain elements for static loads.
Plane stress elements are supported with a simplification.
The :doc:`supported_element_types` section provides a complete list of supported element types.
The theory section describes in detail how configurational forces are computed.
:doc:`examples/example_1_two_phase_bar` provides a detailed walkthrough for the Abaqus plug-in.
:doc:`examples/example_2_ct_specimen_linear_elastic` and
:doc:`examples/example_3_ct_specimen_elastic_plastic`
use the Abaqus plug-in within an automated script to compute the J-integral for a crack.
The library references
(:doc:`reference`,
:doc:`reference_abq`,
:doc:`reference_gen`) contain documentation for each function.
The :doc:`CONTRIBUTING` section is intended for developers who want to add new features or fix bugs.


Table of Contents
=================

.. toctree::
    :caption: Getting started
    :maxdepth: 1

    README.md
    supported_element_types

.. toctree::
    :caption: Theory
    :maxdepth: 1

    theory/Configurational_forces.rst
    theory/Deformation_gradient.rst
    theory/Finite_elements.rst
    theory/First_Piola_Kirchhoff_stress.rst
    theory/Abbreviations.rst

.. toctree::
    :caption: Examples
    :maxdepth: 1

    examples/example_1_two_phase_bar
    examples/example_2_ct_specimen_linear_elastic
    examples/example_3_ct_specimen_elastic_plastic

.. toctree::
    :caption: Library references
    :maxdepth: 1

    reference
    reference_gen
    reference_abq

.. toctree::
    :caption: Appendix
    :maxdepth: 1

    CONTRIBUTING.md
    license


Reference
=========

.. [1] R. Mueller and G. A. Maugin,
    “On material forces and finite element discretizations",
    Computational Mechanics, vol. 29, no. 1, pp. 52–60, Jul. 2002, doi: `10.1007/s00466-002-0322-2 <https://doi.org/10.1007/s00466-002-0322-2>`_.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
