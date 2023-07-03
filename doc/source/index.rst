Welcome to the conforce documentation!
======================================

.. image:: conforce_logo.png
    :width: 400
    :alt: conforce logo

Conforce is available as Python 3 package and as Abaqus Plug-in.
It computes configurational forces from the displacement-, stress-, and energy-density-fields
for finite elements.
The implementation is based on the paper "On material forces and finite element discretizations" by Mueller [1]_.
Configurational forces correspond to the change of the energy density if some geometry measure changes.
A common use case for configurational forces is fracture mechanics
where they correspond the energy release if a crack grows.

How to read the documentation
=============================

The :doc:`README` section
describes how to install and use the Python 3 package and the Abaqus Plug-in.
Conforce supports 3D and 2D plane strain elements in a static load case.
Plane stress elements are supported with a simplification.
The :doc:`supported_element_types` section provides a full list of supported element types.
The :doc:`theory` section describes how configurational forces are computed in detail.
:doc:`examples/example_1_two_phase_bar` provides a detail walk-through for the Abaqus Plug-in.
:doc:`examples/example_2_ct_specimen_linear_elastic` and
:doc:`examples/example_3_ct_specimen_elastic_plastic`
use the Abaqus Plug-in within an automated script to compute the J-integral for a crack.
The library references
(:doc:`reference`,
:doc:`reference_abq`,
:doc:`reference_shared`) contain documentation for every function.
The :doc:`CONTRIBUTING` section is intended for developers who want to add new features or fix bugs.


Table of Contents
=================

.. toctree::
    :caption: Getting started
    :maxdepth: 1
	
    README.md
    supported_element_types
    theory

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
    reference_abq
    reference_shared

.. toctree::
    :caption: Appendix
    :maxdepth: 1

    CONTRIBUTING.md
    license


Reference
=========

.. [1] R. Mueller and G. A. Maugin,
    "On material forces and finite element discretizations"
    Computational Mechanics, vol. 29, no. 1, pp. 52-60, Jul. 2002, doi: 10.1007/s00466-002-0322-2.



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
