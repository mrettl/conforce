Welcome to the conforce documentation!
======================================

.. image:: conforce_logo_icon.png
    :width: 400 px
    :alt: conforce logo

Conforce is a Python package for computing configurational forces
out of the displacement-, stress-, and energy-density-fields obtained by Finite Element Analysis.

Features
========

- Supported methods:
    - deformation based formulation (see :py:func:`conforce.expressions.eval_CS_dbf`)
    - motion based formulation (see :py:func:`conforce.expressions.eval_CS_mbf`)
- Applicable element types:
    - Conforce supports element types defined in :py:attr:`conforce_shared.cf_c.map_type_to_info`.
    - Furthermore, element types not directly supported can be treated like an other element type that is supported.
      This allows to use not supported element types with some assumptions and simplifications.
      The dictionary :py:attr:`conforce_shared.element_type_mapping.map_abaqus_element_type_to_supported_element_type`
      defines which supported element types can imitate not directly supported elements.
- Material orientations are supported:
    - Field outputs for the displacement and stress field in the global coordinate system are computed automatically.
- Static configurational stresses and forces may consider:
    - Elastic energy density (Abaqus field output SENER)
    - Plastic energy density (Abaqus field output PENER)
    - Viscous energy density (Abaqus field output VENER)
    - `Kinematic energies are not yet supported.`

:py:mod:`conforce.expressions` provides a short description of the computation.


Table of Contents
=================

.. toctree::
    :maxdepth: 2
	
    README.md
    CONTRIBUTING.md

    examples

    references

    license


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
