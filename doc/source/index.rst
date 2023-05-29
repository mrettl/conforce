.. image:: conforce_logo_icon.png
    :alt: conforce logo

Features
========

- Supported methods:
    - deformation based formulation (see :py:func:`conforce.expressions.eval_CS_dbf`)
    - motion based formulation (see :py:func:`conforce.expressions.eval_CS_mbf`)
- Applicable element types:
    - Conforce supports element types defined in :py:attr:`conforce_shared.cf_c.map_type_to_info`.
    - Furthermore, some element types are treated like a supported element type.
      This allows to apply not directly supported element types.
      The dictionary defined in :py:attr:`conforce_shared.element_type_mapping.map_abaqus_element_type_to_supported_element_type`
      defines which element types can be treated as if they would be a supported element type.
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
    reference
    reference_abq



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
