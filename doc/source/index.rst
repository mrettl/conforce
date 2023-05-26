.. image:: conforce_logo_icon.png
    :alt: conforce logo

Conforce supports element types defined in :py:attr:`conforce_shared.cf_c.map_type_to_info` directly.

Furthermore, some element types are treated like a supported element type.
This allows to apply also not directly supported element types.
The dictionary defined in :py:attr:`conforce_shared.element_type_mapping.map_abaqus_element_type_to_supported_element_type`
defines which element types can be treated as they were a supported element type.

For a short theory on configurational forces see :py:mod:`conforce.expressions`.


Table of Contents
-----------------

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
