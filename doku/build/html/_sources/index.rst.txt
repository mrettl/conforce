CF3D: Implementing Configurational Forces in FEM using automatic code generation
==========================================

This package provides a toolset to derive the numerical implementation of Configurational Forces on Finite Elements
from a symbolic set of equations. The symbolic expressions are element independent, the implementations for 
various element types are generated as C-code by the package in a fully automated way. This also includes the generation of a
Python wrapper, which simplifies the usage from Python.
For the linear elastic and elastic cases, an implementation is already provided within this package.

For a detailed explanation of the symbolic expression have a look at :ref:`Example`.

Quick start
-----------

To install CF3D, check at first the :ref:`Prerequisites` of your python installation.
Upon meeting all the criteria, the package can be installed with pip, or you can clone or download the repo.
For setting up the Python environment, we recommend to virtual Conda environments. See
`the conda guide <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands>`_
for infos about creating and managing Conda environments.

The generated code and the ctypes wrapper can also be used from other Python environments, which also includes deprecated
Python versions like 2.7 or environments with limited expandability like Abaqus Python.

Prerequisites
-------------

It is recommended to use virtual environments (`anaconda <https://www.anaconda.com/>`_).
This package is written and tested in Python 3.7 and relies on here listed packages.

| `numpy <https://numpy.org/>`_ 1.18.1
| `sympy <https://www.sympy.org/en/index.html>`_ 1.5.1

Contributing
------------

Clone the repository and add changes to it. Test the changes and make a pull request.

Modules
-------

.. autosummary::
    :toctree: generated/
    :template: module_template.rst

    ConF3D.Auxiliary_functions
    ConF3D.ele_def


.. toctree::
    :hidden:
    
    Example


Authors
-------
- Markus Tauscher

License
-------

This project is licensed under the MIT License.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
