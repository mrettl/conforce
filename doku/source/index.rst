CF3D: Implementing Configurational Forces in FEM using automatic code generation
================================================================================

This package provides a toolset to derive the numerical implementation of configurational forces for Finite Elements
from a symbolic set of equations. 
Therefore the analytic equations are discretized on Finite elements by use of the corresponding shape functions.
This is done for static cases, including plasticity.
The symbolic expressions are element independent; the implementations for 
various element types are generated as C-code by the package in a fully automated way. This also includes the generation of a
Python wrapper, which simplifies the usage from Python.

Deatils of the symbolic expression are given in :ref:`Example`.

After the numerical implementation is generated and compiled, 
have a look at the :ref:`Interface_description` and the application examples.

Quick start
-----------

To use CF3D, first check the :ref:`prerequisites <Prerequisites>` of your Python installation.
All files needed for code generation are located in ./cf3D/src. This files can be copied to a user-defined folder,
an installation isn't necessary to use the package.

For setting up the Python environment, we recommend to use a
`virtual Conda environment <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands>`_.

The generated code and the ctypes wrapper can also be used from other Python environments, which includes deprecated
Python versions like 2.7 or environments with limited expandability like Abaqus Python.

Application examples:
---------------------
The two application example including a already generated and acompiled numerical implementation of configurational forces are located in ./cf3D/examples.

| :ref:`Two_phase_bar`
| :ref:`CT_specimen`


.. _Prerequisites

Prerequisites
-------------

It is recommended to use virtual environments (`anaconda <https://www.anaconda.com/>`_).
This package is written and tested in Python 3.7 and relies on the following packages.

| `numpy <https://numpy.org/>`_ 1.18.1
| `sympy <https://www.sympy.org/en/index.html>`_ 1.5.1
| `clang <https://clang.llvm.org/>`_ 11.0.0

Contributing
------------

Clone the repository and add changes to it. Test the changes and make a pull request.

Modules
-------

.. autosummary::
    :toctree: generated/
    :template: module_template.rst

    Auxiliary_functions
    ele_def


.. toctree::
    :hidden:
    
    Example
    Interface_description


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
