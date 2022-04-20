# encoding: utf-8
#
# Description: Sample code to demonstrate the functionality of sphinx and show documentation styles
#
# Author: Matthias Drvoderic
#
# ----------------------------
# (c) 2020, Matthias Drvoderic
# ----------------------------
"""
Headline of module1.py
======================

Documentation of this module can be put here. Of course, links to other objects like :class:`DocExampleClass`
available also work from here. You

Here is a good place to put a overview or maybe a theoretical introduction of this module.

E.g have a look at the source files of scipy and numpy. All the documentation of the modules is done this way.

A module is an object with arbitrary named attributes that can be referred to. Simply, normally one *.py* file is one
module. Classes, functions etc. which work together or form a logical group should be in one module.

The *__init__.py* file marks a directory as a python package directory. The import statement :code:`import
sample_project` will look into the *sample_project* folder and look for the *__init__.py* file and run it. If no
*__init__.py* is found, python will no longer look for submodules in this directory. The *__init__.py* is normally
empty but can be used to import submodules with more convenient names, provide convenience functions, make checks e.g
check if the right python version is installed or if necessary packages are installed etc. """

import matplotlib.pyplot as plt


class DocExampleClass(object):
    """
    This is the first line of docstring which should provide a very short explanation.

    Here, a more detailed explanation of the object should be provided. Maybe
    a short theoretical introduction on what it does and what the limitations ar
    would be fitting here.

    Parameters
    ----------
    x : float
        Some value
    y : float
        Some value

    Examples
    --------
    Examples are a great way to quickly show how it can be used.

    >>> 1 + 2
    3
    """

    def __init__(self, x, y):
        pass

    def method1(self, z):
        """
        Also just for documentation purpose

        Parameters
        ----------
        z : float
            Some value

        Returns
        -------
        result : float
            The result of this complicated computation will be None ;)
        """
        pass

    def method2(self):
        """
        Also does nothing

        Just for demonstration of the automated documentation!

        Returns
        -------
        x : float
            In fact, does not return anything..

        Raises
        ------
        ValueError
            Not really, only for demonstration ;)
        """
        pass


def sample_function(x):
    """
    Example of how to document a function. The first line is, as always a quick description

    Then a block with a more detailed description should follow

    Parameters
    ----------
    x : float
        Here all parameters of the function have to be listed.

    Returns
    -------
    y : float
        The return of the function.

    Raises
    ------
    ValueError
        If something goes wrong in the computation, or one of the inputs is the wrong type, or....

    Notes
    -----
    Links to other pages or sections of your documentation also work from in here. But to link to other *rst* pages,
    an explicit link has to be used. To link to other documented code, just type the name of the class, function etc. in
    quotes e.g. ``:class:`DocExampleClass``` will reference to :class:`DocExampleClass`
    """
    pass


def __function_without_doc(x):
    return x**2
