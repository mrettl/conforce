# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys


# -- Project import -----------------------------------------------------------

# inserts the absolute path if the sampleproject_KKV_slim folder into the path-variables to be able to import the
# project if it is not installed in the local python distribution. If it is installed, just import it.
# ATTENTION!!
# In order to work, every module, class or function which should be included in the documentation has to be accessible!
# If something is not imported, it will not be included in the documentation
# If you do not import submodules in the __init__.py of the module, import it here.

project_dir = os.path.abspath(os.path.join(__file__, "../../../src"))
sys.path.insert(0, project_dir)


# -- General configuration ---------------------------------------------------

needs_sphinx = '1.1'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.


extensions = ['sphinx.ext.autodoc',  # autodocumentation module
              # 'sphinx.ext.imgmath',  # mathematical expressions can be rendered as png images
              'numpydoc',           # documentation in numpy-style
              'sphinx.ext.napoleon',  # numpy docstrings
              'sphinx.ext.intersphinx',
              'sphinx.ext.coverage',
              'sphinx.ext.autosummary',  # make autosummarys
              'sphinx.ext.viewcode',
              ]


# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

source_suffix = '.rst'
master_doc = 'index'
pygments_style = 'sphinx'

# -- Render options ---------

# imgmath_image_format = 'svg'  # This would render mathematical expressions as svg images
imgmath_latex_preamble = r'\usepackage{xcolor}'  # mathematical expressions can be colored

# -- Options for HTML output ---------------------------------------------------

html_theme = 'scipy_theme'
html_theme_path = ['.']
# html_logo = '_static/scipyshiny_small.png'
# html_static_path = ['_static']
html_theme_options = {
    "edit_link": False,
    "sidebar": "right",
    "scipy_org_logo": False,
    "rootlinks": []
    # "rootlinks": [("http://scipy.org/", "Scipy.org"),
    #               ("http://docs.scipy.org/", "Docs")]
}


# -- autosummary settings -------------------------------------------------------
autosummary_generate = True
autosummary_imported_members = False

# numpydoc options
numpydoc_show_inherited_class_members = False

# -- Project information -----------------------------------------------------

project = 'ConF3D'
copyright = '2022, Markus Tauscher'
author = 'Markus Tauscher'
version = '0.1'
release = '0.1'