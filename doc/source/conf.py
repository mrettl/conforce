# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup ---------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import shutil

sys.path.insert(0, os.path.abspath('../..'))
sys.path.append(os.path.abspath("./dummy_packages"))

import conforce_shared

# -- Project information ------------------------------------------------------

project = conforce_shared.project
project_copyright = conforce_shared.project_copyright
author = conforce_shared.author
version = conforce_shared.version
release = conforce_shared.version

# -- General configuration ----------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',  # autodocumentation module
    'sphinx.ext.intersphinx',  # Link to other projects’ documentation
    'sphinx.ext.autosummary',  # make autosummarys
    'sphinx.ext.viewcode',  # Add links to highlighted source code
    'sphinx.ext.todo',  # Support for to do
    'sphinx.ext.mathjax',  # Render math via JavaScript
    'sphinx.ext.doctest',  # Test snippets in the documentation
    'recommonmark',  # For markdown syntas (.md-files)
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
add_module_names = False

# If true, the index is generated twice: once as a single page with all the entries, 
# and once as one page per starting letter. Default is False.
html_split_index = True

# -- Options for HTML output --------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'pyramid'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# add table of content to sidebar
html_sidebars = {
    '**': ['searchbox.html', 'relations.html', 'sourcelink.html', 'globaltoc.html', 'localtoc.html']
}

# add logo
html_logo = "_static/conforce_logo_icon_small.png"

# -- Extension configuration --------------------------------------------------
autoclass_content = "init"

# -- Options for intersphinx extension ----------------------------------------
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'sympy': ('https://docs.sympy.org/latest', None)
}

# -- Options for to do extension ----------------------------------------------
todo_include_todos = True

# -- Options for autosummary extension ----------------------------------------
autosummary_generate = True

# -- Options for doctest extension --------------------------------------------
doctest_global_setup = """
import os
os.chdir('..')
"""
doctest_global_cleanup = """
import os
os.chdir('doc')
"""

# -- Create _static dir -------------------------------------------------------

static_dir_path = os.path.abspath(os.path.join(__file__, "../_static"))
if not os.path.exists(static_dir_path):
    print('make directory ' + static_dir_path)
    os.mkdir(static_dir_path)

# -- Delete generated toc -----------------------------------------------------
generated_toc_path = os.path.abspath(os.path.join(__file__, "generated"))
if os.path.exists(generated_toc_path):
    shutil.rmtree(generated_toc_path)

# -- Copy files for root directory --------------------------------------------

for src, dest in [
    ("../../../README.md", "../README.md"),
    ("../../../conforce_logo_icon.png", "../conforce_logo_icon.png"),
    ("../../../plugin_gui.png", "../plugin_gui.png"),
    ("../../../CONTRIBUTING.md", "../CONTRIBUTING.md"),
    ("../../../LICENSE.txt", "../LICENSE.txt"),
    ("../../../examples/Example_1_Two-phase_bar/README.rst", "../example_1.rst"),
    ("../../../examples/Example_1_Two-phase_bar/example_1_images", "../example_1_images"),
]:
    source = os.path.abspath(os.path.join(__file__, src))
    destination = os.path.abspath(os.path.join(__file__, dest))

    mode = ""

    if not os.path.exists(destination):
        mode = "add"

    else:
        source_time = os.path.getmtime(source)
        destination_time = os.path.getmtime(destination)

        if source_time > destination_time:
            mode = "update"

    if mode != "":
        print(f'{mode} {dest}')
        if os.path.isdir(source):
            if os.path.exists(destination):
                shutil.rmtree(destination)
            shutil.copytree(
                src=source,
                dst=destination,
                dirs_exist_ok=False
            )
        else:
            shutil.copyfile(
                src=source,
                dst=destination
            )
