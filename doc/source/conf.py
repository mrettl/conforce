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

DOC_SOURCE_DIR = os.path.abspath(__file__ + '/..')
DOC_DIR = os.path.abspath(DOC_SOURCE_DIR + '/..')
HOME_DIR = os.path.abspath(DOC_SOURCE_DIR + '/../..')

print(__file__)
print(DOC_SOURCE_DIR)
print(DOC_DIR)
print(HOME_DIR)

sys.path.append(os.path.abspath(DOC_SOURCE_DIR + "/dummy_packages"))
sys.path.insert(0, HOME_DIR + "/examples/Example_1_Two_phase_bar")
sys.path.insert(0, HOME_DIR + "/examples/Example_2_CT_specimen_linear_elastic")
sys.path.insert(0, HOME_DIR + "/examples/Example_3_CT_specimen_elastic_plastic")
sys.path.insert(0, HOME_DIR)

import conforce_shared
from conforce_shared import cf_c
from conforce_shared import element_type_mapping

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
    'recommonmark',  # For markdown syntas (.md-files)
    'sphinx_rtd_theme'  # Read-the-Docs theme
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
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# add table of content to sidebar
html_sidebars = {
    '**': ['searchbox.html', 'relations.html', 'sourcelink.html', 'globaltoc.html', 'localtoc.html']
}

# add logo
html_logo = "_static/conforce_logo_small.png"
html_favicon = "_static/favicon.ico"

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


# -- Generate supported_element_types.rst -------------------------------------

with open(f"{DOC_SOURCE_DIR}/supported_element_types.rst", "w") as fh:
    fh.write("Supported element types\n"
             "=======================\n"
             "\n")

    fh.write("""
Directly supported element types
--------------------------------

This element types can be used by methods defined in :py:mod:`conforce_shared.cf_c`.

""")

    for el_type, el_info in cf_c.map_type_to_info.items():  # type: str, cf_c.ElementInfo
        fh.write(f"- **{el_type}**: "
                 f"{el_info.number_of_dimensions}D element with "
                 f"{el_info.number_of_nodes} nodes and "
                 f"{el_info.number_of_integration_points} integration points\n")

    fh.write("""
Indirectly supported element types 
----------------------------------

This (Abaqus) element types can be replaced by the directly supported element types on the right.
This may involve some simplifications and approximations.
For example the out-of-plane strain of plane stress elements is neglected.

.. seealso:: :py:mod:`conforce_shared.element_type_mapping`

""")

    for el_type, supported_el_type in element_type_mapping.map_abaqus_element_type_to_supported_element_type.items():
        fh.write(f"- **{el_type}** -> **{supported_el_type}**\n")

# -- Copy files for root directory ---------------------------------------------

for src, dest in [
    (f"{HOME_DIR}/README.md", f"{DOC_SOURCE_DIR}/README.md"),
    (f"{HOME_DIR}/conforce_logo.png", f"{DOC_SOURCE_DIR}/conforce_logo.png"),
    (f"{HOME_DIR}/plugin_gui.png", f"{DOC_SOURCE_DIR}/plugin_gui.png"),
    (f"{HOME_DIR}/CONTRIBUTING.md", f"{DOC_SOURCE_DIR}/CONTRIBUTING.md"),
    (f"{HOME_DIR}/LICENSE.txt", f"{DOC_SOURCE_DIR}/LICENSE.txt"),
    (f"{HOME_DIR}/theory/theory_images", f"{DOC_SOURCE_DIR}/generated/theory_images"),
    (f"{HOME_DIR}/examples/Example_1_Two_phase_bar/example_1_images", f"{DOC_SOURCE_DIR}/examples/example_1_images"),
    (f"{HOME_DIR}/examples/Example_2_CT_specimen_linear_elastic/example_2_images", f"{DOC_SOURCE_DIR}/examples/example_2_images"),
    (f"{HOME_DIR}/examples/Example_3_CT_specimen_elastic_plastic/example_3_images", f"{DOC_SOURCE_DIR}/examples/example_3_images"),
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
