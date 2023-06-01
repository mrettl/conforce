Add new examples
================

New examples are not added automatically to the Sphinx documentation.
To update the Sphinx documentation, follow the steps:

1. Append in `conf.py` the path to the new example folder to `sys.path`, 
   such that python scripts can be imported by Sphinx.
2. Create new \*.rst files in the folder `doc/source/examples` and write the documentation.
3. Put a reference into the examples section of `index.rst` to the new \*.rst files.
4. Let `conf.py` copy images, ... to `doc/source/examples`
