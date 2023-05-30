Add new examples
================

New examples are not added automatically to the Sphinx documentation.
To update the Sphinx documentation, follow the steps:

1. Append in `conf.py` the path to the new example folder to `sys.path`, 
   such that python scripts can be imported by Sphinx.
2. Create new \*.rst files in the folder `doc/source/examples` and write the documentation.
3. Add the new \*.rst files to `doc/source/examples.rst`
3. Add the new \*.rst files to `doc/source/examples.rst`
4. Let `conf.py` copy images, ... to `doc/source/examples`
