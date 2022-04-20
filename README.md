# A sample Python project

A sample project that exists as an aid to develop installable python packages.
This sample is as slim as possible and does not meet the standards for public distribution via
GitHub or similar. If you want to include your package in the PyPI follow the [packaging guide][packaging guide].
A blank sample project for public distribution via the PyPI can be downloaded from GitHub.
[public sample project][https://github.com/pypa/sampleproject]

This sample project aims to be installed in your local python distribution.
To install this package open a command window, change into the package folder and execute

`pip install .`

To uninstall open a command window and execute

`pip uninstall sample_project`

Most of the configuration for a Python project is done in the `setup.py` file,
an example of which is included in this project. You should edit this file
accordingly to adapt this sample project to your needs.

----

This is the README file for the project.

The file should use UTF-8 encoding and can be written using
[reStructuredText][rst] or [markdown][md use] with the appropriate [key set][md use].
It will be shown on the front-page of the GitLab project.

Typical contents for this file would include an overview of the project, basic
usage examples, etc. Generally, including the project changelog in here is not a
good idea, although a simple “What's New” section for the most recent version
may be appropriate.

[packaging guide]: https://packaging.python.org
[rst]: http://docutils.sourceforge.net/rst.html
[md]: https://tools.ietf.org/html/rfc7764#section-3.5 "CommonMark variant"
[md use]: https://packaging.python.org/specifications/core-metadata/#description-content-type-optional
