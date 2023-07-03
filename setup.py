"""
This is the setup script for the Python 3 package.
"""
from setuptools import setup

import conforce


def read_file(path: str):
    with open(path, "r") as fh:
        return fh.read()


def read_requirements(path: str):
    return [
        requirement.strip()
        for requirement
        in read_file(path).splitlines()
        if not requirement.startswith("#") or requirement.strip() == ""
    ]


setup(
    name=conforce.project,
    version=conforce.version,
    author=conforce.author,
    author_email='matthias.rettl@unileoben.ac.at',
    description=conforce.description,
    long_description=read_file("README.md"),
    long_description_content_type="text/markdown",
    url='https://conforce.readthedocs.io/',
    project_urls={
        'Documentation': 'https://conforce.readthedocs.io/',
        'Source': 'https://github.com/mrettl/conforce',
        'Tracker': 'https://github.com/mrettl/conforce/issues',
    },
    license_files=["LICENSE.txt"],
    packages=['conforce_3', 'conforce'],
    package_dir={'.': ''},
    package_data={
        'conforce_3': ['*.json', '*.c'],
        'conforce': ['*.dll', '*.so']
    },
    python_requires=">=3.7",
    install_requires=read_requirements("requirements.txt"),
    extras_require={
        "doc": read_requirements("doc/requirements.txt"),
        "examples": read_requirements("examples/requirements.txt")
    },
    classifiers=[
        "Intended Audience :: Science/Research",
        "Operating System :: POSIX :: Linux",
        "Operating System :: Microsoft :: Windows",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Mathematics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3 :: Only",
        (
            "Development Status :: 3 - Alpha" if "-alpha" in conforce.version.lower()
            else "Development Status :: 4 - Beta" if "-beta" in conforce.version.lower()
            else "Development Status :: 5 - Production/Stable"
        )
    ]
)
