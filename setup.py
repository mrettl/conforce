from setuptools import setup

import conforce_shared

setup(
    name=conforce_shared.project,
    version=conforce_shared.version,
    description=conforce_shared.description,
    url=conforce_shared.helpUrl,
    license_files=["LICENSE.txt"],
    author=conforce_shared.author,
    author_email='matthias.rettl@unileoben.ac.at',
    packages=['conforce', 'conforce_shared'],
    package_dir={'.': ''},
    package_data={'conforce_shared': ['*.dll', '*.so']},
    install_requires=[
        "numpy~=1.21.5",
        "sympy~=1.10.1",
    ],
    extras_require={
        "documentation": [
            "Sphinx>=5.0",
            "recommonmark>=0.7",
        ]
    },
)