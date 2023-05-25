from setuptools import setup

from cf_shared import version

setup(
    name='ConfForce',
    version=version,
    description='This packages provides methods to compute configurational forces.',
    url='https://github.com/mrettl/ConfigurationalForcesPlugin',
    license_files=["LICENSE.txt"],
    author='Matthias Rettl, Markus Tauscher, Sigfried Frankl, Martin Pletz',
    author_email='matthias.rettl@unileoben.ac.at',
    packages=['cf', 'cf_shared'],
    package_dir={'.': ''},
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
