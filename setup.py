from setuptools import setup, find_packages
from os import path
from io import open

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='ConF3D',  # Required
    version='0.1.0',  # Required
    description='Configurational forces for FEM postprocessing',  # Optional
    long_description=long_description,  # Optional
    long_description_content_type='text/markdown',  # Optional
    url='',  # Optional
    author='Markus Tauscher',  # Optional
    author_email='markus.tauscher@unileoben.ac.at',  # Optional
    classifiers=[  # Optional
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 1 - In Development',

        'Intended Audience :: Science/Research',
        'Topic :: Documentation :: Tutorial',

        'License :: OSI Approved :: MIT License',

        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],

    keywords='configurational forces, finite element method',  # Optional

    package_dir={'': 'src'},  # Optional
    packages=find_packages(where='src'),  # Required
    # packages=['src/ConF3D'],
    python_requires='>=2.7, >=3.0',
    install_requires=required,  # Optional
    # package_data={  # Optional
    #     '': [''],
    # },
)
