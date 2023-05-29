from setuptools import setup

import conforce_shared


def read_file(path: str):
    with open(path, "r") as fh:
        return fh.read()


setup(
    name=conforce_shared.project,
    version=conforce_shared.version,
    author=conforce_shared.author,
    author_email='matthias.rettl@unileoben.ac.at',
    description=conforce_shared.description,
    long_description=read_file("README.md"),
    long_description_content_type="text/markdown",
    url='https://cf-configurational-forces-plug-in.readthedocs.io/en/main/',
    project_urls={
        'Documentation': 'https://cf-configurational-forces-plug-in.readthedocs.io/en/main/',
        'Source': 'https://github.com/mrettl/conforce',
        'Tracker': 'https://github.com/mrettl/conforce/issues',
    },
    license=read_file("LICENSE.txt"),
    license_files=["LICENSE.txt"],
    packages=['conforce', 'conforce_shared'],
    package_dir={'.': ''},
    package_data={'conforce_shared': ['*.dll', '*.so']},
    install_requires=[
        requirement.strip()
        for requirement
        in read_file("requirements.txt").splitlines()
        if not requirement.startswith("#") or requirement.strip() == ""
    ],
    extras_require={
        "doc": [
            requirement.strip()
            for requirement
            in read_file("doc/requirements.txt").splitlines()
            if not requirement.startswith("#") or requirement.strip() == ""
        ]
    },
    classifiers=[
        "Intended Audience :: Science/Research",
        "Operating System :: POSIX :: Linux",
        "Operating System :: Microsoft :: Windows",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Programming Language :: Python :: 3 :: Only",
        (
            "Development Status :: 3 - Alpha" if "-alpha" in conforce_shared.version.lower()
            else "Development Status :: 4 - Beta" if "-beta" in conforce_shared.version.lower()
            else "Development Status :: 5 - Production/Stable"
        )
    ]
)
