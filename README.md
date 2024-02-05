Readme
======

![conforce](conforce_logo.png)

[![Documentation Status](https://readthedocs.org/projects/conforce/badge/?version=latest)](https://conforce.readthedocs.io/en/latest/?badge=latest)

ConForce is a Python package and Abaqus plug-in for the computation of configurational (material) forces.
The computation is automatically generated and compiled from symbolic calculations into C code.
A C binding allows access to the C functions from Python.
This makes ConForce fast enough to compute an entire finite element model 
with several thousand elements within a few seconds.

ConForce features:

- the computation of the following quantities for common 2D and 3D elements:
  - static configurational forces
  - static configurational stresses
  - first Piola-Kirchhoff stresses
  - deformation gradients
- two formulations for the configurational stresses and forces:
  - motion based formulation (Eshebly's formulation)
  - displacement based formulation
- material models:
  - non-linear elasticity in the large strain framework
  - small-strain plasticity
- material orientations: The stresses and displacements are automatically rotated to the global coordinate system.


Supported Operating Systems:
- Windows (64-bit)
- CentOS (64-bit)
- other Linux distributions might work


ConForce is available as:
- Python 3 package
- Abaqus plug-in

## Python 3 Package

The Python package requires Python 3 and contains:

- methods to compute configurational forces, stresses, ...
- methods for the symbolic computation and code generation

### Installation 

#### PyPi

We recommend to install the latest stable version from PyPi.

````shell
pip install conforce
````

#### Remote Repository

The (unstable) development version of ConForce can be installed 
from our remote repository.

1. Download and unzip the [ConForce repository](https://github.com/mrettl/conforce).
2. Open a shell in the unzipped folder where the *setup.py* file is located.
3. Install ConForce by typing the following commands into the shell
   ````shell
   pip install "."
   ````
   In order to run examples or build the documentation,
   additional packages are required.
   The following options are available.
   - Install requirements for the examples:
     ````shell
     pip install ".[examples]"
     ````
   - Install requirements to build the documentation
     ````shell
     pip install ".[doc]"
     ````
   - Install all requirements
     ````shell
     pip install ".[examples, doc]"
     ````
   
4. To uninstall the package, type:
    ````shell
    pip uninstall conforce
    ````

If git is available, the package can be installed directly from the remote repository by
````shell
pip install "conforce[examples, doc] @ git+https://github.com/mrettl/conforce"
````
The extras `examples` and `doc` are optional.

### Usage

Configurational forces are computed using the `cf_c` module.
The energy density `e`, the undeformed positions `X`,
the displacement `U`, and the symmetric Cauchy stress tensor `S`
are passed to the `compute_CF` function along with the element type and the computation method.

````python
from conforce import cf_c

cf_c.compute_CF(
    e_at_int_points=[[10.]],
    X_at_nodes=[[
        [0., 0.],
        [1., 0.],
        [1., 1.],
        [0., 1.],
    ]],
    U_at_nodes=[[
        [0.0, 0.0],
        [0.1, 0.0],
        [0.1, 0.0],
        [0.0, 0.0],
    ]],
    S_at_int_points=[[[
        [100., 0.0],
        [0.0, 0.0]
    ]]],
    element_type="CPE4R",
    method="dbf"
)
````

This yields the configurational forces at the element nodes for a single CPE4R element
and the displacement based formulation ("dbf"):

````
array([[[ 0., -5.],
        [ 0., -5.],
        [ 0.,  5.],
        [-0.,  5.]]])
````


## Abaqus plug-in

The Abaqus plug-in contains:

- methods to compute configurational forces, stresses, ...
- Abaqus specific code to read and write to *.odb files
- a GUI for the Abaqus plug-in

Supported Abaqus versions:
- \>= Abaqus 2017

### Installation

Open 
[https://github.com/mrettl/conforce/releases](https://github.com/mrettl/conforce/releases)
and download the latest zip file called `conforce_plugin_{version}.zip`.
The zip file has the following structure:

- `conforce_plugin_{version}.zip`
  - `conforce`
    - `conforce_abq_plugin.py`
    - ...

There are several possibilities in which folder the plug-in files can be stored.
Decide for one of the following valid `plugin-folder`s:

- `plugin_central_dir`: This folder is defined in the environment file `abaqus_v6.env`.
  For many installations this folder is located in `C:\\SIMULIA\\CAE\\plugins\\{year}`.
  Plug-ins stored in this folder are accessible for all users.
- `current_dir\abaqus_plugins`: This is the folder in which Abaqus is started.
  Plug-ins are only accessible if Abaqus is stared again in the same folder.

Next, unzip the downloaded files and put them in the `plugin-folder` you selected:
The folder tree should look like:

- `plugin-folder`
  - `conforce`
    - `conforce_abq_plugin.py`
    - ...

Start Abaqus and navigate in the toolbar to `Plug-ins -> Conf. Force`.
The GUI of the plug-in should open.
If no toolbar entry called `Conf. Force` exists, check if the plug-in is stored in the right folder 
for the Abaqus version you are using.

### Usage

Open an *.odb file in Abaqus.
Navigate in the toolbar to `Plug-ins -> Conf. Force`.
The plug-in gui opens:

![plugin gui](plugin_gui.png)

Click Apply to compute the requested field outputs.


## Report Bugs, Ask questions

Use the [github issue tracker](https://github.com/mrettl/conforce/issues) to report
bugs and troubles or to ask questions.
If you like our tool, please leave a star on our Github repository.
