[![Documentation Status](https://readthedocs.org/projects/cf-configurational-forces-plug-in/badge/?version=main)](https://cf-configurational-forces-plug-in.readthedocs.io/en/main/?badge=main)

Readme
======

# ![confore](conforce_logo_icon.png)

This package provides:

- methods to compute for various element types the following quantities:
  - configurational forces (deformation and motion based, quasi-static)
  - configurational stresses (deformation and motion based, quasi-static)
  - First Piola-Kirchhoff stresses
  - deformation gradients
- an Abaqus Plug-in to compute those quantities in the post-processing
- a code generation that generates C-code from symbolic functions
- a C-code binding to use the generated and compiled C-code

Supported OS:
- Windows (64-bit)
- CentOS (64-bit)
- other Linux distributions might work

Supported Abaqus versions:
- \>= Abaqus 2017

## Abaqus Plug-in

The Abaqus Plug-in contains:
- methods to compute configurational force, ...
- Abaqus specific code to read and write to *.odb files
- a GUI for the Abaqus Plug-in

### Installation

Download the zip file called `conforce_plugin_{version}.zip`.
The zip file has the following structure:

- `conforce_plugin_{version}.zip`
  - `conforce`
    - `conforce_abq`
      - ...
    - `conforce_shared`
      - ...
    - `conforce_abq_plugin.py`
    - ...

To install the Abaqus Plug-in, unzip the file in one of the following valid `plugin-folder`:

- `plugin_central_dir`: This is an Abaqus environment parameter.
  A Plug-in in this folder is accessible regardless of the home
  or current directory. For many installations this folder is located in
  `C:\\SIMULIA\\CAE\\plugins\\{year}`.
- `home_dir\abaqus_plugins`: This is the folder in which Abaqus is started.
- `current_dir\abaqus_plugins`: The current work directory can be changed inside Abaqus.
  Navigate to `File -> Set Work Directory ...` to define this folder.

The folder tree should look like:

- `plugin-folder`
  - `conforce`
    - `conforce_abq`
      - ...
    - `conforce_shared`
      - ...
    - `conforce_abq_plugin.py`
    - ...

Start Abaqus and navigate in the toolbar to `Plug-ins -> Conf. Force`.
The GUI of the Plug-in should open.
If no toolbar entry called `Conf. Force` exists, check if the Plug-in is stored in the right folder 
for the Abaqus version you are using.


### Usage

Open an *.odb file in Abaqus.
Navigate in the toolbar to `Plug-ins -> Conf. Force`.
The Plug-in gui opens:

![plugin gui](plugin_gui.png)

Click Apply to compute the requested field outputs.


## Package

The installed package contains:
- methods to compute configurational force, ...
- methods for symbolic computations, code generation, ...

### Installation

1. Clone the repository into a folder
2. Open a shell in this folder where the setup.py file is located. Activate your environment, if you are using conda.
3. Type:
    ````shell
    pip install .
    ````
4. To uninstall the package, type:
    ````shell
    pip uninstall conforce
    ````

Alternatively, if git is installed, the package can be installed directly without cloning the repository by
  ````shell
  pip install git+https://github.com/mrettl/conforce
  ````


### Usage

Configurational forces are computed using the `cf_c` module.
The energy density `e`, the undeformed positions `X`,
the displacement `U` and the symmetric Cauchy stress tensor `S`
are passed to the `compute_CF` function alongside the element type and the computation method.

````python
from conforce_shared import cf_c

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

This yields the configurational forces:

````
array([[[ 0., -5.],
        [ 0., -5.],
        [ 0.,  5.],
        [-0.,  5.]]])
````
