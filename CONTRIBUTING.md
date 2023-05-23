# Contributing

This is a guideline for developers.
If you want to use this package, we refer to the README file. 

## Getting started

Clone or download this repository to your disk.
Make sure all dependencies declared in `requirements.txt` are installed.

Only for Windows users:
Download and install `gcc`.
Binaries for `gcc` are included in [w64devkit](https://github.com/skeeto/w64devkit/releases).
Create a system or user environment variable called `PATH` 
and append the path `w64devkit\bin` where `gcc.exe` is located.

Try to run all tests located in the modules `cf` and `cf_shared`.
Additionally, an Abaqus installation is required to execute code located in `cf_abq`.

## Test Code

Test code is written:
 - as `Doctest` in the docstrings of modules, classes or functions
 - as `unittest.TestCase`. Modules containing Unit Tests start with `test_` and are placed 
   in a subfolder of the module called `tests`.


## Documentation

Sphinx is used to create the documentation.
Delete the folder `doc/source/generated` before running sphinx.
To build the documentation open a shell and type

```
activate
cd doc
make html
```

If sphinx raises a warning, check if the warning can be easily fixed.

## Conventions and Guidelines

The following is a set of guidelines and conventions.

### Packages
    
 - The package
   - `cf` contains modules that require Python 3.* and do not run in Abaqus Python
   - `cf_abq` contains modules that require Abaqus Python and do not run in Python 3.*
   - `cf_shared` contains modules that run in both, Abaqus Python and Python 3.*

### Naming Conventions

The naming conventions follow the [Style Guid for Python Code](https://peps.python.org/pep-0008/).
However, there are a few exceptions.


#### Variables

 - A mathematical naming convention is used for variables.
   - Matrices and vectors are upper case.
   - Scalars are lower case.
   - Matrix components are scalars and thus are lower case. The index is written without an underscore. 
     E.g. `X = [x0, x1, x2]`
   - If (sympy) symbols and concrete values are used simultaneously, the variable for the concrete value ends with an underscore. E.g.:
     `J` (symbolic Jacobian matrix) and `J_` (numpy array of the Jacobian matrix)
   - Variables referring to specific points are named according to `{variable_name}_at_{point_name}`.
   - The derivative of `H` with respect to `R` is written as `dH_dR`.

#### Functions

 - Function names are lower case except they refer to a mathematical symbol that is written upper case 
    (e.g. `eval_J` refers to `J` the Jacobian matrix)

## Versioning
Use [Semantic Versioning](https://semver.org/spec/v2.0.0.html).