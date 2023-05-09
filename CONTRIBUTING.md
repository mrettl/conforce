# Contributing

The following is a set of guidelines and conventions.

## Naming Conventions

The naming conventions follows the [Style Guid for Python Code](https://peps.python.org/pep-0008/).
However, there are a few exceptions.

### Packages
    
 - The package
   - `cf` contains modules that require Python 3.* and do not run in Abaqus Python
   - `cf_abq` contains modules that require Abaqus Python and do not run in Python 3.*
   - `cf_shared` contains modules that run in both, Abaqus Python and Python 3.*

### Variables

 - A mathematical naming convention is used for variables.
   - Matrices are upper case.
   - Scalars are lower case.
   - Matrix components are scalars and thus are lower case. The index is written without an underscore. 
     E.g. `X = [x0, x1, x2]`
   - If (sympy) symbols and concrete values are used simultaneously, the variable for the concrete value ends with an underscore. E.g.:
     `J` (symbolic Jacobian matrix) and `J_` (numpy array of the Jacobian matrix)
   - Variables referring to specific points are named according to `{variable_name}_at_{point_name}`.
   - The derivative of `H` with respect to `R` is written as `dH_dR`.

### Functions

 - Function names are lower case except they refer to a mathematical symbol that is written upper case 
    (e.g. `eval_J` refers to `J` the Jacobian matrix)

### Classes

 - Class names are nouns in mixed case. The first character of the noun is upper case.


## Test Code

Test code is written:
 - as `Doctest` in the docstrings of modules, classes or functions
 - as `unittest.TestCase`. Modules containing Unit Tests start with `test_` and are placed 
   in a subfolder called `tests`.


## Documentation Conventions

Sphinx is used to create the documentation.
The docstrings are written in the reStructuredText format
and should not raise warnings in Sphinx (Delete the folder `doc/source/generated` before running sphinx).



