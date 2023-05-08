# Contributing

The following is a set of guidelines and conventions.

## Naming Conventions

The naming conventions follows the [Style Guid for Python Code](https://peps.python.org/pep-0008/).
However, there are a few exceptions.

### Modules

 - The name of modules requiring Abaqus begins with `abaqus_`.
 - The name of unit test files starts with `test_` or `test_abaqus_` (for Abaqus modules).
 - Modules that can be used in Python 3.* and Abaqus, have a note for the inter-compatibility in the module doc-string.

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


## Documentation Conventions

Sphinx is used to create the documentation.
The docstrings are written in the reStructuredText format
and should not raise warnings in Sphinx (Delete the folder `doc/source/generated` before running sphinx).



