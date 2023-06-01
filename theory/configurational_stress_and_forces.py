r"""
Configurational forces are the derivative of the energy with respect to a change in the geometry.

>>> import numpy as np
>>> import sympy as sy
>>> import conforce.expressions as expr
>>> import conforce.element_definitions as el_def





Configurational stresses `CS`
-----------------------------

mbf = energy momentum tensor going back to Eshelby

>>> energy_density

.. todo:: CS


Configurational forces `CF`
---------------------------

>>> el_type = el_def.CPE4R
>>> R_at_nodes = el_def.R_at_nodes_of_element[el_type]
>>> exponents = el_def.exponents_of_shape_functions_of_element[el_type]
>>> R_at_int_points = el_def.R_at_integration_points_of_element[el_type]
>>> int_weights = el_def.weights_of_integration_points_of_element[el_type]
>>> n, d = R_at_nodes.shape
>>> ips = len(int_weights)

>>> X_at_nodes = sy.IndexedBase("X_at_nodes", shape=(n, d))
>>> U_at_nodes = sy.IndexedBase("U_at_nodes", shape=(n, d))
>>> S_at_int_points = sy.IndexedBase("S_at_int_points", shape=(ips, d, d))
>>> e_at_int_points = sy.IndexedBase("e_at_int_points", shape=(ips,))

>>> def compute_CS_and_CF(is_dbf):
...     comp = expr.Computation(
...         R_at_nodes,
...         exponents,
...         R_at_int_points,
...         int_weights,
...         X_at_nodes,
...         U_at_nodes,
...         S_at_int_points,
...         e_at_int_points,
...         is_dbf=is_dbf
...     )
...     return comp.CS, comp.CF_at_nodes
>>> CS_dbf, CF_dbf_at_nodes = compute_CS_and_CF(is_dbf=True)
>>> CS_mbf, CF_mbf_at_nodes = compute_CS_and_CF(is_dbf=False)



.. todo:: CF

"""
import doctest

if __name__ == '__main__':
    pass # doctest.testmod()

import numpy as np
import sympy as sy
import conforce.expressions as expr
import conforce.element_definitions as el_def
from conforce.symbolic_util import create_replacement_rules

el_type = el_def.CPE4R
R_at_nodes = el_def.R_at_nodes_of_element[el_type]
exponents = el_def.exponents_of_shape_functions_of_element[el_type]
R_at_int_points = el_def.R_at_integration_points_of_element[el_type]
int_weights = el_def.weights_of_integration_points_of_element[el_type]
n, d = R_at_nodes.shape
ips = len(int_weights)

S = sy.MatrixSymbol("S", d, d)
e = sy.symbols("e")

X_at_nodes = R_at_nodes
strain = sy.MatrixSymbol("strain", d, d)
dU_dX = strain

S_at_int_points = sy.IndexedBase("S_at_int_points", shape=(ips, d, d))
e_at_int_points = sy.IndexedBase("e_at_int_points", shape=(ips,))

R = expr.eval_R(d)
H = expr.eval_H(R, R_at_nodes, exponents)
dH_dR = expr.eval_dH_dR(H, R)
dX_dR = expr.eval_dX_dR(X_at_nodes, dH_dR)
dH_dX = expr.eval_dH_dX(dH_dR, dX_dR)
F = expr.eval_F(d, dU_dX)
P = expr.eval_P(F, S)
CS = expr.eval_CS_mbf(d, e, F, P)
CF = expr.eval_CF_at_nodes(dH_dX, CS, dX_dR, int_weights, create_replacement_rules(R, *R_at_int_points))

CF.xreplace({
    strain: sy.Matrix([
        [1, 0],
        [0, 0]
    ]),
    e: 0.5,
    S: sy.Matrix([
        [1, 0],
        [0, 0]
    ])
})

print("ok")
