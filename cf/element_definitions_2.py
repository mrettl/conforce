import numpy as np

from cf.tests.one_element_abaqus_runner import simulate_one_element
from cf import element_definitions as el_def
from cf.expressions import eval_H, eval_R
from cf.symbolic_util import create_replacement_rules, apply_replacement_rules

R_at_nodes_of_element_type = dict()
exponents_of_shape_functions_of_element = dict()


#
R_at_nodes_of_element_type["C3D20"] = np.array([
    [-1., -1., -1.], [1.0, -1., -1.], [1.0, 1.0, -1.], [-1., 1.0, -1.],
    [-1., -1., 1.0], [1.0, -1., 1.0], [1.0, 1.0, 1.0], [-1., 1.0, 1.0],
    [0.0, -1., -1.], [1.0, 0.0, -1.], [0.0, 1.0, -1.], [-1., 0.0, -1.],
    [0.0, -1., 1.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0], [-1., 0.0, 1.0],
    [-1., -1., 0.0], [1.0, -1., 0.0], [1.0, 1.0, 0.0], [-1., 1.0, 0.0]
])
exponents_of_shape_functions_of_element['C3D20'] = np.array([
    [0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0],
    [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1],
    [2, 0, 0], [0, 2, 0], [2, 1, 0], [1, 2, 0],
    [2, 0, 1], [0, 2, 1], [2, 1, 1], [1, 2, 1],
    [0, 0, 2], [1, 0, 2], [0, 1, 2], [1, 1, 2]
])

R_at_nodes_of_element_type["C3D20R"] = R_at_nodes_of_element_type["C3D20"]
exponents_of_shape_functions_of_element['C3D20R'] = exponents_of_shape_functions_of_element['C3D20']

R_at_nodes_of_element_type["C3D8"] = R_at_nodes_of_element_type["C3D20"][:8, :]
exponents_of_shape_functions_of_element['C3D8'] = exponents_of_shape_functions_of_element['C3D20'][:8, :]

R_at_nodes_of_element_type["C3D8R"] = R_at_nodes_of_element_type["C3D8"]
exponents_of_shape_functions_of_element['C3D8R'] = exponents_of_shape_functions_of_element['C3D8']

R_at_nodes_of_element_type["CPE8"] = np.concatenate([
    R_at_nodes_of_element_type["C3D20"][0:4, :2],
    R_at_nodes_of_element_type["C3D20"][8:12, :2],
], axis=0)
exponents_of_shape_functions_of_element["CPE8"] = np.concatenate([
    exponents_of_shape_functions_of_element["C3D20"][0:4, :2],
    exponents_of_shape_functions_of_element["C3D20"][8:12, :2],
], axis=0)

R_at_nodes_of_element_type["CPE8R"] = R_at_nodes_of_element_type["CPE8"]
exponents_of_shape_functions_of_element["CPE8R"] = exponents_of_shape_functions_of_element["CPE8"]

R_at_nodes_of_element_type["CPE4"] = R_at_nodes_of_element_type["CPE8"][:4, :]
exponents_of_shape_functions_of_element["CPE4"] = exponents_of_shape_functions_of_element["CPE8"][:4, :]

R_at_nodes_of_element_type["CPE4R"] = R_at_nodes_of_element_type["CPE4"]
exponents_of_shape_functions_of_element["CPE4R"] = exponents_of_shape_functions_of_element["CPE4"]


#
R_at_nodes_of_element_type["C3D10"] = np.array([
    [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
    [0.5, 0.0, 0.0], [0.5, 0.5, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]
])
exponents_of_shape_functions_of_element['C3D10'] = np.array([
    [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1],
    [2, 0, 0], [1, 1, 0], [0, 2, 0], [0, 0, 2], [1, 0, 1], [0, 1, 1],
])

R_at_nodes_of_element_type["C3D4"] = R_at_nodes_of_element_type["C3D10"][:4, :]
exponents_of_shape_functions_of_element["C3D4"] = exponents_of_shape_functions_of_element["C3D10"][:4, :]


#
R_at_nodes_of_element_type["C3D15"] = np.array([
    [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0],
    [0.5, 0.0, 0.0], [0.5, 0.5, 0.0], [0.0, 0.5, 0.0],
    [0.5, 0.0, 1.0], [0.5, 0.5, 1.0], [0.0, 0.5, 1.0],
    [0.0, 0.0, 0.5], [1.0, 0.0, 0.5], [0.0, 1.0, 0.5],
])
exponents_of_shape_functions_of_element['C3D15'] = np.array([
    [0, 0, 0], [1, 0, 0], [0, 1, 0],
    [0, 0, 1], [1, 0, 1], [0, 1, 1],
    [2, 0, 0], [1, 1, 0], [0, 2, 0],
    [2, 0, 1], [1, 1, 1], [0, 2, 1],
    [0, 0, 2], [1, 0, 2], [0, 1, 2],
])

R_at_nodes_of_element_type["C3D6"] = R_at_nodes_of_element_type["C3D15"][:6, :]
exponents_of_shape_functions_of_element["C3D6"] = exponents_of_shape_functions_of_element["C3D15"][:6, :]

R_at_nodes_of_element_type["CPE6"] = np.concatenate([
    R_at_nodes_of_element_type["C3D15"][0:3, :2],
    R_at_nodes_of_element_type["C3D15"][6:9, :2],
], axis=0)
R_at_nodes_of_element_type["CPE6"] = np.concatenate([
    exponents_of_shape_functions_of_element["C3D15"][0:3, :2],
    exponents_of_shape_functions_of_element["C3D15"][6:9, :2],
], axis=0)

R_at_nodes_of_element_type["CPE3"] = R_at_nodes_of_element_type["CPE6"][:3, :]
exponents_of_shape_functions_of_element["CPE3"] = exponents_of_shape_functions_of_element["CPE6"][:3, :]



def sim():
    a = {
        element_type: simulate_one_element(
            R_at_nodes=R_at_nodes,
            U_at_nodes=1e-7 * R_at_nodes,
            element_type=element_type,
            load_name="test",
            folder="simulations/tests/test_element_definitions"
        )
        for element_type, R_at_nodes in R_at_nodes_of_element_type.items()
    }

    for element_type in el_def.R_at_nodes_of_element.keys():
        np.testing.assert_array_almost_equal(
            np.array(a[element_type]["integration_points"]["IVOL"]).flatten(),
            el_def.weights_of_integration_points_of_element[element_type]
        )

        np.testing.assert_array_almost_equal(
            np.array(a[element_type]["integration_points"]["COORD"]),
            el_def.R_at_integration_points_of_element[element_type]
        )

        np.testing.assert_array_almost_equal(
            np.array(a[element_type]["nodes"]["COORD"]),
            el_def.R_at_nodes_of_element[element_type]
        )


if __name__ == '__main__':
    for element_type in exponents_of_shape_functions_of_element.keys():
        R_at_nodes = R_at_nodes_of_element_type[element_type]
        exponents = exponents_of_shape_functions_of_element[element_type]
        n, d = R_at_nodes.shape
        R = eval_R(d)
        for i in range(1, n):
            R_at_nodes_ = R_at_nodes[:i]
            repl = create_replacement_rules(
                R,
                *R_at_nodes_
            )

            h = eval_H(
                R=R,
                R_at_nodes_=R_at_nodes_,
                exponents_=exponents[:i]
            )

            a = np.array(apply_replacement_rules(
                h,
                *repl
            ), dtype=float)[:, :, 0]
            np.testing.assert_array_almost_equal(
                a,
                np.eye(i)
            )

            print(h)

    print("Ok")
