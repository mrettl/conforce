import numpy as np

from cf.tests.one_element_abaqus_runner import simulate_one_element
from cf import element_definitions as el_def
from cf.expressions import eval_H, eval_R
from cf.symbolic_util import create_replacement_rules, apply_replacement_rules

R_at_nodes_of_element_type = dict()
exponents_of_shape_functions_of_element = dict()
adjacent_nodes_of_element = dict()


#
def compute_adjacent_node_of_quadratic_elements(R_at_nodes_of_quad_element):
    n = len(R_at_nodes_of_quad_element)
    adjacent_nodes = np.zeros([n, n], dtype=int)
    node_ids = np.arange(0, n)

    for i in range(n):
        for j in range(i+1, n):
            w_basis = np.eye(3)

            pi = R_at_nodes_of_quad_element[i, :]
            pj = R_at_nodes_of_quad_element[j, :]

            # Gram-Schmidt process to create a orthonormal vector basis
            v_1 = pj - pi
            v_1 = v_1 / np.linalg.norm(v_1)

            if np.linalg.norm(np.cross(w_basis[0], v_1)) < 1e-6:
                w_basis = w_basis[1:]

            w_2 = w_basis[0]
            v_2 = (
                    w_2
                    - np.dot(v_1, w_2) / np.dot(v_1, v_1) * v_1
            )

            v_3 = np.cross(v_1, v_2)
            v_3 = v_3 / np.linalg.norm(v_3)

            # filter nodes not lying on the line from the nodes i to the node j
            offsets = np.column_stack([
                R_at_nodes_of_quad_element @ v_2,
                R_at_nodes_of_quad_element @ v_3])
            offsets -= offsets[i]
            on_line_mask = np.linalg.norm(offsets, axis=1) < 1e-6
            node_ids_on_line = node_ids[on_line_mask]

            # compute positions on the line with the i-th node as start point
            positions = R_at_nodes_of_quad_element @ v_1
            positions = np.abs(positions - positions[i])
            positions = positions[on_line_mask]

            if len(positions) < 3:
                continue

            remove_node_i_mask = node_ids_on_line != i
            positions_wo_i = positions[remove_node_i_mask]
            node_ids_on_line_wo_i = node_ids_on_line[remove_node_i_mask]

            adjacent_mask = np.abs(positions_wo_i - np.min(positions_wo_i)) < 1e-6
            node_ids_adjacent = node_ids_on_line_wo_i[adjacent_mask]

            adjacent_nodes[i, node_ids_adjacent] = 1
            adjacent_nodes[node_ids_adjacent, i] = 1

    return adjacent_nodes


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
adjacent_nodes_of_element["C3D20"] = np.array([
    [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1],
    [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
])


R_at_nodes_of_element_type["C3D20R"] = R_at_nodes_of_element_type["C3D20"]
exponents_of_shape_functions_of_element['C3D20R'] = exponents_of_shape_functions_of_element['C3D20']
adjacent_nodes_of_element['C3D20R'] = adjacent_nodes_of_element['C3D20']

R_at_nodes_of_element_type["C3D8"] = R_at_nodes_of_element_type["C3D20"][:8, :]
exponents_of_shape_functions_of_element['C3D8'] = exponents_of_shape_functions_of_element['C3D20'][:8, :]
adjacent_nodes_of_element['C3D8'] = adjacent_nodes_of_element['C3D20'][:8, :8]

R_at_nodes_of_element_type["C3D8R"] = R_at_nodes_of_element_type["C3D8"]
exponents_of_shape_functions_of_element['C3D8R'] = exponents_of_shape_functions_of_element['C3D8']
adjacent_nodes_of_element['C3D8R'] = adjacent_nodes_of_element['C3D8']

R_at_nodes_of_element_type["CPE8"] = R_at_nodes_of_element_type["C3D20"][[0, 1, 2, 3, 8, 9, 10, 11], :2]
exponents_of_shape_functions_of_element["CPE8"] = exponents_of_shape_functions_of_element["C3D20"][[0, 1, 2, 3, 8, 9, 10, 11], :2]
adjacent_nodes_of_element['CPE8'] = adjacent_nodes_of_element['CPE8'][[0, 1, 2, 3, 8, 9, 10, 11], [0, 1, 2, 3, 8, 9, 10, 11]]

R_at_nodes_of_element_type["CPE8R"] = R_at_nodes_of_element_type["CPE8"]
exponents_of_shape_functions_of_element["CPE8R"] = exponents_of_shape_functions_of_element["CPE8"]
adjacent_nodes_of_element["CPE8R"] = adjacent_nodes_of_element["CPE8"]

R_at_nodes_of_element_type["CPE4"] = R_at_nodes_of_element_type["CPE8"][:4, :]
exponents_of_shape_functions_of_element["CPE4"] = exponents_of_shape_functions_of_element["CPE8"][:4, :]
adjacent_nodes_of_element["CPE4"] = adjacent_nodes_of_element["CPE8"][:4, :4]

R_at_nodes_of_element_type["CPE4R"] = R_at_nodes_of_element_type["CPE4"]
exponents_of_shape_functions_of_element["CPE4R"] = exponents_of_shape_functions_of_element["CPE4"]
adjacent_nodes_of_element["CPE4R"] = adjacent_nodes_of_element["CPE4"]

#
R_at_nodes_of_element_type["C3D10"] = np.array([
    [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
    [0.5, 0.0, 0.0], [0.5, 0.5, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]
])
exponents_of_shape_functions_of_element['C3D10'] = np.array([
    [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1],
    [2, 0, 0], [1, 1, 0], [0, 2, 0], [0, 0, 2], [1, 0, 1], [0, 1, 1],
])
adjacent_nodes_of_element["C3D10"] = np.array([
    [0, 0, 0, 0, 1, 0, 1, 1, 0, 0],
    [0, 0, 0, 0, 1, 1, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 1, 1, 0, 0, 1],
    [0, 0, 0, 0, 0, 0, 0, 1, 1, 1],
    [1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 1, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 1, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 0, 1, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 1, 0, 0, 0, 0, 0, 0]
])

R_at_nodes_of_element_type["C3D4"] = R_at_nodes_of_element_type["C3D10"][:4, :]
exponents_of_shape_functions_of_element["C3D4"] = exponents_of_shape_functions_of_element["C3D10"][:4, :]
adjacent_nodes_of_element["C3D4"] = adjacent_nodes_of_element["C3D10"][:4, :4]


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
adjacent_nodes_of_element["C3D15"] = np.array([
    [0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1],
    [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
])

R_at_nodes_of_element_type["C3D6"] = R_at_nodes_of_element_type["C3D15"][:6, :]
exponents_of_shape_functions_of_element["C3D6"] = exponents_of_shape_functions_of_element["C3D15"][:6, :]
adjacent_nodes_of_element["C3D6"] = adjacent_nodes_of_element["C3D15"][:6, :6]

R_at_nodes_of_element_type["CPE6"] = R_at_nodes_of_element_type["C3D15"][[0, 1, 2, 6, 7, 8], :2]
exponents_of_shape_functions_of_element["CPE6"] = exponents_of_shape_functions_of_element["C3D15"][[0, 1, 2, 6, 7, 8], :2]
adjacent_nodes_of_element["CPE6"] = adjacent_nodes_of_element["C3D15"][[0, 1, 2, 6, 7, 8], [0, 1, 2, 6, 7, 8]]

R_at_nodes_of_element_type["CPE3"] = R_at_nodes_of_element_type["CPE6"][:3, :]
exponents_of_shape_functions_of_element["CPE3"] = exponents_of_shape_functions_of_element["CPE6"][:3, :]
adjacent_nodes_of_element["CPE3"] = adjacent_nodes_of_element["CPE6"][:3, :3]



def test():
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


def sim():
    a = {
        element_type: simulate_one_element(
            R_at_nodes=R_at_nodes,
            U_at_nodes=1e-7 * R_at_nodes,
            element_type=element_type,
            load_name="test",
            folder="res/tests/test_element_definitions"
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
    print("Ok")
