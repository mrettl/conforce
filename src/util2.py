from typing import Union, List, Tuple

import numpy as np
import sympy as sy

ARRAY_LIKE = Union[List, Tuple, np.ndarray]


# reference space coordinates
r, s, t = sy.symbols("r s t", real=True)

# real space coordinates
x, y, z = sy.symbols("x y z", real=True)


class Element(object):
    def __init__(
            self,
            coordinates: sy.MatrixBase,
            reference_coordinates: sy.MatrixBase,
            reference_shapes: sy.MatrixBase,
            nodes: ARRAY_LIKE,
    ):
        self.nodes = nodes
        self.coordinates = coordinates
        self.reference_coordinates = reference_coordinates

        # jac[i, j] = d(coordinates[i])/d(ref_coordinates[j])
        self.jacobian = self.coordinates.jacobian(self.reference_coordinates)
        self.inv_jacobian = self.jacobian.inv()
        self.det = sy.det(self.jacobian)

        self.reference_shapes = reference_shapes

        self.d_shapes_d_ref_coordinates = self.reference_shapes.jacobian(self.reference_coordinates)
        self.d_shapes_d_coordinates = self.d_shapes_d_ref_coordinates * self.inv_jacobian

    def count_dimensions(self):
        return len(self.coordinates)

    def count_nodes(self):
        return len(self.nodes)

    @staticmethod
    def _assure_matrix(values: ARRAY_LIKE) -> sy.MatrixBase:
        if not isinstance(values, sy.MatrixBase):
            values = np.asarray(values)
            if isinstance(values, np.ndarray) and values.ndim == 1:
                values = values.reshape((-1, 1))

            values = sy.Matrix(values)

        return values

    def __call__(self, values_at_nodes: ARRAY_LIKE = None):
        if values_at_nodes is None:
            values_at_nodes = self.nodes

        values_at_nodes_mat = self._assure_matrix(values_at_nodes)
        return values_at_nodes_mat.T * self.reference_shapes

    def diff(self, axis: int, values_at_nodes: ARRAY_LIKE = None):
        if values_at_nodes is None:
            values_at_nodes = self.nodes

        values_at_nodes_mat = self._assure_matrix(values_at_nodes)
        return values_at_nodes_mat.T * self.d_shapes_d_coordinates[:, axis]

    def create_real_space_element(self, nodes: ARRAY_LIKE):
        nodes = nodes[:, :self.count_dimensions()]
        return Element(
            self(nodes),
            self.reference_coordinates,
            self.reference_shapes,
            nodes
        )

    @classmethod
    def create_ref_space_element(
            cls,
            nodes: ARRAY_LIKE,
            shape_powers: ARRAY_LIKE
    ):
        nodes = np.array(nodes, dtype=float)
        shape_powers = np.array(shape_powers, dtype=int)
        coordinates = sy.Matrix((r, s, t))

        is_2d = shape_powers.shape[1] == 2 or np.all(0 == shape_powers[:, 2])
        if is_2d:
            nodes = nodes[:, :2]
            shape_powers = shape_powers[:, :2]
            coordinates = sy.Matrix(coordinates[:2])

        num_points = nodes.shape[0]
        num_powers = shape_powers.shape[0]
        assert num_points == num_powers

        # symbolic powers of shape functions
        symbolic_shape_powers = sy.Matrix([
            sy.Mul(*[
                sy.Pow(rst, power)
                for rst, power in zip(coordinates, powers)
            ]).nsimplify()
            for powers in shape_powers
        ])

        # a[i, j] = i-th shape power evaluated at j-th point
        a = sy.Matrix([
            symbolic_shape_powers.T.subs({
                rst: xyz
                for rst, xyz
                in zip(coordinates, point)
            })
            for point in nodes
        ]).T

        # solve coef_matrix * a = I,
        # such that every shape function is one at exactly one point and zero at the other points
        coef_matrix = a.inv()

        # combine shapes powers to shape functions
        shape_functions = coef_matrix * symbolic_shape_powers

        return cls(coordinates, coordinates, shape_functions, nodes)


class IntegrationScheme(object):
    symbolic_int_weight = sy.symbols("w", real=True)

    def __init__(self, reference_coordinates: sy.MatrixBase, integration_points: ARRAY_LIKE, integration_weight: ARRAY_LIKE):
        self._reference_coordinates = reference_coordinates
        n_dim = len(self._reference_coordinates)
        self._integration_points = np.array(integration_points, dtype=float)[:, :n_dim]
        self._integration_weights = np.array(integration_weight, dtype=float)

    def integrate(self, det, expression=1.):
        function = sy.lambdify(
            list(self._reference_coordinates) + [self.symbolic_int_weight],
            expression * det * self.symbolic_int_weight
        )
        result = np.sum(
            function(*self._integration_points.T, self._integration_weights),
            axis=-1
        )
        return result
