from functools import lru_cache
from typing import Union, List, Optional, Tuple

import numpy as np
import sympy as sy

ARRAY_LIKE = Union[List, Tuple, np.ndarray, sy.MatrixBase]


# reference space coordinates
r, s, t = sy.symbols("r s t", real=True)

# real space coordinates
x, y, z = sy.symbols("x y z", real=True)

# symbolic stress tensor
S = sy.MatrixSymbol("S", 3, 3)


class IntegrationScheme(object):
    def __init__(
            self,
            integration_points: ARRAY_LIKE,
            integration_weight: ARRAY_LIKE
    ):
        self._integration_points = np.array(integration_points, dtype=float)
        self._integration_weights = np.array(integration_weight, dtype=float)

    def integrate(
            self,
            symbolic_coordinates,
            det: Union[sy.Expr, float],
            expression: Union[sy.Expr, sy.MatrixBase, float] = 1.,
            value_replacements_at_integration_points: List[dict] = None
    ):
        if value_replacements_at_integration_points is None:
            value_replacements_at_integration_points = [dict()] * len(self._integration_points)

        assert len(value_replacements_at_integration_points) == len(self._integration_points)

        terms = list()
        for integration_point, integration_weight, value_replacement_at_integration_point \
                in zip(
                    self._integration_points,
                    self._integration_weights,
                    value_replacements_at_integration_points
                ):

            replacement_rules = {
                symbol: value
                for symbol, value in zip(symbolic_coordinates, integration_point)
            }
            replacement_rules.update(value_replacement_at_integration_point)

            terms.append(
                sy.Mul(integration_weight, (expression * det).xreplace(replacement_rules))
            )

        result = sy.Add(*terms)
        return result


class Element(object):
    def __init__(
            self,
            coordinates: sy.MatrixBase,
            reference_coordinates: sy.MatrixBase,
            reference_shapes: sy.MatrixBase,
            nodes: ARRAY_LIKE,
            integration_scheme: Optional[IntegrationScheme]
    ):
        self.coordinates = coordinates
        self.reference_coordinates = reference_coordinates
        self.reference_shapes = reference_shapes
        self.nodes = nodes
        self.integration_scheme = integration_scheme

        # jac[i, j] = d(coordinates[i])/d(ref_coordinates[j])
        self.jacobian = self.coordinates.jacobian(self.reference_coordinates)
        self.inv_jacobian = self.jacobian.inv()
        self.det = self.jacobian.det()

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

    def diff(self, axis: Union[int, slice, range] = None, values_at_nodes: ARRAY_LIKE = None):
        if values_at_nodes is None:
            values_at_nodes = self.nodes

        values_at_nodes_mat = self._assure_matrix(values_at_nodes)
        return values_at_nodes_mat.T * self.d_shapes_d_coordinates[:, axis]

    @classmethod
    def create_ref_space_element(
            cls,
            nodes: ARRAY_LIKE,
            shape_powers: ARRAY_LIKE,
            integration_scheme: Optional[IntegrationScheme]
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

        return cls(coordinates, coordinates, shape_functions, nodes, integration_scheme)

    def create_real_space_element(self, nodes: ARRAY_LIKE):
        nodes = nodes[:, :self.count_dimensions()]
        return Element(
            self(nodes),
            self.reference_coordinates,
            self.reference_shapes,
            nodes,
            self.integration_scheme
        )


class Deformation(object):
    def __init__(
            self,
            element: Element,
            internal_energy_density: Union[sy.Symbol, float],
            displacements_at_nodes: ARRAY_LIKE,
            stress_tensors_at_integration_points: ARRAY_LIKE,
    ):

        self.num_dim = element.count_dimensions()
        self.element = element
        self.internal_energy_density = internal_energy_density
        self.displacements_at_nodes = displacements_at_nodes
        self.stress_tensors_at_integration_points = stress_tensors_at_integration_points

        sym_stress_tensor = S[:self.num_dim, :self.num_dim]
        self.d_displacement_d_coordinate = displacements_at_nodes.T * element.d_shapes_d_coordinates
        self.deformation_gradient = sy.eye(self.num_dim) + self.d_displacement_d_coordinate
        self.piola_stress_tensor = self.deformation_gradient.det() * sym_stress_tensor * self.deformation_gradient.inv().T

    @property
    @lru_cache()
    def configurational_stress_tensor_mbf(self):
        configurational_stress_tensor = (
                self.internal_energy_density * sy.eye(self.num_dim)
                - self.deformation_gradient.T * self.piola_stress_tensor
        )

        return configurational_stress_tensor

    @property
    @lru_cache()
    def configurational_stress_tensor_dbf(self):
        configurational_stress_tensor = (
                self.internal_energy_density * sy.eye(self.num_dim)
                - self.d_displacement_d_coordinate.T * self.piola_stress_tensor
        )
        return configurational_stress_tensor

    def _configurational_force(self, configurational_stress_tensor):
        return self.element.integration_scheme.integrate(
            symbolic_coordinates=self.element.reference_coordinates,
            det=self.element.det,
            expression=(
                    self.element.d_shapes_d_coordinates
                    * configurational_stress_tensor.T
                    * self.element.det
            ),
            value_replacements_at_integration_points=[
                {S: stress_tensor}
                for stress_tensor in self.stress_tensors_at_integration_points
            ]
        )

    @property
    @lru_cache()
    def configurational_force_mbf(self):
        return self._configurational_force(self.configurational_stress_tensor_mbf)

    @property
    @lru_cache()
    def configurational_force_dbf(self):
        return self._configurational_force(self.configurational_stress_tensor_dbf)
