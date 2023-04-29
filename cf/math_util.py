from itertools import product, count
from typing import Union, Iterable

import numpy as np
import sympy as sy
from sympy.codegen import ast


def create_symbolic_matrix(name: str, row_spec, col_spec, *variables: sy.Symbol, is_symmetric: bool = False):
    if isinstance(row_spec, int):
        row_names = list(range(row_spec))
    else:
        row_names = row_spec
    rows = len(row_names)

    if isinstance(col_spec, int):
        col_names = list(range(col_spec))
    else:
        col_names = col_spec
    cols = len(col_names)

    M = sy.zeros(rows, cols)
    cell_kwargs = dict(real=True)
    for row in range(rows):
        for col in range(cols):
            if is_symmetric and col < row:
                M[row, col] = M[col, row]
            else:
                cell_name = name.format(row=row_names[row], col=col_names[col])
                if len(variables) == 0:
                    cell = sy.Symbol(
                        cell_name,
                        **cell_kwargs
                    )
                else:
                    cell = sy.Function(
                        cell_name,
                        **cell_kwargs
                    )(*variables)
                M[row, col] = cell

    return M.as_immutable()


def inverse(matrix: sy.MatrixBase) -> sy.Matrix:
    """
    Inverts a symbolic matrix 2x2 or 3x3 matrix.
    This method is faster than the corresponding sympy method.

    >>> m22 = sy.MatrixSymbol('m', 2, 2).as_explicit()
    >>> inverse(m22) == m22.inv(method='ADJ')
    True

    >>> m33 = sy.MatrixSymbol('m', 3, 3).as_explicit()
    >>> inverse(m33) == m33.inv(method='ADJ')
    True

    :param matrix: SymPy matrix of shape (3, 3) or (2, 2)
    :return: inverted SymPy matrix of shape (3, 3) or (2, 2)
    """
    m = matrix
    shape = m.shape

    if shape == (2, 2):
        inv = sy.zeros(2, 2)
        determinant = m[0, 0] * m[1, 1] - m[0, 1] * m[1, 0]
        determinant_inv = 1 / determinant
        inv[0, 0] = m[1, 1] * determinant_inv
        inv[0, 1] = -1 * m[0, 1] * determinant_inv
        inv[1, 0] = -1 * m[1, 0] * determinant_inv
        inv[1, 1] = m[0, 0] * determinant_inv

    elif shape == (3, 3):
        inv = sy.zeros(3, 3)
        determinant = (
            m[0, 0]*m[1, 1]*m[2, 2]
            + m[0, 1]*m[1, 2]*m[2, 0]
            + m[0, 2]*m[1, 0]*m[2, 1]
            - m[0, 0]*m[1, 2]*m[2, 1]
            - m[0, 1]*m[1, 0]*m[2, 2]
            - m[0, 2]*m[1, 1]*m[2, 0]
        )

        determinant_inv = 1 / determinant
        inv[0, 0] = (m[1, 1]*m[2, 2] - m[1, 2]*m[2, 1]) * determinant_inv
        inv[0, 1] = (m[0, 2]*m[2, 1] - m[0, 1]*m[2, 2]) * determinant_inv
        inv[0, 2] = (m[0, 1]*m[1, 2] - m[0, 2]*m[1, 1]) * determinant_inv
        inv[1, 0] = (m[1, 2]*m[2, 0] - m[1, 0]*m[2, 2]) * determinant_inv
        inv[1, 1] = (m[0, 0]*m[2, 2] - m[0, 2]*m[2, 0]) * determinant_inv
        inv[1, 2] = (m[0, 2]*m[1, 0] - m[0, 0]*m[1, 2]) * determinant_inv
        inv[2, 0] = (m[1, 0]*m[2, 1] - m[1, 1]*m[2, 0]) * determinant_inv
        inv[2, 1] = (m[0, 1]*m[2, 0] - m[0, 0]*m[2, 1]) * determinant_inv
        inv[2, 2] = (m[0, 0]*m[1, 1] - m[0, 1]*m[1, 0]) * determinant_inv

    else:
        raise NotImplementedError()

    return inv


def create_replacement_rules(R, list_of_R_):
    return [
        {
                r_i: r_i_
                for r_i, r_i_
                in zip(R, R_)
        }
        for R_ in np.array(list_of_R_)
    ]


def apply_replacement_rules(expr_: Union[sy.MatrixExpr, sy.MatrixBase], replacement_rules):
    num_rules = len(replacement_rules)
    if hasattr(expr_, "shape"):
        expr_shape = [int(dim) for dim in expr_.shape]
        expr_slice = [slice(None)] * len(expr_shape)  # type: list
    else:
        expr_shape = []
        expr_slice = []  # type: list

    new_shape = [num_rules] + expr_shape

    replaced_expressions = sy.MutableDenseNDimArray(np.zeros(new_shape))
    for i in range(num_rules):
        new_slice = tuple([i] + expr_slice)
        replaced_expressions[new_slice] = expr_.xreplace(replacement_rules[i])

    return replaced_expressions.as_immutable()


class TermCollector:
    def __init__(self, symbols_to_expressions: dict, R: sy.MatrixBase, points: Iterable, point_names: Iterable[str]):
        self.symbols_to_expressions = symbols_to_expressions
        self.R = R

        self.used_symbols_to_names = dict()
        self.names_to_used_expressions = dict()
        self.point_name_map = {
            tuple(point): name
            for point, name
            in zip(np.array(points, dtype=float), point_names)
        }

    def reset(self):
        self.used_symbols_to_names.clear()
        self.names_to_used_expressions.clear()

    def collect_expressions_for_array(self, result_array: sy.MatrixBase):
        for idx in product(*[range(dim) for dim in result_array.shape]):
            self.collect_expressions_for_symbol(result_array[idx])

    def collect_expressions_for_symbol(self, result: sy.Expr):
        if result in self.used_symbols_to_names:
            return self.used_symbols_to_names[result]

        if result.is_Function:
            point = tuple(np.array(result.args, dtype=float))
            point_name = self.point_name_map[point]
            new_result_symbol = sy.Symbol(
                result.name + "_at_" + point_name,
                real=True
            )

            if result in self.symbols_to_expressions:
                used_expression_ = self.symbols_to_expressions[result]

            else:
                original_function = result.func(*self.R)
                original_expression_ = self.symbols_to_expressions[original_function]
                evaluation_rules = create_replacement_rules(self.R, [result.args])[0]
                used_expression_ = original_expression_.xreplace(evaluation_rules)
        else:
            if result in self.symbols_to_expressions:
                used_expression_ = self.symbols_to_expressions[result]

            else:
                used_expression_ = None

            new_result_symbol = result

        if used_expression_ == result:
            used_expression_ = None

        elif hasattr(used_expression_, "atoms"):
            replacements = {
                next_symbol: self.collect_expressions_for_symbol(next_symbol)
                for next_symbol in used_expression_.atoms(sy.Symbol)
            }
            replacements.update({
                next_function: self.collect_expressions_for_symbol(next_function)
                for next_function in used_expression_.atoms(sy.Function)
            })

            used_expression_ = used_expression_.xreplace(replacements)

        self.used_symbols_to_names[result] = new_result_symbol
        self.names_to_used_expressions[new_result_symbol] = used_expression_

        return new_result_symbol

    def to_assignments(self, input_symbols, cse: bool = True):
        input_symbols = {
            str(symbol)
            for symbol in input_symbols
        }
        assignments = ast.CodeBlock(*[
            ast.Assignment(lhs, rhs)
            for lhs, rhs in self.names_to_used_expressions.items()
            if str(lhs) not in input_symbols
        ])

        if cse:
            assignments = assignments.cse(symbols=(sy.Symbol(f"tmp{i}", real=True) for i in count()))

        return assignments

    def collect_assignments(self, input_symbols, result_array: sy.MatrixBase, cse: bool = True):
        self.reset()
        self.collect_expressions_for_array(result_array)
        return self.to_assignments(input_symbols, cse)
