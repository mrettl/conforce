"""
This module contains methods to work with symbolic expressions.
"""

from itertools import product, count, chain
from typing import Union, Iterable

import numpy as np
import sympy as sy
from sympy.codegen import ast


def create_symbolic_matrix(
        name_format: str, row_spec, col_spec,
        *variables: sy.Symbol, is_symmetric: bool = False):
    """
    Create a sympy matrix with symbolic components.

    **Examples**

    Create a column vector:

    >>> create_symbolic_matrix("A{row}", 2, 1)
    Matrix([
    [A0],
    [A1]])

    Create a symmetric matrix whose components are functions of a.

    >>> create_symbolic_matrix("B{row}{col}", ["x", "y"], ["r", "s"], sy.symbols("a"), is_symmetric=True)
    Matrix([
    [Bxr(a), Bxs(a)],
    [Bxs(a), Bys(a)]])

    :param name_format: A string defining the component names. The occurrence of `{row}` and `{col}`
        is replaced with the row and column name.
    :param row_spec: An integer or a sequence of strings definig the number of rows and the row names.
        An integer defines the number of rows and the row names as "0", "1", ..., `number of rows - 1`.
        If a sequence of string is given, the number of rows is the length of the sequence
        and the row names are the entries in the sequence.
    :param col_spec: An integer or a sequence of strings definig the number of columns and the column names.
        An integer defines the number of columns and the column names as "0", "1", ..., `number of columns - 1`.
        If a sequence of string is given, the number of columns is the length of the sequence
        and the column names are the entries in the sequence.
    :param variables: The matrix components are functions of these variables or if variables is an empty sequence,
        the matrix components are independent symbols.
    :param is_symmetric: Create a symmetric matrix if True.
    :return: immutable sympy matrix
    """
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
                cell_name = name_format.format(row=row_names[row], col=col_names[col])
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
    Invert a symbolic matrix 2x2 or 3x3 matrix.
    This method is faster than the corresponding sympy method.

    **Examples**

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


def create_replacement_rules(symbols, *symbol_values):
    """
    Create a list of dictionaries.
    Each dictionary maps symbols to values.

    .. seealso:: :py:func:`apply_replacement_rules`

    :param symbols: A sequence of length (s, ) containing sympy symbols
    :param symbol_values: array of shape (n, s) containing values that should replace sympy symbols.
    :return: List of length n containing dictionaries that map s symbols to values.
    """
    return [
        {
                r_i: r_i_
                for r_i, r_i_
                in zip(symbols, R_)
        }
        for R_ in np.array(symbol_values)
    ]


def apply_replacement_rules(expression: Union[sy.MatrixExpr, sy.MatrixBase], *replacement_rules):
    """
    A Replacement rule define symbols and values.
    This method searches in the expression for symbols and replaces them by values
    as defined by a replacement rule.
    Each replacement rule is applied independently.

    **Examples**

    In the following expression some symbols should be replaced by values.

    >>> r, s, t = sy.symbols("r s t")
    >>> expr = 1*r + 2*s + 3*t

    Three replacement rules map `r` and `s` to values.

    >>> repl = create_replacement_rules(
    ...     [r, s],
    ...     [1, 0],
    ...     [0, 1],
    ...     [1, 1]
    ... )
    >>> repl
    [{r: 1, s: 0}, {r: 0, s: 1}, {r: 1, s: 1}]

    Each of the three replacement rules is applied independently.
    The `t` symbol is not replaceed.

    >>> apply_replacement_rules(expr, *repl)
    [3*t + 1, 3*t + 2, 3*t + 3]

    .. seealso:: :py:func:`create_replacement_rules`

    :param expression: sympy expression
    :param replacement_rules: A replacement rule is a dictionary mapping symbols to values.
        (see :py:func:`create_replacement_rules`).
    :return: List containing replaced expression. The i-th entry in the list corresponds
        to the i-th replacement rule.
    """
    num_rules = len(replacement_rules)
    if hasattr(expression, "shape"):
        expr_shape = [int(dim) for dim in expression.shape]
        expr_slice = [slice(None)] * len(expr_shape)  # type: list
    else:
        expr_shape = []
        expr_slice = []  # type: list

    new_shape = [num_rules] + expr_shape

    replaced_expressions = sy.MutableDenseNDimArray(np.zeros(new_shape))
    for i in range(num_rules):
        new_slice = tuple([i] + expr_slice)
        replaced_expressions[new_slice] = expression.xreplace(replacement_rules[i])

    return replaced_expressions.as_immutable()


def expand_matrices_in_symbols_to_expressions(symbols_to_expressions):
    """
    A mapping of symbols to expression might contain matrix symbols.
    This method generates a mapping for each matrix component.

    **Examples**

    The matrix `M` is mapped to the `2*T`.
    This method returns additional mappings
    for all matrix components of `M`.

    >>> M = sy.MatrixSymbol("M", 2, 1)
    >>> T = sy.MatrixSymbol("T", 2, 1)
    >>> expand_matrices_in_symbols_to_expressions({
    ...     M: 2 * T,
    ...     T[0]: 1,
    ...     T[1]: 1,
    ... })
    {M[0, 0]: 2*T[0, 0], M[1, 0]: 2*T[1, 0]}


    :param symbols_to_expressions: Dictionary mapping symbols to expressions
    :return: Dictionary with additional mappings of matrix components to expressions
    """
    return dict(chain(
        *[
            {
                symbols[idx]: expressions[idx]
                for idx in product(*[range(int(dim)) for dim in symbols.shape])
            }.items()
            for symbols, expressions in symbols_to_expressions.items()
            if hasattr(symbols, "shape") and not hasattr(symbols, "indices")
        ]
    ))


class TermCollector:
    _MATRIX_ELEMENT_TYPE = type(sy.MatrixSymbol("TEST", 1, 1)[0, 0])

    def __init__(self, symbols_to_expressions: dict, R: sy.MatrixBase, points: Iterable, point_names: Iterable[str]):
        """
        This class takes a mapping of symbols to expressions
        and provides methods to create an abstract :py:class:`ast.CodeBlock`
        that contains a list of :py:class:`ast.Assignment` objects.
        Those objects can be easily converted to real code.

        **symbols_to_expressions**

        `symbols_to_expressions` maps symbols to sympy expressions.
        Allowed symbols are:

            - :py:class:`sy.Symbol`
            - :py:class:`sy.Function` of `R`
            - :py:class:`sy.MatrixSymbol`
            - :py:class:`sy.Indexed`
            - :py:class:`sympy.matrices.expressions.matexpr.MatrixElement`

        Any sympy expression is allowed that has the same shape as the
        coresponding symbol.

        Each point at which a :py:class:`sy.Function` is evaluated,
        must be declared in the `points` argument.

        **Examples**

        The matrix `M` is expressed in terms of matrix `T`.
        For each matrix component an :py:class`ast.Assignment` is created.
        The assignments can be executed from the first one to the last one,
        because it is guaranteed that a dependent variable occurs after
        all its dependencies are computed.

        >>> M = sy.MatrixSymbol("M", 2, 1)
        >>> T = sy.MatrixSymbol("T", 2, 1)
        >>> collector = TermCollector(
        ...     symbols_to_expressions={
        ...         M: 2 * T,
        ...         T[0]: 1,
        ...         T[1]: 1,
        ...     },
        ...     R=sy.Matrix([]),
        ...     points=[],
        ...     point_names=[]
        ... )
        >>> collector.collect_assignments(
        ...     input_symbols=[],
        ...     result_symbol=M
        ... )
        CodeBlock(
        Assignment(T[0, 0], 1),
        Assignment(M[0, 0], 2*T[0, 0]),
        Assignment(T[1, 0], 1),
        Assignment(M[1, 0], 2*T[1, 0])
        )

        Another example demonstrates the use of a :py:class:`sy.Function`.
        The function is evaluated at `r=0` and `r=1`.
        For this reason two points and point names are defined.
        The point names are used to define variable names for the function value at a specific point.
        E.g. `g_at_p0` is the variable name of `g(0)`.

        >>> r, f, k, o, d = sy.symbols("r f k o d")
        >>> g = sy.Function("g")
        >>> collector = TermCollector(
        ...     symbols_to_expressions={
        ...         f: g(0) + k*(g(1) - g(0)),
        ...         k: g(1) - g(0),
        ...         g(r): o + d * r
        ...     },
        ...     R=sy.Matrix([r]),
        ...     points=[[0], [1]],
        ...     point_names=["p0", "p1"]
        ... )

        Three input variables are defined.
        No assignments will be made to input variables, even if it would be possible.
        E.g. `g_at_p0` is `g(0)` and could be expressed as `g_at_p0=o`.

        >>> collector.collect_assignments(
        ...     input_symbols=[o, d, "g_at_p0"],
        ...     result_symbol=f,
        ...     cse=False
        ... )
        CodeBlock(
        Assignment(g_at_p1, d + o),
        Assignment(k, -g_at_p0 + g_at_p1),
        Assignment(f, g_at_p0 + k*(-g_at_p0 + g_at_p1))
        )

        Furthermore, the common subexpression elimination (cse) of sympy
        elimiates duplicate expressions.
        E.g.: `-g_at_p0 + g_at_p1` is computed twice in the above example.
        Cse introduces a new variable (tmp0), instead.

        >>> collector.to_assignments(input_symbols=[o, d, "g_at_p0", "g_at_p1"], cse=True)
        CodeBlock(
        Assignment(tmp0, -g_at_p0 + g_at_p1),
        Assignment(k, tmp0),
        Assignment(f, g_at_p0 + k*tmp0)
        )


        :param symbols_to_expressions: Dictionary mapping symbols to sympy expressions
        :param R: sympy Matrix of shape (d,) containing symbolic reference space coordinates
        :param points: Array of shape (p, d) containing points. Each point that is evalutated
            by a :py:class:`sy.Function` in `symbols_to_expressions` object has to be declared.
        :param point_names: Array of shape (p,) containing names of points defined in `points`.
        """
        symbols_to_expressions.update(
            expand_matrices_in_symbols_to_expressions(symbols_to_expressions)
        )
        self.symbols_to_expressions = symbols_to_expressions
        self.R = R

        self.used_symbols_to_names = dict()
        self.names_to_used_expressions = dict()
        self.point_name_map = {
            tuple(point): name
            for point, name
            in zip(np.array(points, dtype=float), point_names)
        }

    __doc__ = __init__.__doc__

    def reset(self):
        """
        Resets all collected assignments
        """
        self.used_symbols_to_names.clear()
        self.names_to_used_expressions.clear()

    def collect_expressions_for_array(self, result_array: sy.MatrixBase):
        """
        Collect all assignments needed to compute the `result_array`.

        :param result_array: symbolic matrix or symbolic array
        """
        for idx in product(*[range(dim) for dim in result_array.shape]):
            self.collect_expressions_for_symbol(result_array[idx])

    def collect_expressions_for_symbol(self, result: sy.Expr):
        """
        Collect all assignments needed to compute the `result`.

        :param result: symbolic expression (not a matrix or array)
        :return: A symbol with a valid variable name. This is relevant for functions.
            E.g. f(0) is not a valid variable name, but f_at_point0 is valid.
        """
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
                evaluation_rules = create_replacement_rules(self.R, result.args)[0]
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
                next_array: self.collect_expressions_for_symbol(next_array)
                for next_array in used_expression_.atoms(self._MATRIX_ELEMENT_TYPE)
            })
            replacements.update({
                next_function: self.collect_expressions_for_symbol(next_function)
                for next_function in used_expression_.atoms(sy.Function)
            })

            used_expression_ = used_expression_.xreplace(replacements)

        self.used_symbols_to_names[result] = new_result_symbol
        self.names_to_used_expressions[new_result_symbol] = used_expression_

        return new_result_symbol

    def to_assignments(self, input_symbols, cse: bool = True):
        """
        Create a code block with assignments.

        :param input_symbols: given symbols to which no assignments will be made
        :param cse: perform a common subexpression elimination
        :return: A code block containing previously collected assignments.
        """
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

    def collect_assignments(self, input_symbols, result_symbol: Union[sy.MatrixBase, sy.Expr], cse: bool = True):
        """
        Collect all assignments needed to compute the `result` with the given input symbols.

        :param input_symbols: given symbols to which no assignments will be made
        :param result_symbol: symbol that should be computed
        :param cse: perform a common subexpression elimination
        :return: A code block containing previously collected assignments.
        """
        self.reset()
        if hasattr(result_symbol, "shape") and not hasattr(result_symbol, "indices"):
            self.collect_expressions_for_array(result_symbol)
        else:
            self.collect_expressions_for_symbol(result_symbol)

        return self.to_assignments(input_symbols, cse)
