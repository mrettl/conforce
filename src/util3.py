from typing import Union, List, Dict, Any, Iterable, TextIO, Optional
from itertools import product, chain, count
import subprocess

import numpy as np
import sympy as sy

from sympy.codegen import ast
from sympy.printing.c import C99CodePrinter

import ele_def


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


def eval_H(R: sy.MatrixBase, nodes_: np.ndarray, powers_: np.ndarray):
    nodes_ = np.array(nodes_, dtype=float)
    powers_ = np.array(powers_, dtype=int)

    num_points_ = nodes_.shape[0]
    num_powers_ = powers_.shape[0]
    assert num_points_ == num_powers_

    POWERS_ = sy.Matrix([
        sy.Mul(*[
            sy.Pow(rst, power)
            for rst, power in zip(R.T.tolist()[0], powers)
        ]).nsimplify()
        for powers in powers_
    ])
    """symbolic powers of shape functions"""

    A_ = sy.Matrix([
        POWERS_.T.subs({
            rst: xyz
            for rst, xyz
            in zip(R.T.tolist()[0], point)
        })
        for point in nodes_
    ]).T
    """A_[i, j] = i-th shape power evaluated at j-th point"""

    COEFFICIENTS_ = A_.inv()
    """
    coefficients solve the system COEFFICIENTS_ * A_ = I,
    such that every shape function is one at exactly one point and zero at the other points
    """

    H_ = COEFFICIENTS_ * POWERS_
    """shape functions are the combination of coefficients and shapes powers"""

    return H_


def eval_dH_dR(H: sy.MatrixBase, R: sy.MatrixBase):
    n_ = H.shape[0]
    d_ = R.shape[0]

    dH_dR_ = sy.zeros(n_, d_)
    for i in range(n_):
        for j in range(d_):
            dH_dR_[i, j] = H[i].diff(R[j])

    return dH_dR_


def eval_J(X_at_nodes: Union[sy.MatrixBase, np.ndarray], dH_dR: sy.MatrixBase):
    return X_at_nodes.T * dH_dR


def eval_dH_dX(dH_dR: sy.MatrixBase, J: sy.MatrixBase):
    return dH_dR * inverse(J)


def eval_dU_dX(U_at_nodes: Union[sy.MatrixBase, np.ndarray], dH_dX: sy.MatrixBase):
    return sy.Matrix(U_at_nodes.T) * dH_dX


def eval_F(d_: int, dU_dX: sy.MatrixBase):
    return sy.eye(d_) + dU_dX


def eval_P(F: sy.MatrixBase, S: sy.MatrixBase):
    return S * inverse(F).T * F.det()


def eval_CS_mbf(d_: int, e: sy.Expr, F: sy.MatrixBase, P: sy.MatrixBase):
    return e * sy.eye(d_) - F.T * P


def eval_CS_dbf(d_: int, e: sy.Expr, dU_dX: sy.MatrixBase, P: sy.MatrixBase):
    return e * sy.eye(d_) - dU_dX.T * P


def eval_CF(
        dH_dX: sy.MatrixBase,
        CS: sy.MatrixBase,
        J: sy.MatrixBase,
        int_weights: np.ndarray,
        int_points_replacements: List[Dict[sy.Expr, Any]]):

    CF_contributions = list()
    for w, int_point_replacement in zip(int_weights, int_points_replacements):
        expr_ = dH_dX * (CS.T * J.det() * w)

        expr_ = expr_.xreplace(int_point_replacement)
        CF_contributions.append(expr_)

    CF_at_nodes_ = sy.MatAdd(*CF_contributions).doit()
    return CF_at_nodes_


def create_replacement_rules(R, list_of_R_):
    return [
        {
                r_i: r_i_
                for r_i, r_i_
                in zip(R, R_)
        }
        for R_ in np.array(list_of_R_)
    ]


def apply_replacement_rules(expr_: sy.MatrixExpr, replacement_rules):
    num_rules = len(replacement_rules)
    expr_shape = [int(dim) for dim in expr_.shape]
    new_shape = [num_rules] + expr_shape
    expr_slice = [slice(None)] * len(expr_shape)  # type: list

    replaced_expressions = sy.MutableDenseNDimArray(np.zeros(new_shape))
    for i in range(num_rules):
        new_slice = tuple([i] + expr_slice)
        replaced_expressions[new_slice] = expr_.xreplace(replacement_rules[i])

    return replaced_expressions.as_immutable()


def map_symbols_to_expressions(array, array_):
    return {
        array[idx]: array_[idx]
        for idx in product(*[range(dim) for dim in array.shape])
    }


def compute_CF(
        R_at_nodes_, powers_,
        R_at_int_points, int_weights,
        X_at_nodes_, U_at_nodes_, S_at_int_points_, e_,
        is_dbf: bool = True):

    mat_symbols_to_expression = dict()
    symbols_to_expressions = dict()
    n_, d_ = R_at_nodes_.shape
    dim_n = list(range(n_))

    R = create_symbolic_matrix("{row}", ["r", "s", "t"][:d_], 1)
    dim_r = [r.name for r in R]

    X = create_symbolic_matrix("{row}", ["x", "y", "z"][:d_], 1, *R)
    dim_x = [x.name for x in X]

    U = create_symbolic_matrix("U{row}", dim_x, 1, *R)
    S = create_symbolic_matrix("S{row}{col}", dim_x, dim_x, *R, is_symmetric=True)

    e = sy.Symbol("e", real=True)
    symbols_to_expressions[e] = e_

    #
    nodes_replacements = create_replacement_rules(R, R_at_nodes_)
    int_points_replacements = create_replacement_rules(R, R_at_int_points)

    X_at_nodes = sy.Matrix(apply_replacement_rules(X, nodes_replacements)[:, :, 0]).as_immutable()
    mat_symbols_to_expression[X_at_nodes] = X_at_nodes_

    U_at_nodes = sy.Matrix(apply_replacement_rules(U, nodes_replacements)[:, :, 0]).as_immutable()
    mat_symbols_to_expression[U_at_nodes] = U_at_nodes_

    S_at_int_points = apply_replacement_rules(S, int_points_replacements)
    mat_symbols_to_expression[S_at_int_points] = S_at_int_points_

    #
    H = create_symbolic_matrix("H{row}", dim_n, 1, *R)
    H_ = eval_H(R, R_at_nodes_, powers_).doit()
    mat_symbols_to_expression[H] = H_

    dH_dR = create_symbolic_matrix("dH{row}_d{col}", dim_n, dim_r, *R)
    dH_dR_ = eval_dH_dR(H_, R).doit()
    mat_symbols_to_expression[dH_dR] = dH_dR_

    J = create_symbolic_matrix("J{row}{col}", dim_x, dim_r, *R)
    J_ = eval_J(X_at_nodes, dH_dR_).doit()
    mat_symbols_to_expression[J] = J_

    dH_dX = create_symbolic_matrix("dH{row}_d{col}", dim_n, dim_x, *R)
    dH_dX_ = eval_dH_dX(dH_dR_, J).doit()
    mat_symbols_to_expression[dH_dX] = dH_dX_

    dU_dX = create_symbolic_matrix("dU{row}_d{col}", dim_x, dim_x, *R)
    dU_dX_ = eval_dU_dX(U_at_nodes, dH_dX).doit()
    mat_symbols_to_expression[dU_dX] = dU_dX_

    #
    F = create_symbolic_matrix("F{row}{col}", dim_x, dim_x, *R)
    F_ = eval_F(d_, dU_dX).doit()
    mat_symbols_to_expression[F] = F_

    P = create_symbolic_matrix("P{row}{col}", dim_x, dim_x, *R)
    P_ = eval_P(F, S).doit()
    mat_symbols_to_expression[P] = P_

    CS = create_symbolic_matrix("CS{row}{col}", dim_x, dim_x, *R)
    if is_dbf:
        CS_ = eval_CS_dbf(d_, e, dU_dX, P).doit()
    else:
        CS_ = eval_CS_mbf(d_, e, F, P).doit()

    mat_symbols_to_expression[CS] = CS_

    #
    CF_at_nodes = sy.MatrixSymbol("CF", n_, d_)
    CF_at_nodes_ = eval_CF(
        dH_dX,
        CS,
        J,
        int_weights,
        int_points_replacements
    ).doit()
    mat_symbols_to_expression[CF_at_nodes] = CF_at_nodes_

    symbols_to_expressions.update(chain(
        *[
            map_symbols_to_expressions(symbols, expressions).items()
            for symbols, expressions in mat_symbols_to_expression.items()
        ]
    ))

    return R, CF_at_nodes, symbols_to_expressions


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

    def collect_expressions_for_array(self, result_array: sy.MatrixBase):
        for i, j in product(*[range(dim) for dim in result_array.shape]):
            self.collect_expressions_for_symbol(result_array[i, j])

    def collect_expressions_for_symbol(self, result: sy.Expr):
        if result in self.used_symbols_to_names:
            return self.used_symbols_to_names[result]

        print(result)

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

    def doit(self, input_symbols, result_array: sy.MatrixBase, cse: bool = True):
        self.collect_expressions_for_array(result_array)
        return self.to_assignments(input_symbols, cse)


class CCodePrinter(C99CodePrinter):
    def _print_MatrixElement(self, expr: sy.Expr):
        return f"{self._print(str(expr.args[0]))}[{']['.join([str(i) for i in expr.indices])}]"

    def _print_Indexed(self, expr: sy.Indexed):
        return self._print_MatrixElement(expr)


class CPyCodeGenerator(object):
    def __init__(self, name: str, compile_at_exit: bool):
        self._name = str(name)
        self._compile_at_exit = bool(compile_at_exit)

        self._py_file_handle = None  # type: Optional[TextIO]
        self._c_file_handle = None  # type: Optional[TextIO]

    def __enter__(self):
        self._py_file_handle = open(self._name + ".py", "w")
        self._c_file_handle = open(self._name + ".c", "w")

        return self

    def __exit__(self, *_):
        self._py_file_handle.close()
        self._c_file_handle.close()

        if self._compile_at_exit:
            self._compile_c_code()

    def _compile_c_code(self):
        subprocess.call(f"gcc -shared -O3 -o {self._name}.dll {self._name}.c", shell=True)

    def write_headers(self):
        self._c_file_handle.write(
            "#ifdef _WIN32\n"
            "#    define CF_API __declspec(dllexport)\n"
            "#else\n"
            "#    define CF_API\n"
            "#endif\n"
            "\n"
            "#include <math.h>\n"  # contains power
            "#include <stdint.h>\n"  # contains size_t
            "#include <stddef.h>\n"  # contains size_t
            "#include <stdio.h>\n"  # contains print
        )

        self._py_file_handle.write(
            f"import numpy as np\n"
            f"import ctypes\n"
            f"import os\n"
            f"\n"
            f"# load c library\n"
            f"if os.name == 'nt':\n"
            f"    # windows\n"
            f"    lib = ctypes.cdll.LoadLibrary('{self._name}.dll')\n"
            f"else:\n"
            f"    # linux\n"
            f"    lib = ctypes.cdll.LoadLibrary('./{self._name}.so')\n"
            f"\n"
        )

    def write_function(
            self,
            assignments: List[ast.Assignment],
            element_typ: str,
            e: sy.Symbol,
            X_at_nodes: sy.IndexedBase,
            U_at_nodes: sy.IndexedBase,
            S_at_int_points: sy.IndexedBase,
            CF_at_nodes: sy.IndexedBase,
            is_dbf: bool
    ):

        c_printer = CCodePrinter()
        n_, d_ = X_at_nodes.shape
        ips_ = S_at_int_points.shape[0]

        # c code for one element
        self._c_file_handle.write(
            f"\n"
            f"/**\n"
            f" * Computes the configurational forces for one element of typ {element_typ}.\n"
            f" * The element has n={n_} nodes and ips={ips_} integration points.\n"
            f" * The computation is {'displacement' if is_dbf else 'motion'} based.\n"
            f" * A static load case is assumed.\n"
            f" * \n"
            f" * @param[in] {e} internal energy density\n"
            f" * @param[in] {X_at_nodes} Coordinates at n nodes of the element.\n"
            f" * @param[in] {U_at_nodes} Displacements at n nodes of the element.\n"
            f" * @param[in] {S_at_int_points} Symmetric stress tensors at ips integration points.\n"
            f" * @param[out] {CF_at_nodes} The resulting configurational forces for n nodes are written \n"
            f" *             to this array.\n"
            f" */\n"
            f"void compute_static_{'dbf' if is_dbf else 'mbf'}_CF_for_{element_typ}(\n"
            f"        double {e}, \n"
            f"        double {X_at_nodes}[{n_}][{d_}], \n"
            f"        double {U_at_nodes}[{n_}][{d_}], \n"
            f"        double {S_at_int_points}[{ips_}][{d_}][{d_}], \n"
            f"        double {CF_at_nodes}[{n_}][{d_}] \n"
            f") {{\n"
            f"    // COMPUTER GENERATED CODE:\n"
        )

        for assignment in assignments:
            lhs = assignment.lhs
            rhs = assignment.rhs

            lhs_code = c_printer.doprint(lhs)
            rhs_code = c_printer.doprint(rhs)
            if isinstance(lhs, sy.matrices.expressions.matexpr.MatrixElement):
                self._c_file_handle.write(f"    {lhs_code} = {rhs_code};\n")
            else:
                self._c_file_handle.write(f"    double {lhs_code} = {rhs_code};\n")

        self._c_file_handle.write("}\n")

        # c code for multiple elements
        self._c_file_handle.write(
            f"\n"
            f"/**\n"
            f" * Computes the configurational forces for num_elem elements of typ {element_typ}.\n"
            f" * Each element has n={n_} nodes and ips={ips_} integration points.\n"
            f" * The computation is {'displacement' if is_dbf else 'motion'} based.\n"
            f" * A static load case is assumed.\n"
            f" * \n"
            f" * @param[in] num_elem number of elements\n"
            f" * @param[in] {e} internal energy density\n"
            f" * @param[in] {X_at_nodes} Coordinates at n nodes of num_elem elements.\n"
            f" * @param[in] {U_at_nodes} Displacements at n nodes of num_elem elements.\n"
            f" * @param[in] {S_at_int_points} Symmetric stress tensors at ips integration points \n"
            f" *            for num_elem elements.\n"
            f" * @param[out] {CF_at_nodes} The resulting configurational forces for n nodes and \n"
            f" *             num_elem elements are written to this array.\n"
            f" */\n"
            f"CF_API void compute_static_{'dbf' if is_dbf else 'mbf'}_CF_for_multiple_{element_typ}(\n"
            f"        int num_elem, \n"
            f"        double {e}[], \n"
            f"        double {X_at_nodes}[][{n_}][{d_}], \n"
            f"        double {U_at_nodes}[][{n_}][{d_}], \n"
            f"        double {S_at_int_points}[][{ips_}][{d_}][{d_}], \n"
            f"        double {CF_at_nodes}[][{n_}][{d_}] \n"
            f") {{\n"
            f"    // COMPUTER GENERATED CODE:\n"
            f"    #pragma omp parallel for\n"
            f"    for (size_t i=0; i<num_elem; i++) {{\n"
            f"        compute_static_{'dbf' if is_dbf else 'mbf'}_CF_for_{element_typ}(\n"
            f"            {e}[i], \n"
            f"            {X_at_nodes}[i], \n"
            f"            {U_at_nodes}[i], \n"
            f"            {S_at_int_points}[i], \n"
            f"            {CF_at_nodes}[i]\n"
            f"        );\n"
            f"    }}\n"
            f"}}\n"
        )

        # python code for multiple elements
        self._py_file_handle.write(f"""
def compute_static_{'dbf' if is_dbf else 'mbf'}_CF_for_multiple_{element_typ}(
        {e},
        {X_at_nodes},
        {U_at_nodes},
        {S_at_int_points}):

    {e} = np.ascontiguousarray({e}, dtype=np.float64)
    {X_at_nodes} = np.ascontiguousarray({X_at_nodes}, dtype=np.float64)
    {U_at_nodes} = np.ascontiguousarray({U_at_nodes}, dtype=np.float64)
    {S_at_int_points} = np.ascontiguousarray({S_at_int_points}, dtype=np.float64)

    num_elem = {X_at_nodes}.shape[0]

    {CF_at_nodes} = np.zeros((num_elem, {n_}, {d_}), dtype=np.float64)

    assert {e}.shape == (num_elem, )
    assert {U_at_nodes}.shape == (num_elem, {n_}, {d_})
    assert {U_at_nodes}.shape == (num_elem, {n_}, {d_})
    assert {S_at_int_points}.shape == (num_elem, {ips_}, {d_}, {d_})

    lib.compute_static_{'dbf' if is_dbf else 'mbf'}_CF_for_multiple_{element_typ}(
        ctypes.c_size_t(num_elem),
        {e}.ctypes.data_as(ctypes.c_void_p),
        {X_at_nodes}.ctypes.data_as(ctypes.c_void_p),
        {U_at_nodes}.ctypes.data_as(ctypes.c_void_p),
        {S_at_int_points}.ctypes.data_as(ctypes.c_void_p),
        {CF_at_nodes}.ctypes.data_as(ctypes.c_void_p)
    )

    return {CF_at_nodes}
""")


def main():
    is_dbf = True
    typ = "CPE4"  # "C3D20"  # "CPE4"
    R_at_nodes_ = ele_def.ref_nodes[typ]
    powers_ = ele_def.poly_power[typ]
    R_at_int_points_ = ele_def.int_points[typ]
    int_weights_ = ele_def.int_weights[typ]

    n_, d_ = R_at_nodes_.shape
    ips_ = len(int_weights_)
    int_point_names = [f"ip{idx}" for idx in range(ips_)]
    node_names = [f"node{idx}" for idx in range(n_)]

    X_at_nodes = sy.IndexedBase("X_at_nodes", shape=(n_, d_))
    U_at_nodes = sy.IndexedBase("U_at_nodes", shape=(n_, d_))
    S_at_int_points = sy.IndexedBase("S_at_int_points", shape=(ips_, d_, d_))
    e = sy.Symbol("e", real=True)

    # compute CF
    R, CF_at_nodes, symbols_to_expressions = compute_CF(
        R_at_nodes_,
        powers_,
        R_at_int_points_,
        int_weights_,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points,
        e,
        is_dbf=is_dbf
    )

    # create abstract code assignments
    assignments = TermCollector(
        symbols_to_expressions=symbols_to_expressions,
        R=R,
        points=np.concatenate([R_at_int_points_, R_at_nodes_], dtype=float),
        point_names=[*int_point_names, *node_names]
    ).doit(
        input_symbols=[e, X_at_nodes, U_at_nodes, S_at_int_points],
        result_array=CF_at_nodes,
        cse=True
    )

    # generate and compile c and python code
    with CPyCodeGenerator("test3", compile_at_exit=True) as compiler:
        compiler.write_headers()

        compiler.write_function(
            assignments=assignments,
            element_typ=typ,
            e=e,
            X_at_nodes=X_at_nodes,
            U_at_nodes=U_at_nodes,
            S_at_int_points=S_at_int_points,
            CF_at_nodes=CF_at_nodes,
            is_dbf=is_dbf
        )

    print("ok")


if __name__ == '__main__':
    main()
