"""
This module contains methods to generate and compile C-code.
"""

from __future__ import annotations
from typing import List, TextIO, Optional
import subprocess
import os

import numpy as np
import sympy as sy

from sympy.codegen import ast
from sympy.printing.c import C99CodePrinter

from conforce_3 import element_definitions
from conforce_3.expressions import Computation
from conforce_3.symbolic_util import TermCollector, expand_matrices_in_symbols_to_expressions, apply_replacement_rules


def write_code_for_all_element_types():
    """
    Generate and compile C- and Python-code for all element types defined in :py:mod:`conforce_3.element_definitions`.
    Places the generated code to :py:mod:`conforce.cf_c`.
    """
    types = element_definitions.R_at_nodes_of_element.keys()

    # generate and compile c and python code
    with CPyCodeCompiler(
            name="cf_c", folder="conforce",
            compile_at_exit=False, write_header_at_enter=True
    ) as compiler:
        for element_type in types:
            for is_dbf in [True, False]:
                print(f"element_type={element_type}, is_dbf={is_dbf}")
                write_code_for_element_type(
                    element_type=element_type,
                    is_dbf=is_dbf,
                    write_F=is_dbf,  # F is independent of is_dbf and is written only once for (is_dbf=True)
                    write_P=is_dbf,  # P is independent of is_dbf and is written only once for (is_dbf=True)
                    write_CS=True,
                    write_CF=True,
                    compiler=compiler
                )


def write_code_for_element_type(
        element_type: str,
        is_dbf: bool,
        write_F: bool,
        write_P: bool,
        write_CS: bool,
        write_CF: bool,
        compiler: CPyCodeCompiler
):
    """
    Generate Code for an element type defined in :py:mod:`conforce_3.element_definitions`.

    :param element_type: Name of the element type
    :param is_dbf: Deformation based or motion based formulation
    :param write_F: use compiler to generate code for the deformation gradient
    :param write_P: use compiler to generate code for the first Piola-Kirchhoff stress tensor
    :param write_CS: use compiler to generate code for the configurational stress tensor
    :param write_CF: use compiler to generate code for the configurational forces
    :param compiler: provides methods to generate code
    """
    R_at_nodes_ = element_definitions.R_at_nodes_of_element[element_type]
    exponents_ = element_definitions.exponents_of_shape_functions_of_element[element_type]
    R_at_int_points_ = element_definitions.R_at_integration_points_of_element[element_type]
    int_weights_ = element_definitions.weights_of_integration_points_of_element[element_type]

    n_, d_ = R_at_nodes_.shape
    ips_ = len(int_weights_)
    int_point_names = [f"ip{idx}" for idx in range(ips_)]
    node_names = [f"node{idx}" for idx in range(n_)]

    X_at_nodes = sy.IndexedBase("X_at_nodes", shape=(n_, d_))
    U_at_nodes = sy.IndexedBase("U_at_nodes", shape=(n_, d_))
    S_at_int_points = sy.IndexedBase("S_at_int_points", shape=(ips_, d_, d_))
    e_at_int_points = sy.IndexedBase("e_at_int_points", shape=(ips_, ))

    # compute all until CF
    computation = Computation(
        R_at_nodes_,
        exponents_,
        R_at_int_points_,
        int_weights_,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points,
        e_at_int_points,
        is_dbf=is_dbf
    )

    # R
    R = computation.R

    # F
    F_at_int_points_ = apply_replacement_rules(
        computation.F, *computation.replacements_by_int_points)
    F_at_int_points = sy.IndexedBase("F_at_int_points", shape=(ips_, d_, d_))
    computation.map_symbolic_to_expression(F_at_int_points, F_at_int_points_)

    # P
    P_at_int_points_ = apply_replacement_rules(
        computation.P, *computation.replacements_by_int_points)
    P_at_int_points = sy.IndexedBase("P_at_int_points", shape=(ips_, d_, d_))
    computation.map_symbolic_to_expression(P_at_int_points, P_at_int_points_)

    # CS
    CS_at_int_points_ = apply_replacement_rules(
        computation.CS, *computation.replacements_by_int_points)
    CS_at_int_points = sy.IndexedBase("CS_at_int_points", shape=(ips_, d_, d_))
    computation.map_symbolic_to_expression(CS_at_int_points, CS_at_int_points_)

    # CF
    CF_at_nodes = computation.CF_at_nodes

    #
    symbols_to_expressions = computation.symbols_to_expressions
    symbols_to_expressions.update(
        expand_matrices_in_symbols_to_expressions(symbols_to_expressions))

    # create abstract code assignments
    term_collector = TermCollector(
        symbols_to_expressions=symbols_to_expressions,
        R=R,
        points=np.concatenate([R_at_int_points_, R_at_nodes_], dtype=float),
        point_names=[*int_point_names, *node_names]
    )

    # write element info
    compiler.write_element_info(
        element_type=element_type,
        d=d_,
        n=n_,
        ips=ips_
    )

    # F
    if write_F:
        assignments = term_collector.collect_assignments(
            input_symbols=[X_at_nodes, U_at_nodes],
            result_symbol=F_at_int_points,
            cse=True
        )

        compiler.write_function_for_F(
            assignments=assignments,
            element_typ=element_type,
            X_at_nodes=X_at_nodes,
            U_at_nodes=U_at_nodes,
            ips=ips_,
            F_at_int_points=F_at_int_points
        )

    # P
    if write_P:
        assignments = term_collector.collect_assignments(
            input_symbols=[X_at_nodes, U_at_nodes, S_at_int_points],
            result_symbol=P_at_int_points,
            cse=True
        )

        compiler.write_function_for_P(
            assignments=assignments,
            element_typ=element_type,
            X_at_nodes=X_at_nodes,
            U_at_nodes=U_at_nodes,
            S_at_int_points=S_at_int_points,
            P_at_int_points=P_at_int_points
        )

    # CS
    if write_CS:
        assignments = term_collector.collect_assignments(
            input_symbols=[e_at_int_points, X_at_nodes, U_at_nodes, S_at_int_points],
            result_symbol=CS_at_int_points,
            cse=True
        )

        compiler.write_function_for_CS(
            assignments=assignments,
            element_typ=element_type,
            e=e_at_int_points,
            X_at_nodes=X_at_nodes,
            U_at_nodes=U_at_nodes,
            S_at_int_points=S_at_int_points,
            CS_at_int_points=CS_at_int_points,
            is_dbf=is_dbf
        )

    # CF
    if write_CF:
        assignments = term_collector.collect_assignments(
            input_symbols=[e_at_int_points, X_at_nodes, U_at_nodes, S_at_int_points],
            result_symbol=CF_at_nodes,
            cse=True
        )

        compiler.write_function_for_CF(
            assignments=assignments,
            element_typ=element_type,
            e=e_at_int_points,
            X_at_nodes=X_at_nodes,
            U_at_nodes=U_at_nodes,
            S_at_int_points=S_at_int_points,
            CF_at_nodes=CF_at_nodes,
            is_dbf=is_dbf
        )


class CCodePrinter(C99CodePrinter):
    def _print_MatrixElement(self, expr: sy.Expr):
        return f"{self._print(str(expr.args[0]))}[{']['.join([str(i) for i in expr.indices])}]"

    def _print_Indexed(self, expr: sy.Indexed):
        return self._print_MatrixElement(expr)


class CPyCodeCompiler(object):
    def __init__(self, name: str, compile_at_exit: bool, write_header_at_enter: bool, folder: str = None):
        """
        Generate the following files:

         - Python file containing a binding to the generated library
         - C-file containing the source code for the library
         - The Library is the compiled C-file with the file extension \\*.dll on Windows or \\*.so on Linux

        **Examples**

        The following code generates

         - test_library.py,
         - _test_library.c,
         - _test_library.dll or _test_library.so

        in the folder `res/tests/codegen` and writes an element information
        for the element type `CPE4`

        >>> os.makedirs("res/tests/codegen", exist_ok=True)
        >>> with CPyCodeCompiler(
        ...         "test_library",
        ...         compile_at_exit=True,
        ...         write_header_at_enter=True,
        ...         folder="res/tests/codegen"
        ... ) as compiler:
        ...     compiler.write_element_info("CPE4", 2, 4, 4)  # doctest:+SKIP

        The generated package can be imported.

        >>> from res.tests.codegen import test_library
        >>> test_library.map_type_to_info["CPE4"]
        ElementInfo(number_of_dimensions=2, number_of_nodes=4, number_of_integration_points=4)

        :param name: File name (without extension) of the generated files
        :param compile_at_exit: Compile C-code to a library after the end of the with statement
        :param write_header_at_enter: Import and declare all necessary packages and attributes at the
            top of the Python and C-file.
        :param folder: Folder in which the generated code and the compiled library is located
        """
        self._name = str(name)
        self._compile_at_exit = bool(compile_at_exit)
        self._write_header_at_entry = bool(write_header_at_enter)
        if folder is None:
            folder = os.path.abspath(os.path.join(__file__, os.pardir))
        self._folder = os.path.abspath(folder)

        assert os.path.isdir(self._folder)
        assert os.path.basename(self._name) == self._name

        self._c_printer = CCodePrinter()
        self._py_file_handle = None  # type: Optional[TextIO]
        self._c_file_handle = None  # type: Optional[TextIO]

    def __enter__(self):
        self._py_file_handle = open(os.path.join(self._folder, self._name + ".py"), "w", encoding="utf-8")
        self._c_file_handle = open(os.path.join(self._folder, "_" + self._name + ".c"), "w", encoding="utf-8")

        if self._write_header_at_entry:
            self.write_headers()

        return self

    def __exit__(self, *_):
        self._py_file_handle.close()
        self._c_file_handle.close()

        if self._compile_at_exit:
            self._compile_c_code()

    def _compile_c_code(self):
        working_directory = os.path.abspath(".")
        try:
            os.chdir(self._folder)

            if os.name == 'nt':
                command = [
                    "gcc", "-shared", "-fPIC", "-std=c99",
                    "-o", f"_{self._name}.dll", f"_{self._name}.c"
                ]
            else:
                command = [
                    "gcc", "-shared", "-fPIC", "-std=c99",
                    "-o", f"_{self._name}.so", f"_{self._name}.c"
                ]

            print(f"{self._folder}>{command}")
            subprocess.call(command, shell=True)
        finally:
            os.chdir(working_directory)

    def _write_assignments(self, assignments):
        for assignment in assignments:
            lhs = assignment.lhs
            rhs = assignment.rhs

            lhs_code = self._c_printer.doprint(lhs)
            rhs_code = self._c_printer.doprint(rhs)
            if "[" in lhs_code:
                self._c_file_handle.write(f"    {lhs_code} = {rhs_code};\n")
            else:
                self._c_file_handle.write(f"    double {lhs_code} = {rhs_code};\n")

    def write_headers(self):
        """
        Write all necessary declarations and import statements
        into the C- and Python-file.

        .. note::

            This function must not be called more than once.

        """
        c_header_file_name = os.path.join(__file__, os.pardir, "_codegen_c_header.c")
        with open(c_header_file_name, "r", encoding="utf-8") as fh:
            self._c_file_handle.write(
                fh.read()
                .replace(
                    " * DO NOT INCLUDE OR RUN THIS FILE. THIS IS JUST A TEMPLATE USED BY THE CODE GENERATION!\n",
                    "")
            )

        python_header_file_name = os.path.join(__file__, os.pardir, "_codegen_python_header.py")
        with open(python_header_file_name, "r", encoding="utf-8") as fh:
            self._py_file_handle.write(
                fh.read()
                .replace("REPLACE_THIS_BY_LIBRARY_FILE_NAME", "_" + self._name)
                .replace("# DO NOT IMPORT OR RUN THIS FILE. "
                         "THIS IS JUST A TEMPLATE USED BY THE CODE GENERATION!\n", "")
            )

    def write_element_info(
            self,
            element_type: str,
            d: int,
            n: int,
            ips: int
    ):
        """
        Put information about an element into a dictionary in the generated Python file.

        :param element_type: name of the element type
        :param d: number of dimensions
        :param n: number of nodes
        :param ips: number of integration points
        """
        self._py_file_handle.write(f'''
map_type_to_info['{element_type}'] = ElementInfo(
    number_of_dimensions={d},
    number_of_nodes={n},
    number_of_integration_points={ips}
)
''')

    def write_function_for_F(
            self,
            assignments: List[ast.Assignment],
            element_typ: str,
            X_at_nodes: sy.IndexedBase,
            U_at_nodes: sy.IndexedBase,
            ips: int,
            F_at_int_points: sy.IndexedBase
    ):
        """
        Write a C-function that computes the deformation gradient
        and write a Python-function that binds to the C-function.

        :param assignments: list of assignments.
            Components of X_at_nodes, U_at_nodes may only occur on the right hand side.
            Each component of F_at_int_points must occur on the left hand side of one assignment.
        :param element_typ: name of the element type
        :param X_at_nodes: name of the matrix that contains coordinates at nodes
        :param U_at_nodes: name of the matrix that contains displacements at nodes
        :param ips: number of integration points
        :param F_at_int_points: name of the matrix that contains the deformation gradient at the integration points
        """
        n_, d_ = X_at_nodes.shape
        ips_ = int(ips)

        function_name_one_element = f"compute_F_for_one_element_of_type_{element_typ}"
        function_name_n_elements = f"compute_F_for_{element_typ}"

        # c code for one element
        self._c_file_handle.write(f'''
/**
 * Computes the deformation gradients for one element of typ {element_typ}.
 * The element has n={n_} nodes and ips={ips_} integration points.
 * 
 * @param[in] {X_at_nodes} Coordinates at n nodes of the element.
 * @param[in] {U_at_nodes} Displacements at n nodes of the element.
 * @param[out] {F_at_int_points} Empty array into which the function write the deformation gradients.
 *             The deformation gradients are evaluated on ips integration points.
 */
void {function_name_one_element}(
        double {X_at_nodes}[{n_}][{d_}], 
        double {U_at_nodes}[{n_}][{d_}], 
        double {F_at_int_points}[{ips_}][{d_}][{d_}]
) {{
    // COMPUTER GENERATED CODE:
''')
        self._write_assignments(assignments)
        self._c_file_handle.write(f'''
}}
''')

        # c code for multiple elements
        self._c_file_handle.write(f'''
/**
 * Compute the deformation gradients for num_elem elements of typ {element_typ}.
 * Each element has n={n_} nodes and ips={ips_} integration points.
 * 
 * @param[in] num_elem number of elements
 * @param[in] {X_at_nodes} Array of shape (num_elem, n, {d_}) 
 *            containing the coordinates at n nodes of num_elem elements.
 * @param[in] {U_at_nodes} Array of shape (num_elem, n, {d_}) 
 *            containing the displacements at n nodes of num_elem elements.
 * @param[out] {F_at_int_points} Empty array of shape (num_elem, ips, {d_}, {d_}).
 *             The function writes the deformation gradients into this array.
 *             The deformation gradients are evaluated on ips integration points for num_elem elements.
 */
CF_API void {function_name_n_elements}(
        int num_elem, 
        double {X_at_nodes}[][{n_}][{d_}], 
        double {U_at_nodes}[][{n_}][{d_}], 
        double {F_at_int_points}[][{ips_}][{d_}][{d_}] 
) {{
    #pragma omp parallel for
    for (size_t i=0; i<num_elem; i++) {{
        {function_name_one_element}(
            {X_at_nodes}[i], 
            {U_at_nodes}[i], 
            {F_at_int_points}[i]
        );
    }}
}}
''')

        # python code for multiple elements
        self._py_file_handle.write(f'''

def {function_name_n_elements}(
        {X_at_nodes},
        {U_at_nodes}):
    """
    Compute the deformation gradients for num_elem elements of typ {element_typ}.
    Each element has n={n_} nodes and ips={ips_} integration points.

    :param {X_at_nodes}: Array of shape (num_elem, n, {d_}) containing the coordinates at n nodes of num_elem elements.
    :param {U_at_nodes}: Array of shape (num_elem, n, {d_}) containing the displacements at n nodes of num_elem elements.
    :return: {F_at_int_points}: Array of shape (num_elem, ips, {d_}, {d_}) containing the deformation gradients
        evaluated on ips integration points for num_elem element.
    """

    {X_at_nodes} = np.ascontiguousarray({X_at_nodes}, dtype=np.float64)
    {U_at_nodes} = np.ascontiguousarray({U_at_nodes}, dtype=np.float64)

    num_elem = {X_at_nodes}.shape[0]

    {F_at_int_points} = np.zeros((num_elem, {ips_}, {d_}, {d_}), dtype=np.float64)

    assert {X_at_nodes}.shape == (num_elem, {n_}, {d_})
    assert {U_at_nodes}.shape == (num_elem, {n_}, {d_})

    lib.{function_name_n_elements}(
        ctypes.c_size_t(num_elem),
        {X_at_nodes}.ctypes.data_as(ctypes.c_void_p),
        {U_at_nodes}.ctypes.data_as(ctypes.c_void_p),
        {F_at_int_points}.ctypes.data_as(ctypes.c_void_p)
    )

    return {F_at_int_points}


map_type_to_F_function['{element_typ}'] = {function_name_n_elements}
''')

    def write_function_for_P(
            self,
            assignments: List[ast.Assignment],
            element_typ: str,
            X_at_nodes: sy.IndexedBase,
            U_at_nodes: sy.IndexedBase,
            S_at_int_points: sy.IndexedBase,
            P_at_int_points: sy.IndexedBase
    ):
        """
        Write a C-function that computes the first Piola-Kirchhoff stress tensor
        and write a Python-function that binds to the C-function.

        :param assignments: list of assignments.
            Components of X_at_nodes, U_at_nodes, S_at_int_points may only occur on the right hand side.
            Each component of P_at_int_points must occur on the left hand side of one assignment.
        :param element_typ: name of the element type
        :param X_at_nodes: name of the matrix that contains coordinates at nodes
        :param U_at_nodes: name of the matrix that contains displacements at nodes
        :param S_at_int_points: name of the matrix that contains stress tensors at the integration points
        :param P_at_int_points: name of the matrix that contains the
            first Piola-Kirchhoff stress tensors at the integration points
        """
        n_, d_ = X_at_nodes.shape
        ips_ = S_at_int_points.shape[0]

        function_name_one_element = f"compute_P_for_one_element_of_type_{element_typ}"
        function_name_n_elements = f"compute_P_for_{element_typ}"

        # c code for one element
        self._c_file_handle.write(f'''
/**
 * Compute the first Piola-Kirchhoff stress tensors for one element of typ {element_typ}.
 * The element has n={n_} nodes and ips={ips_} integration points.
 * 
 * @param[in] {X_at_nodes} Coordinates at n nodes of the element.
 * @param[in] {U_at_nodes} Displacements at n nodes of the element.
 * @param[in] {S_at_int_points} Symmetric stress tensors at ips integration points.
 * @param[out] {P_at_int_points} Empty array into which the 1. Piola-Kirchhoff stress tensors are written.
 *             The 1. Piola-Kirchhoff stress tensors are evaluated on ips integration points.
 */
void {function_name_one_element}(
        double {X_at_nodes}[{n_}][{d_}], 
        double {U_at_nodes}[{n_}][{d_}], 
        double {S_at_int_points}[{ips_}][{d_}][{d_}], 
        double {P_at_int_points}[{ips_}][{d_}][{d_}]
) {{
    // COMPUTER GENERATED CODE:
''')
        self._write_assignments(assignments)
        self._c_file_handle.write(f'''
}}
''')

        # c code for multiple elements
        self._c_file_handle.write(f'''
/**
 * Compute the first Piola-Kirchhoff stress tensors for num_elem elements of typ {element_typ}.
 * Each element has n={n_} nodes and ips={ips_} integration points.
 * 
 * @param[in] num_elem number of elements
 * @param[in] {X_at_nodes} Array of shape (num_elem, n, {d_}) 
 *            containing the coordinates at n nodes of num_elem elements.
 * @param[in] {U_at_nodes} Array of shape (num_elem, n, {d_}) 
 *            containing the displacements at n nodes of num_elem elements.
 * @param[in] {S_at_int_points} Array of shape (num_elem, ips, {d_}, {d_}) 
 *            containing the symmetric stress tensors at ips integration points for num_elem elements.
 * @param[out] {P_at_int_points} Empty array of shape (num_elem, ips, {d_}, {d_}).
 *             The function writes the 1. Piola-Kirchhoff stress tensors into this array.
 *             The 1. Piola-Kirchhoff stress tensors are evaluated on ips integration points for num_elem elements.
 */
CF_API void {function_name_n_elements}(
        int num_elem, 
        double {X_at_nodes}[][{n_}][{d_}], 
        double {U_at_nodes}[][{n_}][{d_}], 
        double {S_at_int_points}[][{ips_}][{d_}][{d_}], 
        double {P_at_int_points}[][{ips_}][{d_}][{d_}] 
) {{
    // COMPUTER GENERATED CODE:
    #pragma omp parallel for
    for (size_t i=0; i<num_elem; i++) {{
        {function_name_one_element}(
            {X_at_nodes}[i], 
            {U_at_nodes}[i], 
            {S_at_int_points}[i], 
            {P_at_int_points}[i]
        );
    }}
}}
''')

        # python code for multiple elements
        self._py_file_handle.write(f'''

def {function_name_n_elements}(
        {X_at_nodes},
        {U_at_nodes},
        {S_at_int_points}):
    """
    Compute the first Piola-Kirchhoff stress tensors for num_elem elements of typ {element_typ}.
    Each element has n={n_} nodes and ips={ips_} integration points.

    :param {X_at_nodes}: Array of shape (num_elem, n, {d_}) containing the coordinates at n nodes of num_elem elements.
    :param {U_at_nodes}: Array of shape (num_elem, n, {d_}) containing the displacements at n nodes of num_elem elements.
    :param {S_at_int_points}: Array of shape (num_elem, ips, {d_}, {d_})
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: {P_at_int_points}: Array of shape (num_elem, ips, {d_}, {d_}) containing the 1. Piola-Kirchhoff stress tensors
        evaluated on ips integration points for num_elem element.
    """

    {X_at_nodes} = np.ascontiguousarray({X_at_nodes}, dtype=np.float64)
    {U_at_nodes} = np.ascontiguousarray({U_at_nodes}, dtype=np.float64)
    {S_at_int_points} = np.ascontiguousarray({S_at_int_points}, dtype=np.float64)

    num_elem = {X_at_nodes}.shape[0]

    {P_at_int_points} = np.zeros((num_elem, {ips_}, {d_}, {d_}), dtype=np.float64)

    assert {X_at_nodes}.shape == (num_elem, {n_}, {d_})
    assert {U_at_nodes}.shape == (num_elem, {n_}, {d_})
    assert {S_at_int_points}.shape == (num_elem, {ips_}, {d_}, {d_})

    lib.{function_name_n_elements}(
        ctypes.c_size_t(num_elem),
        {X_at_nodes}.ctypes.data_as(ctypes.c_void_p),
        {U_at_nodes}.ctypes.data_as(ctypes.c_void_p),
        {S_at_int_points}.ctypes.data_as(ctypes.c_void_p),
        {P_at_int_points}.ctypes.data_as(ctypes.c_void_p)
    )

    return {P_at_int_points}


map_type_to_P_function['{element_typ}'] = {function_name_n_elements}
''')

    def write_function_for_CS(
            self,
            assignments: List[ast.Assignment],
            element_typ: str,
            e: sy.Symbol,
            X_at_nodes: sy.IndexedBase,
            U_at_nodes: sy.IndexedBase,
            S_at_int_points: sy.IndexedBase,
            CS_at_int_points: sy.IndexedBase,
            is_dbf: bool
    ):
        """
        Write a C-function that computes the configurational stresses
        and write a Python-function that binds to the C-function.

        :param assignments: list of assignments.
            Components of e, X_at_nodes, U_at_nodes, S_at_int_points may only occur on the right hand side.
            Each component of CS_at_int_points must occur on the left hand side of one assignment.
        :param element_typ: name of the element type
        :param e: name of the energy density
        :param X_at_nodes: name of the matrix that contains coordinates at nodes
        :param U_at_nodes: name of the matrix that contains displacements at nodes
        :param S_at_int_points: name of the matrix that contains stress tensors at the integration points
        :param CS_at_int_points: name of the matrix that contains the
            configurational stress tensors at the integration points
        :param is_dbf: displacement based or motion based formulation
        """
        n_, d_ = X_at_nodes.shape
        ips_ = S_at_int_points.shape[0]
        method = 'dbf' if is_dbf else 'mbf'
        method_name = f"{'deformation' if is_dbf else 'motion'} based"

        function_name_one_element = f"compute_CS_for_one_element_of_type_{element_typ}_using_{method}"
        function_name_n_elements = f"compute_CS_for_{element_typ}_using_{method}"

        # c code for one element
        self._c_file_handle.write(f'''
/**
 * Compute the configurational stresses for one element of typ {element_typ}.
 * The element has n={n_} nodes and ips={ips_} integration points.
 * The computation is {method_name}.
 * 
 * @param[in] {e} Helmholtz free energy density
 * @param[in] {X_at_nodes} Coordinates at n nodes of the element.
 * @param[in] {U_at_nodes} Displacements at n nodes of the element.
 * @param[in] {S_at_int_points} Symmetric stress tensors at ips integration points.
 * @param[out] {CS_at_int_points} Empty array into which the configurational stresses are written.
 *             The configurational stresses are evaluated on ips integration points.
 */
void {function_name_one_element}(
        double {e}[{ips_}], 
        double {X_at_nodes}[{n_}][{d_}], 
        double {U_at_nodes}[{n_}][{d_}], 
        double {S_at_int_points}[{ips_}][{d_}][{d_}], 
        double {CS_at_int_points}[{ips_}][{d_}][{d_}]
) {{
    // COMPUTER GENERATED CODE:
''')
        self._write_assignments(assignments)
        self._c_file_handle.write(f'''
}}\n''')

        # c code for multiple elements
        self._c_file_handle.write(f'''
/**
 * Compute the configurational stresses for num_elem elements of typ {element_typ}.
 * Each element has n={n_} nodes and ips={ips_} integration points.
 * The computation is {method_name}.
 * 
 * @param[in] num_elem number of elements
 * @param[in] {e} Array of shape (num_elem, ips) containing the Helmholtz free energy densities
 * @param[in] {X_at_nodes} Array of shape (num_elem, n, {d_})
 *            containing the coordinates at n nodes of num_elem elements.
 * @param[in] {U_at_nodes} Array of shape (num_elem, n, {d_})
 *            containing the displacements at n nodes of num_elem elements.
 * @param[in] {S_at_int_points} Array of shape (num_elem, ips, {d_}, {d_})
 *            containing the symmetric stress tensors at ips integration points for num_elem elements.
 * @param[out] {CS_at_int_points} Empty array of shape (num_elem, ips, {d_}, {d_}).
 *             The function writes the configurational stresses into this array.
 *             The configurational stresses are evaluated on ips integration points for num_elem elements.
 */
CF_API void {function_name_n_elements}(
        int num_elem, 
        double {e}[][{ips_}], 
        double {X_at_nodes}[][{n_}][{d_}], 
        double {U_at_nodes}[][{n_}][{d_}], 
        double {S_at_int_points}[][{ips_}][{d_}][{d_}], 
        double {CS_at_int_points}[][{ips_}][{d_}][{d_}] 
) {{
    // COMPUTER GENERATED CODE:
    #pragma omp parallel for
    for (size_t i=0; i<num_elem; i++) {{
        {function_name_one_element}(
            {e}[i], 
            {X_at_nodes}[i], 
            {U_at_nodes}[i], 
            {S_at_int_points}[i], 
            {CS_at_int_points}[i]
        );
    }}
}}
''')

        # python code for multiple elements
        self._py_file_handle.write(f'''

def {function_name_n_elements}(
        {e},
        {X_at_nodes},
        {U_at_nodes},
        {S_at_int_points}):
    """
    Compute the configurational stresses for num_elem elements of typ {element_typ}.
    Each element has n={n_} nodes and ips={ips_} integration points.

    :param {e}: Array of shape (num_elem, ips) containing the Helmholtz free energy densities of num_elem elements.
    :param {X_at_nodes}: Array of shape (num_elem, n, {d_}) containing the coordinates at n nodes of num_elem elements.
    :param {U_at_nodes}: Array of shape (num_elem, n, {d_}) containing the displacements at n nodes of num_elem elements.
    :param {S_at_int_points}: Array of shape (num_elem, ips, {d_}, {d_})
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: {CS_at_int_points}: Array of shape (num_elem, ips, {d_}, {d_}) containing the configurational stresses
        evaluated on ips integration points for num_elem element.
    """

    {e} = np.ascontiguousarray({e}, dtype=np.float64)
    {X_at_nodes} = np.ascontiguousarray({X_at_nodes}, dtype=np.float64)
    {U_at_nodes} = np.ascontiguousarray({U_at_nodes}, dtype=np.float64)
    {S_at_int_points} = np.ascontiguousarray({S_at_int_points}, dtype=np.float64)

    num_elem = {X_at_nodes}.shape[0]

    {CS_at_int_points} = np.zeros((num_elem, {ips_}, {d_}, {d_}), dtype=np.float64)

    assert {e}.shape == (num_elem, {ips_})
    assert {X_at_nodes}.shape == (num_elem, {n_}, {d_})
    assert {U_at_nodes}.shape == (num_elem, {n_}, {d_})
    assert {S_at_int_points}.shape == (num_elem, {ips_}, {d_}, {d_})

    lib.{function_name_n_elements}(
        ctypes.c_size_t(num_elem),
        {e}.ctypes.data_as(ctypes.c_void_p),
        {X_at_nodes}.ctypes.data_as(ctypes.c_void_p),
        {U_at_nodes}.ctypes.data_as(ctypes.c_void_p),
        {S_at_int_points}.ctypes.data_as(ctypes.c_void_p),
        {CS_at_int_points}.ctypes.data_as(ctypes.c_void_p)
    )

    return {CS_at_int_points}


map_type_and_method_to_CS_function[('{element_typ}', '{method}')] = {function_name_n_elements}
''')

    def write_function_for_CF(
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
        """
        Write a C-function that computes the configurational forces
        and write a Python-function that binds to the C-function.

        :param assignments: list of assignments.
            Components of e, X_at_nodes, U_at_nodes, S_at_int_points may only occur on the right hand side.
            Each component of CF_at_nodes must occur on the left hand side of one assignment.
        :param element_typ: name of the element type
        :param e: name of the energy density
        :param X_at_nodes: name of the matrix that contains coordinates at nodes
        :param U_at_nodes: name of the matrix that contains displacements at nodes
        :param S_at_int_points: name of the matrix that contains stress tensors at the integration points
        :param CF_at_nodes: name of the matrix that contains the
            configurational forces at the element nodes
        :param is_dbf: displacement based or motion based formulation
        """
        n_, d_ = X_at_nodes.shape
        ips_ = S_at_int_points.shape[0]
        method = 'dbf' if is_dbf else 'mbf'
        method_name = f"{'deformation' if is_dbf else 'motion'} based"

        function_name_one_element = f"compute_CF_for_one_element_of_type_{element_typ}_using_{method}"
        function_name_n_elements = f"compute_CF_for_{element_typ}_using_{method}"

        # c code for one element
        self._c_file_handle.write(f'''
/**
 * Compute the configurational forces for one element of typ {element_typ}.
 * The element has n={n_} nodes and ips={ips_} integration points.
 * The computation is {method_name}.
 * 
 * @param[in] {e} Helmholtz free energy density
 * @param[in] {X_at_nodes} Coordinates at n nodes of the element.
 * @param[in] {U_at_nodes} Displacements at n nodes of the element.
 * @param[in] {S_at_int_points} Symmetric stress tensors at ips integration points.
 * @param[out] {CF_at_nodes} Empty array into which the configurational forces are written.
 *             The configurational forces are evaluated on n nodes of the element.
 */
void {function_name_one_element}(
        double {e}[{ips_}], 
        double {X_at_nodes}[{n_}][{d_}], 
        double {U_at_nodes}[{n_}][{d_}], 
        double {S_at_int_points}[{ips_}][{d_}][{d_}], 
        double {CF_at_nodes}[{n_}][{d_}]
) {{
    // COMPUTER GENERATED CODE:
''')
        self._write_assignments(assignments)
        self._c_file_handle.write(f'''
}}\n''')

        # c code for multiple elements
        self._c_file_handle.write(f'''
/**
 * Compute the configurational forces for num_elem elements of typ {element_typ}.
 * Each element has n={n_} nodes and ips={ips_} integration points.
 * The computation is {method_name}.
 * 
 * @param[in] num_elem number of elements
 * @param[in] {e} Array of shape (num_elem, ips) containing the Helmholtz free energy densities
 * @param[in] {X_at_nodes} Array of shape (num_elem, n, {d_}) 
 *            containing the coordinates at n nodes of num_elem elements.
 * @param[in] {U_at_nodes} Array of shape (num_elem, n, {d_}) 
 *            containing the displacements at n nodes of num_elem elements.
 * @param[in] {S_at_int_points} Array of shape (num_elem, ips, {d_}, {d_}) 
 *            containing the symmetric stress tensors at ips integration points for num_elem elements.
 * @param[out] {CF_at_nodes} Empty array of shape (num_elem, n, {d_}).
 *             The function writes the configurational forces into this array.
 *             The configurational forces are evaluated on n nodes for num_elem elements.
 */
CF_API void {function_name_n_elements}(
        int num_elem, 
        double {e}[][{ips_}], 
        double {X_at_nodes}[][{n_}][{d_}], 
        double {U_at_nodes}[][{n_}][{d_}], 
        double {S_at_int_points}[][{ips_}][{d_}][{d_}], 
        double {CF_at_nodes}[][{n_}][{d_}] 
) {{
    // COMPUTER GENERATED CODE:
    #pragma omp parallel for
    for (size_t i=0; i<num_elem; i++) {{
        {function_name_one_element}(
            {e}[i], 
            {X_at_nodes}[i], 
            {U_at_nodes}[i], 
            {S_at_int_points}[i], 
            {CF_at_nodes}[i]
        );
    }}
}}
''')

        # python code for multiple elements
        self._py_file_handle.write(f'''

def {function_name_n_elements}(
        {e},
        {X_at_nodes},
        {U_at_nodes},
        {S_at_int_points}):
    """
    Compute the configurational forces for num_elem elements of typ {element_typ}.
    Each element has n={n_} nodes and ips={ips_} integration points.

    :param {e}: Array of shape (num_elem, ips) containing the Helmholtz free energy densities of num_elem elements.
    :param {X_at_nodes}: Array of shape (num_elem, n, {d_}) containing the coordinates at n nodes of num_elem elements.
    :param {U_at_nodes}: Array of shape (num_elem, n, {d_}) containing the displacements at n nodes of num_elem elements.
    :param {S_at_int_points}: Array of shape (num_elem, ips, {d_}, {d_})
        containing the symmetric stress tensors at ips integration points for num_elem elements.
    :return: {CF_at_nodes}: Array of shape (num_elem, n, {d_}) containing the configurational forces
        evaluated on n nodes for num_elem element.
    """

    {e} = np.ascontiguousarray({e}, dtype=np.float64)
    {X_at_nodes} = np.ascontiguousarray({X_at_nodes}, dtype=np.float64)
    {U_at_nodes} = np.ascontiguousarray({U_at_nodes}, dtype=np.float64)
    {S_at_int_points} = np.ascontiguousarray({S_at_int_points}, dtype=np.float64)

    num_elem = {X_at_nodes}.shape[0]

    {CF_at_nodes} = np.zeros((num_elem, {n_}, {d_}), dtype=np.float64)

    assert {e}.shape == (num_elem, {ips_})
    assert {X_at_nodes}.shape == (num_elem, {n_}, {d_})
    assert {U_at_nodes}.shape == (num_elem, {n_}, {d_})
    assert {S_at_int_points}.shape == (num_elem, {ips_}, {d_}, {d_})

    lib.{function_name_n_elements}(
        ctypes.c_size_t(num_elem),
        {e}.ctypes.data_as(ctypes.c_void_p),
        {X_at_nodes}.ctypes.data_as(ctypes.c_void_p),
        {U_at_nodes}.ctypes.data_as(ctypes.c_void_p),
        {S_at_int_points}.ctypes.data_as(ctypes.c_void_p),
        {CF_at_nodes}.ctypes.data_as(ctypes.c_void_p)
    )

    return {CF_at_nodes}


map_type_and_method_to_CF_function[('{element_typ}', '{method}')] = {function_name_n_elements}
''')


if __name__ == '__main__':
    write_code_for_all_element_types()
