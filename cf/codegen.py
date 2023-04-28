from typing import List, TextIO, Optional
import subprocess
import os
from datetime import date
from itertools import product

import numpy as np
import sympy as sy

from sympy.codegen import ast
from sympy.printing.c import C99CodePrinter

from cf import element_definitions
from cf.expressions import compute_CF, TermCollector
import cf.expressions as expr


def write_code_for_all_element_types(*types):
    if len(types) == 0:
        types = element_definitions.R_at_nodes_of_element.keys()

    # generate and compile c and python code
    with CPyCodeCompiler("cf_c", compile_at_exit=True, write_header_at_enter=True) as compiler:
        for element_type in types:
            for is_dbf in [True, False]:
                print(f"element_type={element_type}, is_dbf={is_dbf}")
                write_code_for_element_type(
                    element_type=element_type,
                    is_dbf=is_dbf,
                    compiler=compiler
                )


def write_code_for_element_type(element_type: str, is_dbf: bool, compiler):
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

    # compute CF
    R, CF_at_nodes, symbols_to_expressions = compute_CF(
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

    #
    int_repl_rules = expr.create_replacement_rules(
        expr.eval_R(d_),
        R_at_int_points_
    )

    # F
    F_ = expr.create_symbolic_matrix(
        "F{row}{col}",
        ["x", "y", "z"][:d_],
        ["x", "y", "z"][:d_],
        *expr.eval_R(d_)
    )
    F_at_int_points_ = expr.apply_replacement_rules(
        F_,
        int_repl_rules
    )
    F_at_int_points = sy.IndexedBase("F_at_int_points", shape=(ips_, d_, d_))
    symbols_to_expressions.update({
        F_at_int_points[idx]: F_at_int_points_[idx]
        for idx in product(*[range(int(dim)) for dim in F_at_int_points.shape])
    })

    # P
    P_ = expr.create_symbolic_matrix(
        "P{row}{col}",
        ["x", "y", "z"][:d_],
        ["x", "y", "z"][:d_],
        *expr.eval_R(d_)
    )
    P_at_int_points_ = expr.apply_replacement_rules(
        P_,
        int_repl_rules
    )
    P_at_int_points = sy.IndexedBase("P_at_int_points", shape=(ips_, d_, d_))
    symbols_to_expressions.update({
        P_at_int_points[idx]: P_at_int_points_[idx]
        for idx in product(*[range(int(dim)) for dim in P_at_int_points.shape])
    })

    # CS
    CS_ = expr.create_symbolic_matrix(
        "CS{row}{col}",
        ["x", "y", "z"][:d_],
        ["x", "y", "z"][:d_],
        *expr.eval_R(d_)
    )
    CS_at_int_points_ = expr.apply_replacement_rules(
        CS_,
        int_repl_rules
    )
    CS_at_int_points = sy.IndexedBase("CS_at_int_points", shape=(ips_, d_, d_))
    symbols_to_expressions.update({
        CS_at_int_points[idx]: CS_at_int_points_[idx]
        for idx in product(*[range(int(dim)) for dim in CS_at_int_points.shape])
    })

    # create abstract code assignments
    term_collector = TermCollector(
        symbols_to_expressions=symbols_to_expressions,
        R=R,
        points=np.concatenate([R_at_int_points_, R_at_nodes_], dtype=float),
        point_names=[*int_point_names, *node_names]
    )

    # F
    if not is_dbf:
        assignments = term_collector.collect_assignments(
            input_symbols=[X_at_nodes, U_at_nodes],
            result_array=F_at_int_points,
            cse=True
        )

        compiler.write_function_for_F(
            assignments=assignments,
            element_typ=element_type,
            X_at_nodes=X_at_nodes,
            U_at_nodes=U_at_nodes,
            number_of_integration_points=ips_,
            F_at_int_points=F_at_int_points
        )

    # P
    if not is_dbf:
        assignments = term_collector.collect_assignments(
            input_symbols=[X_at_nodes, U_at_nodes, S_at_int_points],
            result_array=P_at_int_points,
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
    assignments = term_collector.collect_assignments(
        input_symbols=[e_at_int_points, X_at_nodes, U_at_nodes, S_at_int_points],
        result_array=CS_at_int_points,
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
    assignments = term_collector.collect_assignments(
        input_symbols=[e_at_int_points, X_at_nodes, U_at_nodes, S_at_int_points],
        result_array=CF_at_nodes,
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
        self._name = str(name)
        self._compile_at_exit = bool(compile_at_exit)
        self._write_header_at_entry = bool(write_header_at_enter)
        if folder is None:
            folder = os.path.abspath(os.path.join(__file__, os.pardir))
        self._folder = os.path.abspath(folder)

        assert os.path.isdir(self._folder)
        assert os.path.basename(self._name) == self._name

        self._py_file_handle = None  # type: Optional[TextIO]
        self._c_file_handle = None  # type: Optional[TextIO]

    def __enter__(self):
        self._py_file_handle = open(os.path.join(self._folder, self._name + ".py"), "w")
        self._c_file_handle = open(os.path.join(self._folder, self._name + ".c"), "w")

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
            command = f"gcc -shared -O3 -o {self._name}.dll {self._name}.c"
            print(f"{self._folder}>{command}")
            subprocess.call(command, shell=True)
        finally:
            os.chdir(working_directory)

    def write_headers(self):
        self._c_file_handle.write(f'''// generated by cf on {date.today().strftime("%d.%m.%Y")} (dd.mm.yyyy)
#ifdef _WIN32
#    define CF_API __declspec(dllexport)
#else
#    define CF_API
#endif

#include <math.h>  // contains power
#include <stddef.h>  // contains size_t
''')

        self._py_file_handle.write(f'''"""
# TODO: doc
"""
import numpy as np
import ctypes
import os

# load c library
if os.name == 'nt':
    # windows
    lib = ctypes.cdll.LoadLibrary(os.path.abspath(os.path.join(
        __file__, os.path.pardir, '{self._name}.dll' 
    )))
else:
    # linux
    lib = ctypes.cdll.LoadLibrary(os.path.abspath(os.path.join(
        __file__, os.path.pardir, '{self._name}.so'
    )))

# function lookup dictionaries
_map_typ_to_F_function = dict()
_map_typ_to_P_function = dict()
_map_typ_and_method_to_CS_function = dict()
_map_typ_and_method_to_CF_function = dict()


def compute_F(
        X_at_nodes,
        U_at_nodes,
        element_type
):
    fun = _map_typ_to_F_function[str(element_type)]
    return fun(
        X_at_nodes,
        U_at_nodes
    )
 

def compute_P(
        X_at_nodes,
        U_at_nodes,
        S_at_int_points,
        element_type
):
    fun = _map_typ_to_P_function[str(element_type)]
    return fun(
        X_at_nodes,
        U_at_nodes,
        S_at_int_points
    )


def compute_CS(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points,
        element_type,
        method
):
    fun = _map_typ_and_method_to_CS_function[(str(element_type), str(method))]
    return fun(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points
    )


def compute_CF(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points,
        element_type,
        method
):
    fun = _map_typ_and_method_to_CF_function[(str(element_type), str(method))]
    return fun(
        e_at_int_points,
        X_at_nodes,
        U_at_nodes,
        S_at_int_points
    )
''')

    def write_function_for_F(
            self,
            assignments: List[ast.Assignment],
            element_typ: str,
            X_at_nodes: sy.IndexedBase,
            U_at_nodes: sy.IndexedBase,
            number_of_integration_points: int,
            F_at_int_points: sy.IndexedBase
    ):

        c_printer = CCodePrinter()
        n_, d_ = X_at_nodes.shape
        ips_ = int(number_of_integration_points)

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
 * @param[out] {F_at_int_points} The resulting deformation gradients for ips integration points are written 
 *             to this array.
 */
void {function_name_one_element}(
        double {X_at_nodes}[{n_}][{d_}], 
        double {U_at_nodes}[{n_}][{d_}], 
        double {F_at_int_points}[{ips_}][{d_}][{d_}]
) {{
    // COMPUTER GENERATED CODE:
''')

        for assignment in assignments:
            lhs = assignment.lhs
            rhs = assignment.rhs

            lhs_code = c_printer.doprint(lhs)
            rhs_code = c_printer.doprint(rhs)
            if "[" in lhs_code:
                self._c_file_handle.write(f"    {lhs_code} = {rhs_code};\n")
            else:
                self._c_file_handle.write(f"    double {lhs_code} = {rhs_code};\n")

        self._c_file_handle.write("}\n")

        # c code for multiple elements
        self._c_file_handle.write(f'''
/**
 * Computes the deformation gradients for num_elem elements of typ {element_typ}.
 * Each element has n={n_} nodes and ips={ips_} integration points.
 * A static load case is assumed.
 * 
 * @param[in] num_elem number of elements
 * @param[in] {X_at_nodes} Coordinates at n nodes of num_elem elements.
 * @param[in] {U_at_nodes} Displacements at n nodes of num_elem elements.
 * @param[out] {F_at_int_points} The resulting deformation gradients for ips integration points are written 
 *             to this array.
 */
CF_API void {function_name_n_elements}(
        int num_elem, 
        double {X_at_nodes}[][{n_}][{d_}], 
        double {U_at_nodes}[][{n_}][{d_}], 
        double {F_at_int_points}[][{ips_}][{d_}][{d_}] 
) {{
    // COMPUTER GENERATED CODE:
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
        # TODO: doc
        self._py_file_handle.write(f"""

def {function_name_n_elements}(
        {X_at_nodes},
        {U_at_nodes}):

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


_map_typ_to_F_function['{element_typ}'] = {function_name_n_elements}
""")

    def write_function_for_P(
            self,
            assignments: List[ast.Assignment],
            element_typ: str,
            X_at_nodes: sy.IndexedBase,
            U_at_nodes: sy.IndexedBase,
            S_at_int_points: sy.IndexedBase,
            P_at_int_points: sy.IndexedBase
    ):
        c_printer = CCodePrinter()
        n_, d_ = X_at_nodes.shape
        ips_ = S_at_int_points.shape[0]

        function_name_one_element = f"compute_P_for_one_element_of_type_{element_typ}"
        function_name_n_elements = f"compute_P_for_{element_typ}"

        # c code for one element
        self._c_file_handle.write(f'''
/**
 * Computes the first Piola-Kirchhoff stress tensor for one element of typ {element_typ}.
 * The element has n={n_} nodes and ips={ips_} integration points.
 * 
 * @param[in] {X_at_nodes} Coordinates at n nodes of the element.
 * @param[in] {U_at_nodes} Displacements at n nodes of the element.
 * @param[in] {S_at_int_points} Symmetric stress tensors at ips integration points.
 * @param[out] {P_at_int_points} The resulting Piola-Kirchhoff stress tensors for ips integration points are written 
 *             to this array.
 */
void {function_name_one_element}(
        double {X_at_nodes}[{n_}][{d_}], 
        double {U_at_nodes}[{n_}][{d_}], 
        double {S_at_int_points}[{ips_}][{d_}][{d_}], 
        double {P_at_int_points}[{ips_}][{d_}][{d_}]
) {{
    // COMPUTER GENERATED CODE:
''')

        for assignment in assignments:
            lhs = assignment.lhs
            rhs = assignment.rhs

            lhs_code = c_printer.doprint(lhs)
            rhs_code = c_printer.doprint(rhs)
            if "[" in lhs_code:
                self._c_file_handle.write(f"    {lhs_code} = {rhs_code};\n")
            else:
                self._c_file_handle.write(f"    double {lhs_code} = {rhs_code};\n")

        self._c_file_handle.write("}\n")

        # c code for multiple elements
        self._c_file_handle.write(f'''
/**
 * Computes the first Piola-Kirchhoff stress tensor for num_elem elements of typ {element_typ}.
 * Each element has n={n_} nodes and ips={ips_} integration points.
 * A static load case is assumed.
 * 
 * @param[in] num_elem number of elements
 * @param[in] {X_at_nodes} Coordinates at n nodes of num_elem elements.
 * @param[in] {U_at_nodes} Displacements at n nodes of num_elem elements.
 * @param[in] {S_at_int_points} Symmetric stress tensors at ips integration points 
 *            for num_elem elements.
 * @param[out] {P_at_int_points} The resulting Piola-Kirchhoff stress tensors for ips integration points are written 
 *             to this array.
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
        # TODO: doc
        self._py_file_handle.write(f"""

def {function_name_n_elements}(
        {X_at_nodes},
        {U_at_nodes},
        {S_at_int_points}):

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


_map_typ_to_P_function['{element_typ}'] = {function_name_n_elements}
""")

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

        c_printer = CCodePrinter()
        n_, d_ = X_at_nodes.shape
        ips_ = S_at_int_points.shape[0]
        method = 'dbf' if is_dbf else 'mbf'
        method_name = f"{'displacement' if is_dbf else 'motion'} based"

        function_name_one_element = f"compute_CS_for_one_element_of_type_{element_typ}_using_{method}"
        function_name_n_elements = f"compute_CS_for_{element_typ}_using_{method}"

        # c code for one element
        self._c_file_handle.write(f'''
/**
 * Computes the configurational stresses for one element of typ {element_typ}.
 * The element has n={n_} nodes and ips={ips_} integration points.
 * The computation is {method_name}.
 * A static load case is assumed.
 * 
 * @param[in] {e} internal energy density
 * @param[in] {X_at_nodes} Coordinates at n nodes of the element.
 * @param[in] {U_at_nodes} Displacements at n nodes of the element.
 * @param[in] {S_at_int_points} Symmetric stress tensors at ips integration points.
 * @param[out] {CS_at_int_points} The resulting configurational stresses for ips integration points are written 
 *             to this array.
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

        for assignment in assignments:
            lhs = assignment.lhs
            rhs = assignment.rhs

            lhs_code = c_printer.doprint(lhs)
            rhs_code = c_printer.doprint(rhs)
            if "[" in lhs_code:
                self._c_file_handle.write(f"    {lhs_code} = {rhs_code};\n")
            else:
                self._c_file_handle.write(f"    double {lhs_code} = {rhs_code};\n")

        self._c_file_handle.write("}\n")

        # c code for multiple elements
        self._c_file_handle.write(f'''
/**
 * Computes the configurational stresses for num_elem elements of typ {element_typ}.
 * Each element has n={n_} nodes and ips={ips_} integration points.
 * The computation is {method_name}.
 * A static load case is assumed.
 * 
 * @param[in] num_elem number of elements
 * @param[in] {e} internal energy density
 * @param[in] {X_at_nodes} Coordinates at n nodes of num_elem elements.
 * @param[in] {U_at_nodes} Displacements at n nodes of num_elem elements.
 * @param[in] {S_at_int_points} Symmetric stress tensors at ips integration points 
 *            for num_elem elements.
 * @param[out] {CS_at_int_points} The resulting configurational stresses for ips integration points are written 
 *             to this array.
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
        # TODO: doc
        self._py_file_handle.write(f"""

def {function_name_n_elements}(
        {e},
        {X_at_nodes},
        {U_at_nodes},
        {S_at_int_points}):

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


_map_typ_and_method_to_CS_function[('{element_typ}', '{method}')] = {function_name_n_elements}
""")

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

        c_printer = CCodePrinter()
        n_, d_ = X_at_nodes.shape
        ips_ = S_at_int_points.shape[0]
        method = 'dbf' if is_dbf else 'mbf'
        method_name = f"{'displacement' if is_dbf else 'motion'} based"

        function_name_one_element = f"compute_CF_for_one_element_of_type_{element_typ}_using_{method}"
        function_name_n_elements = f"compute_CF_for_{element_typ}_using_{method}"

        # c code for one element
        self._c_file_handle.write(f'''
/**
 * Computes the configurational forces for one element of typ {element_typ}.
 * The element has n={n_} nodes and ips={ips_} integration points.
 * The computation is {method_name}.
 * A static load case is assumed.
 * 
 * @param[in] {e} internal energy density
 * @param[in] {X_at_nodes} Coordinates at n nodes of the element.
 * @param[in] {U_at_nodes} Displacements at n nodes of the element.
 * @param[in] {S_at_int_points} Symmetric stress tensors at ips integration points.
 * @param[out] {CF_at_nodes} The resulting configurational forces for n nodes are written 
 *             to this array.
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

        for assignment in assignments:
            lhs = assignment.lhs
            rhs = assignment.rhs

            lhs_code = c_printer.doprint(lhs)
            rhs_code = c_printer.doprint(rhs)
            if "[" in lhs_code:
                self._c_file_handle.write(f"    {lhs_code} = {rhs_code};\n")
            else:
                self._c_file_handle.write(f"    double {lhs_code} = {rhs_code};\n")

        self._c_file_handle.write("}\n")

        # c code for multiple elements
        self._c_file_handle.write(f'''
/**
 * Computes the configurational forces for num_elem elements of typ {element_typ}.
 * Each element has n={n_} nodes and ips={ips_} integration points.
 * The computation is {method_name}.
 * A static load case is assumed.
 * 
 * @param[in] num_elem number of elements
 * @param[in] {e} internal energy density
 * @param[in] {X_at_nodes} Coordinates at n nodes of num_elem elements.
 * @param[in] {U_at_nodes} Displacements at n nodes of num_elem elements.
 * @param[in] {S_at_int_points} Symmetric stress tensors at ips integration points 
 *            for num_elem elements.
 * @param[out] {CF_at_nodes} The resulting configurational forces for n nodes and 
 *             num_elem elements are written to this array.
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
        # TODO: doc
        self._py_file_handle.write(f"""

def {function_name_n_elements}(
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


_map_typ_and_method_to_CF_function[('{element_typ}', '{method}')] = {function_name_n_elements}
""")


if __name__ == '__main__':
    write_code_for_all_element_types("CPE4")
    print("ok")
