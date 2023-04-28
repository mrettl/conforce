from typing import List, TextIO, Optional
import subprocess
import os

import numpy as np
import sympy as sy

from sympy.codegen import ast
from sympy.printing.c import C99CodePrinter

from cf import element_definitions
from cf.expressions import compute_CF, TermCollector


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

    # create abstract code assignments
    assignments = TermCollector(
        symbols_to_expressions=symbols_to_expressions,
        R=R,
        points=np.concatenate([R_at_int_points_, R_at_nodes_], dtype=float),
        point_names=[*int_point_names, *node_names]
    ).doit(
        input_symbols=[e_at_int_points, X_at_nodes, U_at_nodes, S_at_int_points],
        result_array=CF_at_nodes,
        cse=True
    )

    compiler.write_function(
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
        self._write_header_at_enter = bool(write_header_at_enter)
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
        self._c_file_handle.write(
            "#ifdef _WIN32\n"
            "#    define CF_API __declspec(dllexport)\n"
            "#else\n"
            "#    define CF_API\n"
            "#endif\n"
            "\n"
            "#include <math.h>\n"  # contains power
            "#include <stddef.h>\n"  # contains size_t
        )

        self._py_file_handle.write(
            f"import numpy as np\n"
            f"import ctypes\n"
            f"import os\n"
            f"\n"
            f"# load c library\n"
            f"if os.name == 'nt':\n"
            f"    # windows\n"
            f"    lib = ctypes.cdll.LoadLibrary(os.path.abspath(os.path.join(\n"
            f"        __file__, os.path.pardir, '{self._name}.dll' \n"
            f"    )))\n"
            f"else:\n"
            f"    # linux\n"
            f"    lib = ctypes.cdll.LoadLibrary(os.path.abspath(os.path.join(\n"
            f"        __file__, os.path.pardir, '{self._name}.so'\n"
            f"    )))\n"
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
            f"        double {e}[{ips_}], \n"
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
            f"        double {e}[][{ips_}], \n"
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
        # TODO: doc
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

    assert {e}.shape == (num_elem, {ips_})
    assert {X_at_nodes}.shape == (num_elem, {n_}, {d_})
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
    # generate and compile c and python code
    with CPyCodeCompiler("cf_c", compile_at_exit=True) as compiler:
        compiler.write_headers()
        write_code_for_element_type("CPE4", True, compiler)


if __name__ == '__main__':
    write_code_for_all_element_types("CPE4", "C3D8")
    print("ok")
