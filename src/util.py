r"""
In this module all auxiliary functions are defined. The basic functionalities are:

- Generating symbolic shape functions
- Generating numeric matrices and Tensors
- Analytical helper functions
- Code generation of analytical expressions
- Generating wrappers and headers
- C-code generation including integration
"""
from typing import Union, List, Tuple
import numpy as np
import sympy as sy

ARRAY_LIKE = Union[List, Tuple, np.ndarray]


class Shapes(object):
    """
    Assembly of shape functions.
    The shapes can be created from element points (nodes)
    and the shape powers (e.g. `x**p0 * y**p1 * z**p2`)

    >>> points=np.array([
    ...         [0, 0, 0],
    ...         [1, 0, 0],
    ...         [0, 1, 0]
    ...     ])
    >>> shapes = Shapes.from_points(
    ...     nodes=points,
    ...     shape_powers=[
    ...         [0, 0, 0],
    ...         [1, 0, 0],
    ...         [0, 1, 0]
    ...     ]
    ... )

    This creates a coefficient matrix `coef` that guarantees

    >>> all([
    ...     shapes.eval_shape(i, points[i]) == 1
    ...     for i in range(len(points))
    ... ])
    True

    and

    >>> all([
    ...     shapes.eval_shape(j, points[i]) == 0
    ...     for j in range(len(points))
    ...     for i in range(len(points))
    ...     if j != i
    ... ])
    True

    """

    def __init__(self, coef: ARRAY_LIKE, shape_powers: ARRAY_LIKE):
        self.coef = np.array(coef)
        self.shape_powers = np.array(shape_powers)

    @classmethod
    def from_points(cls, points: ARRAY_LIKE, shape_powers: ARRAY_LIKE):
        """
        creates shape functions, such that each shape function is one at one point and zero at the other points

        Parameters
        ----------
        points : (number of points,3) nd-array
            Coordinates of points
        shape_powers : (number of points,3) nd-array
            Powers of shape functions in 3 dimensional space

        Returns
        -------
        coef : (number of points, r coef, s coef, t coef) nd-array
            Coefficients of shape functions;
        """
        points = np.asarray(points, dtype=float)
        shape_powers = np.asarray(shape_powers, dtype=int)

        num_points = points.shape[0]
        num_powers = shape_powers.shape[0]
        assert num_points == num_powers

        # a[i, j] = i-th shape power evaluated at j-th point
        a = np.prod(
                points.reshape((num_points, 1, 3))
                ** shape_powers.reshape((1, num_powers, 3)),
                axis=-1
            ).T

        # solve coef_matrix * a = I,
        # such that every shape function is one at exactly one point and zero at the other points
        coef_matrix = np.linalg.inv(a)

        # shape_i(x,y,z) = sum([coef[i, j, k, l] * x**j * y**k * z**l for i, j, k in ...])
        coef = np.zeros([num_points] + [np.max(power) + 1 for power in shape_powers.T], dtype=float)
        for i, shape_power in enumerate(shape_powers):
            coef[:, shape_power[0], shape_power[1], shape_power[2]] = coef_matrix[:, i]

        return cls(coef, shape_powers)

    @property
    def num_shapes(self):
        return len(self.coef)

    def eval_shape(self, j: int, point: ARRAY_LIKE):
        """
        Evaluates the j-th shape function at the given point

        Parameters
        ----------
        j: int
            index of shape function
        point: array like of shape (3,)
            point on which to evaluate the j-shape function

        Returns
        -------
        value: float
        """
        coef = self.coef
        shape_powers = self.shape_powers

        return sum([
            coef[j, shape_powers[i, 0], shape_powers[i, 1], shape_powers[i, 2]]
            * point[0] ** shape_powers[i, 0]
            * point[1] ** shape_powers[i, 1]
            * point[2] ** shape_powers[i, 2]
            for i in range(self.num_shapes)
        ])

    def to_symbolic(self, r: sy.Symbol, s: sy.Symbol, t: sy.Symbol):
        """
        Parameters
        ----------
        r: first coordinate
        s: second coordinate
        t: third coordinate

        Returns
        -------
        symbolic shape functions of (r, s, t).
        """
        return sy.lambdify((r, s, t), sy.Matrix([
            sy.nsimplify(self.eval_shape(j, (r, s, t)))
            for j in range(self.num_shapes)
        ]))


def tensor_from_vector_notation(vec) -> sy.MutableMatrix:
    """
    Create a symmetric second order Tensor from vector notation.

    This function is Abaqus specific. 
    eg. vector notation of the stress tensor:
    (S11,S22,S33,S12,S13,S23)

    Parameters
    ----------
    vec : (6, ) array_like
     consisting of SymPy symbols or SymPy matrix
     Second order tensor in vector notation

    Returns
    -------
    tensor : SymPy matrix of shape (3, 3)
        Second order tensor (sympy matrix object)
    """
    return sy.Matrix([
            [vec[0], vec[3], vec[4]],
            [vec[3], vec[1], vec[5]],
            [vec[4], vec[5], vec[2]],
        ])


def vector_from_tensor_notation(tensor: sy.Matrix):
    """
    returns the abaqus vector notation of a tensor.

    Parameters
    ----------
    tensor: array_like (3, 3) matrix

    Returns
    -------
    vector of shape (3, 1)

    """
    return sy.Matrix([
        tensor[0, 0], tensor[1, 1], tensor[2, 2],
        tensor[0, 1], tensor[0, 2], tensor[1, 2],
    ])


def get_jac(shapeFunc, points_xyz, r, s, t):
    """
    Calculate the jacobi matrix symbolically for a element

    Parameters
    ----------
    shapeFunc : (number of nodes, 1) SymPy matrix 
        Symbolic definition of the shape functions
    points_xyz : (number of nodes,3)
        Symbolic definition of the nodal coordinates
    r : SymPy symbol
        Natural coordinate r
    s : SymPy symbol
        Natural coordinate s
    t : SymPy symbol
        Natural coordinate t

    Returns
    -------
    jac : 
        Symbolic definition of the jacobi matrix of an element
        depending on the natural coordinates r,s,t
    """
    jac = np.zeros((3, 3), dtype=np.object)
    for node in range(shapeFunc.shape[0]):
        jac[0, 0] += sy.diff(shapeFunc[node], r) * points_xyz[node, 0]
        jac[0, 1] += sy.diff(shapeFunc[node], r) * points_xyz[node, 1]
        jac[0, 2] += sy.diff(shapeFunc[node], r) * points_xyz[node, 2]
        jac[1, 0] += sy.diff(shapeFunc[node], s) * points_xyz[node, 0]
        jac[1, 1] += sy.diff(shapeFunc[node], s) * points_xyz[node, 1]
        jac[1, 2] += sy.diff(shapeFunc[node], s) * points_xyz[node, 2]
        jac[2, 0] += sy.diff(shapeFunc[node], t) * points_xyz[node, 0]
        jac[2, 1] += sy.diff(shapeFunc[node], t) * points_xyz[node, 1]
        jac[2, 2] += sy.diff(shapeFunc[node], t) * points_xyz[node, 2]
    return sy.Matrix(jac)


def get_val_at_int_Pkt(shapeFunc, value):
    """
    Calculates a value at the integration point of an element

    Parameters
    ----------
    shapeFunc : (number of nodes, 1) SymPy matrix or nd-array
        Symbolic definition of the shape functions
    value : (number of nodes,3) SymPy matrix or nd-array
        Nodal value 

    Returns
    -------
    val : (3,) SymPy matrix 
        Value of Element Nodal results inside the element
    """
    val = np.zeros((3), dtype=np.object)
    for node in range(shapeFunc.shape[0]):
        val[0] += shapeFunc[node] * value[node, 0]
        val[1] += shapeFunc[node] * value[node, 1]
        val[2] += shapeFunc[node] * value[node, 2]
    return sy.Matrix(val)


def inversion_and_det_3x3(m_in):
    """
    Inverts a symbolic 3x3 matrix and calculates the determinant.
    This can be also accomplished with SymPy, the use of this function is for performance reasons only.

    Parameters
    ----------
    m_in : (3,3) SymPy matrix
        Symbolic 3x3 matrix to be inverted

    Returns
    -------
    minv : (3,3) SymPy matrix
        Inverted symbolic matrix
    determinant : SymPy expression
        Determinant of the matrix
    """
    shape = m_in.shape
    m = m_in.reshape(1, 9)

    minv = sy.Matrix(np.zeros((1, 9), dtype=np.object))
    determinant = m[0, 0] * m[0, 4] * m[0, 8] + m[0, 3] * m[0, 7] * m[0, 2] + m[0, 6] * m[0, 1] * m[0, 5] - m[0, 0] * m[
        0, 5] * m[0, 7] - m[0, 2] * m[0, 4] * m[0, 6] - m[0, 1] * m[0, 3] * m[0, 8]
    determinant_inv = 1 / determinant
    minv[0, 0] = (m[0, 4] * m[0, 8] - m[0, 5] * m[0, 7]) * determinant_inv
    minv[0, 1] = (m[0, 2] * m[0, 7] - m[0, 1] * m[0, 8]) * determinant_inv
    minv[0, 2] = (m[0, 1] * m[0, 5] - m[0, 2] * m[0, 4]) * determinant_inv
    minv[0, 3] = (m[0, 5] * m[0, 6] - m[0, 3] * m[0, 8]) * determinant_inv
    minv[0, 4] = (m[0, 0] * m[0, 8] - m[0, 2] * m[0, 6]) * determinant_inv
    minv[0, 5] = (m[0, 2] * m[0, 3] - m[0, 0] * m[0, 5]) * determinant_inv
    minv[0, 6] = (m[0, 3] * m[0, 7] - m[0, 4] * m[0, 6]) * determinant_inv
    minv[0, 7] = (m[0, 1] * m[0, 6] - m[0, 0] * m[0, 7]) * determinant_inv
    minv[0, 8] = (m[0, 0] * m[0, 4] - m[0, 1] * m[0, 3]) * determinant_inv
    return minv.reshape(shape[0], shape[0]), determinant


def inversion_and_det_2x2(m):
    """
    Inverts a symbolic 2x2 matrix and calculates the determinant.

    This can be also accomplished with SymPy, the use of this function is for performance reasons only.

    Parameters
    ----------
    m : (2,2) SymPy matrix
        Symbolic 2x2 matrix to be inverted
    
    Returns
    -------
    minv : (2,2) SymPy matrix
        Inverted symbolic matrix
    determinant : SymPy expression
        Determinant of the matrix
    """
    minv = sy.Matrix(np.zeros((2, 2), dtype=np.object))
    determinant = m[0, 0] * m[1, 1] - m[0, 1] * m[1, 0]
    determinant_inv = 1 / determinant
    minv[0, 0] = m[1, 1] * determinant_inv
    minv[0, 1] = -1 * m[0, 1] * determinant_inv
    minv[1, 0] = -1 * m[1, 0] * determinant_inv
    minv[1, 1] = m[0, 0] * determinant_inv
    return minv, determinant


def get_det_and_inv(jac):
    """
    Inverts a symbolic matrix and calculates the determinant.

    This functions works for 3d-elements and plane strain 2d-elements.
    If the determinant of the 3x3 matrix is zero, the 2d- plane strain case 
    is assumed.

    Parameters
    ----------
    jac : (2 or 3,2 or 3) SymPy matrix
        Symbolic 3x3 matrix to be inverted

    Returns
    -------
    inv_jac : (2 or 3,2 or 3) SymPy matrix
        Inverted symbolic matrix
    det : SymPy expression
        Determinant of the matrix
    """
    inv_jac,det = inversion_and_det_3x3(jac)
    
    #If determinant is zero assume 2x2 sub-matrix (eg. CPE4 Jacobian)
    if det == 0:
        inv_jac = np.zeros((3,3),dtype=np.object)
        inv_jac[:2,:2],det = inversion_and_det_2x2(jac[:2,:2])
    return sy.Matrix(inv_jac),det


def evaluate_rst(sym_Mat, Points):
    New_Mat = np.zeros((Points.shape[0]), dtype=np.object)
    for num_Point in range(Points.shape[0]):
        subs = (("r", Points[num_Point, 0]), ("s", Points[num_Point, 1]), ("t", Points[num_Point, 2]))
        New_Mat[num_Point] = sym_Mat.subs(subs)
    return New_Mat


def calc_dN_rst(shapeFunc, r, s, t):
    """
    Calculates the derivative of the shape functions with respect to
    the natural coordinates r,s,t

    Parameters
    ----------
    shapeFunc : (number of nodes, 1) SymPy matrix 
        Symbolic definition of the shape functions
    r : SymPy symbol
        Natural coordinate r
    s : SymPy symbol
        Natural coordinate s
    t : SymPy symbol
        Natural coordinate t

    Returns
    -------
    dN_rst : (number of nodes,3) SymPy matrix
        Derived shape functions
    """
    dN_rst = sy.Matrix(np.zeros((shapeFunc.shape[0], 3), dtype=np.object))
    for node in range(shapeFunc.shape[0]):
        dN_rst[node, 0] += sy.diff(shapeFunc[node], r)
        dN_rst[node, 1] += sy.diff(shapeFunc[node], s)
        dN_rst[node, 2] += sy.diff(shapeFunc[node], t)
    return dN_rst


def lambdify_numba(args, arg_names, expr, expr_name="_lambdified"):
    """
    Generates Numba code for a given symbolic expression. 

    This function is called by gen_Integration_Numba
    This function isn't tested and is an alternative to the function which generates C-code.

    Parameters
    ----------
    args : Tuple(symbolic objects)
        Tuple of symbolically defined arguments
    arg_names : Tuple(str,)
        Tuple of the symbol names
    expr : SymPy expression
        Symbolic expression for which the numerical implementation should be generated
    expr_name : string
        Name of the generated function.

    Returns
    -------
    s : string
        Numba function as string 
    """
    import sympy
    sy.printing.fcode
    from sympy.printing.pycode import PythonCodePrinter as Printer
    printer = Printer({'fully_qualified_modules': False, 'inline': True, 'allow_unknown_functions': True})

    arg_names_str = ""
    for i in range(len(arg_names)):
        arg_names_str += arg_names[i]
        if i < len(args) - 1:
            arg_names_str += ","

    ###############
    # generate outputs
    outputs = ""
    if isinstance(expr, (list, tuple)):
        for i in range(len(expr)):
            outputs += ",res_" + str(i)
    else:
        outputs += ",res_0"

    s = ""
    s = ""
    s += "from numpy import *\n"
    s += "import numba as nb\n"
    s += "\n"
    s += "@nb.njit(fastmath=True,error_model='numpy')\n"
    s += "def " + expr_name + "(" + arg_names_str + outputs + "):\n"

    # Expand Matrix args
    for ii in range(len(args)):
        arg = args[ii]
        if isinstance(arg, (sympy.Matrix, sympy.ImmutableMatrix)):
            for i in range(arg.shape[0]):
                for j in range(arg.shape[1]):
                    s += "    " + str(arg[i, j]) + " = " + arg_names[ii] + "[" + str(i) + "," + str(j) + "]\n"
        if isinstance(arg, (tuple, list)):
            for i in range(len(arg)):
                s += "    " + str(arg[i]) + " = " + arg_names[ii] + "[" + str(i) + "]\n"
        if hasattr(arg, 'name'):
            if arg.name != arg_names[ii]:
                s += "    " + str(arg) + " = " + arg_names[ii] + "\n"
    s += "    \n"

    # Simplify expr
    sub_exprs, simplified_rhs = sympy.cse(expr, order='none')

    # Write Subexpressions
    for sub_expr in sub_exprs:
        s += "    " + str(sub_expr[0]) + " = " + printer.doprint(sub_expr[1]) + "\n"
    s += "    \n"

    # Write rhs_expressions and generate return
    res_counter = 0
    return_s = "    return "
    for rhs in simplified_rhs:
        # If the expr is a matrix
        if rhs.is_Matrix:
            rhs_result = sympy.MatrixSymbol('rhs_result', rhs.shape[0], rhs.shape[1])
            for i in range(rhs.shape[0]):
                for j in range(rhs.shape[1]):
                    s += "    res_" + str(res_counter) + "[" + str(i) + "," + str(j) + "]" + " = " + \
                         printer.doprint(rhs[i * rhs.shape[1] + j]) + "\n"
        else:
            s += "    res_" + str(res_counter) + " = " + printer.doprint(rhs) + "\n"
        return_s += "res_" + str(res_counter) + " "
        res_counter += 1

    # Return output
    s += "    \n"
    s += return_s
    s += "\n"
    return s


def gen_Integration_Numba(rst, coord, Element_U, S, SENER, PENER, int_points, weights, expr, expr_name="C3D20"):
    """
    Generates Numba code for the static implementation of configurational forces at element nodal position.

    The generation of Numba code is not officially supported by this package, therfore this function is left as an example.
    This function performs the numeric integration of function values determined at the integration points.
    The numeric representation at the integration point is generated by lambdify_numba.

    Parameters
    ----------
    rst : (3) SymPy matrix
        Symbolically defined coordinates
    coord : (Number of Nodes,3) SymPy matrix
        Symbolically defined nodal coordinates
    Element_U : (Number of Nodes,3) SymPy matrix
        Symbolically defined nodal displacements
    S : (6,1) or (3,3) SymPy matrix
        Symbolically defined stress tensor at integration point
    SENER : SymPy symbol
        Symbolically defined strain energy at integration point
    PENER : SymPy symbol
        Symbolically defined plastic energy at integration point
    int_points : (Number of integration points,3) nd-array
        Natural coordinates of the integration point of an element
    weights : (Number of integration points) nd-array
        Integration weights
    expr : SymPy expression
        Symbolic expression for which the numerical implementation should be generated
    expr_name : string
        Name of the generated function.

    Returns
    -------
    s : string
        Numba function as string
    """
    s = ""
    s += lambdify_numba((rst, coord, Element_U, S, SENER, PENER), ("rst", "coord", "Element_U", "S", "SENER", "PENER"),
                        expr, expr_name=expr_name)
    s += "\n"
    s += "@nb.njit(fastmath=True,parallel=True,error_model='numpy',cache=True)\n"
    s += "def Integration_" + expr_name + "(Coords,Element_U,S,PENER,SENER):\n"
    s += '    """\n'
    s += '    Inputs:\n'
    s += '        Coords.........Element_Nodal_Coordinates [num_elem,numNodes,3]\n'
    s += '        Element_U......Element_Nodal_Displacements [num_elem,numNodes,3]\n'
    s += '        S..............Stress at integration points [num_elem,num_int_points,6]\n'
    s += '        PENER..........Plastic strain energy [num_elem,num_int_points]\n'
    s += '        SENER..........Elastic strain energy [num_elem,num_int_points]\n'
    s += '    """\n'
    # print constants
    s += '    weights=np.array(('
    for i in range(weights.shape[0]):
        s += str(weights[i]) + ","
    s += '))\n'
    s += '    int_points=np.array(\n'
    for i in range(int_points.shape[0]):
        s += '        ('
        for j in range(int_points.shape[1]):
            s += str(int_points[i, j]) + ","
        if i == int_points.shape[0] - 1:
            s += ')'
        else:
            s += ')\n'
    s += ')\n'
    s += '\n'
    s += '    Conf_Force=np.zeros((Coords.shape[0],Coords.shape[1],3))\n'
    s += '    for i in nb.prange(Coords.shape[0]):\n'
    s += '        res_0=np.empty((Coords.shape[1],3))\n'
    s += '        for j in range(int_points.shape[0]):\n'
    s += '            res_0=C_Force_' + expr_name + '(int_points[j],Coords[i],Element_U[i],S_vec[i,j],PENER[i,j],SENER[i,j],res_0)\n'
    s += '            for num_node in range(Coords.shape[1]):\n'
    s += '                for ii in range(3):\n'
    s += '                    Conf_Force[i,num_node,ii]+=res_0[num_node,ii]*weights[j]\n'
    s += '    return Conf_Force\n'
    s += '\n'
    return s


def gen_C_header():
    """
    Generates the header for the implementation in C.

    Additionally a function to calculate nodal unique values from 
    element nodal values is generated.

    Returns
    -------
    s : string
        C-header and a function to calculate nodal unique values as string
    """
    s = '#ifdef _WIN32\n'
    s += '#    define Conf_Forces_API __declspec(dllexport)\n'
    s += '#else\n'
    s += '#    define Conf_Forces_API\n'
    s += '#endif\n'
    s += '#include <math.h>\n'
    s += '#include <stdint.h>\n'
    s += '#include <stddef.h>\n\n'
    s += 'Conf_Forces_API void calc_Nodal_C(size_t num_elem_nodal,double CF_out[][3],double CF[][3],int64_t inverse[]){\n'
    s += '    for (size_t i=0;i<num_elem_nodal;i++){\n'
    s += '        for (size_t j=0;j<3;j++){\n'
    s += '            CF_out[inverse[i]][j]+=CF[i][j];\n'
    s += '        }\n'
    s += '    }\n'
    s += '}\n\n'
    return s


def lambdify_C(args, arg_names, expr, expr_name="_lambdified"):
    """
    Generates C99 code for a given symbolic expression.

    This function is called by gen_Integration_C_static or gen_Integration_C_dynamic. In this package this function is used to 
    generate the numeric representation of configurational forces at the integration point, or more specifically the inner part of the integral.

    Parameters
    ----------
    args : Tuple(symbolic objects)
        Tuple of symbolically defined arguments
    arg_names : Tuple(str,)
        Tuple of the symbol names
    expr : SymPy expression
        Symbolic expression for which the numerical implementation should be generated
    expr_name : string
        Name of the generated function.

    Returns
    -------
    s : string
        C-function as string 
    """

    import sympy
    #Code Printer is sympy version dependent
    try:
        from sympy.printing.c import C99CodePrinter as Printer
    except:
        from sympy.printing.ccode import C99CodePrinter as Printer
    from sympy.codegen.rewriting import create_expand_pow_optimization
    printer = Printer()
    expand_opt = create_expand_pow_optimization(3)

    arg_names_str = ""
    for i in range(len(arg_names)):
        if isinstance(args[i], (sympy.Matrix, sympy.ImmutableMatrix, sympy.MatrixSymbol, list, tuple)):
            arg_names_str += "double *" + arg_names[i]
        else:
            arg_names_str += "double " + arg_names[i]
        if i < len(args) - 1:
            arg_names_str += ","

    ###############
    # generate outputs
    outputs = ""
    if isinstance(expr, (list, tuple)):
        for i in range(len(expr)):
            outputs += ",double *res_" + str(i)
    else:
        outputs += ",double *res_0"

    s = ""
    s = ""
    s += "\n"
    s += "void " + expr_name + "(" + arg_names_str + outputs + "){\n"

    # Expand Matrix args
    for ii in range(len(args)):
        arg = args[ii]
        if isinstance(arg, (sympy.Matrix, sympy.ImmutableMatrix)):
            for i in range(arg.shape[0]):
                for j in range(arg.shape[1]):
                    ind = i * arg.shape[1] + j
                    s += "    double " + str(arg[i, j]) + " = " + arg_names[ii] + "[" + str(ind) + "];\n"
        if isinstance(arg, (tuple, list)):
            for i in range(len(arg)):
                s += "    double " + str(arg[i]) + " = " + arg_names[ii] + "[" + str(i) + "];\n"
        if hasattr(arg, 'name'):
            if arg.name != arg_names[ii]:
                s += "    double *" + str(arg) + " = " + arg_names[ii] + ";\n"
    s += "    \n"

    # Simplify expr
    sub_exprs, simplified_rhs = sympy.cse(expr, order='none')

    # Write Subexpressions
    for sub_expr in sub_exprs:
        if isinstance(sub_expr[0], (sympy.Matrix, sympy.ImmutableMatrix, sympy.MatrixSymbol, list, tuple)):
            s += "    double *" + str(sub_expr[0]) + " = " + printer.doprint(expand_opt(sub_expr[1])) + ";\n"
        else:
            s += "    double " + str(sub_expr[0]) + " = " + printer.doprint(expand_opt(sub_expr[1])) + ";\n"
    s += "    \n"

    # Write rhs_expressions and generate return
    res_counter = 0
    return_s = "    return "
    for rhs in simplified_rhs:
        # If the expr is a matrix
        if rhs.is_Matrix:
            rhs_result = sympy.MatrixSymbol('rhs_result', rhs.shape[0], rhs.shape[1])
            for i in range(rhs.shape[0]):
                for j in range(rhs.shape[1]):
                    ind = i * rhs.shape[1] + j
                    s += "    res_" + str(res_counter) + "[" + str(ind) + "]" + " = " + \
                         printer.doprint(rhs[i * rhs.shape[1] + j]) + ";\n"
        else:
            s += "    res_" + str(res_counter) + " = " + printer.doprint(rhs) + ";\n"
        return_s += "res_" + str(res_counter) + " "
        res_counter += 1

    # Return output
    s += "}\n"

    return s


def gen_Integration_C_dyn(rst, coord, rho, Element_U, Element_V, Element_A, S, PENER, SENER, int_points, weights, expr,
                          expr_name):
    """
    Generates C99 code for the dynamic implementation of configurational forces at element nodal position.

    This function is called by :func:`gen_Configurational_Forces_Dynamic` and performs the numeric integration of
    function values determined at the integration points.
    Also have a look at this function too see how rst, coord, rho, Element_U, Element_V, Element_A ,S ,SENER ,PENER are defined.
    The numeric representation at the integration point is generated by :func:`lambdify_C`.
    
    Parameters
    ----------
    rst : (3) SymPy matrix
        Symbolically defined coordinates
    coord : (Number of Nodes,3) SymPy matrix
        Symbolically defined nodal coordinates
    rho : SymPy symbol
        Symbolically defined material density
    Element_U : (Number of Nodes,3) SymPy matrix
        Symbolically defined nodal displacements
    Element_V : (Number of Nodes,3) SymPy matrix
        Symbolically defined nodal velocity
    Element_A : (Number of Nodes,3) SymPy matrix
        Symbolically defined nodal acceleration
    S : (6,1) or (3,3) SymPy matrix
        Symbolically defined stress tensor at integration point
    SENER : SymPy symbol
        Symbolically defined strain energy at integration point
    int_points : (Number of integration points,3) nd-array
        Natural coordinates of the integration point of an element
    weights : (Number of integration points) nd-array
        Integration weights
    expr : SymPy expression
        Symbolic expression for which the numerical implementation should be generated
    expr_name : string
        Name of the generated function.

    Returns
    -------
    s : string
        C-function as string
    """
    numNodes = str(coord.shape[0])
    num_int_points = str(int_points.shape[0])
    s = ""
    s += lambdify_C((rst, coord, rho, Element_U, Element_V, Element_A, S, PENER, SENER),
                    ("rst", "coord", "rho", "Element_U", "Element_V", "Element_A", "S", "PENER", "SENER"), expr,
                    expr_name=expr_name)
    s += "\n"
    s += "Conf_Forces_API void Integration_" + expr_name + "(size_t num_elem,double Coords[][" + numNodes + "][3],double *rho," + \
         "double Element_U[][" + numNodes + "][3],double Element_V[][" + numNodes + "][3],double Element_A[][" + numNodes + "][3],double S[][" + num_int_points + "][6],\n    " + \
         "double PENER[][" + num_int_points + "],double SENER[][" + num_int_points + "],double Conf_Force[][" + numNodes + "][3]){\n"
    s += '    \n'
    s += '    //Inputs:\n'
    s += '    //    Coords.........Element_Nodal_Coordinates  [num_elem][numNodes,3]\n'
    s += '    //    Element_U......Element_Nodal_Displacements[num_elem][numNodes,3]\n'
    s += '    //    Element_V......Element_Nodal_Velocity     [num_elem][numNodes,3]\n'
    s += '    //    Element_A......Element_Nodal_Acceleration [num_elem][numNodes,3]\n'
    s += '    //    S..............Stress at integration points [num_elem][num_int_points,6]\n'
    s += '    //    PENER..........Plastic strain energy [num_elem][num_int_points]\n'
    s += '    //    SENER..........Elastic strain energy [num_elem][num_int_points]\n'
    s += '    //    \n'
    # print constants
    s += '    double weights[' + num_int_points + ']={'
    for i in range(weights.shape[0]):
        s += str(weights[i]) + ","
    s += '};\n'
    s += '    double int_points[' + num_int_points + '][3]={\n'
    for i in range(int_points.shape[0]):
        s += '        {'
        for j in range(int_points.shape[1]):
            s += str(int_points[i, j]) + ","
        if i == int_points.shape[0] - 1:
            s += '}'
        else:
            s += '},\n'
    s += '};\n'
    s += '\n'
    s += '    #pragma omp parallel for\n'
    s += '    for (size_t i=0;i<num_elem;i++){\n'
    s += '        double TMP[' + numNodes + '][3];\n'
    s += '        for (size_t j=0;j<' + num_int_points + ';j++){\n'
    s += '            ' + expr_name + '((double*) int_points[j],(double*) Coords[i],rho[i],(double*) Element_U[i],(double*) Element_V[i],(double*) Element_A[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);\n'
    s += '            for (size_t num_node=0;num_node<' + numNodes + ';num_node++){\n'
    s += '                for (size_t ii=0;ii<3;ii++){\n'
    s += '                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];\n'
    s += '                }\n'
    s += '            }\n'
    s += '        }\n'
    s += '    }\n'
    s += '}\n'
    return s


def gen_Integration_C_static(rst, coord, Element_U, S, PENER, SENER, int_points, weights, expr, expr_name):
    """
    Generates C99 code for the static implementation of configurational forces at element nodal position.

    This function is called by :func:`gen_Configurational_Forces_Static` and performs the numeric integration of
    function values determined at the integration points.
    Have a look at this function too see how rst, coord, Element_U, S, SENER, PENER are defined.
    The numeric representation at the integration point is generated by :func:`lambdify_C`.
    
    Parameters
    ----------
    rst : (3) SymPy matrix
        Symbolically defined coordinates
    coord : (Number of Nodes,3) SymPy matrix
        Symbolically defined nodal coordinates
    Element_U : (Number of Nodes,3) SymPy matrix
        Symbolically defined nodal displacements
    S : (6,1) or (3,3) SymPy matrix
        Symbolically defined stress tensor at integration point
    SENER : SymPy symbol
        Symbolically defined strain energy at integration point
    PENER : SymPy symbol
        Symbolically defined plastic energy at integration point
    int_points : (Number of integration points,3) nd-array
        Natural coordinates of the integration point of an element
    weights : (Number of integration points) nd-array
        Integration weights
    expr : SymPy expression
        Symbolic expression for which the numerical implementation should be generated
    expr_name : string
        Name of the generated function.

    Returns
    -------
    s : string
        C-function as string
    """
    numNodes = str(coord.shape[0])
    num_int_points = str(int_points.shape[0])
    s = ""
    s += lambdify_C((rst, coord, Element_U, S, PENER, SENER), ("rst", "coord", "Element_U", "S", "PENER", "SENER"),
                    expr, expr_name=expr_name)
    s += "\n"
    s += "Conf_Forces_API void Integration_" + expr_name + "(size_t num_elem,double Coords[][" + numNodes + "][3]," + \
         "double Element_U[][" + numNodes + "][3],double S[][" + num_int_points + "][6],\n    " + \
         "double PENER[][" + num_int_points + "],double SENER[][" + num_int_points + "],double Conf_Force[][" + numNodes + "][3]){\n"
    s += '    \n'
    s += '    //Inputs:\n'
    s += '    //    Coords.........Element_Nodal_Coordinates [num_elem][numNodes,3]\n'
    s += '    //    Element_U......Element_Nodal_Displacements [num_elem][numNodes,3]\n'
    s += '    //    S..............Stress at integration points [num_elem][num_int_points,6]\n'
    s += '    //    PENER..........Plastic strain energy [num_elem][num_int_points]\n'
    s += '    //    SENER..........Elastic strain energy [num_elem][num_int_points]\n'
    s += '    //    \n'
    # print constants
    s += '    double weights[' + num_int_points + ']={'
    for i in range(weights.shape[0]):
        s += str(weights[i]) + ","
    s += '};\n'
    s += '    double int_points[' + num_int_points + '][3]={\n'
    for i in range(int_points.shape[0]):
        s += '        {'
        for j in range(int_points.shape[1]):
            s += str(int_points[i, j]) + ","
        if i == int_points.shape[0] - 1:
            s += '}'
        else:
            s += '},\n'
    s += '};\n'
    s += '\n'
    s += '    #pragma omp parallel for\n'
    s += '    for (size_t i=0;i<num_elem;i++){\n'
    s += '        double TMP[' + numNodes + '][3];\n'
    s += '        for (size_t j=0;j<' + num_int_points + ';j++){\n'
    s += '            ' + expr_name + '((double*) int_points[j],(double*) Coords[i],(double*) Element_U[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);\n'
    s += '            for (size_t num_node=0;num_node<' + numNodes + ';num_node++){\n'
    s += '                for (size_t ii=0;ii<3;ii++){\n'
    s += '                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];\n'
    s += '                }\n'
    s += '            }\n'
    s += '        }\n'
    s += '    }\n'
    s += '}\n'
    return s


def gen_Python_Wrapper_Header(lib_name):
    """
    Generates the header of the Python wrapper used to call the C-implementations from Python.

    This function also defines the function calc_Nodal, which is used to calculate unique nodal results
    of configurational forces calculated on element nodal position.

    Parameters
    ----------
    lib_name : string
        Name of the compiled library

    Returns
    -------
    s : string
        Python code as string
    """
    s = 'import numpy as np\n'
    s += 'import ctypes\n'
    s += 'import os\n\n'
    s += '# Determine Operating System\n'
    s += "if os.name == 'nt':\n"
    s += '    lib = ctypes.cdll.LoadLibrary("' + lib_name + '.dll")\n'
    s += 'else:\n'
    s += '    lib = ctypes.cdll.LoadLibrary("./' + lib_name + '.so")\n\n'
    s += 'def calc_Nodal(Element_Connectivity,CF):\n'
    s += '    """\n'
    s += '    Input:\n'
    s += '    Element_Connectivity.......List of np.arrays or np.array\n'
    s += '    CF_in......................List of np.arrays or np.array\n'
    s += '    Output:\n'
    s += '    Node_labels................Array of unique Nodelabels where CF is calcualted\n'
    s += '    CF_out.....................Array of unique Nodal Configurational Forces [num_nodes,3]\n'
    s += '    """\n'
    s += '    if not isinstance(Element_Connectivity,list):\n'
    s += '        Element_Connectivity=[Element_Connectivity]\n'
    s += '        CF=[CF]\n\n'
    s += '    Element_Connectivity=np.concatenate([elem_con.reshape(-1) for elem_con in Element_Connectivity])\n'
    s += '    CF=np.concatenate([cf.reshape(-1,3) for cf in CF])\n'
    s += '    assert Element_Connectivity.shape[0]==CF.shape[0]\n\n'
    s += '    Node_labels,inverse=np.unique(Element_Connectivity,return_inverse=True)\n\n'
    s += '    CF=np.ascontiguousarray(CF,dtype=np.float64)\n'
    s += '    inverse=np.ascontiguousarray(inverse,dtype=np.int64)\n'
    s += '    CF_out=np.zeros((Node_labels.shape[0],3),dtype=np.float64)\n\n'
    s += '    num_elem_nodal=ctypes.c_size_t(inverse.shape[0])\n'
    s += '    CF_out_p=CF_out.ctypes.data_as(ctypes.c_void_p)\n'
    s += '    CF_p=CF.ctypes.data_as(ctypes.c_void_p)\n'
    s += '    inverse_p=inverse.ctypes.data_as(ctypes.c_void_p)\n\n'
    s += '    lib.calc_Nodal_C(num_elem_nodal,CF_out_p,CF_p,inverse_p)\n'
    s += '    return Node_labels,CF_out\n\n'
    return s


def gen_Python_Wrapper_static(expr_name, numNodes, num_int_points):
    """
    Generates the Python wrapper used to call the static configurational force implementation from Python.

    This wrapper also includes basic checks for array types and datatypes.

    Parameters
    ----------
    expr_name : string
        Name of the generated function.
    
    numNodes : int
        Number of nodes
    num_int_points : int
        Number of integration points

    Returns
    -------
    s : string
        Python code as string
    """
    s = "def calc_Conf_Force_" + expr_name + "(Coords,Element_U,S_vec,PENER,SENER,method='mbf'):\n"
    s += '    """\n'
    s += '    Inputs:\n'
    s += '        Coords.........Element_Nodal_Coordinates [num_elem,numNodes,3]\n'
    s += '        Element_U......Element_Nodal_Displacements [num_elem,numNodes,3]\n'
    s += '        S_vec......... Stress at integration points [num_elem,num_int_points,6]\n'
    s += '        PENER..........Plastic strain energy [num_elem,num_int_points]\n'
    s += '        SENER..........Elastic strain energy [num_elem,num_int_points]\n'
    s += '    """\n'
    s += '    #This is for safety, the C code will crash without warnings on wrong inputs...\n'
    s += '    ################################################\n'
    s += '    num_elem=Coords.shape[0]\n'
    s += '    numNodes=' + str(numNodes) + '\n'
    s += '    num_int_points=' + str(num_int_points) + '\n\n'
    s += '    assert Coords.shape==(num_elem,numNodes,3)\n'
    s += '    assert Element_U.shape==(num_elem,numNodes,3)\n'
    s += '    assert S_vec.shape==(num_elem,num_int_points,6)\n'
    s += '    assert PENER.shape==(num_elem,num_int_points)\n'
    s += '    assert SENER.shape==(num_elem,num_int_points)\n\n'
    s += '    #Ensure that all arrays are c contiguous\n'
    s += '    Coords=np.ascontiguousarray(Coords,dtype=np.float64)\n'
    s += '    Element_U=np.ascontiguousarray(Element_U,dtype=np.float64)\n'
    s += '    S_vec=np.ascontiguousarray(S_vec,dtype=np.float64)\n'
    s += '    PENER=np.ascontiguousarray(PENER,dtype=np.float64)\n'
    s += '    SENER=np.ascontiguousarray(SENER,dtype=np.float64)\n'
    s += '    ################################################\n\n'
    s += '    num_elem_=ctypes.c_size_t(Coords.shape[0])\n'
    s += '    Coords_p=Coords.ctypes.data_as(ctypes.c_void_p)\n'
    s += '    Element_U_p=Element_U.ctypes.data_as(ctypes.c_void_p)\n'
    s += '    S_vec_p=S_vec.ctypes.data_as(ctypes.c_void_p)\n'
    s += '    PENER_p=PENER.ctypes.data_as(ctypes.c_void_p)\n'
    s += '    SENER_p=SENER.ctypes.data_as(ctypes.c_void_p)\n\n'
    s += '    Conf_Force=np.zeros((Coords.shape[0],Coords.shape[1],3))\n'
    s += '    Conf_Force_p=Conf_Force.ctypes.data_as(ctypes.c_void_p)\n\n'
    s += "    if method=='mbf':\n"
    s += '        lib.Integration_' + expr_name + '_mbf(num_elem_,Coords_p,Element_U_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)\n'
    s += "    elif method=='dbf':\n"
    s += '        lib.Integration_' + expr_name + '_dbf(num_elem_,Coords_p,Element_U_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)\n'
    s += '    else:\n'
    s += '        print("Method not found!")\n'
    s += '    return Conf_Force\n\n'
    return s


def gen_Python_Wrapper_dynamic(expr_name, numNodes, num_int_points):
    """
    Generates the Python wrapper used to call the dynamic configurational force implementation from Python.

    This wrapper also includes basic checks for array types and datatypes.
    
    Parameters
    ----------
    expr_name : string
        Name of the generated function.
    
    numNodes : int
        Number of nodes
    num_int_points : int
        Number of integration points

    Returns
    -------
    s : string
        Python code as string
    """
    s = "def calc_Conf_Force_" + expr_name + "(Coords,rho,Element_U,Element_V,Element_A,S_vec,PENER,SENER):\n"
    s += '    """\n'
    s += '    Inputs:\n'
    s += '        Coords.........Element_Nodal_Coordinates    [num_elem,numNodesPerElem,3]\n'
    s += '        rho............Density                      [num_elem]\n'
    s += '        Element_U......Element_Nodal_Displacements  [num_elem,numNodesPerElem,3]\n'
    s += '        Element_V......Element_Velocity             [num_elem,numNodesPerElem,3]\n'
    s += '        Element_A......Element_Acceleration         [num_elem,numNodesPerElem,3]\n'
    s += '        S_vec......... Stress at integration points [num_elem,num_int_points,6]\n'
    s += '        PENER..........Plastic strain energy        [num_elem,num_int_points]\n'
    s += '        SENER..........Elastic strain energy        [num_elem,num_int_points]\n'
    s += '    """\n'
    s += '    #This is for safety, the C code will crash without warnings on wrong inputs...\n'
    s += '    ################################################\n'
    s += '    num_elem=Coords.shape[0]\n'
    s += '    numNodes=' + str(numNodes) + '\n'
    s += '    num_int_points=' + str(num_int_points) + '\n\n'
    s += '    assert Coords.shape==(num_elem,numNodes,3)\n'
    s += '    assert rho.shape==(num_elem,)\n'
    s += '    assert Element_U.shape==(num_elem,numNodes,3)\n'
    s += '    assert Element_V.shape==(num_elem,numNodes,3)\n'
    s += '    assert Element_A.shape==(num_elem,numNodes,3)\n'
    s += '    assert S_vec.shape==(num_elem,num_int_points,6)\n'
    s += '    assert PENER.shape==(num_elem,num_int_points)\n'
    s += '    assert SENER.shape==(num_elem,num_int_points)\n\n'
    s += '    #Ensure that all arrays are c contiguous\n'
    s += '    Coords=np.ascontiguousarray(Coords,dtype=np.float64)\n'
    s += '    rho=np.ascontiguousarray(rho,dtype=np.float64)\n'
    s += '    Element_U=np.ascontiguousarray(Element_U,dtype=np.float64)\n'
    s += '    Element_V=np.ascontiguousarray(Element_V,dtype=np.float64)\n'
    s += '    Element_A=np.ascontiguousarray(Element_A,dtype=np.float64)\n'
    s += '    S_vec=np.ascontiguousarray(S_vec,dtype=np.float64)\n'
    s += '    PENER=np.ascontiguousarray(PENER,dtype=np.float64)\n'
    s += '    SENER=np.ascontiguousarray(SENER,dtype=np.float64)\n'
    s += '    ################################################\n\n'
    s += '    num_elem_=ctypes.c_size_t(Coords.shape[0])\n'
    s += '    Coords_p=Coords.ctypes.data_as(ctypes.c_void_p)\n'
    s += '    rho_p   =rho.ctypes.data_as(ctypes.c_void_p)\n'
    s += '    Element_U_p=Element_U.ctypes.data_as(ctypes.c_void_p)\n'
    s += '    Element_V_p=Element_V.ctypes.data_as(ctypes.c_void_p)\n'
    s += '    Element_A_p=Element_A.ctypes.data_as(ctypes.c_void_p)\n'
    s += '    S_vec_p=S_vec.ctypes.data_as(ctypes.c_void_p)\n'
    s += '    PENER_p=PENER.ctypes.data_as(ctypes.c_void_p)\n'
    s += '    SENER_p=SENER.ctypes.data_as(ctypes.c_void_p)\n\n'
    s += '    Conf_Force=np.zeros((Coords.shape[0],Coords.shape[1],3))\n'
    s += '    Conf_Force_p=Conf_Force.ctypes.data_as(ctypes.c_void_p)\n\n'
    s += '    lib.Integration_' + expr_name + '(num_elem_,Coords_p,rho_p,Element_U_p,Element_V_p,Element_A_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)\n'
    s += '    return Conf_Force\n\n'
    return s


####### Static Implementation #################################
def gen_Configurational_Forces_Static(poly_power, bild_points, int_points, int_weights, typ, method='mbf'):
    """
    Generates the static configurational force implementation.
    This is a entry function and uses functions defined above to generate the numerical implementation 
    of configurational forces for the static and elasto-plastic case.
    All Inputs are defined in ele_def.py

    Parameters
    ----------
    poly_power : (at least number of nodes,3) nd-array
        Polynomial coefficients of shape functions in 3 dimensional space
    bild_points : (number of nodes,3) nd-array
        Natural coordinates of nodes
    int_points : (number of points, 3) nd-array
        Polynomial coefficients of shape functions in 3 dimensional space
    int_weights : (number of points) nd-array
        Integration weights
    typ : str
        Element type
    method : str (mbf/dbf)
        Motion based or displacement based formulation

    Returns
    -------
    s : string
        Python code as string
    """
    # generate shape functions
    shapes = Shapes.from_points(bild_points, poly_power)
    shapeFuncCoef_bild = shapes.coef

    # Generate some symbols
    num_nodes = shapeFuncCoef_bild.shape[0]
    num_int_points = int_points.shape[0]
    r, s, t = sy.symbols("r s t")
    rst = sy.Matrix((r, s, t))
    coord = sy.MatrixSymbol('coord', num_nodes, 3)
    Element_U = sy.MatrixSymbol('U', num_nodes, 3)
    SENER, PENER = sy.symbols("SENER PENER")
    S_Ten = sy.symbols("S11 S22 S33 S12 S13 S23")
    
    S = tensor_from_vector_notation(S_Ten)
    
    # Generate Shape functions
    shapeFunc = sy.Matrix(shapes.eval_shapes((r, s, t)))

    # Calculate Jacobi matrix
    jac = get_jac(shapeFunc, coord, r, s, t)

    # Calculate inverse Jacobi matrix
    jac_inv, jacobi_det = get_det_and_inv(jac)

    # calculate du/dX
    Jacobi_Element_U = get_jac(shapeFunc, Element_U, r, s, t)
    dU_dx = (jac_inv * Jacobi_Element_U).T

    # Calculate deformation gradient
    Def_grad = sy.Matrix(np.eye(3)) + dU_dx

    # Calculate dN_dxyz
    dN_rst = calc_dN_rst(shapeFunc, r, s, t)
    dN_dxyz = dN_rst * jac_inv.T

    # calculate 1st Piola Kirchhoff stress
    Def_grad_inv, Def_grad_det = get_det_and_inv(Def_grad)
    Piola_1 = Def_grad_det * S * Def_grad_inv.T

    # Calculate Configurational stress
    ENER = SENER + PENER
    if method == 'mbf':
        C = ENER * sy.Matrix(np.eye(3)) - Def_grad.T * Piola_1
    else:
        C = ENER * sy.Matrix(np.eye(3)) - dU_dx.T * Piola_1

    # Calculate the inner part of the integral
    C_Force = dN_dxyz * C.T * jacobi_det

    # generate_Code
    expr_name = ""
    if method == "mbf":
        expr_name = typ + "_mbf"
    else:
        expr_name = typ + "_dbf"

    impl = gen_Integration_C_static(rst, coord, Element_U, S_Ten, PENER, SENER, int_points, int_weights, C_Force,
                                    expr_name)
    return impl


####### Dynamic Implementation ################################
def gen_Configurational_Forces_Dynamic(poly_power, bild_points, int_points, int_weights, typ):
    """
    Generates the dynamic configurational force implementation.
    This is a entry function and uses functions defined above to generate the numerical implementation 
    of configurational forces for the dynamic case.
    All Inputs are defined in ele_def.py

    Parameters
    ----------
    poly_power : (at least number of nodes,3) nd-array
        Polynomial coefficients of shape functions in 3 dimensional space
    bild_points : (number of nodes,3) nd-array
        Natural coordinates of nodes
    int_points : (number of points, 3) nd-array
        Polynomial coefficients of shape functions in 3 dimensional space
    int_weights : (number of points) nd-array
        Integration weights
    typ : str
        Element-type

    Returns
    -------
    s : string
        Python code as string
    """
    # generate shape functions
    shapes = Shapes.from_points(bild_points, poly_power)
    shapeFuncCoef_bild = shapes.coef

    # Generate some symbols
    num_nodes = shapeFuncCoef_bild.shape[0]
    num_int_points = int_points.shape[0]
    r, s, t = sy.symbols("r s t")
    rst = sy.Matrix((r, s, t))
    coord = sy.MatrixSymbol('coord', num_nodes, 3)
    Element_U = sy.MatrixSymbol('U', num_nodes, 3)
    SENER, PENER = sy.symbols("SENER PENER")
    S_Ten = sy.symbols("S11 S22 S33 S12 S13 S23")

    rho = sy.symbols("rho")
    Element_V = sy.MatrixSymbol('V', num_nodes, 3)
    Element_A = sy.MatrixSymbol('A', num_nodes, 3)

    # Generate tensor, vaild for Abaqus notation
    S = tensor_from_vector_notation(S_Ten)

    # Generate Shape functions
    shapeFunc = sy.Matrix(shapes.eval_shapes((r, s, t)))

    # Calculate Jacobi matrix
    jac = get_jac(shapeFunc, coord, r, s, t)

    # Calculate inverse Jacobi matrix
    jac_inv, jacobi_det = get_det_and_inv(jac)

    # calculate du/DX
    Jacobi_Element_U = get_jac(shapeFunc, Element_U, r, s, t)
    dU_dx = (jac_inv * Jacobi_Element_U).T

    # calculate dv/dx
    Jacobi_Element_V = get_jac(shapeFunc, Element_V, r, s, t)
    dV_dx = (jac_inv * Jacobi_Element_V).T

    # Calculate deformation gradient and time derivative of deformation gradient
    Def_grad = dU_dx + sy.Matrix(np.eye(3))
    Def_grad_V = dV_dx

    ##################################################
    # Def_grad_mbf_dbf=Def_grad #mbf
    Def_grad_mbf_dbf = dU_dx  # dbf
    ##################################################
    # Acceleration and Velocity with respect to r,s,t
    acc = get_val_at_int_Pkt(shapeFunc, Element_A)
    velocity = get_val_at_int_Pkt(shapeFunc, Element_V)

    # Calculate Pseudo momentum vector
    P_m = -rho * (Def_grad_V.T * velocity + Def_grad_mbf_dbf.T * acc)

    # Calculate dN_dxyz
    dN_rst = calc_dN_rst(shapeFunc, r, s, t)
    dN_dxyz = dN_rst * jac_inv.T

    # calculate 1st Piola Kirchhoff stress
    Def_grad_inv, Def_grad_det = get_det_and_inv(Def_grad)
    Piola_1 = Def_grad_det * S * Def_grad_inv.T

    # calculate kinetic energy density
    T = (1. / 2.) * rho * (velocity[0] ** 2 + velocity[1] ** 2 + velocity[2] ** 2)

    # Configurational Stress ermitteln
    ENER = SENER + PENER
    C = (ENER - T) * sy.Matrix(np.eye(3)) - Def_grad_mbf_dbf.T * Piola_1

    # Calculate the inner part of the integral
    C_Force = (shapeFunc * P_m.T + dN_dxyz * C.T) * jacobi_det

    # generate_Code
    impl = gen_Integration_C_dyn(rst, coord, rho, Element_U, Element_V, Element_A, S_Ten, PENER, SENER, int_points,
                                 int_weights, C_Force, expr_name=typ)
    return impl
