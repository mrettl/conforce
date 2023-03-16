import numpy as np
import ctypes
import os

# Determine Operating System
if os.name == 'nt':
    lib = ctypes.cdll.LoadLibrary("Conf_Forces.dll")
else:
    lib = ctypes.cdll.LoadLibrary("./Conf_Forces.so")


def resulting_nodal_forces(element_connectivity, element_node_forces):
    """
    Computes the resulting forces for all nodes, that are used inside element_connectivity.
    For nodes connected to more elements, each element has a contributing force vector.
    The resulting force vector is the sum of all those contributions.

    Parameters
    ----------
    element_connectivity: list of length E with an entry for each element.
        Each entry is a list containing integer node labels.
    element_node_forces: forces computed for each element node.
        Like element_connectivity, element_node_forces is a list of E entries.
        Each entry is a list containing force vectors of length 3.
        The force vector `element_node_forces[i][j]` corresponds to the node label `element_connectivity[i][j]`.

    Returns
    -------
    node_labels: np.ndarray of shape (N,) containing unique node labels
    resulting_node_forces: np.ndarray of shape (N, 3) containing resulting force vectors for each node label
    """
    element_connectivity = np.concatenate(element_connectivity)
    element_node_forces = np.concatenate(element_node_forces)
    assert element_connectivity.shape[0] == element_node_forces.shape[0]

    node_labels, inverse = np.unique(element_connectivity, return_inverse=True)

    element_node_forces = np.ascontiguousarray(element_node_forces, dtype=np.float64)
    inverse = np.ascontiguousarray(inverse, dtype=np.int64)
    resulting_node_forces = np.zeros((node_labels.shape[0], 3), dtype=np.float64)

    num_elem_nodal = ctypes.c_size_t(inverse.shape[0])
    resulting_node_forces_out_p = resulting_node_forces.ctypes.data_as(ctypes.c_void_p)
    element_nodal_forces_p = element_node_forces.ctypes.data_as(ctypes.c_void_p)
    inverse_p = inverse.ctypes.data_as(ctypes.c_void_p)

    lib.calc_Nodal_C(num_elem_nodal, resulting_node_forces_out_p, element_nodal_forces_p, inverse_p)
    return node_labels, resulting_node_forces


def calc_Conf_Force_CPE4_static(Coords, Element_U, S_vec, PENER, SENER, method='mbf'):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates [num_elem,numNodes,3]
        Element_U......Element_Nodal_Displacements [num_elem,numNodes,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy [num_elem,num_int_points]
        SENER..........Elastic strain energy [num_elem,num_int_points]
    """
    # This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem = Coords.shape[0]
    numNodes = 4
    num_int_points = 4

    assert Coords.shape == (num_elem, numNodes, 3)
    assert Element_U.shape == (num_elem, numNodes, 3)
    assert S_vec.shape == (num_elem, num_int_points, 6)
    assert PENER.shape == (num_elem, num_int_points)
    assert SENER.shape == (num_elem, num_int_points)

    # Ensure that all arrays are c contiguous
    Coords = np.ascontiguousarray(Coords, dtype=np.float64)
    Element_U = np.ascontiguousarray(Element_U, dtype=np.float64)
    S_vec = np.ascontiguousarray(S_vec, dtype=np.float64)
    PENER = np.ascontiguousarray(PENER, dtype=np.float64)
    SENER = np.ascontiguousarray(SENER, dtype=np.float64)
    ################################################

    num_elem_ = ctypes.c_size_t(Coords.shape[0])
    Coords_p = Coords.ctypes.data_as(ctypes.c_void_p)
    Element_U_p = Element_U.ctypes.data_as(ctypes.c_void_p)
    S_vec_p = S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p = PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p = SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force = np.zeros((Coords.shape[0], Coords.shape[1], 3))
    Conf_Force_p = Conf_Force.ctypes.data_as(ctypes.c_void_p)

    if method == 'mbf':
        lib.Integration_CPE4_static_mbf(num_elem_, Coords_p, Element_U_p, S_vec_p, PENER_p, SENER_p, Conf_Force_p)
    elif method == 'dbf':
        lib.Integration_CPE4_static_dbf(num_elem_, Coords_p, Element_U_p, S_vec_p, PENER_p, SENER_p, Conf_Force_p)
    else:
        print("Method not found!")
    return Conf_Force


def calc_Conf_Force_CPE4_dynamic(Coords, rho, Element_U, Element_V, Element_A, S_vec, PENER, SENER):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates    [num_elem,numNodesPerElem,3]
        rho............Density                      [num_elem]
        Element_U......Element_Nodal_Displacements  [num_elem,numNodesPerElem,3]
        Element_V......Element_Velocity             [num_elem,numNodesPerElem,3]
        Element_A......Element_Acceleration         [num_elem,numNodesPerElem,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy        [num_elem,num_int_points]
        SENER..........Elastic strain energy        [num_elem,num_int_points]
    """
    # This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem = Coords.shape[0]
    numNodes = 4
    num_int_points = 4

    assert Coords.shape == (num_elem, numNodes, 3)
    assert rho.shape == (num_elem,)
    assert Element_U.shape == (num_elem, numNodes, 3)
    assert Element_V.shape == (num_elem, numNodes, 3)
    assert Element_A.shape == (num_elem, numNodes, 3)
    assert S_vec.shape == (num_elem, num_int_points, 6)
    assert PENER.shape == (num_elem, num_int_points)
    assert SENER.shape == (num_elem, num_int_points)

    # Ensure that all arrays are c contiguous
    Coords = np.ascontiguousarray(Coords, dtype=np.float64)
    rho = np.ascontiguousarray(rho, dtype=np.float64)
    Element_U = np.ascontiguousarray(Element_U, dtype=np.float64)
    Element_V = np.ascontiguousarray(Element_V, dtype=np.float64)
    Element_A = np.ascontiguousarray(Element_A, dtype=np.float64)
    S_vec = np.ascontiguousarray(S_vec, dtype=np.float64)
    PENER = np.ascontiguousarray(PENER, dtype=np.float64)
    SENER = np.ascontiguousarray(SENER, dtype=np.float64)
    ################################################

    num_elem_ = ctypes.c_size_t(Coords.shape[0])
    Coords_p = Coords.ctypes.data_as(ctypes.c_void_p)
    rho_p = rho.ctypes.data_as(ctypes.c_void_p)
    Element_U_p = Element_U.ctypes.data_as(ctypes.c_void_p)
    Element_V_p = Element_V.ctypes.data_as(ctypes.c_void_p)
    Element_A_p = Element_A.ctypes.data_as(ctypes.c_void_p)
    S_vec_p = S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p = PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p = SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force = np.zeros((Coords.shape[0], Coords.shape[1], 3))
    Conf_Force_p = Conf_Force.ctypes.data_as(ctypes.c_void_p)

    lib.Integration_CPE4_dynamic(num_elem_, Coords_p, rho_p, Element_U_p, Element_V_p, Element_A_p, S_vec_p, PENER_p,
                                 SENER_p, Conf_Force_p)
    return Conf_Force


def calc_Conf_Force_CPE4R_static(Coords, Element_U, S_vec, PENER, SENER, method='mbf'):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates [num_elem,numNodes,3]
        Element_U......Element_Nodal_Displacements [num_elem,numNodes,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy [num_elem,num_int_points]
        SENER..........Elastic strain energy [num_elem,num_int_points]
    """
    # This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem = Coords.shape[0]
    numNodes = 4
    num_int_points = 1

    assert Coords.shape == (num_elem, numNodes, 3)
    assert Element_U.shape == (num_elem, numNodes, 3)
    assert S_vec.shape == (num_elem, num_int_points, 6)
    assert PENER.shape == (num_elem, num_int_points)
    assert SENER.shape == (num_elem, num_int_points)

    # Ensure that all arrays are c contiguous
    Coords = np.ascontiguousarray(Coords, dtype=np.float64)
    Element_U = np.ascontiguousarray(Element_U, dtype=np.float64)
    S_vec = np.ascontiguousarray(S_vec, dtype=np.float64)
    PENER = np.ascontiguousarray(PENER, dtype=np.float64)
    SENER = np.ascontiguousarray(SENER, dtype=np.float64)
    ################################################

    num_elem_ = ctypes.c_size_t(Coords.shape[0])
    Coords_p = Coords.ctypes.data_as(ctypes.c_void_p)
    Element_U_p = Element_U.ctypes.data_as(ctypes.c_void_p)
    S_vec_p = S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p = PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p = SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force = np.zeros((Coords.shape[0], Coords.shape[1], 3))
    Conf_Force_p = Conf_Force.ctypes.data_as(ctypes.c_void_p)

    if method == 'mbf':
        lib.Integration_CPE4R_static_mbf(num_elem_, Coords_p, Element_U_p, S_vec_p, PENER_p, SENER_p, Conf_Force_p)
    elif method == 'dbf':
        lib.Integration_CPE4R_static_dbf(num_elem_, Coords_p, Element_U_p, S_vec_p, PENER_p, SENER_p, Conf_Force_p)
    else:
        print("Method not found!")
    return Conf_Force


def calc_Conf_Force_CPE4R_dynamic(Coords, rho, Element_U, Element_V, Element_A, S_vec, PENER, SENER):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates    [num_elem,numNodesPerElem,3]
        rho............Density                      [num_elem]
        Element_U......Element_Nodal_Displacements  [num_elem,numNodesPerElem,3]
        Element_V......Element_Velocity             [num_elem,numNodesPerElem,3]
        Element_A......Element_Acceleration         [num_elem,numNodesPerElem,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy        [num_elem,num_int_points]
        SENER..........Elastic strain energy        [num_elem,num_int_points]
    """
    # This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem = Coords.shape[0]
    numNodes = 4
    num_int_points = 1

    assert Coords.shape == (num_elem, numNodes, 3)
    assert rho.shape == (num_elem,)
    assert Element_U.shape == (num_elem, numNodes, 3)
    assert Element_V.shape == (num_elem, numNodes, 3)
    assert Element_A.shape == (num_elem, numNodes, 3)
    assert S_vec.shape == (num_elem, num_int_points, 6)
    assert PENER.shape == (num_elem, num_int_points)
    assert SENER.shape == (num_elem, num_int_points)

    # Ensure that all arrays are c contiguous
    Coords = np.ascontiguousarray(Coords, dtype=np.float64)
    rho = np.ascontiguousarray(rho, dtype=np.float64)
    Element_U = np.ascontiguousarray(Element_U, dtype=np.float64)
    Element_V = np.ascontiguousarray(Element_V, dtype=np.float64)
    Element_A = np.ascontiguousarray(Element_A, dtype=np.float64)
    S_vec = np.ascontiguousarray(S_vec, dtype=np.float64)
    PENER = np.ascontiguousarray(PENER, dtype=np.float64)
    SENER = np.ascontiguousarray(SENER, dtype=np.float64)
    ################################################

    num_elem_ = ctypes.c_size_t(Coords.shape[0])
    Coords_p = Coords.ctypes.data_as(ctypes.c_void_p)
    rho_p = rho.ctypes.data_as(ctypes.c_void_p)
    Element_U_p = Element_U.ctypes.data_as(ctypes.c_void_p)
    Element_V_p = Element_V.ctypes.data_as(ctypes.c_void_p)
    Element_A_p = Element_A.ctypes.data_as(ctypes.c_void_p)
    S_vec_p = S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p = PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p = SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force = np.zeros((Coords.shape[0], Coords.shape[1], 3))
    Conf_Force_p = Conf_Force.ctypes.data_as(ctypes.c_void_p)

    lib.Integration_CPE4R_dynamic(num_elem_, Coords_p, rho_p, Element_U_p, Element_V_p, Element_A_p, S_vec_p, PENER_p,
                                  SENER_p, Conf_Force_p)
    return Conf_Force


def calc_Conf_Force_CPE8_static(Coords, Element_U, S_vec, PENER, SENER, method='mbf'):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates [num_elem,numNodes,3]
        Element_U......Element_Nodal_Displacements [num_elem,numNodes,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy [num_elem,num_int_points]
        SENER..........Elastic strain energy [num_elem,num_int_points]
    """
    # This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem = Coords.shape[0]
    numNodes = 8
    num_int_points = 9

    assert Coords.shape == (num_elem, numNodes, 3)
    assert Element_U.shape == (num_elem, numNodes, 3)
    assert S_vec.shape == (num_elem, num_int_points, 6)
    assert PENER.shape == (num_elem, num_int_points)
    assert SENER.shape == (num_elem, num_int_points)

    # Ensure that all arrays are c contiguous
    Coords = np.ascontiguousarray(Coords, dtype=np.float64)
    Element_U = np.ascontiguousarray(Element_U, dtype=np.float64)
    S_vec = np.ascontiguousarray(S_vec, dtype=np.float64)
    PENER = np.ascontiguousarray(PENER, dtype=np.float64)
    SENER = np.ascontiguousarray(SENER, dtype=np.float64)
    ################################################

    num_elem_ = ctypes.c_size_t(Coords.shape[0])
    Coords_p = Coords.ctypes.data_as(ctypes.c_void_p)
    Element_U_p = Element_U.ctypes.data_as(ctypes.c_void_p)
    S_vec_p = S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p = PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p = SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force = np.zeros((Coords.shape[0], Coords.shape[1], 3))
    Conf_Force_p = Conf_Force.ctypes.data_as(ctypes.c_void_p)

    if method == 'mbf':
        lib.Integration_CPE8_static_mbf(num_elem_, Coords_p, Element_U_p, S_vec_p, PENER_p, SENER_p, Conf_Force_p)
    elif method == 'dbf':
        lib.Integration_CPE8_static_dbf(num_elem_, Coords_p, Element_U_p, S_vec_p, PENER_p, SENER_p, Conf_Force_p)
    else:
        print("Method not found!")
    return Conf_Force


def calc_Conf_Force_CPE8_dynamic(Coords, rho, Element_U, Element_V, Element_A, S_vec, PENER, SENER):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates    [num_elem,numNodesPerElem,3]
        rho............Density                      [num_elem]
        Element_U......Element_Nodal_Displacements  [num_elem,numNodesPerElem,3]
        Element_V......Element_Velocity             [num_elem,numNodesPerElem,3]
        Element_A......Element_Acceleration         [num_elem,numNodesPerElem,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy        [num_elem,num_int_points]
        SENER..........Elastic strain energy        [num_elem,num_int_points]
    """
    # This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem = Coords.shape[0]
    numNodes = 8
    num_int_points = 9

    assert Coords.shape == (num_elem, numNodes, 3)
    assert rho.shape == (num_elem,)
    assert Element_U.shape == (num_elem, numNodes, 3)
    assert Element_V.shape == (num_elem, numNodes, 3)
    assert Element_A.shape == (num_elem, numNodes, 3)
    assert S_vec.shape == (num_elem, num_int_points, 6)
    assert PENER.shape == (num_elem, num_int_points)
    assert SENER.shape == (num_elem, num_int_points)

    # Ensure that all arrays are c contiguous
    Coords = np.ascontiguousarray(Coords, dtype=np.float64)
    rho = np.ascontiguousarray(rho, dtype=np.float64)
    Element_U = np.ascontiguousarray(Element_U, dtype=np.float64)
    Element_V = np.ascontiguousarray(Element_V, dtype=np.float64)
    Element_A = np.ascontiguousarray(Element_A, dtype=np.float64)
    S_vec = np.ascontiguousarray(S_vec, dtype=np.float64)
    PENER = np.ascontiguousarray(PENER, dtype=np.float64)
    SENER = np.ascontiguousarray(SENER, dtype=np.float64)
    ################################################

    num_elem_ = ctypes.c_size_t(Coords.shape[0])
    Coords_p = Coords.ctypes.data_as(ctypes.c_void_p)
    rho_p = rho.ctypes.data_as(ctypes.c_void_p)
    Element_U_p = Element_U.ctypes.data_as(ctypes.c_void_p)
    Element_V_p = Element_V.ctypes.data_as(ctypes.c_void_p)
    Element_A_p = Element_A.ctypes.data_as(ctypes.c_void_p)
    S_vec_p = S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p = PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p = SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force = np.zeros((Coords.shape[0], Coords.shape[1], 3))
    Conf_Force_p = Conf_Force.ctypes.data_as(ctypes.c_void_p)

    lib.Integration_CPE8_dynamic(num_elem_, Coords_p, rho_p, Element_U_p, Element_V_p, Element_A_p, S_vec_p, PENER_p,
                                 SENER_p, Conf_Force_p)
    return Conf_Force


def calc_Conf_Force_CPE8R_static(Coords, Element_U, S_vec, PENER, SENER, method='mbf'):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates [num_elem,numNodes,3]
        Element_U......Element_Nodal_Displacements [num_elem,numNodes,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy [num_elem,num_int_points]
        SENER..........Elastic strain energy [num_elem,num_int_points]
    """
    # This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem = Coords.shape[0]
    numNodes = 8
    num_int_points = 4

    assert Coords.shape == (num_elem, numNodes, 3)
    assert Element_U.shape == (num_elem, numNodes, 3)
    assert S_vec.shape == (num_elem, num_int_points, 6)
    assert PENER.shape == (num_elem, num_int_points)
    assert SENER.shape == (num_elem, num_int_points)

    # Ensure that all arrays are c contiguous
    Coords = np.ascontiguousarray(Coords, dtype=np.float64)
    Element_U = np.ascontiguousarray(Element_U, dtype=np.float64)
    S_vec = np.ascontiguousarray(S_vec, dtype=np.float64)
    PENER = np.ascontiguousarray(PENER, dtype=np.float64)
    SENER = np.ascontiguousarray(SENER, dtype=np.float64)
    ################################################

    num_elem_ = ctypes.c_size_t(Coords.shape[0])
    Coords_p = Coords.ctypes.data_as(ctypes.c_void_p)
    Element_U_p = Element_U.ctypes.data_as(ctypes.c_void_p)
    S_vec_p = S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p = PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p = SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force = np.zeros((Coords.shape[0], Coords.shape[1], 3))
    Conf_Force_p = Conf_Force.ctypes.data_as(ctypes.c_void_p)

    if method == 'mbf':
        lib.Integration_CPE8R_static_mbf(num_elem_, Coords_p, Element_U_p, S_vec_p, PENER_p, SENER_p, Conf_Force_p)
    elif method == 'dbf':
        lib.Integration_CPE8R_static_dbf(num_elem_, Coords_p, Element_U_p, S_vec_p, PENER_p, SENER_p, Conf_Force_p)
    else:
        print("Method not found!")
    return Conf_Force


def calc_Conf_Force_CPE8R_dynamic(Coords, rho, Element_U, Element_V, Element_A, S_vec, PENER, SENER):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates    [num_elem,numNodesPerElem,3]
        rho............Density                      [num_elem]
        Element_U......Element_Nodal_Displacements  [num_elem,numNodesPerElem,3]
        Element_V......Element_Velocity             [num_elem,numNodesPerElem,3]
        Element_A......Element_Acceleration         [num_elem,numNodesPerElem,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy        [num_elem,num_int_points]
        SENER..........Elastic strain energy        [num_elem,num_int_points]
    """
    # This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem = Coords.shape[0]
    numNodes = 8
    num_int_points = 4

    assert Coords.shape == (num_elem, numNodes, 3)
    assert rho.shape == (num_elem,)
    assert Element_U.shape == (num_elem, numNodes, 3)
    assert Element_V.shape == (num_elem, numNodes, 3)
    assert Element_A.shape == (num_elem, numNodes, 3)
    assert S_vec.shape == (num_elem, num_int_points, 6)
    assert PENER.shape == (num_elem, num_int_points)
    assert SENER.shape == (num_elem, num_int_points)

    # Ensure that all arrays are c contiguous
    Coords = np.ascontiguousarray(Coords, dtype=np.float64)
    rho = np.ascontiguousarray(rho, dtype=np.float64)
    Element_U = np.ascontiguousarray(Element_U, dtype=np.float64)
    Element_V = np.ascontiguousarray(Element_V, dtype=np.float64)
    Element_A = np.ascontiguousarray(Element_A, dtype=np.float64)
    S_vec = np.ascontiguousarray(S_vec, dtype=np.float64)
    PENER = np.ascontiguousarray(PENER, dtype=np.float64)
    SENER = np.ascontiguousarray(SENER, dtype=np.float64)
    ################################################

    num_elem_ = ctypes.c_size_t(Coords.shape[0])
    Coords_p = Coords.ctypes.data_as(ctypes.c_void_p)
    rho_p = rho.ctypes.data_as(ctypes.c_void_p)
    Element_U_p = Element_U.ctypes.data_as(ctypes.c_void_p)
    Element_V_p = Element_V.ctypes.data_as(ctypes.c_void_p)
    Element_A_p = Element_A.ctypes.data_as(ctypes.c_void_p)
    S_vec_p = S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p = PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p = SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force = np.zeros((Coords.shape[0], Coords.shape[1], 3))
    Conf_Force_p = Conf_Force.ctypes.data_as(ctypes.c_void_p)

    lib.Integration_CPE8R_dynamic(num_elem_, Coords_p, rho_p, Element_U_p, Element_V_p, Element_A_p, S_vec_p, PENER_p,
                                  SENER_p, Conf_Force_p)
    return Conf_Force


def calc_Conf_Force_CPE3_static(Coords, Element_U, S_vec, PENER, SENER, method='mbf'):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates [num_elem,numNodes,3]
        Element_U......Element_Nodal_Displacements [num_elem,numNodes,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy [num_elem,num_int_points]
        SENER..........Elastic strain energy [num_elem,num_int_points]
    """
    # This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem = Coords.shape[0]
    numNodes = 3
    num_int_points = 1

    assert Coords.shape == (num_elem, numNodes, 3)
    assert Element_U.shape == (num_elem, numNodes, 3)
    assert S_vec.shape == (num_elem, num_int_points, 6)
    assert PENER.shape == (num_elem, num_int_points)
    assert SENER.shape == (num_elem, num_int_points)

    # Ensure that all arrays are c contiguous
    Coords = np.ascontiguousarray(Coords, dtype=np.float64)
    Element_U = np.ascontiguousarray(Element_U, dtype=np.float64)
    S_vec = np.ascontiguousarray(S_vec, dtype=np.float64)
    PENER = np.ascontiguousarray(PENER, dtype=np.float64)
    SENER = np.ascontiguousarray(SENER, dtype=np.float64)
    ################################################

    num_elem_ = ctypes.c_size_t(Coords.shape[0])
    Coords_p = Coords.ctypes.data_as(ctypes.c_void_p)
    Element_U_p = Element_U.ctypes.data_as(ctypes.c_void_p)
    S_vec_p = S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p = PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p = SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force = np.zeros((Coords.shape[0], Coords.shape[1], 3))
    Conf_Force_p = Conf_Force.ctypes.data_as(ctypes.c_void_p)

    if method == 'mbf':
        lib.Integration_CPE3_static_mbf(num_elem_, Coords_p, Element_U_p, S_vec_p, PENER_p, SENER_p, Conf_Force_p)
    elif method == 'dbf':
        lib.Integration_CPE3_static_dbf(num_elem_, Coords_p, Element_U_p, S_vec_p, PENER_p, SENER_p, Conf_Force_p)
    else:
        print("Method not found!")
    return Conf_Force


def calc_Conf_Force_CPE3_dynamic(Coords, rho, Element_U, Element_V, Element_A, S_vec, PENER, SENER):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates    [num_elem,numNodesPerElem,3]
        rho............Density                      [num_elem]
        Element_U......Element_Nodal_Displacements  [num_elem,numNodesPerElem,3]
        Element_V......Element_Velocity             [num_elem,numNodesPerElem,3]
        Element_A......Element_Acceleration         [num_elem,numNodesPerElem,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy        [num_elem,num_int_points]
        SENER..........Elastic strain energy        [num_elem,num_int_points]
    """
    # This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem = Coords.shape[0]
    numNodes = 3
    num_int_points = 1

    assert Coords.shape == (num_elem, numNodes, 3)
    assert rho.shape == (num_elem,)
    assert Element_U.shape == (num_elem, numNodes, 3)
    assert Element_V.shape == (num_elem, numNodes, 3)
    assert Element_A.shape == (num_elem, numNodes, 3)
    assert S_vec.shape == (num_elem, num_int_points, 6)
    assert PENER.shape == (num_elem, num_int_points)
    assert SENER.shape == (num_elem, num_int_points)

    # Ensure that all arrays are c contiguous
    Coords = np.ascontiguousarray(Coords, dtype=np.float64)
    rho = np.ascontiguousarray(rho, dtype=np.float64)
    Element_U = np.ascontiguousarray(Element_U, dtype=np.float64)
    Element_V = np.ascontiguousarray(Element_V, dtype=np.float64)
    Element_A = np.ascontiguousarray(Element_A, dtype=np.float64)
    S_vec = np.ascontiguousarray(S_vec, dtype=np.float64)
    PENER = np.ascontiguousarray(PENER, dtype=np.float64)
    SENER = np.ascontiguousarray(SENER, dtype=np.float64)
    ################################################

    num_elem_ = ctypes.c_size_t(Coords.shape[0])
    Coords_p = Coords.ctypes.data_as(ctypes.c_void_p)
    rho_p = rho.ctypes.data_as(ctypes.c_void_p)
    Element_U_p = Element_U.ctypes.data_as(ctypes.c_void_p)
    Element_V_p = Element_V.ctypes.data_as(ctypes.c_void_p)
    Element_A_p = Element_A.ctypes.data_as(ctypes.c_void_p)
    S_vec_p = S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p = PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p = SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force = np.zeros((Coords.shape[0], Coords.shape[1], 3))
    Conf_Force_p = Conf_Force.ctypes.data_as(ctypes.c_void_p)

    lib.Integration_CPE3_dynamic(num_elem_, Coords_p, rho_p, Element_U_p, Element_V_p, Element_A_p, S_vec_p, PENER_p,
                                 SENER_p, Conf_Force_p)
    return Conf_Force


def calc_Conf_Force_CPE6_static(Coords, Element_U, S_vec, PENER, SENER, method='mbf'):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates [num_elem,numNodes,3]
        Element_U......Element_Nodal_Displacements [num_elem,numNodes,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy [num_elem,num_int_points]
        SENER..........Elastic strain energy [num_elem,num_int_points]
    """
    # This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem = Coords.shape[0]
    numNodes = 6
    num_int_points = 3

    assert Coords.shape == (num_elem, numNodes, 3)
    assert Element_U.shape == (num_elem, numNodes, 3)
    assert S_vec.shape == (num_elem, num_int_points, 6)
    assert PENER.shape == (num_elem, num_int_points)
    assert SENER.shape == (num_elem, num_int_points)

    # Ensure that all arrays are c contiguous
    Coords = np.ascontiguousarray(Coords, dtype=np.float64)
    Element_U = np.ascontiguousarray(Element_U, dtype=np.float64)
    S_vec = np.ascontiguousarray(S_vec, dtype=np.float64)
    PENER = np.ascontiguousarray(PENER, dtype=np.float64)
    SENER = np.ascontiguousarray(SENER, dtype=np.float64)
    ################################################

    num_elem_ = ctypes.c_size_t(Coords.shape[0])
    Coords_p = Coords.ctypes.data_as(ctypes.c_void_p)
    Element_U_p = Element_U.ctypes.data_as(ctypes.c_void_p)
    S_vec_p = S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p = PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p = SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force = np.zeros((Coords.shape[0], Coords.shape[1], 3))
    Conf_Force_p = Conf_Force.ctypes.data_as(ctypes.c_void_p)

    if method == 'mbf':
        lib.Integration_CPE6_static_mbf(num_elem_, Coords_p, Element_U_p, S_vec_p, PENER_p, SENER_p, Conf_Force_p)
    elif method == 'dbf':
        lib.Integration_CPE6_static_dbf(num_elem_, Coords_p, Element_U_p, S_vec_p, PENER_p, SENER_p, Conf_Force_p)
    else:
        print("Method not found!")
    return Conf_Force


def calc_Conf_Force_CPE6_dynamic(Coords, rho, Element_U, Element_V, Element_A, S_vec, PENER, SENER):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates    [num_elem,numNodesPerElem,3]
        rho............Density                      [num_elem]
        Element_U......Element_Nodal_Displacements  [num_elem,numNodesPerElem,3]
        Element_V......Element_Velocity             [num_elem,numNodesPerElem,3]
        Element_A......Element_Acceleration         [num_elem,numNodesPerElem,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy        [num_elem,num_int_points]
        SENER..........Elastic strain energy        [num_elem,num_int_points]
    """
    # This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem = Coords.shape[0]
    numNodes = 6
    num_int_points = 3

    assert Coords.shape == (num_elem, numNodes, 3)
    assert rho.shape == (num_elem,)
    assert Element_U.shape == (num_elem, numNodes, 3)
    assert Element_V.shape == (num_elem, numNodes, 3)
    assert Element_A.shape == (num_elem, numNodes, 3)
    assert S_vec.shape == (num_elem, num_int_points, 6)
    assert PENER.shape == (num_elem, num_int_points)
    assert SENER.shape == (num_elem, num_int_points)

    # Ensure that all arrays are c contiguous
    Coords = np.ascontiguousarray(Coords, dtype=np.float64)
    rho = np.ascontiguousarray(rho, dtype=np.float64)
    Element_U = np.ascontiguousarray(Element_U, dtype=np.float64)
    Element_V = np.ascontiguousarray(Element_V, dtype=np.float64)
    Element_A = np.ascontiguousarray(Element_A, dtype=np.float64)
    S_vec = np.ascontiguousarray(S_vec, dtype=np.float64)
    PENER = np.ascontiguousarray(PENER, dtype=np.float64)
    SENER = np.ascontiguousarray(SENER, dtype=np.float64)
    ################################################

    num_elem_ = ctypes.c_size_t(Coords.shape[0])
    Coords_p = Coords.ctypes.data_as(ctypes.c_void_p)
    rho_p = rho.ctypes.data_as(ctypes.c_void_p)
    Element_U_p = Element_U.ctypes.data_as(ctypes.c_void_p)
    Element_V_p = Element_V.ctypes.data_as(ctypes.c_void_p)
    Element_A_p = Element_A.ctypes.data_as(ctypes.c_void_p)
    S_vec_p = S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p = PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p = SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force = np.zeros((Coords.shape[0], Coords.shape[1], 3))
    Conf_Force_p = Conf_Force.ctypes.data_as(ctypes.c_void_p)

    lib.Integration_CPE6_dynamic(num_elem_, Coords_p, rho_p, Element_U_p, Element_V_p, Element_A_p, S_vec_p, PENER_p,
                                 SENER_p, Conf_Force_p)
    return Conf_Force


def calc_Conf_Force_C3D8_static(Coords, Element_U, S_vec, PENER, SENER, method='mbf'):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates [num_elem,numNodes,3]
        Element_U......Element_Nodal_Displacements [num_elem,numNodes,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy [num_elem,num_int_points]
        SENER..........Elastic strain energy [num_elem,num_int_points]
    """
    # This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem = Coords.shape[0]
    numNodes = 8
    num_int_points = 8

    assert Coords.shape == (num_elem, numNodes, 3)
    assert Element_U.shape == (num_elem, numNodes, 3)
    assert S_vec.shape == (num_elem, num_int_points, 6)
    assert PENER.shape == (num_elem, num_int_points)
    assert SENER.shape == (num_elem, num_int_points)

    # Ensure that all arrays are c contiguous
    Coords = np.ascontiguousarray(Coords, dtype=np.float64)
    Element_U = np.ascontiguousarray(Element_U, dtype=np.float64)
    S_vec = np.ascontiguousarray(S_vec, dtype=np.float64)
    PENER = np.ascontiguousarray(PENER, dtype=np.float64)
    SENER = np.ascontiguousarray(SENER, dtype=np.float64)
    ################################################

    num_elem_ = ctypes.c_size_t(Coords.shape[0])
    Coords_p = Coords.ctypes.data_as(ctypes.c_void_p)
    Element_U_p = Element_U.ctypes.data_as(ctypes.c_void_p)
    S_vec_p = S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p = PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p = SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force = np.zeros((Coords.shape[0], Coords.shape[1], 3))
    Conf_Force_p = Conf_Force.ctypes.data_as(ctypes.c_void_p)

    if method == 'mbf':
        lib.Integration_C3D8_static_mbf(num_elem_, Coords_p, Element_U_p, S_vec_p, PENER_p, SENER_p, Conf_Force_p)
    elif method == 'dbf':
        lib.Integration_C3D8_static_dbf(num_elem_, Coords_p, Element_U_p, S_vec_p, PENER_p, SENER_p, Conf_Force_p)
    else:
        print("Method not found!")
    return Conf_Force


def calc_Conf_Force_C3D8_dynamic(Coords, rho, Element_U, Element_V, Element_A, S_vec, PENER, SENER):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates    [num_elem,numNodesPerElem,3]
        rho............Density                      [num_elem]
        Element_U......Element_Nodal_Displacements  [num_elem,numNodesPerElem,3]
        Element_V......Element_Velocity             [num_elem,numNodesPerElem,3]
        Element_A......Element_Acceleration         [num_elem,numNodesPerElem,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy        [num_elem,num_int_points]
        SENER..........Elastic strain energy        [num_elem,num_int_points]
    """
    # This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem = Coords.shape[0]
    numNodes = 8
    num_int_points = 8

    assert Coords.shape == (num_elem, numNodes, 3)
    assert rho.shape == (num_elem,)
    assert Element_U.shape == (num_elem, numNodes, 3)
    assert Element_V.shape == (num_elem, numNodes, 3)
    assert Element_A.shape == (num_elem, numNodes, 3)
    assert S_vec.shape == (num_elem, num_int_points, 6)
    assert PENER.shape == (num_elem, num_int_points)
    assert SENER.shape == (num_elem, num_int_points)

    # Ensure that all arrays are c contiguous
    Coords = np.ascontiguousarray(Coords, dtype=np.float64)
    rho = np.ascontiguousarray(rho, dtype=np.float64)
    Element_U = np.ascontiguousarray(Element_U, dtype=np.float64)
    Element_V = np.ascontiguousarray(Element_V, dtype=np.float64)
    Element_A = np.ascontiguousarray(Element_A, dtype=np.float64)
    S_vec = np.ascontiguousarray(S_vec, dtype=np.float64)
    PENER = np.ascontiguousarray(PENER, dtype=np.float64)
    SENER = np.ascontiguousarray(SENER, dtype=np.float64)
    ################################################

    num_elem_ = ctypes.c_size_t(Coords.shape[0])
    Coords_p = Coords.ctypes.data_as(ctypes.c_void_p)
    rho_p = rho.ctypes.data_as(ctypes.c_void_p)
    Element_U_p = Element_U.ctypes.data_as(ctypes.c_void_p)
    Element_V_p = Element_V.ctypes.data_as(ctypes.c_void_p)
    Element_A_p = Element_A.ctypes.data_as(ctypes.c_void_p)
    S_vec_p = S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p = PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p = SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force = np.zeros((Coords.shape[0], Coords.shape[1], 3))
    Conf_Force_p = Conf_Force.ctypes.data_as(ctypes.c_void_p)

    lib.Integration_C3D8_dynamic(num_elem_, Coords_p, rho_p, Element_U_p, Element_V_p, Element_A_p, S_vec_p, PENER_p,
                                 SENER_p, Conf_Force_p)
    return Conf_Force


def calc_Conf_Force_C3D8R_static(Coords, Element_U, S_vec, PENER, SENER, method='mbf'):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates [num_elem,numNodes,3]
        Element_U......Element_Nodal_Displacements [num_elem,numNodes,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy [num_elem,num_int_points]
        SENER..........Elastic strain energy [num_elem,num_int_points]
    """
    # This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem = Coords.shape[0]
    numNodes = 8
    num_int_points = 1

    assert Coords.shape == (num_elem, numNodes, 3)
    assert Element_U.shape == (num_elem, numNodes, 3)
    assert S_vec.shape == (num_elem, num_int_points, 6)
    assert PENER.shape == (num_elem, num_int_points)
    assert SENER.shape == (num_elem, num_int_points)

    # Ensure that all arrays are c contiguous
    Coords = np.ascontiguousarray(Coords, dtype=np.float64)
    Element_U = np.ascontiguousarray(Element_U, dtype=np.float64)
    S_vec = np.ascontiguousarray(S_vec, dtype=np.float64)
    PENER = np.ascontiguousarray(PENER, dtype=np.float64)
    SENER = np.ascontiguousarray(SENER, dtype=np.float64)
    ################################################

    num_elem_ = ctypes.c_size_t(Coords.shape[0])
    Coords_p = Coords.ctypes.data_as(ctypes.c_void_p)
    Element_U_p = Element_U.ctypes.data_as(ctypes.c_void_p)
    S_vec_p = S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p = PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p = SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force = np.zeros((Coords.shape[0], Coords.shape[1], 3))
    Conf_Force_p = Conf_Force.ctypes.data_as(ctypes.c_void_p)

    if method == 'mbf':
        lib.Integration_C3D8R_static_mbf(num_elem_, Coords_p, Element_U_p, S_vec_p, PENER_p, SENER_p, Conf_Force_p)
    elif method == 'dbf':
        lib.Integration_C3D8R_static_dbf(num_elem_, Coords_p, Element_U_p, S_vec_p, PENER_p, SENER_p, Conf_Force_p)
    else:
        print("Method not found!")
    return Conf_Force


def calc_Conf_Force_C3D8R_dynamic(Coords, rho, Element_U, Element_V, Element_A, S_vec, PENER, SENER):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates    [num_elem,numNodesPerElem,3]
        rho............Density                      [num_elem]
        Element_U......Element_Nodal_Displacements  [num_elem,numNodesPerElem,3]
        Element_V......Element_Velocity             [num_elem,numNodesPerElem,3]
        Element_A......Element_Acceleration         [num_elem,numNodesPerElem,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy        [num_elem,num_int_points]
        SENER..........Elastic strain energy        [num_elem,num_int_points]
    """
    # This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem = Coords.shape[0]
    numNodes = 8
    num_int_points = 1

    assert Coords.shape == (num_elem, numNodes, 3)
    assert rho.shape == (num_elem,)
    assert Element_U.shape == (num_elem, numNodes, 3)
    assert Element_V.shape == (num_elem, numNodes, 3)
    assert Element_A.shape == (num_elem, numNodes, 3)
    assert S_vec.shape == (num_elem, num_int_points, 6)
    assert PENER.shape == (num_elem, num_int_points)
    assert SENER.shape == (num_elem, num_int_points)

    # Ensure that all arrays are c contiguous
    Coords = np.ascontiguousarray(Coords, dtype=np.float64)
    rho = np.ascontiguousarray(rho, dtype=np.float64)
    Element_U = np.ascontiguousarray(Element_U, dtype=np.float64)
    Element_V = np.ascontiguousarray(Element_V, dtype=np.float64)
    Element_A = np.ascontiguousarray(Element_A, dtype=np.float64)
    S_vec = np.ascontiguousarray(S_vec, dtype=np.float64)
    PENER = np.ascontiguousarray(PENER, dtype=np.float64)
    SENER = np.ascontiguousarray(SENER, dtype=np.float64)
    ################################################

    num_elem_ = ctypes.c_size_t(Coords.shape[0])
    Coords_p = Coords.ctypes.data_as(ctypes.c_void_p)
    rho_p = rho.ctypes.data_as(ctypes.c_void_p)
    Element_U_p = Element_U.ctypes.data_as(ctypes.c_void_p)
    Element_V_p = Element_V.ctypes.data_as(ctypes.c_void_p)
    Element_A_p = Element_A.ctypes.data_as(ctypes.c_void_p)
    S_vec_p = S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p = PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p = SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force = np.zeros((Coords.shape[0], Coords.shape[1], 3))
    Conf_Force_p = Conf_Force.ctypes.data_as(ctypes.c_void_p)

    lib.Integration_C3D8R_dynamic(num_elem_, Coords_p, rho_p, Element_U_p, Element_V_p, Element_A_p, S_vec_p, PENER_p,
                                  SENER_p, Conf_Force_p)
    return Conf_Force


def calc_Conf_Force_C3D20_static(Coords, Element_U, S_vec, PENER, SENER, method='mbf'):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates [num_elem,numNodes,3]
        Element_U......Element_Nodal_Displacements [num_elem,numNodes,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy [num_elem,num_int_points]
        SENER..........Elastic strain energy [num_elem,num_int_points]
    """
    # This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem = Coords.shape[0]
    numNodes = 20
    num_int_points = 27

    assert Coords.shape == (num_elem, numNodes, 3)
    assert Element_U.shape == (num_elem, numNodes, 3)
    assert S_vec.shape == (num_elem, num_int_points, 6)
    assert PENER.shape == (num_elem, num_int_points)
    assert SENER.shape == (num_elem, num_int_points)

    # Ensure that all arrays are c contiguous
    Coords = np.ascontiguousarray(Coords, dtype=np.float64)
    Element_U = np.ascontiguousarray(Element_U, dtype=np.float64)
    S_vec = np.ascontiguousarray(S_vec, dtype=np.float64)
    PENER = np.ascontiguousarray(PENER, dtype=np.float64)
    SENER = np.ascontiguousarray(SENER, dtype=np.float64)
    ################################################

    num_elem_ = ctypes.c_size_t(Coords.shape[0])
    Coords_p = Coords.ctypes.data_as(ctypes.c_void_p)
    Element_U_p = Element_U.ctypes.data_as(ctypes.c_void_p)
    S_vec_p = S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p = PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p = SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force = np.zeros((Coords.shape[0], Coords.shape[1], 3))
    Conf_Force_p = Conf_Force.ctypes.data_as(ctypes.c_void_p)

    if method == 'mbf':
        lib.Integration_C3D20_static_mbf(num_elem_, Coords_p, Element_U_p, S_vec_p, PENER_p, SENER_p, Conf_Force_p)
    elif method == 'dbf':
        lib.Integration_C3D20_static_dbf(num_elem_, Coords_p, Element_U_p, S_vec_p, PENER_p, SENER_p, Conf_Force_p)
    else:
        print("Method not found!")
    return Conf_Force


def calc_Conf_Force_C3D20_dynamic(Coords, rho, Element_U, Element_V, Element_A, S_vec, PENER, SENER):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates    [num_elem,numNodesPerElem,3]
        rho............Density                      [num_elem]
        Element_U......Element_Nodal_Displacements  [num_elem,numNodesPerElem,3]
        Element_V......Element_Velocity             [num_elem,numNodesPerElem,3]
        Element_A......Element_Acceleration         [num_elem,numNodesPerElem,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy        [num_elem,num_int_points]
        SENER..........Elastic strain energy        [num_elem,num_int_points]
    """
    # This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem = Coords.shape[0]
    numNodes = 20
    num_int_points = 27

    assert Coords.shape == (num_elem, numNodes, 3)
    assert rho.shape == (num_elem,)
    assert Element_U.shape == (num_elem, numNodes, 3)
    assert Element_V.shape == (num_elem, numNodes, 3)
    assert Element_A.shape == (num_elem, numNodes, 3)
    assert S_vec.shape == (num_elem, num_int_points, 6)
    assert PENER.shape == (num_elem, num_int_points)
    assert SENER.shape == (num_elem, num_int_points)

    # Ensure that all arrays are c contiguous
    Coords = np.ascontiguousarray(Coords, dtype=np.float64)
    rho = np.ascontiguousarray(rho, dtype=np.float64)
    Element_U = np.ascontiguousarray(Element_U, dtype=np.float64)
    Element_V = np.ascontiguousarray(Element_V, dtype=np.float64)
    Element_A = np.ascontiguousarray(Element_A, dtype=np.float64)
    S_vec = np.ascontiguousarray(S_vec, dtype=np.float64)
    PENER = np.ascontiguousarray(PENER, dtype=np.float64)
    SENER = np.ascontiguousarray(SENER, dtype=np.float64)
    ################################################

    num_elem_ = ctypes.c_size_t(Coords.shape[0])
    Coords_p = Coords.ctypes.data_as(ctypes.c_void_p)
    rho_p = rho.ctypes.data_as(ctypes.c_void_p)
    Element_U_p = Element_U.ctypes.data_as(ctypes.c_void_p)
    Element_V_p = Element_V.ctypes.data_as(ctypes.c_void_p)
    Element_A_p = Element_A.ctypes.data_as(ctypes.c_void_p)
    S_vec_p = S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p = PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p = SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force = np.zeros((Coords.shape[0], Coords.shape[1], 3))
    Conf_Force_p = Conf_Force.ctypes.data_as(ctypes.c_void_p)

    lib.Integration_C3D20_dynamic(num_elem_, Coords_p, rho_p, Element_U_p, Element_V_p, Element_A_p, S_vec_p, PENER_p,
                                  SENER_p, Conf_Force_p)
    return Conf_Force


def calc_Conf_Force_C3D20R_static(Coords, Element_U, S_vec, PENER, SENER, method='mbf'):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates [num_elem,numNodes,3]
        Element_U......Element_Nodal_Displacements [num_elem,numNodes,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy [num_elem,num_int_points]
        SENER..........Elastic strain energy [num_elem,num_int_points]
    """
    # This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem = Coords.shape[0]
    numNodes = 20
    num_int_points = 8

    assert Coords.shape == (num_elem, numNodes, 3)
    assert Element_U.shape == (num_elem, numNodes, 3)
    assert S_vec.shape == (num_elem, num_int_points, 6)
    assert PENER.shape == (num_elem, num_int_points)
    assert SENER.shape == (num_elem, num_int_points)

    # Ensure that all arrays are c contiguous
    Coords = np.ascontiguousarray(Coords, dtype=np.float64)
    Element_U = np.ascontiguousarray(Element_U, dtype=np.float64)
    S_vec = np.ascontiguousarray(S_vec, dtype=np.float64)
    PENER = np.ascontiguousarray(PENER, dtype=np.float64)
    SENER = np.ascontiguousarray(SENER, dtype=np.float64)
    ################################################

    num_elem_ = ctypes.c_size_t(Coords.shape[0])
    Coords_p = Coords.ctypes.data_as(ctypes.c_void_p)
    Element_U_p = Element_U.ctypes.data_as(ctypes.c_void_p)
    S_vec_p = S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p = PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p = SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force = np.zeros((Coords.shape[0], Coords.shape[1], 3))
    Conf_Force_p = Conf_Force.ctypes.data_as(ctypes.c_void_p)

    if method == 'mbf':
        lib.Integration_C3D20R_static_mbf(num_elem_, Coords_p, Element_U_p, S_vec_p, PENER_p, SENER_p, Conf_Force_p)
    elif method == 'dbf':
        lib.Integration_C3D20R_static_dbf(num_elem_, Coords_p, Element_U_p, S_vec_p, PENER_p, SENER_p, Conf_Force_p)
    else:
        print("Method not found!")
    return Conf_Force


def calc_Conf_Force_C3D20R_dynamic(Coords, rho, Element_U, Element_V, Element_A, S_vec, PENER, SENER):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates    [num_elem,numNodesPerElem,3]
        rho............Density                      [num_elem]
        Element_U......Element_Nodal_Displacements  [num_elem,numNodesPerElem,3]
        Element_V......Element_Velocity             [num_elem,numNodesPerElem,3]
        Element_A......Element_Acceleration         [num_elem,numNodesPerElem,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy        [num_elem,num_int_points]
        SENER..........Elastic strain energy        [num_elem,num_int_points]
    """
    # This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem = Coords.shape[0]
    numNodes = 20
    num_int_points = 8

    assert Coords.shape == (num_elem, numNodes, 3)
    assert rho.shape == (num_elem,)
    assert Element_U.shape == (num_elem, numNodes, 3)
    assert Element_V.shape == (num_elem, numNodes, 3)
    assert Element_A.shape == (num_elem, numNodes, 3)
    assert S_vec.shape == (num_elem, num_int_points, 6)
    assert PENER.shape == (num_elem, num_int_points)
    assert SENER.shape == (num_elem, num_int_points)

    # Ensure that all arrays are c contiguous
    Coords = np.ascontiguousarray(Coords, dtype=np.float64)
    rho = np.ascontiguousarray(rho, dtype=np.float64)
    Element_U = np.ascontiguousarray(Element_U, dtype=np.float64)
    Element_V = np.ascontiguousarray(Element_V, dtype=np.float64)
    Element_A = np.ascontiguousarray(Element_A, dtype=np.float64)
    S_vec = np.ascontiguousarray(S_vec, dtype=np.float64)
    PENER = np.ascontiguousarray(PENER, dtype=np.float64)
    SENER = np.ascontiguousarray(SENER, dtype=np.float64)
    ################################################

    num_elem_ = ctypes.c_size_t(Coords.shape[0])
    Coords_p = Coords.ctypes.data_as(ctypes.c_void_p)
    rho_p = rho.ctypes.data_as(ctypes.c_void_p)
    Element_U_p = Element_U.ctypes.data_as(ctypes.c_void_p)
    Element_V_p = Element_V.ctypes.data_as(ctypes.c_void_p)
    Element_A_p = Element_A.ctypes.data_as(ctypes.c_void_p)
    S_vec_p = S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p = PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p = SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force = np.zeros((Coords.shape[0], Coords.shape[1], 3))
    Conf_Force_p = Conf_Force.ctypes.data_as(ctypes.c_void_p)

    lib.Integration_C3D20R_dynamic(num_elem_, Coords_p, rho_p, Element_U_p, Element_V_p, Element_A_p, S_vec_p, PENER_p,
                                   SENER_p, Conf_Force_p)
    return Conf_Force


def calc_Conf_Force_C3D4_static(Coords, Element_U, S_vec, PENER, SENER, method='mbf'):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates [num_elem,numNodes,3]
        Element_U......Element_Nodal_Displacements [num_elem,numNodes,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy [num_elem,num_int_points]
        SENER..........Elastic strain energy [num_elem,num_int_points]
    """
    # This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem = Coords.shape[0]
    numNodes = 4
    num_int_points = 1

    assert Coords.shape == (num_elem, numNodes, 3)
    assert Element_U.shape == (num_elem, numNodes, 3)
    assert S_vec.shape == (num_elem, num_int_points, 6)
    assert PENER.shape == (num_elem, num_int_points)
    assert SENER.shape == (num_elem, num_int_points)

    # Ensure that all arrays are c contiguous
    Coords = np.ascontiguousarray(Coords, dtype=np.float64)
    Element_U = np.ascontiguousarray(Element_U, dtype=np.float64)
    S_vec = np.ascontiguousarray(S_vec, dtype=np.float64)
    PENER = np.ascontiguousarray(PENER, dtype=np.float64)
    SENER = np.ascontiguousarray(SENER, dtype=np.float64)
    ################################################

    num_elem_ = ctypes.c_size_t(Coords.shape[0])
    Coords_p = Coords.ctypes.data_as(ctypes.c_void_p)
    Element_U_p = Element_U.ctypes.data_as(ctypes.c_void_p)
    S_vec_p = S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p = PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p = SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force = np.zeros((Coords.shape[0], Coords.shape[1], 3))
    Conf_Force_p = Conf_Force.ctypes.data_as(ctypes.c_void_p)

    if method == 'mbf':
        lib.Integration_C3D4_static_mbf(num_elem_, Coords_p, Element_U_p, S_vec_p, PENER_p, SENER_p, Conf_Force_p)
    elif method == 'dbf':
        lib.Integration_C3D4_static_dbf(num_elem_, Coords_p, Element_U_p, S_vec_p, PENER_p, SENER_p, Conf_Force_p)
    else:
        print("Method not found!")
    return Conf_Force


def calc_Conf_Force_C3D4_dynamic(Coords, rho, Element_U, Element_V, Element_A, S_vec, PENER, SENER):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates    [num_elem,numNodesPerElem,3]
        rho............Density                      [num_elem]
        Element_U......Element_Nodal_Displacements  [num_elem,numNodesPerElem,3]
        Element_V......Element_Velocity             [num_elem,numNodesPerElem,3]
        Element_A......Element_Acceleration         [num_elem,numNodesPerElem,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy        [num_elem,num_int_points]
        SENER..........Elastic strain energy        [num_elem,num_int_points]
    """
    # This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem = Coords.shape[0]
    numNodes = 4
    num_int_points = 1

    assert Coords.shape == (num_elem, numNodes, 3)
    assert rho.shape == (num_elem,)
    assert Element_U.shape == (num_elem, numNodes, 3)
    assert Element_V.shape == (num_elem, numNodes, 3)
    assert Element_A.shape == (num_elem, numNodes, 3)
    assert S_vec.shape == (num_elem, num_int_points, 6)
    assert PENER.shape == (num_elem, num_int_points)
    assert SENER.shape == (num_elem, num_int_points)

    # Ensure that all arrays are c contiguous
    Coords = np.ascontiguousarray(Coords, dtype=np.float64)
    rho = np.ascontiguousarray(rho, dtype=np.float64)
    Element_U = np.ascontiguousarray(Element_U, dtype=np.float64)
    Element_V = np.ascontiguousarray(Element_V, dtype=np.float64)
    Element_A = np.ascontiguousarray(Element_A, dtype=np.float64)
    S_vec = np.ascontiguousarray(S_vec, dtype=np.float64)
    PENER = np.ascontiguousarray(PENER, dtype=np.float64)
    SENER = np.ascontiguousarray(SENER, dtype=np.float64)
    ################################################

    num_elem_ = ctypes.c_size_t(Coords.shape[0])
    Coords_p = Coords.ctypes.data_as(ctypes.c_void_p)
    rho_p = rho.ctypes.data_as(ctypes.c_void_p)
    Element_U_p = Element_U.ctypes.data_as(ctypes.c_void_p)
    Element_V_p = Element_V.ctypes.data_as(ctypes.c_void_p)
    Element_A_p = Element_A.ctypes.data_as(ctypes.c_void_p)
    S_vec_p = S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p = PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p = SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force = np.zeros((Coords.shape[0], Coords.shape[1], 3))
    Conf_Force_p = Conf_Force.ctypes.data_as(ctypes.c_void_p)

    lib.Integration_C3D4_dynamic(num_elem_, Coords_p, rho_p, Element_U_p, Element_V_p, Element_A_p, S_vec_p, PENER_p,
                                 SENER_p, Conf_Force_p)
    return Conf_Force


def calc_Conf_Force_C3D10_static(Coords, Element_U, S_vec, PENER, SENER, method='mbf'):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates [num_elem,numNodes,3]
        Element_U......Element_Nodal_Displacements [num_elem,numNodes,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy [num_elem,num_int_points]
        SENER..........Elastic strain energy [num_elem,num_int_points]
    """
    # This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem = Coords.shape[0]
    numNodes = 10
    num_int_points = 4

    assert Coords.shape == (num_elem, numNodes, 3)
    assert Element_U.shape == (num_elem, numNodes, 3)
    assert S_vec.shape == (num_elem, num_int_points, 6)
    assert PENER.shape == (num_elem, num_int_points)
    assert SENER.shape == (num_elem, num_int_points)

    # Ensure that all arrays are c contiguous
    Coords = np.ascontiguousarray(Coords, dtype=np.float64)
    Element_U = np.ascontiguousarray(Element_U, dtype=np.float64)
    S_vec = np.ascontiguousarray(S_vec, dtype=np.float64)
    PENER = np.ascontiguousarray(PENER, dtype=np.float64)
    SENER = np.ascontiguousarray(SENER, dtype=np.float64)
    ################################################

    num_elem_ = ctypes.c_size_t(Coords.shape[0])
    Coords_p = Coords.ctypes.data_as(ctypes.c_void_p)
    Element_U_p = Element_U.ctypes.data_as(ctypes.c_void_p)
    S_vec_p = S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p = PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p = SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force = np.zeros((Coords.shape[0], Coords.shape[1], 3))
    Conf_Force_p = Conf_Force.ctypes.data_as(ctypes.c_void_p)

    if method == 'mbf':
        lib.Integration_C3D10_static_mbf(num_elem_, Coords_p, Element_U_p, S_vec_p, PENER_p, SENER_p, Conf_Force_p)
    elif method == 'dbf':
        lib.Integration_C3D10_static_dbf(num_elem_, Coords_p, Element_U_p, S_vec_p, PENER_p, SENER_p, Conf_Force_p)
    else:
        print("Method not found!")
    return Conf_Force


def calc_Conf_Force_C3D10_dynamic(Coords, rho, Element_U, Element_V, Element_A, S_vec, PENER, SENER):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates    [num_elem,numNodesPerElem,3]
        rho............Density                      [num_elem]
        Element_U......Element_Nodal_Displacements  [num_elem,numNodesPerElem,3]
        Element_V......Element_Velocity             [num_elem,numNodesPerElem,3]
        Element_A......Element_Acceleration         [num_elem,numNodesPerElem,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy        [num_elem,num_int_points]
        SENER..........Elastic strain energy        [num_elem,num_int_points]
    """
    # This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem = Coords.shape[0]
    numNodes = 10
    num_int_points = 4

    assert Coords.shape == (num_elem, numNodes, 3)
    assert rho.shape == (num_elem,)
    assert Element_U.shape == (num_elem, numNodes, 3)
    assert Element_V.shape == (num_elem, numNodes, 3)
    assert Element_A.shape == (num_elem, numNodes, 3)
    assert S_vec.shape == (num_elem, num_int_points, 6)
    assert PENER.shape == (num_elem, num_int_points)
    assert SENER.shape == (num_elem, num_int_points)

    # Ensure that all arrays are c contiguous
    Coords = np.ascontiguousarray(Coords, dtype=np.float64)
    rho = np.ascontiguousarray(rho, dtype=np.float64)
    Element_U = np.ascontiguousarray(Element_U, dtype=np.float64)
    Element_V = np.ascontiguousarray(Element_V, dtype=np.float64)
    Element_A = np.ascontiguousarray(Element_A, dtype=np.float64)
    S_vec = np.ascontiguousarray(S_vec, dtype=np.float64)
    PENER = np.ascontiguousarray(PENER, dtype=np.float64)
    SENER = np.ascontiguousarray(SENER, dtype=np.float64)
    ################################################

    num_elem_ = ctypes.c_size_t(Coords.shape[0])
    Coords_p = Coords.ctypes.data_as(ctypes.c_void_p)
    rho_p = rho.ctypes.data_as(ctypes.c_void_p)
    Element_U_p = Element_U.ctypes.data_as(ctypes.c_void_p)
    Element_V_p = Element_V.ctypes.data_as(ctypes.c_void_p)
    Element_A_p = Element_A.ctypes.data_as(ctypes.c_void_p)
    S_vec_p = S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p = PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p = SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force = np.zeros((Coords.shape[0], Coords.shape[1], 3))
    Conf_Force_p = Conf_Force.ctypes.data_as(ctypes.c_void_p)

    lib.Integration_C3D10_dynamic(num_elem_, Coords_p, rho_p, Element_U_p, Element_V_p, Element_A_p, S_vec_p, PENER_p,
                                  SENER_p, Conf_Force_p)
    return Conf_Force


def calc_Conf_Force_C3D10R_static(Coords, Element_U, S_vec, PENER, SENER, method='mbf'):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates [num_elem,numNodes,3]
        Element_U......Element_Nodal_Displacements [num_elem,numNodes,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy [num_elem,num_int_points]
        SENER..........Elastic strain energy [num_elem,num_int_points]
    """
    # This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem = Coords.shape[0]
    numNodes = 10
    num_int_points = 1

    assert Coords.shape == (num_elem, numNodes, 3)
    assert Element_U.shape == (num_elem, numNodes, 3)
    assert S_vec.shape == (num_elem, num_int_points, 6)
    assert PENER.shape == (num_elem, num_int_points)
    assert SENER.shape == (num_elem, num_int_points)

    # Ensure that all arrays are c contiguous
    Coords = np.ascontiguousarray(Coords, dtype=np.float64)
    Element_U = np.ascontiguousarray(Element_U, dtype=np.float64)
    S_vec = np.ascontiguousarray(S_vec, dtype=np.float64)
    PENER = np.ascontiguousarray(PENER, dtype=np.float64)
    SENER = np.ascontiguousarray(SENER, dtype=np.float64)
    ################################################

    num_elem_ = ctypes.c_size_t(Coords.shape[0])
    Coords_p = Coords.ctypes.data_as(ctypes.c_void_p)
    Element_U_p = Element_U.ctypes.data_as(ctypes.c_void_p)
    S_vec_p = S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p = PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p = SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force = np.zeros((Coords.shape[0], Coords.shape[1], 3))
    Conf_Force_p = Conf_Force.ctypes.data_as(ctypes.c_void_p)

    if method == 'mbf':
        lib.Integration_C3D10R_static_mbf(num_elem_, Coords_p, Element_U_p, S_vec_p, PENER_p, SENER_p, Conf_Force_p)
    elif method == 'dbf':
        lib.Integration_C3D10R_static_dbf(num_elem_, Coords_p, Element_U_p, S_vec_p, PENER_p, SENER_p, Conf_Force_p)
    else:
        print("Method not found!")
    return Conf_Force


def calc_Conf_Force_C3D10R_dynamic(Coords, rho, Element_U, Element_V, Element_A, S_vec, PENER, SENER):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates    [num_elem,numNodesPerElem,3]
        rho............Density                      [num_elem]
        Element_U......Element_Nodal_Displacements  [num_elem,numNodesPerElem,3]
        Element_V......Element_Velocity             [num_elem,numNodesPerElem,3]
        Element_A......Element_Acceleration         [num_elem,numNodesPerElem,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy        [num_elem,num_int_points]
        SENER..........Elastic strain energy        [num_elem,num_int_points]
    """
    # This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem = Coords.shape[0]
    numNodes = 10
    num_int_points = 1

    assert Coords.shape == (num_elem, numNodes, 3)
    assert rho.shape == (num_elem,)
    assert Element_U.shape == (num_elem, numNodes, 3)
    assert Element_V.shape == (num_elem, numNodes, 3)
    assert Element_A.shape == (num_elem, numNodes, 3)
    assert S_vec.shape == (num_elem, num_int_points, 6)
    assert PENER.shape == (num_elem, num_int_points)
    assert SENER.shape == (num_elem, num_int_points)

    # Ensure that all arrays are c contiguous
    Coords = np.ascontiguousarray(Coords, dtype=np.float64)
    rho = np.ascontiguousarray(rho, dtype=np.float64)
    Element_U = np.ascontiguousarray(Element_U, dtype=np.float64)
    Element_V = np.ascontiguousarray(Element_V, dtype=np.float64)
    Element_A = np.ascontiguousarray(Element_A, dtype=np.float64)
    S_vec = np.ascontiguousarray(S_vec, dtype=np.float64)
    PENER = np.ascontiguousarray(PENER, dtype=np.float64)
    SENER = np.ascontiguousarray(SENER, dtype=np.float64)
    ################################################

    num_elem_ = ctypes.c_size_t(Coords.shape[0])
    Coords_p = Coords.ctypes.data_as(ctypes.c_void_p)
    rho_p = rho.ctypes.data_as(ctypes.c_void_p)
    Element_U_p = Element_U.ctypes.data_as(ctypes.c_void_p)
    Element_V_p = Element_V.ctypes.data_as(ctypes.c_void_p)
    Element_A_p = Element_A.ctypes.data_as(ctypes.c_void_p)
    S_vec_p = S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p = PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p = SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force = np.zeros((Coords.shape[0], Coords.shape[1], 3))
    Conf_Force_p = Conf_Force.ctypes.data_as(ctypes.c_void_p)

    lib.Integration_C3D10R_dynamic(num_elem_, Coords_p, rho_p, Element_U_p, Element_V_p, Element_A_p, S_vec_p, PENER_p,
                                   SENER_p, Conf_Force_p)
    return Conf_Force
