import numpy as np
import ctypes
import os

# Determine Operating System
if os.name == 'nt':
    lib = ctypes.cdll.LoadLibrary("Conf_Forces.dll")
else:
    lib = ctypes.cdll.LoadLibrary("./Conf_Forces.so")

def calc_Nodal(Element_Connectivity,CF):
    """
    Input:
    Element_Connectivity.......List of np.arrays or np.array
    CF_in......................List of np.arrays or np.array
    Output:
    Node_labels................Array of unique Nodelabels where CF is calcualted
    CF_out.....................Array of unique Nodal Configurational Forces [num_nodes,3]
    """
    if not isinstance(Element_Connectivity,list):
        Element_Connectivity=[Element_Connectivity]
        CF=[CF]

    Element_Connectivity=np.concatenate([elem_con.reshape(-1) for elem_con in Element_Connectivity])
    CF=np.concatenate([cf.reshape(-1,3) for cf in CF])
    assert Element_Connectivity.shape[0]==CF.shape[0]

    Node_labels,inverse=np.unique(Element_Connectivity,return_inverse=True)

    CF=np.ascontiguousarray(CF,dtype=np.float64)
    inverse=np.ascontiguousarray(inverse,dtype=np.int64)
    CF_out=np.zeros((Node_labels.shape[0],3),dtype=np.float64)

    num_elem_nodal=ctypes.c_size_t(inverse.shape[0])
    CF_out_p=CF_out.ctypes.data_as(ctypes.c_void_p)
    CF_p=CF.ctypes.data_as(ctypes.c_void_p)
    inverse_p=inverse.ctypes.data_as(ctypes.c_void_p)

    lib.calc_Nodal_C(num_elem_nodal,CF_out_p,CF_p,inverse_p)
    return Node_labels,CF_out

def calc_Conf_Force_CPE4_static(Coords,Element_U,S_vec,PENER,SENER,method='mbf'):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates [num_elem,numNodes,3]
        Element_U......Element_Nodal_Displacements [num_elem,numNodes,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy [num_elem,num_int_points]
        SENER..........Elastic strain energy [num_elem,num_int_points]
    """
    #This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem=Coords.shape[0]
    numNodes=4
    num_int_points=4

    assert Coords.shape==(num_elem,numNodes,3)
    assert Element_U.shape==(num_elem,numNodes,3)
    assert S_vec.shape==(num_elem,num_int_points,6)
    assert PENER.shape==(num_elem,num_int_points)
    assert SENER.shape==(num_elem,num_int_points)

    #Ensure that all arrays are c contiguous
    Coords=np.ascontiguousarray(Coords,dtype=np.float64)
    Element_U=np.ascontiguousarray(Element_U,dtype=np.float64)
    S_vec=np.ascontiguousarray(S_vec,dtype=np.float64)
    PENER=np.ascontiguousarray(PENER,dtype=np.float64)
    SENER=np.ascontiguousarray(SENER,dtype=np.float64)
    ################################################

    num_elem_=ctypes.c_size_t(Coords.shape[0])
    Coords_p=Coords.ctypes.data_as(ctypes.c_void_p)
    Element_U_p=Element_U.ctypes.data_as(ctypes.c_void_p)
    S_vec_p=S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p=PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p=SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force=np.zeros((Coords.shape[0],Coords.shape[1],3))
    Conf_Force_p=Conf_Force.ctypes.data_as(ctypes.c_void_p)

    if method=='mbf':
        lib.Integration_CPE4_static_mbf(num_elem_,Coords_p,Element_U_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    elif method=='dbf':
        lib.Integration_CPE4_static_dbf(num_elem_,Coords_p,Element_U_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    else:
        print("Method not found!")
    return Conf_Force

def calc_Conf_Force_CPE4_dynamic(Coords,rho,Element_U,Element_V,Element_A,S_vec,PENER,SENER):
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
    #This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem=Coords.shape[0]
    numNodes=4
    num_int_points=4

    assert Coords.shape==(num_elem,numNodes,3)
    assert rho.shape==(num_elem,)
    assert Element_U.shape==(num_elem,numNodes,3)
    assert Element_V.shape==(num_elem,numNodes,3)
    assert Element_A.shape==(num_elem,numNodes,3)
    assert S_vec.shape==(num_elem,num_int_points,6)
    assert PENER.shape==(num_elem,num_int_points)
    assert SENER.shape==(num_elem,num_int_points)

    #Ensure that all arrays are c contiguous
    Coords=np.ascontiguousarray(Coords,dtype=np.float64)
    rho=np.ascontiguousarray(rho,dtype=np.float64)
    Element_U=np.ascontiguousarray(Element_U,dtype=np.float64)
    Element_V=np.ascontiguousarray(Element_V,dtype=np.float64)
    Element_A=np.ascontiguousarray(Element_A,dtype=np.float64)
    S_vec=np.ascontiguousarray(S_vec,dtype=np.float64)
    PENER=np.ascontiguousarray(PENER,dtype=np.float64)
    SENER=np.ascontiguousarray(SENER,dtype=np.float64)
    ################################################

    num_elem_=ctypes.c_size_t(Coords.shape[0])
    Coords_p=Coords.ctypes.data_as(ctypes.c_void_p)
    rho_p   =rho.ctypes.data_as(ctypes.c_void_p)
    Element_U_p=Element_U.ctypes.data_as(ctypes.c_void_p)
    Element_V_p=Element_V.ctypes.data_as(ctypes.c_void_p)
    Element_A_p=Element_A.ctypes.data_as(ctypes.c_void_p)
    S_vec_p=S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p=PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p=SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force=np.zeros((Coords.shape[0],Coords.shape[1],3))
    Conf_Force_p=Conf_Force.ctypes.data_as(ctypes.c_void_p)

    lib.Integration_CPE4_dynamic(num_elem_,Coords_p,rho_p,Element_U_p,Element_V_p,Element_A_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    return Conf_Force

def calc_Conf_Force_CPE4R_static(Coords,Element_U,S_vec,PENER,SENER,method='mbf'):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates [num_elem,numNodes,3]
        Element_U......Element_Nodal_Displacements [num_elem,numNodes,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy [num_elem,num_int_points]
        SENER..........Elastic strain energy [num_elem,num_int_points]
    """
    #This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem=Coords.shape[0]
    numNodes=4
    num_int_points=1

    assert Coords.shape==(num_elem,numNodes,3)
    assert Element_U.shape==(num_elem,numNodes,3)
    assert S_vec.shape==(num_elem,num_int_points,6)
    assert PENER.shape==(num_elem,num_int_points)
    assert SENER.shape==(num_elem,num_int_points)

    #Ensure that all arrays are c contiguous
    Coords=np.ascontiguousarray(Coords,dtype=np.float64)
    Element_U=np.ascontiguousarray(Element_U,dtype=np.float64)
    S_vec=np.ascontiguousarray(S_vec,dtype=np.float64)
    PENER=np.ascontiguousarray(PENER,dtype=np.float64)
    SENER=np.ascontiguousarray(SENER,dtype=np.float64)
    ################################################

    num_elem_=ctypes.c_size_t(Coords.shape[0])
    Coords_p=Coords.ctypes.data_as(ctypes.c_void_p)
    Element_U_p=Element_U.ctypes.data_as(ctypes.c_void_p)
    S_vec_p=S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p=PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p=SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force=np.zeros((Coords.shape[0],Coords.shape[1],3))
    Conf_Force_p=Conf_Force.ctypes.data_as(ctypes.c_void_p)

    if method=='mbf':
        lib.Integration_CPE4R_static_mbf(num_elem_,Coords_p,Element_U_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    elif method=='dbf':
        lib.Integration_CPE4R_static_dbf(num_elem_,Coords_p,Element_U_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    else:
        print("Method not found!")
    return Conf_Force

def calc_Conf_Force_CPE4R_dynamic(Coords,rho,Element_U,Element_V,Element_A,S_vec,PENER,SENER):
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
    #This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem=Coords.shape[0]
    numNodes=4
    num_int_points=1

    assert Coords.shape==(num_elem,numNodes,3)
    assert rho.shape==(num_elem,)
    assert Element_U.shape==(num_elem,numNodes,3)
    assert Element_V.shape==(num_elem,numNodes,3)
    assert Element_A.shape==(num_elem,numNodes,3)
    assert S_vec.shape==(num_elem,num_int_points,6)
    assert PENER.shape==(num_elem,num_int_points)
    assert SENER.shape==(num_elem,num_int_points)

    #Ensure that all arrays are c contiguous
    Coords=np.ascontiguousarray(Coords,dtype=np.float64)
    rho=np.ascontiguousarray(rho,dtype=np.float64)
    Element_U=np.ascontiguousarray(Element_U,dtype=np.float64)
    Element_V=np.ascontiguousarray(Element_V,dtype=np.float64)
    Element_A=np.ascontiguousarray(Element_A,dtype=np.float64)
    S_vec=np.ascontiguousarray(S_vec,dtype=np.float64)
    PENER=np.ascontiguousarray(PENER,dtype=np.float64)
    SENER=np.ascontiguousarray(SENER,dtype=np.float64)
    ################################################

    num_elem_=ctypes.c_size_t(Coords.shape[0])
    Coords_p=Coords.ctypes.data_as(ctypes.c_void_p)
    rho_p   =rho.ctypes.data_as(ctypes.c_void_p)
    Element_U_p=Element_U.ctypes.data_as(ctypes.c_void_p)
    Element_V_p=Element_V.ctypes.data_as(ctypes.c_void_p)
    Element_A_p=Element_A.ctypes.data_as(ctypes.c_void_p)
    S_vec_p=S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p=PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p=SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force=np.zeros((Coords.shape[0],Coords.shape[1],3))
    Conf_Force_p=Conf_Force.ctypes.data_as(ctypes.c_void_p)

    lib.Integration_CPE4R_dynamic(num_elem_,Coords_p,rho_p,Element_U_p,Element_V_p,Element_A_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    return Conf_Force

def calc_Conf_Force_CPE8_static(Coords,Element_U,S_vec,PENER,SENER,method='mbf'):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates [num_elem,numNodes,3]
        Element_U......Element_Nodal_Displacements [num_elem,numNodes,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy [num_elem,num_int_points]
        SENER..........Elastic strain energy [num_elem,num_int_points]
    """
    #This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem=Coords.shape[0]
    numNodes=8
    num_int_points=9

    assert Coords.shape==(num_elem,numNodes,3)
    assert Element_U.shape==(num_elem,numNodes,3)
    assert S_vec.shape==(num_elem,num_int_points,6)
    assert PENER.shape==(num_elem,num_int_points)
    assert SENER.shape==(num_elem,num_int_points)

    #Ensure that all arrays are c contiguous
    Coords=np.ascontiguousarray(Coords,dtype=np.float64)
    Element_U=np.ascontiguousarray(Element_U,dtype=np.float64)
    S_vec=np.ascontiguousarray(S_vec,dtype=np.float64)
    PENER=np.ascontiguousarray(PENER,dtype=np.float64)
    SENER=np.ascontiguousarray(SENER,dtype=np.float64)
    ################################################

    num_elem_=ctypes.c_size_t(Coords.shape[0])
    Coords_p=Coords.ctypes.data_as(ctypes.c_void_p)
    Element_U_p=Element_U.ctypes.data_as(ctypes.c_void_p)
    S_vec_p=S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p=PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p=SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force=np.zeros((Coords.shape[0],Coords.shape[1],3))
    Conf_Force_p=Conf_Force.ctypes.data_as(ctypes.c_void_p)

    if method=='mbf':
        lib.Integration_CPE8_static_mbf(num_elem_,Coords_p,Element_U_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    elif method=='dbf':
        lib.Integration_CPE8_static_dbf(num_elem_,Coords_p,Element_U_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    else:
        print("Method not found!")
    return Conf_Force

def calc_Conf_Force_CPE8_dynamic(Coords,rho,Element_U,Element_V,Element_A,S_vec,PENER,SENER):
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
    #This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem=Coords.shape[0]
    numNodes=8
    num_int_points=9

    assert Coords.shape==(num_elem,numNodes,3)
    assert rho.shape==(num_elem,)
    assert Element_U.shape==(num_elem,numNodes,3)
    assert Element_V.shape==(num_elem,numNodes,3)
    assert Element_A.shape==(num_elem,numNodes,3)
    assert S_vec.shape==(num_elem,num_int_points,6)
    assert PENER.shape==(num_elem,num_int_points)
    assert SENER.shape==(num_elem,num_int_points)

    #Ensure that all arrays are c contiguous
    Coords=np.ascontiguousarray(Coords,dtype=np.float64)
    rho=np.ascontiguousarray(rho,dtype=np.float64)
    Element_U=np.ascontiguousarray(Element_U,dtype=np.float64)
    Element_V=np.ascontiguousarray(Element_V,dtype=np.float64)
    Element_A=np.ascontiguousarray(Element_A,dtype=np.float64)
    S_vec=np.ascontiguousarray(S_vec,dtype=np.float64)
    PENER=np.ascontiguousarray(PENER,dtype=np.float64)
    SENER=np.ascontiguousarray(SENER,dtype=np.float64)
    ################################################

    num_elem_=ctypes.c_size_t(Coords.shape[0])
    Coords_p=Coords.ctypes.data_as(ctypes.c_void_p)
    rho_p   =rho.ctypes.data_as(ctypes.c_void_p)
    Element_U_p=Element_U.ctypes.data_as(ctypes.c_void_p)
    Element_V_p=Element_V.ctypes.data_as(ctypes.c_void_p)
    Element_A_p=Element_A.ctypes.data_as(ctypes.c_void_p)
    S_vec_p=S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p=PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p=SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force=np.zeros((Coords.shape[0],Coords.shape[1],3))
    Conf_Force_p=Conf_Force.ctypes.data_as(ctypes.c_void_p)

    lib.Integration_CPE8_dynamic(num_elem_,Coords_p,rho_p,Element_U_p,Element_V_p,Element_A_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    return Conf_Force

def calc_Conf_Force_CPE8R_static(Coords,Element_U,S_vec,PENER,SENER,method='mbf'):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates [num_elem,numNodes,3]
        Element_U......Element_Nodal_Displacements [num_elem,numNodes,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy [num_elem,num_int_points]
        SENER..........Elastic strain energy [num_elem,num_int_points]
    """
    #This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem=Coords.shape[0]
    numNodes=8
    num_int_points=4

    assert Coords.shape==(num_elem,numNodes,3)
    assert Element_U.shape==(num_elem,numNodes,3)
    assert S_vec.shape==(num_elem,num_int_points,6)
    assert PENER.shape==(num_elem,num_int_points)
    assert SENER.shape==(num_elem,num_int_points)

    #Ensure that all arrays are c contiguous
    Coords=np.ascontiguousarray(Coords,dtype=np.float64)
    Element_U=np.ascontiguousarray(Element_U,dtype=np.float64)
    S_vec=np.ascontiguousarray(S_vec,dtype=np.float64)
    PENER=np.ascontiguousarray(PENER,dtype=np.float64)
    SENER=np.ascontiguousarray(SENER,dtype=np.float64)
    ################################################

    num_elem_=ctypes.c_size_t(Coords.shape[0])
    Coords_p=Coords.ctypes.data_as(ctypes.c_void_p)
    Element_U_p=Element_U.ctypes.data_as(ctypes.c_void_p)
    S_vec_p=S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p=PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p=SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force=np.zeros((Coords.shape[0],Coords.shape[1],3))
    Conf_Force_p=Conf_Force.ctypes.data_as(ctypes.c_void_p)

    if method=='mbf':
        lib.Integration_CPE8R_static_mbf(num_elem_,Coords_p,Element_U_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    elif method=='dbf':
        lib.Integration_CPE8R_static_dbf(num_elem_,Coords_p,Element_U_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    else:
        print("Method not found!")
    return Conf_Force

def calc_Conf_Force_CPE8R_dynamic(Coords,rho,Element_U,Element_V,Element_A,S_vec,PENER,SENER):
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
    #This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem=Coords.shape[0]
    numNodes=8
    num_int_points=4

    assert Coords.shape==(num_elem,numNodes,3)
    assert rho.shape==(num_elem,)
    assert Element_U.shape==(num_elem,numNodes,3)
    assert Element_V.shape==(num_elem,numNodes,3)
    assert Element_A.shape==(num_elem,numNodes,3)
    assert S_vec.shape==(num_elem,num_int_points,6)
    assert PENER.shape==(num_elem,num_int_points)
    assert SENER.shape==(num_elem,num_int_points)

    #Ensure that all arrays are c contiguous
    Coords=np.ascontiguousarray(Coords,dtype=np.float64)
    rho=np.ascontiguousarray(rho,dtype=np.float64)
    Element_U=np.ascontiguousarray(Element_U,dtype=np.float64)
    Element_V=np.ascontiguousarray(Element_V,dtype=np.float64)
    Element_A=np.ascontiguousarray(Element_A,dtype=np.float64)
    S_vec=np.ascontiguousarray(S_vec,dtype=np.float64)
    PENER=np.ascontiguousarray(PENER,dtype=np.float64)
    SENER=np.ascontiguousarray(SENER,dtype=np.float64)
    ################################################

    num_elem_=ctypes.c_size_t(Coords.shape[0])
    Coords_p=Coords.ctypes.data_as(ctypes.c_void_p)
    rho_p   =rho.ctypes.data_as(ctypes.c_void_p)
    Element_U_p=Element_U.ctypes.data_as(ctypes.c_void_p)
    Element_V_p=Element_V.ctypes.data_as(ctypes.c_void_p)
    Element_A_p=Element_A.ctypes.data_as(ctypes.c_void_p)
    S_vec_p=S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p=PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p=SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force=np.zeros((Coords.shape[0],Coords.shape[1],3))
    Conf_Force_p=Conf_Force.ctypes.data_as(ctypes.c_void_p)

    lib.Integration_CPE8R_dynamic(num_elem_,Coords_p,rho_p,Element_U_p,Element_V_p,Element_A_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    return Conf_Force

def calc_Conf_Force_CPE3_static(Coords,Element_U,S_vec,PENER,SENER,method='mbf'):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates [num_elem,numNodes,3]
        Element_U......Element_Nodal_Displacements [num_elem,numNodes,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy [num_elem,num_int_points]
        SENER..........Elastic strain energy [num_elem,num_int_points]
    """
    #This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem=Coords.shape[0]
    numNodes=3
    num_int_points=1

    assert Coords.shape==(num_elem,numNodes,3)
    assert Element_U.shape==(num_elem,numNodes,3)
    assert S_vec.shape==(num_elem,num_int_points,6)
    assert PENER.shape==(num_elem,num_int_points)
    assert SENER.shape==(num_elem,num_int_points)

    #Ensure that all arrays are c contiguous
    Coords=np.ascontiguousarray(Coords,dtype=np.float64)
    Element_U=np.ascontiguousarray(Element_U,dtype=np.float64)
    S_vec=np.ascontiguousarray(S_vec,dtype=np.float64)
    PENER=np.ascontiguousarray(PENER,dtype=np.float64)
    SENER=np.ascontiguousarray(SENER,dtype=np.float64)
    ################################################

    num_elem_=ctypes.c_size_t(Coords.shape[0])
    Coords_p=Coords.ctypes.data_as(ctypes.c_void_p)
    Element_U_p=Element_U.ctypes.data_as(ctypes.c_void_p)
    S_vec_p=S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p=PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p=SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force=np.zeros((Coords.shape[0],Coords.shape[1],3))
    Conf_Force_p=Conf_Force.ctypes.data_as(ctypes.c_void_p)

    if method=='mbf':
        lib.Integration_CPE3_static_mbf(num_elem_,Coords_p,Element_U_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    elif method=='dbf':
        lib.Integration_CPE3_static_dbf(num_elem_,Coords_p,Element_U_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    else:
        print("Method not found!")
    return Conf_Force

def calc_Conf_Force_CPE3_dynamic(Coords,rho,Element_U,Element_V,Element_A,S_vec,PENER,SENER):
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
    #This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem=Coords.shape[0]
    numNodes=3
    num_int_points=1

    assert Coords.shape==(num_elem,numNodes,3)
    assert rho.shape==(num_elem,)
    assert Element_U.shape==(num_elem,numNodes,3)
    assert Element_V.shape==(num_elem,numNodes,3)
    assert Element_A.shape==(num_elem,numNodes,3)
    assert S_vec.shape==(num_elem,num_int_points,6)
    assert PENER.shape==(num_elem,num_int_points)
    assert SENER.shape==(num_elem,num_int_points)

    #Ensure that all arrays are c contiguous
    Coords=np.ascontiguousarray(Coords,dtype=np.float64)
    rho=np.ascontiguousarray(rho,dtype=np.float64)
    Element_U=np.ascontiguousarray(Element_U,dtype=np.float64)
    Element_V=np.ascontiguousarray(Element_V,dtype=np.float64)
    Element_A=np.ascontiguousarray(Element_A,dtype=np.float64)
    S_vec=np.ascontiguousarray(S_vec,dtype=np.float64)
    PENER=np.ascontiguousarray(PENER,dtype=np.float64)
    SENER=np.ascontiguousarray(SENER,dtype=np.float64)
    ################################################

    num_elem_=ctypes.c_size_t(Coords.shape[0])
    Coords_p=Coords.ctypes.data_as(ctypes.c_void_p)
    rho_p   =rho.ctypes.data_as(ctypes.c_void_p)
    Element_U_p=Element_U.ctypes.data_as(ctypes.c_void_p)
    Element_V_p=Element_V.ctypes.data_as(ctypes.c_void_p)
    Element_A_p=Element_A.ctypes.data_as(ctypes.c_void_p)
    S_vec_p=S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p=PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p=SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force=np.zeros((Coords.shape[0],Coords.shape[1],3))
    Conf_Force_p=Conf_Force.ctypes.data_as(ctypes.c_void_p)

    lib.Integration_CPE3_dynamic(num_elem_,Coords_p,rho_p,Element_U_p,Element_V_p,Element_A_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    return Conf_Force

def calc_Conf_Force_CPE6_static(Coords,Element_U,S_vec,PENER,SENER,method='mbf'):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates [num_elem,numNodes,3]
        Element_U......Element_Nodal_Displacements [num_elem,numNodes,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy [num_elem,num_int_points]
        SENER..........Elastic strain energy [num_elem,num_int_points]
    """
    #This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem=Coords.shape[0]
    numNodes=6
    num_int_points=3

    assert Coords.shape==(num_elem,numNodes,3)
    assert Element_U.shape==(num_elem,numNodes,3)
    assert S_vec.shape==(num_elem,num_int_points,6)
    assert PENER.shape==(num_elem,num_int_points)
    assert SENER.shape==(num_elem,num_int_points)

    #Ensure that all arrays are c contiguous
    Coords=np.ascontiguousarray(Coords,dtype=np.float64)
    Element_U=np.ascontiguousarray(Element_U,dtype=np.float64)
    S_vec=np.ascontiguousarray(S_vec,dtype=np.float64)
    PENER=np.ascontiguousarray(PENER,dtype=np.float64)
    SENER=np.ascontiguousarray(SENER,dtype=np.float64)
    ################################################

    num_elem_=ctypes.c_size_t(Coords.shape[0])
    Coords_p=Coords.ctypes.data_as(ctypes.c_void_p)
    Element_U_p=Element_U.ctypes.data_as(ctypes.c_void_p)
    S_vec_p=S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p=PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p=SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force=np.zeros((Coords.shape[0],Coords.shape[1],3))
    Conf_Force_p=Conf_Force.ctypes.data_as(ctypes.c_void_p)

    if method=='mbf':
        lib.Integration_CPE6_static_mbf(num_elem_,Coords_p,Element_U_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    elif method=='dbf':
        lib.Integration_CPE6_static_dbf(num_elem_,Coords_p,Element_U_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    else:
        print("Method not found!")
    return Conf_Force

def calc_Conf_Force_CPE6_dynamic(Coords,rho,Element_U,Element_V,Element_A,S_vec,PENER,SENER):
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
    #This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem=Coords.shape[0]
    numNodes=6
    num_int_points=3

    assert Coords.shape==(num_elem,numNodes,3)
    assert rho.shape==(num_elem,)
    assert Element_U.shape==(num_elem,numNodes,3)
    assert Element_V.shape==(num_elem,numNodes,3)
    assert Element_A.shape==(num_elem,numNodes,3)
    assert S_vec.shape==(num_elem,num_int_points,6)
    assert PENER.shape==(num_elem,num_int_points)
    assert SENER.shape==(num_elem,num_int_points)

    #Ensure that all arrays are c contiguous
    Coords=np.ascontiguousarray(Coords,dtype=np.float64)
    rho=np.ascontiguousarray(rho,dtype=np.float64)
    Element_U=np.ascontiguousarray(Element_U,dtype=np.float64)
    Element_V=np.ascontiguousarray(Element_V,dtype=np.float64)
    Element_A=np.ascontiguousarray(Element_A,dtype=np.float64)
    S_vec=np.ascontiguousarray(S_vec,dtype=np.float64)
    PENER=np.ascontiguousarray(PENER,dtype=np.float64)
    SENER=np.ascontiguousarray(SENER,dtype=np.float64)
    ################################################

    num_elem_=ctypes.c_size_t(Coords.shape[0])
    Coords_p=Coords.ctypes.data_as(ctypes.c_void_p)
    rho_p   =rho.ctypes.data_as(ctypes.c_void_p)
    Element_U_p=Element_U.ctypes.data_as(ctypes.c_void_p)
    Element_V_p=Element_V.ctypes.data_as(ctypes.c_void_p)
    Element_A_p=Element_A.ctypes.data_as(ctypes.c_void_p)
    S_vec_p=S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p=PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p=SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force=np.zeros((Coords.shape[0],Coords.shape[1],3))
    Conf_Force_p=Conf_Force.ctypes.data_as(ctypes.c_void_p)

    lib.Integration_CPE6_dynamic(num_elem_,Coords_p,rho_p,Element_U_p,Element_V_p,Element_A_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    return Conf_Force

def calc_Conf_Force_C3D8_static(Coords,Element_U,S_vec,PENER,SENER,method='mbf'):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates [num_elem,numNodes,3]
        Element_U......Element_Nodal_Displacements [num_elem,numNodes,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy [num_elem,num_int_points]
        SENER..........Elastic strain energy [num_elem,num_int_points]
    """
    #This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem=Coords.shape[0]
    numNodes=8
    num_int_points=8

    assert Coords.shape==(num_elem,numNodes,3)
    assert Element_U.shape==(num_elem,numNodes,3)
    assert S_vec.shape==(num_elem,num_int_points,6)
    assert PENER.shape==(num_elem,num_int_points)
    assert SENER.shape==(num_elem,num_int_points)

    #Ensure that all arrays are c contiguous
    Coords=np.ascontiguousarray(Coords,dtype=np.float64)
    Element_U=np.ascontiguousarray(Element_U,dtype=np.float64)
    S_vec=np.ascontiguousarray(S_vec,dtype=np.float64)
    PENER=np.ascontiguousarray(PENER,dtype=np.float64)
    SENER=np.ascontiguousarray(SENER,dtype=np.float64)
    ################################################

    num_elem_=ctypes.c_size_t(Coords.shape[0])
    Coords_p=Coords.ctypes.data_as(ctypes.c_void_p)
    Element_U_p=Element_U.ctypes.data_as(ctypes.c_void_p)
    S_vec_p=S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p=PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p=SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force=np.zeros((Coords.shape[0],Coords.shape[1],3))
    Conf_Force_p=Conf_Force.ctypes.data_as(ctypes.c_void_p)

    if method=='mbf':
        lib.Integration_C3D8_static_mbf(num_elem_,Coords_p,Element_U_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    elif method=='dbf':
        lib.Integration_C3D8_static_dbf(num_elem_,Coords_p,Element_U_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    else:
        print("Method not found!")
    return Conf_Force

def calc_Conf_Force_C3D8_dynamic(Coords,rho,Element_U,Element_V,Element_A,S_vec,PENER,SENER):
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
    #This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem=Coords.shape[0]
    numNodes=8
    num_int_points=8

    assert Coords.shape==(num_elem,numNodes,3)
    assert rho.shape==(num_elem,)
    assert Element_U.shape==(num_elem,numNodes,3)
    assert Element_V.shape==(num_elem,numNodes,3)
    assert Element_A.shape==(num_elem,numNodes,3)
    assert S_vec.shape==(num_elem,num_int_points,6)
    assert PENER.shape==(num_elem,num_int_points)
    assert SENER.shape==(num_elem,num_int_points)

    #Ensure that all arrays are c contiguous
    Coords=np.ascontiguousarray(Coords,dtype=np.float64)
    rho=np.ascontiguousarray(rho,dtype=np.float64)
    Element_U=np.ascontiguousarray(Element_U,dtype=np.float64)
    Element_V=np.ascontiguousarray(Element_V,dtype=np.float64)
    Element_A=np.ascontiguousarray(Element_A,dtype=np.float64)
    S_vec=np.ascontiguousarray(S_vec,dtype=np.float64)
    PENER=np.ascontiguousarray(PENER,dtype=np.float64)
    SENER=np.ascontiguousarray(SENER,dtype=np.float64)
    ################################################

    num_elem_=ctypes.c_size_t(Coords.shape[0])
    Coords_p=Coords.ctypes.data_as(ctypes.c_void_p)
    rho_p   =rho.ctypes.data_as(ctypes.c_void_p)
    Element_U_p=Element_U.ctypes.data_as(ctypes.c_void_p)
    Element_V_p=Element_V.ctypes.data_as(ctypes.c_void_p)
    Element_A_p=Element_A.ctypes.data_as(ctypes.c_void_p)
    S_vec_p=S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p=PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p=SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force=np.zeros((Coords.shape[0],Coords.shape[1],3))
    Conf_Force_p=Conf_Force.ctypes.data_as(ctypes.c_void_p)

    lib.Integration_C3D8_dynamic(num_elem_,Coords_p,rho_p,Element_U_p,Element_V_p,Element_A_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    return Conf_Force

def calc_Conf_Force_C3D8R_static(Coords,Element_U,S_vec,PENER,SENER,method='mbf'):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates [num_elem,numNodes,3]
        Element_U......Element_Nodal_Displacements [num_elem,numNodes,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy [num_elem,num_int_points]
        SENER..........Elastic strain energy [num_elem,num_int_points]
    """
    #This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem=Coords.shape[0]
    numNodes=8
    num_int_points=1

    assert Coords.shape==(num_elem,numNodes,3)
    assert Element_U.shape==(num_elem,numNodes,3)
    assert S_vec.shape==(num_elem,num_int_points,6)
    assert PENER.shape==(num_elem,num_int_points)
    assert SENER.shape==(num_elem,num_int_points)

    #Ensure that all arrays are c contiguous
    Coords=np.ascontiguousarray(Coords,dtype=np.float64)
    Element_U=np.ascontiguousarray(Element_U,dtype=np.float64)
    S_vec=np.ascontiguousarray(S_vec,dtype=np.float64)
    PENER=np.ascontiguousarray(PENER,dtype=np.float64)
    SENER=np.ascontiguousarray(SENER,dtype=np.float64)
    ################################################

    num_elem_=ctypes.c_size_t(Coords.shape[0])
    Coords_p=Coords.ctypes.data_as(ctypes.c_void_p)
    Element_U_p=Element_U.ctypes.data_as(ctypes.c_void_p)
    S_vec_p=S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p=PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p=SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force=np.zeros((Coords.shape[0],Coords.shape[1],3))
    Conf_Force_p=Conf_Force.ctypes.data_as(ctypes.c_void_p)

    if method=='mbf':
        lib.Integration_C3D8R_static_mbf(num_elem_,Coords_p,Element_U_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    elif method=='dbf':
        lib.Integration_C3D8R_static_dbf(num_elem_,Coords_p,Element_U_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    else:
        print("Method not found!")
    return Conf_Force

def calc_Conf_Force_C3D8R_dynamic(Coords,rho,Element_U,Element_V,Element_A,S_vec,PENER,SENER):
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
    #This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem=Coords.shape[0]
    numNodes=8
    num_int_points=1

    assert Coords.shape==(num_elem,numNodes,3)
    assert rho.shape==(num_elem,)
    assert Element_U.shape==(num_elem,numNodes,3)
    assert Element_V.shape==(num_elem,numNodes,3)
    assert Element_A.shape==(num_elem,numNodes,3)
    assert S_vec.shape==(num_elem,num_int_points,6)
    assert PENER.shape==(num_elem,num_int_points)
    assert SENER.shape==(num_elem,num_int_points)

    #Ensure that all arrays are c contiguous
    Coords=np.ascontiguousarray(Coords,dtype=np.float64)
    rho=np.ascontiguousarray(rho,dtype=np.float64)
    Element_U=np.ascontiguousarray(Element_U,dtype=np.float64)
    Element_V=np.ascontiguousarray(Element_V,dtype=np.float64)
    Element_A=np.ascontiguousarray(Element_A,dtype=np.float64)
    S_vec=np.ascontiguousarray(S_vec,dtype=np.float64)
    PENER=np.ascontiguousarray(PENER,dtype=np.float64)
    SENER=np.ascontiguousarray(SENER,dtype=np.float64)
    ################################################

    num_elem_=ctypes.c_size_t(Coords.shape[0])
    Coords_p=Coords.ctypes.data_as(ctypes.c_void_p)
    rho_p   =rho.ctypes.data_as(ctypes.c_void_p)
    Element_U_p=Element_U.ctypes.data_as(ctypes.c_void_p)
    Element_V_p=Element_V.ctypes.data_as(ctypes.c_void_p)
    Element_A_p=Element_A.ctypes.data_as(ctypes.c_void_p)
    S_vec_p=S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p=PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p=SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force=np.zeros((Coords.shape[0],Coords.shape[1],3))
    Conf_Force_p=Conf_Force.ctypes.data_as(ctypes.c_void_p)

    lib.Integration_C3D8R_dynamic(num_elem_,Coords_p,rho_p,Element_U_p,Element_V_p,Element_A_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    return Conf_Force

def calc_Conf_Force_C3D20_static(Coords,Element_U,S_vec,PENER,SENER,method='mbf'):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates [num_elem,numNodes,3]
        Element_U......Element_Nodal_Displacements [num_elem,numNodes,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy [num_elem,num_int_points]
        SENER..........Elastic strain energy [num_elem,num_int_points]
    """
    #This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem=Coords.shape[0]
    numNodes=20
    num_int_points=27

    assert Coords.shape==(num_elem,numNodes,3)
    assert Element_U.shape==(num_elem,numNodes,3)
    assert S_vec.shape==(num_elem,num_int_points,6)
    assert PENER.shape==(num_elem,num_int_points)
    assert SENER.shape==(num_elem,num_int_points)

    #Ensure that all arrays are c contiguous
    Coords=np.ascontiguousarray(Coords,dtype=np.float64)
    Element_U=np.ascontiguousarray(Element_U,dtype=np.float64)
    S_vec=np.ascontiguousarray(S_vec,dtype=np.float64)
    PENER=np.ascontiguousarray(PENER,dtype=np.float64)
    SENER=np.ascontiguousarray(SENER,dtype=np.float64)
    ################################################

    num_elem_=ctypes.c_size_t(Coords.shape[0])
    Coords_p=Coords.ctypes.data_as(ctypes.c_void_p)
    Element_U_p=Element_U.ctypes.data_as(ctypes.c_void_p)
    S_vec_p=S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p=PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p=SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force=np.zeros((Coords.shape[0],Coords.shape[1],3))
    Conf_Force_p=Conf_Force.ctypes.data_as(ctypes.c_void_p)

    if method=='mbf':
        lib.Integration_C3D20_static_mbf(num_elem_,Coords_p,Element_U_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    elif method=='dbf':
        lib.Integration_C3D20_static_dbf(num_elem_,Coords_p,Element_U_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    else:
        print("Method not found!")
    return Conf_Force

def calc_Conf_Force_C3D20_dynamic(Coords,rho,Element_U,Element_V,Element_A,S_vec,PENER,SENER):
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
    #This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem=Coords.shape[0]
    numNodes=20
    num_int_points=27

    assert Coords.shape==(num_elem,numNodes,3)
    assert rho.shape==(num_elem,)
    assert Element_U.shape==(num_elem,numNodes,3)
    assert Element_V.shape==(num_elem,numNodes,3)
    assert Element_A.shape==(num_elem,numNodes,3)
    assert S_vec.shape==(num_elem,num_int_points,6)
    assert PENER.shape==(num_elem,num_int_points)
    assert SENER.shape==(num_elem,num_int_points)

    #Ensure that all arrays are c contiguous
    Coords=np.ascontiguousarray(Coords,dtype=np.float64)
    rho=np.ascontiguousarray(rho,dtype=np.float64)
    Element_U=np.ascontiguousarray(Element_U,dtype=np.float64)
    Element_V=np.ascontiguousarray(Element_V,dtype=np.float64)
    Element_A=np.ascontiguousarray(Element_A,dtype=np.float64)
    S_vec=np.ascontiguousarray(S_vec,dtype=np.float64)
    PENER=np.ascontiguousarray(PENER,dtype=np.float64)
    SENER=np.ascontiguousarray(SENER,dtype=np.float64)
    ################################################

    num_elem_=ctypes.c_size_t(Coords.shape[0])
    Coords_p=Coords.ctypes.data_as(ctypes.c_void_p)
    rho_p   =rho.ctypes.data_as(ctypes.c_void_p)
    Element_U_p=Element_U.ctypes.data_as(ctypes.c_void_p)
    Element_V_p=Element_V.ctypes.data_as(ctypes.c_void_p)
    Element_A_p=Element_A.ctypes.data_as(ctypes.c_void_p)
    S_vec_p=S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p=PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p=SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force=np.zeros((Coords.shape[0],Coords.shape[1],3))
    Conf_Force_p=Conf_Force.ctypes.data_as(ctypes.c_void_p)

    lib.Integration_C3D20_dynamic(num_elem_,Coords_p,rho_p,Element_U_p,Element_V_p,Element_A_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    return Conf_Force

def calc_Conf_Force_C3D20R_static(Coords,Element_U,S_vec,PENER,SENER,method='mbf'):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates [num_elem,numNodes,3]
        Element_U......Element_Nodal_Displacements [num_elem,numNodes,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy [num_elem,num_int_points]
        SENER..........Elastic strain energy [num_elem,num_int_points]
    """
    #This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem=Coords.shape[0]
    numNodes=20
    num_int_points=8

    assert Coords.shape==(num_elem,numNodes,3)
    assert Element_U.shape==(num_elem,numNodes,3)
    assert S_vec.shape==(num_elem,num_int_points,6)
    assert PENER.shape==(num_elem,num_int_points)
    assert SENER.shape==(num_elem,num_int_points)

    #Ensure that all arrays are c contiguous
    Coords=np.ascontiguousarray(Coords,dtype=np.float64)
    Element_U=np.ascontiguousarray(Element_U,dtype=np.float64)
    S_vec=np.ascontiguousarray(S_vec,dtype=np.float64)
    PENER=np.ascontiguousarray(PENER,dtype=np.float64)
    SENER=np.ascontiguousarray(SENER,dtype=np.float64)
    ################################################

    num_elem_=ctypes.c_size_t(Coords.shape[0])
    Coords_p=Coords.ctypes.data_as(ctypes.c_void_p)
    Element_U_p=Element_U.ctypes.data_as(ctypes.c_void_p)
    S_vec_p=S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p=PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p=SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force=np.zeros((Coords.shape[0],Coords.shape[1],3))
    Conf_Force_p=Conf_Force.ctypes.data_as(ctypes.c_void_p)

    if method=='mbf':
        lib.Integration_C3D20R_static_mbf(num_elem_,Coords_p,Element_U_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    elif method=='dbf':
        lib.Integration_C3D20R_static_dbf(num_elem_,Coords_p,Element_U_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    else:
        print("Method not found!")
    return Conf_Force

def calc_Conf_Force_C3D20R_dynamic(Coords,rho,Element_U,Element_V,Element_A,S_vec,PENER,SENER):
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
    #This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem=Coords.shape[0]
    numNodes=20
    num_int_points=8

    assert Coords.shape==(num_elem,numNodes,3)
    assert rho.shape==(num_elem,)
    assert Element_U.shape==(num_elem,numNodes,3)
    assert Element_V.shape==(num_elem,numNodes,3)
    assert Element_A.shape==(num_elem,numNodes,3)
    assert S_vec.shape==(num_elem,num_int_points,6)
    assert PENER.shape==(num_elem,num_int_points)
    assert SENER.shape==(num_elem,num_int_points)

    #Ensure that all arrays are c contiguous
    Coords=np.ascontiguousarray(Coords,dtype=np.float64)
    rho=np.ascontiguousarray(rho,dtype=np.float64)
    Element_U=np.ascontiguousarray(Element_U,dtype=np.float64)
    Element_V=np.ascontiguousarray(Element_V,dtype=np.float64)
    Element_A=np.ascontiguousarray(Element_A,dtype=np.float64)
    S_vec=np.ascontiguousarray(S_vec,dtype=np.float64)
    PENER=np.ascontiguousarray(PENER,dtype=np.float64)
    SENER=np.ascontiguousarray(SENER,dtype=np.float64)
    ################################################

    num_elem_=ctypes.c_size_t(Coords.shape[0])
    Coords_p=Coords.ctypes.data_as(ctypes.c_void_p)
    rho_p   =rho.ctypes.data_as(ctypes.c_void_p)
    Element_U_p=Element_U.ctypes.data_as(ctypes.c_void_p)
    Element_V_p=Element_V.ctypes.data_as(ctypes.c_void_p)
    Element_A_p=Element_A.ctypes.data_as(ctypes.c_void_p)
    S_vec_p=S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p=PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p=SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force=np.zeros((Coords.shape[0],Coords.shape[1],3))
    Conf_Force_p=Conf_Force.ctypes.data_as(ctypes.c_void_p)

    lib.Integration_C3D20R_dynamic(num_elem_,Coords_p,rho_p,Element_U_p,Element_V_p,Element_A_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    return Conf_Force

def calc_Conf_Force_C3D4_static(Coords,Element_U,S_vec,PENER,SENER,method='mbf'):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates [num_elem,numNodes,3]
        Element_U......Element_Nodal_Displacements [num_elem,numNodes,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy [num_elem,num_int_points]
        SENER..........Elastic strain energy [num_elem,num_int_points]
    """
    #This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem=Coords.shape[0]
    numNodes=4
    num_int_points=1

    assert Coords.shape==(num_elem,numNodes,3)
    assert Element_U.shape==(num_elem,numNodes,3)
    assert S_vec.shape==(num_elem,num_int_points,6)
    assert PENER.shape==(num_elem,num_int_points)
    assert SENER.shape==(num_elem,num_int_points)

    #Ensure that all arrays are c contiguous
    Coords=np.ascontiguousarray(Coords,dtype=np.float64)
    Element_U=np.ascontiguousarray(Element_U,dtype=np.float64)
    S_vec=np.ascontiguousarray(S_vec,dtype=np.float64)
    PENER=np.ascontiguousarray(PENER,dtype=np.float64)
    SENER=np.ascontiguousarray(SENER,dtype=np.float64)
    ################################################

    num_elem_=ctypes.c_size_t(Coords.shape[0])
    Coords_p=Coords.ctypes.data_as(ctypes.c_void_p)
    Element_U_p=Element_U.ctypes.data_as(ctypes.c_void_p)
    S_vec_p=S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p=PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p=SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force=np.zeros((Coords.shape[0],Coords.shape[1],3))
    Conf_Force_p=Conf_Force.ctypes.data_as(ctypes.c_void_p)

    if method=='mbf':
        lib.Integration_C3D4_static_mbf(num_elem_,Coords_p,Element_U_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    elif method=='dbf':
        lib.Integration_C3D4_static_dbf(num_elem_,Coords_p,Element_U_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    else:
        print("Method not found!")
    return Conf_Force

def calc_Conf_Force_C3D4_dynamic(Coords,rho,Element_U,Element_V,Element_A,S_vec,PENER,SENER):
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
    #This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem=Coords.shape[0]
    numNodes=4
    num_int_points=1

    assert Coords.shape==(num_elem,numNodes,3)
    assert rho.shape==(num_elem,)
    assert Element_U.shape==(num_elem,numNodes,3)
    assert Element_V.shape==(num_elem,numNodes,3)
    assert Element_A.shape==(num_elem,numNodes,3)
    assert S_vec.shape==(num_elem,num_int_points,6)
    assert PENER.shape==(num_elem,num_int_points)
    assert SENER.shape==(num_elem,num_int_points)

    #Ensure that all arrays are c contiguous
    Coords=np.ascontiguousarray(Coords,dtype=np.float64)
    rho=np.ascontiguousarray(rho,dtype=np.float64)
    Element_U=np.ascontiguousarray(Element_U,dtype=np.float64)
    Element_V=np.ascontiguousarray(Element_V,dtype=np.float64)
    Element_A=np.ascontiguousarray(Element_A,dtype=np.float64)
    S_vec=np.ascontiguousarray(S_vec,dtype=np.float64)
    PENER=np.ascontiguousarray(PENER,dtype=np.float64)
    SENER=np.ascontiguousarray(SENER,dtype=np.float64)
    ################################################

    num_elem_=ctypes.c_size_t(Coords.shape[0])
    Coords_p=Coords.ctypes.data_as(ctypes.c_void_p)
    rho_p   =rho.ctypes.data_as(ctypes.c_void_p)
    Element_U_p=Element_U.ctypes.data_as(ctypes.c_void_p)
    Element_V_p=Element_V.ctypes.data_as(ctypes.c_void_p)
    Element_A_p=Element_A.ctypes.data_as(ctypes.c_void_p)
    S_vec_p=S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p=PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p=SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force=np.zeros((Coords.shape[0],Coords.shape[1],3))
    Conf_Force_p=Conf_Force.ctypes.data_as(ctypes.c_void_p)

    lib.Integration_C3D4_dynamic(num_elem_,Coords_p,rho_p,Element_U_p,Element_V_p,Element_A_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    return Conf_Force

def calc_Conf_Force_C3D10_static(Coords,Element_U,S_vec,PENER,SENER,method='mbf'):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates [num_elem,numNodes,3]
        Element_U......Element_Nodal_Displacements [num_elem,numNodes,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy [num_elem,num_int_points]
        SENER..........Elastic strain energy [num_elem,num_int_points]
    """
    #This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem=Coords.shape[0]
    numNodes=10
    num_int_points=4

    assert Coords.shape==(num_elem,numNodes,3)
    assert Element_U.shape==(num_elem,numNodes,3)
    assert S_vec.shape==(num_elem,num_int_points,6)
    assert PENER.shape==(num_elem,num_int_points)
    assert SENER.shape==(num_elem,num_int_points)

    #Ensure that all arrays are c contiguous
    Coords=np.ascontiguousarray(Coords,dtype=np.float64)
    Element_U=np.ascontiguousarray(Element_U,dtype=np.float64)
    S_vec=np.ascontiguousarray(S_vec,dtype=np.float64)
    PENER=np.ascontiguousarray(PENER,dtype=np.float64)
    SENER=np.ascontiguousarray(SENER,dtype=np.float64)
    ################################################

    num_elem_=ctypes.c_size_t(Coords.shape[0])
    Coords_p=Coords.ctypes.data_as(ctypes.c_void_p)
    Element_U_p=Element_U.ctypes.data_as(ctypes.c_void_p)
    S_vec_p=S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p=PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p=SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force=np.zeros((Coords.shape[0],Coords.shape[1],3))
    Conf_Force_p=Conf_Force.ctypes.data_as(ctypes.c_void_p)

    if method=='mbf':
        lib.Integration_C3D10_static_mbf(num_elem_,Coords_p,Element_U_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    elif method=='dbf':
        lib.Integration_C3D10_static_dbf(num_elem_,Coords_p,Element_U_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    else:
        print("Method not found!")
    return Conf_Force

def calc_Conf_Force_C3D10_dynamic(Coords,rho,Element_U,Element_V,Element_A,S_vec,PENER,SENER):
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
    #This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem=Coords.shape[0]
    numNodes=10
    num_int_points=4

    assert Coords.shape==(num_elem,numNodes,3)
    assert rho.shape==(num_elem,)
    assert Element_U.shape==(num_elem,numNodes,3)
    assert Element_V.shape==(num_elem,numNodes,3)
    assert Element_A.shape==(num_elem,numNodes,3)
    assert S_vec.shape==(num_elem,num_int_points,6)
    assert PENER.shape==(num_elem,num_int_points)
    assert SENER.shape==(num_elem,num_int_points)

    #Ensure that all arrays are c contiguous
    Coords=np.ascontiguousarray(Coords,dtype=np.float64)
    rho=np.ascontiguousarray(rho,dtype=np.float64)
    Element_U=np.ascontiguousarray(Element_U,dtype=np.float64)
    Element_V=np.ascontiguousarray(Element_V,dtype=np.float64)
    Element_A=np.ascontiguousarray(Element_A,dtype=np.float64)
    S_vec=np.ascontiguousarray(S_vec,dtype=np.float64)
    PENER=np.ascontiguousarray(PENER,dtype=np.float64)
    SENER=np.ascontiguousarray(SENER,dtype=np.float64)
    ################################################

    num_elem_=ctypes.c_size_t(Coords.shape[0])
    Coords_p=Coords.ctypes.data_as(ctypes.c_void_p)
    rho_p   =rho.ctypes.data_as(ctypes.c_void_p)
    Element_U_p=Element_U.ctypes.data_as(ctypes.c_void_p)
    Element_V_p=Element_V.ctypes.data_as(ctypes.c_void_p)
    Element_A_p=Element_A.ctypes.data_as(ctypes.c_void_p)
    S_vec_p=S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p=PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p=SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force=np.zeros((Coords.shape[0],Coords.shape[1],3))
    Conf_Force_p=Conf_Force.ctypes.data_as(ctypes.c_void_p)

    lib.Integration_C3D10_dynamic(num_elem_,Coords_p,rho_p,Element_U_p,Element_V_p,Element_A_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    return Conf_Force

def calc_Conf_Force_C3D10R_static(Coords,Element_U,S_vec,PENER,SENER,method='mbf'):
    """
    Inputs:
        Coords.........Element_Nodal_Coordinates [num_elem,numNodes,3]
        Element_U......Element_Nodal_Displacements [num_elem,numNodes,3]
        S_vec......... Stress at integration points [num_elem,num_int_points,6]
        PENER..........Plastic strain energy [num_elem,num_int_points]
        SENER..........Elastic strain energy [num_elem,num_int_points]
    """
    #This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem=Coords.shape[0]
    numNodes=10
    num_int_points=1

    assert Coords.shape==(num_elem,numNodes,3)
    assert Element_U.shape==(num_elem,numNodes,3)
    assert S_vec.shape==(num_elem,num_int_points,6)
    assert PENER.shape==(num_elem,num_int_points)
    assert SENER.shape==(num_elem,num_int_points)

    #Ensure that all arrays are c contiguous
    Coords=np.ascontiguousarray(Coords,dtype=np.float64)
    Element_U=np.ascontiguousarray(Element_U,dtype=np.float64)
    S_vec=np.ascontiguousarray(S_vec,dtype=np.float64)
    PENER=np.ascontiguousarray(PENER,dtype=np.float64)
    SENER=np.ascontiguousarray(SENER,dtype=np.float64)
    ################################################

    num_elem_=ctypes.c_size_t(Coords.shape[0])
    Coords_p=Coords.ctypes.data_as(ctypes.c_void_p)
    Element_U_p=Element_U.ctypes.data_as(ctypes.c_void_p)
    S_vec_p=S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p=PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p=SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force=np.zeros((Coords.shape[0],Coords.shape[1],3))
    Conf_Force_p=Conf_Force.ctypes.data_as(ctypes.c_void_p)

    if method=='mbf':
        lib.Integration_C3D10R_static_mbf(num_elem_,Coords_p,Element_U_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    elif method=='dbf':
        lib.Integration_C3D10R_static_dbf(num_elem_,Coords_p,Element_U_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    else:
        print("Method not found!")
    return Conf_Force

def calc_Conf_Force_C3D10R_dynamic(Coords,rho,Element_U,Element_V,Element_A,S_vec,PENER,SENER):
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
    #This is for safety, the C code will crash without warnings on wrong inputs...
    ################################################
    num_elem=Coords.shape[0]
    numNodes=10
    num_int_points=1

    assert Coords.shape==(num_elem,numNodes,3)
    assert rho.shape==(num_elem,)
    assert Element_U.shape==(num_elem,numNodes,3)
    assert Element_V.shape==(num_elem,numNodes,3)
    assert Element_A.shape==(num_elem,numNodes,3)
    assert S_vec.shape==(num_elem,num_int_points,6)
    assert PENER.shape==(num_elem,num_int_points)
    assert SENER.shape==(num_elem,num_int_points)

    #Ensure that all arrays are c contiguous
    Coords=np.ascontiguousarray(Coords,dtype=np.float64)
    rho=np.ascontiguousarray(rho,dtype=np.float64)
    Element_U=np.ascontiguousarray(Element_U,dtype=np.float64)
    Element_V=np.ascontiguousarray(Element_V,dtype=np.float64)
    Element_A=np.ascontiguousarray(Element_A,dtype=np.float64)
    S_vec=np.ascontiguousarray(S_vec,dtype=np.float64)
    PENER=np.ascontiguousarray(PENER,dtype=np.float64)
    SENER=np.ascontiguousarray(SENER,dtype=np.float64)
    ################################################

    num_elem_=ctypes.c_size_t(Coords.shape[0])
    Coords_p=Coords.ctypes.data_as(ctypes.c_void_p)
    rho_p   =rho.ctypes.data_as(ctypes.c_void_p)
    Element_U_p=Element_U.ctypes.data_as(ctypes.c_void_p)
    Element_V_p=Element_V.ctypes.data_as(ctypes.c_void_p)
    Element_A_p=Element_A.ctypes.data_as(ctypes.c_void_p)
    S_vec_p=S_vec.ctypes.data_as(ctypes.c_void_p)
    PENER_p=PENER.ctypes.data_as(ctypes.c_void_p)
    SENER_p=SENER.ctypes.data_as(ctypes.c_void_p)

    Conf_Force=np.zeros((Coords.shape[0],Coords.shape[1],3))
    Conf_Force_p=Conf_Force.ctypes.data_as(ctypes.c_void_p)

    lib.Integration_C3D10R_dynamic(num_elem_,Coords_p,rho_p,Element_U_p,Element_V_p,Element_A_p,S_vec_p,PENER_p,SENER_p,Conf_Force_p)
    return Conf_Force

