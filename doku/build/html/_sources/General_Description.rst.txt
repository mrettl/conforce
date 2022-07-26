.. _Interface_description:

API description
===============

:math:`\newcommand{\mytensor}[1] {\boldsymbol{\mathrm{#1}}}`
:math:`\newcommand{\myjaci}[2]   {\displaystyle \sum^i\frac{\partial \mathrm{N}^{\,i}}{\partial #2} #1^{\,i}}`
:math:`\newcommand{\mynderiv}[2] {\displaystyle \frac{\partial\mathrm{N}^{#1}}{\partial\mathrm{r}_{#2}}}`

The derivation of configurational forces (CF) is defined in a symbolic way for all supported element types from a small
set of governing equations shown bellow. For a detailed description of the derivation have a look at :ref:`Example`.

:math:`\mytensor{g}=\int_{V}{\frac{\partial\mytensor{N}}{\partial\mytensor{x}}{\mytensor{C}}^\mathrm{T}}\det\left(\mytensor{J}\right)dV`

:math:`\mytensor{C} = \phi\mytensor{I}-\mytensor{F}^\mathrm{T}\mytensor{S}`

Therefore, the users does not need to implement all details which are necessary to calculate CF from a 
discretized FE- model themselves. A detailed description of the symbolic derivation, including all necessary steps and used auxiliary functions 
is described in :ref:`Example`.


This section deals with the interface and describes the input and output arrays. 
In contrary to the symbolic formulation, the numerical implementations are element-dependent. 
All element-specific functions are depending on the respective shape functions and expression derived from it, like the Jacobian.

Therefore, the post-processing routine has to be called for every element type separately. 
The C-function :code:`Conf_Forces_API void Integration_{Element Type}_static_{mbf/dbf}` as well as the wrapped 
Python function :code:`void calc_Conf_Force_{Element Type}_static(Coords,Element_U,S_vec,PENER,SENER,method='mbf')` are evaluating CFs on element nodal position.


To combine the nodal CFs from all adjacent elements, nodal unique values can be calculated by calling the C-function :code:`void calc_Nodal_C` or 
the Python function :code:`calc_Nodal`.
Note that node labels within the whole model must be unique. A part/assembly structure is therefore only possible if a node label only occurs once in a model. 
In Abaqus CAE, the option “Do not use parts and assemblies in input files” is recommended to avoid this issue entirely.

For each element type, two functions in C - language are generated.
To make the implementations also available in Python, a wrapper for all public functions is automatically generated.
This wrapper also ensures to pass the inputs with the appropriate data type and shape to the C-functions.
First, the Python-API is explained below, followed by the C API.


Python - Interface
------------------


All Python functions are just wrappers for the corresponding C-functions. The only dependencies are ctypes, which is a standard Python package,
and NumPy. This also ensures that this code is usable with a wide range of Python versions. The interface has been tested on Linux and Windows operating systems.
As with the C-functions below, there is a separate function for all supported element types.
All input arrays are checked for consistancy (array-shape) and the appropriate data type is ensured before calling the C-API.

**Calculation of configurational forces**

.. code-block:: Python


    def calc_Conf_Force_{Element Type}_static(Coords,Element_U,S_vec,PENER,SENER,method='mbf'):
        """
        Parameters
        ----------
            Coords : (num_elem,numNodes,3) nd-array
                global nodal coordinates of the elements 
            Element_U : (num_elem,numNodes,3) nd-array
                displacement at the element nodes
            S_vec : (num_elem,num_int_points,6) nd-array
                chauchy stress at integration points 
                stress tensor in the order of(S11,S22,S33,S12,S13,S23)
            PENER : (num_elem,num_int_points) nd-array
                plastic energy density at integration points
            SENER : (num_elem,num_int_points) nd-array
                elastic energy density at integration points
        
        Returns
        -------
            Conf_Force : (number_of_elements,number_of_nodes_per_element,3) nd-array
              configurational forces at element nodal position
        """


**Calculation of nodal unique values**

As already mentioned, the function call :code:`calc_Conf_Force_{Element type}` is element-dependent. 
If there are multiple element types in the model, the corresponding function for each element type has to be called seperately.
To combine the CFs on the nodes from all adjecent elements, nodal unique values can be calculated by calling the function :code:`calc_Nodal`. 
This function accepts both Numpy nd-arrays if only one element type is present in the model or a list of Numpy nd-arrays for multiple element types. Note that node labels within the model must be unique. 
A part/assembly structure is therefore only possible if a node label only occurs once in a model. 
In Abaqus CAE, the option "Do not use parts and assemblies in input files" is recommended to avoid this issue entirely.


.. code-block:: Python


    def calc_Nodal(Element_Connectivity,CF):
        """
        Parameters
        ----------
        Element_Connectivity : list of nd-arrays or nd-array
            connectivity of the element nodal position as nd-array, if multiple elements 
            should be evaluated a list of nd-arrays must be passed
        CF_in : list of np.arrays or np.array
            configurational forces at element nodal position as nd-array, 
            if multiple elements should be evaluated
            a list of nd-arrays must be passed
        Returns
        -------
        Node_labels : (number_of_unique_nodes) nd-array
            node labels coresponding to the CF_out array bellow
        CF_out : (number_of_unique_nodes,3) nd-array
            configurational forces at nodes
    """


C - Interface
-------------


**Inner Part of the configurational forces integral**

This function calculates the inner part of the configurational forces integral :math:`\mytensor{f}` 
using the motion or displacement-based formulation within an element depending on the natural coordinates :math:`r , s , t`.

:math:`\mytensor{f} = \displaystyle \frac{\partial\mytensor{N}}{\partial\mytensor{x}} \, \mytensor{C}^\mathrm{T} \, \mathrm{det}(\mytensor{J})`

It is not a public function and is only meant to be called by a numerical gauss integration function. 

.. code-block:: c

   void {Element Type}_static_{mbf/dbf}(double *rst,double *coord,double *Element_U,double *S,double PENER,
      double SENER, double *res_0)
   
   /*
   double *rst :        array of size 3 
                        natural coordinates, where the function is evaluated.
   
   double *coord :      array of size number_of_nodes*3
                        global nodal coordinates of the element in the order(N1_x,N1_y,N1_z,N2_x,...)
   
   double *Element_U :  array of size number_of_nodes*3
                        displacement of the element nodes in the order (U1_x,U1_y,U1_z,U2_x,...)
   
   double *S :          array of size number_of_integration_points*6
                        chauchy stress at integration points (S1_11,S1_22,S1_33,S1_12,S1_13,S1_23,S2_11,...)
   
   double PENER :       value
                        plastic energy density at an integration point
   
   double SENER :       value
                        elastic energy density at an integration point
   
   double *res_0 :      array of size 3
                        result of the inner part of the integral
                        evaluated at given natural coordinates
   */


**Calculation of configurational forces**

This element-specic function calculates the configurational forces at element nodal position by using Gauss integration. The integration points and weights 
are hard-coded, therefore it isn't necessary to pass them to the function. 
This is a public function and exported from the shared library.

.. code-block:: c

   Conf_Forces_API void Integration_{Element Type}_static_{mbf/dbf}(size_t num_elem,double Coords[][4][3],
      double Element_U[][_][_], double S[][_][6],double PENER[][_],double SENER[][_],double Conf_Force[][_][3])
   
   /*
   size_t num_elem :          value
                              number of elements
   
   double Coords[][_][_] :    array with shape [number_of_elements][number_of_nodes_per_element][3]
                              global nodal coordinates of the elements 
   
   double Element_U[][_][_] : array with shape [number_of_elements][number_of_nodes_per_element][3]
                              displacement at the element nodes
   
   double S[][_][_] :         array of shape [number_of_elements][number_of_integration_points][6]
                              chauchy stress at integration points 
                              stress tensor in the order of(S1_11,S1_22,S1_33,S1_12,S1_13,S1_23)
   
   double PENER[][4] :        array of shape [number_of_elements][number_of_integration_points]
                              plastic energy density at integration points
   
   double SENER[][4] :        array of shape [number_of_elements][number_of_integration_points]
                              elastic energy density at integration points
   
   double Conf_Force[][_][3]: array of shape [number_of_elements][number_of_nodes_per_element][3]
                              configurational forces at element nodal position
   */


**Calculation of unique nodal values**

As already mentioned, the function call :code:`void Integration_{Element Type}_static_{mbf/dbf}` is element-dependent. 
If there are multiple element types in the model, the corresponding function for each element type has to be called seperately.
To combine the CFs on the nodes from all adjecent elements, nodal unique values can be calculated by calling the function :code:`void calc_Nodal_C`. 

The array :code:`int64_t inverse[]` defines how to reconstruct an element-connectivity array from an array of unique node labels.
If there are multiple element-types in the model, the element connectivity arrays and the element nodal configurational force vectors of all elements have to be concatenated. 
In NumPy the inverse array can be created by :code:`unique_arr,inverse = np.unique(ar,return_inverse=True)`.

A part/assembly structure is therefore only possible if a node label only occurs once in a model. 
In Abaqus CAE, the option "Do not use parts and assemblies in input files" is recommended to avoid this issue entirely.

.. code-block:: c

   Conf_Forces_API void calc_Nodal_C(size_t num_elem_nodal,double CF_out[][3],double CF[][3],int64_t inverse[])
   
   /*
   size_t num_elem :          value
                              number of elements
   
   double CF_out[][3] :       array with shape [number_of_unique_nodelabels][3]
                              nodal unique values (Output)
   
   double CF[][3] :           array with shape [number_of_elements*number_of_nodes_per_element][3]
                              values at element nodal position (Input)
   
   int64_t inverse[] :        array of shape [number_of_elements*number_of_nodes_per_element]
                              array of indices which reconstructs the element nodal node labels from unique
                              node labels
   */