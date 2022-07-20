.. _Two_phase_bar:


Two-phase bar
=============

In this example, configurational forces (CF) are calulated at an material interface. CF are pointing towards the material with the lower stiffness. 
The model sketch is shown in Figure 1.

.. figure:: 102_Two_phase_bar.svg

The model consists of two quadrilateral plane strain elements with a thickness of 1 mm. The left and right side have a Young's modulus of :math:`E_{1}=210\,\mathrm{GPa}` 
and :math:`E_{2}=E_{1}/2 = 105\,\mathrm{GPa}` respectively. 
The Poisson's ratio is set to zero to match the 1D-case where an analytical solution is available [1]_.
The bar is fixed on the left side, on the right side a displacement of :math:`u=0.1\,\mathrm{mm}` is applied.


:math:`G=\frac{u^2 E_1 E_2 A [E_2 - E_1]}{l^2[E_1+E_2]^{2}} = 11.67\,\mathrm{N}`

Working example
-------------

This example demonstrates the evaluation of configurational forces from an Abaqus output database. Alternatively, the results of the FE-calculation
(nodal coordinates, element connectivity, displacement, stress, plastic energy, elastic energy)
are provided within this documentation. Additionally a Abaqus cae file and a Abaqus input file is provided with this example.
Therefore, the extraction of the FE-results is optional. 

**Extraction of necessary input data**

The Abaqus Python script :code:`get_Data_from_abq.py` extracts all necessary data from the Abaqus output database (odb file).
The compiled functions for configurational force evaluation are element-dependent, therefore the results have to be extracted on a per-element basis.
In this example, only linear quadrilateral plane-strain elements with reduced integration (CPE4R) are used.

The provided script outputs the data in a Numpy .npz file format. If some results are not available, e.g. plastic strain for a linear elastic calculation, 
an array filled with zeros with the shape described below have to be generated. The same applies for other inputs, e.g. the stress vector must always be of shape 6, regardless of the calculation.

- Coords
    Numpy nd-array of shape (number of elements, number of nodes per element,3)

- Element Connectivity 
    Numpy nd-array of shape (number of elements, number of nodes per element)

- Displacements
    Numpy nd-array of shape (number of time-steps, number of elements, number of nodes per element,3)

- Stress
    Numpy nd-array of shape (number of time-steps, number of elements, number of integration points per element,6)

- Plastic energy
    Numpy nd-array of shape (number of time-steps, number of elements, number of integration points per element)

- Strain energy
    Numpy nd-array of shape (number of time-steps, number of elements, number of integration points per element)


**Calculation of configurational forces**

In the following step, the script :code:`Calculate_CF.py` evaluates the configurational forces for all available time steps. The function call :code:`calc_Conf_Force_[Element type]` is element-dependent. 
If there are multiple element types in the model, the corresponding function for each element type has to be called seperately.
First, CFs are calculated at element nodal position. In a subsequent step, nodal unique values can be calculated by calling the function :code:`calc_Nodal`. This function accepts both Numpy nd-arrays
if only one element type is present in the model or a list of Numpy nd-arrays for multiple element types. Note that node labels within the model must be unique. A part/assembly structure is therefore only 
possible if a node label only occurs once in a model. In Abaqus CAE, the option "Do not use parts and assemblies in input files" is recommended to avoid this issue entirely.

    >>> Node_Labels=np.unique(Element_Connectivity)
    >>> #calculate configurational forces
    >>> CF_Nodal=np.empty((Element_U.shape[0],Node_Labels.shape[0],3))
    >>> for i in range(Element_U.shape[0]):
    >>>     # Calculate on element nodal position
    >>>     CF_Element_Nodal=cf.calc_Conf_Force_CPE4R_static(Coords,Element_U[i],S_vec[i],PENER[i],SENER[i],method='dbf')
    >>>     # Get nodal unique value
    >>>     Node_labels,CF_Nodal[i]=cf.calc_Nodal(Element_Connectivity,CF_Element_Nodal)


.. [1] `Kolling S, Mueller R. On configurational forces in short-time dynamics and their computation with an explicit solver. Comput Mech 2005;35(5):392–9. <https://doi.org/10.1007/s00466-004-0627-4>`_
