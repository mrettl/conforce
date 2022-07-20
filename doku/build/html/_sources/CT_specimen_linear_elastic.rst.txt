.. _CT_specimen_linear_elastic:


Linear elastic CT-specimen
===========================

This example uses Configurational Forces to determine the crack driving force. The sum of the CFs over a volume in crack direction represents the well-known J integral [1]_.
The dimensions of the CT-specimen are shown in the subsequent figure.
The model is provided as an Abaqus Input file and an Abaqus cae file and consists of fully-integrated bilinear quadrilateral plane strain elements with reduced integration (CPE4R). 


.. figure:: 101_ct_probe_002.svg 

For the partitions F and G, a mapped mesh is used. In area F, the mesh size is 0.15 mm. In area G the mesh size is 0.165 mm on the crack line. 
All other areas are meshed using a free meshing algorithm with a mesh size of 1 mm. 

A linear elastic material with a Young's modulus of :math:`E_{1}=200\,\mathrm{GPa}` and a Poisson's ratio :math:`\nu` of 0.3. is used for the whole specimen.
For this example a displacement of :math:`u_{y}=0.5\,\mathrm{mm}` was used, which leads to a reaction force of :math:`P=52.24\,\mathrm{kN}`.
According to ASTM [2]_, the J integral for a CT-specimen can be calculated in the following way:

:math:`J=\cfrac{K^2 (1 - \nu^2)}{E}`

Where :math:`K` represents the fracture thoughness depending on the thickness :math:`B=B_N=25\,\mathrm{mm}` and the geometry factor :math:`f(a/W)`.

:math:`K=\cfrac{P}{(B B_N W)^{1/2}} f(a/W)`

The geometry factor :math:`f(a/W)` is a function of the crack length :math:`a` and the width :math:`W` of the specimen.

:math:`f(a/W)=\cfrac{(2+a/W)(0.886+4.64 a/W -13.32 (a/W)^2+14.72 (a/W)^3 -5.6 (a/W)^4)}{1-(a/W)^{3/2}}`

Summing up the configurational forces in crack direction over three countours around the crack tip yields the J integral :math:`J_{3}=46.57\,\mathrm{kJ/m^2}`.
This can be considered as converged since 30 contours only leads to a minimal additional change of 0.5% in the J integral. Evaluating the analytic formula for the choosen geometry gives 
:math:`J=47.95\,\mathrm{kJ/m^2}` which shows good agreement to the numerical result.


Working example
--------------

This example demonstrates the evaluation of Configurational Forces from an Abaqus output database. Alternatively, the results of the last increment of the FE-calculation are provided within this documentation.
Therefore, the extraction of the FE results is optional.

**Extraction of necessary input data**

The Abaqus Python script :code:`get_Data_from_abq.py` extracts all necessary data from the Abaqus output database. 
Additionally to the :ref:`Two_phase_bar` example, a node set which defines the area where the J-integral is evaluated is exported.
The compiled functions for Configurational Force evaluation are element-dependent, therefore the results have to be extracted on a per-element basis.
In this example, only linear quadrilateral plane-strain elements with reduced integration are used.

The provided script outputs the data in a numpy :code:`.npz` file format. If some results are not available, e.g. plastic strain for a linear elastic calculation, 
an array filled with zeros with the shape described below have to be generated. The same applies for other inputs too, e.g. the stress vector must always be of shape 6, regardless of the calculation.

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

- Node set
    Numpy nd-array of shape (number of nodes)


**Calculation of Configurational Forces**

In the following step, the script :code:`J_Integral.py` evaluates the configurational forces for all available time steps. The function call :code:`calc_Conf_Force_[Element type]` is element-dependent. 
If there are multiple element types in the model, the corresponding function for each element type has to be called seperately.
First, CFs are calculated at element nodal position. In a subsequent step, nodal unique values can be calculated by calling the function :code:`calc_Nodal`. This function accepts both Numpy nd-arrays
if only one element type is present in the model or a list of Numpy nd-arrays for multiple element types. Note that node labels within the model must be unique. A part/assembly structure is therefore only 
possible if a node label only occurs once in a model. In Abaqus CAE, the option "Do not use parts and assemblies in input files" is recommended to avoid this issue entirely.


    >>> Node_Labels=np.unique(Element_Connectivity)
    >>> # Calculate Configurational Forces
    >>> CF_Nodal=np.empty((Element_U.shape[0],Node_Labels.shape[0],3))
    >>> for i in range(Element_U.shape[0]):
    >>>     # Calculate on Element Nodal Position
    >>>     CF_Element_Nodal=cf.calc_Conf_Force_CPE4R_static(Coords,Element_U[i],S_vec[i],PENER[i],SENER[i],method='dbf')
    >>>     # Get Nodal unique value
    >>>     Node_labels,CF_Nodal[i]=cf.calc_Nodal(Element_Connectivity,CF_Element_Nodal)


**Evaluation of the J-integral**

In this step the configurational forces in crack direction are summed up to calculate the J-integral. 

    >>> # Select the node labels which are in a given node set
    >>> idx=np.isin(Node_labels,eval_Node_Labels,assume_unique=True)
    >>> # Sum configurational forces of all selected nodes
    >>> J_dbf_3 = CF_Nodal[:,idx].sum(axis=1)
    >>>
    >>> # Output the result in x-direction
    >>> print("J-integral over 3 Contours: "+ str(J_dbf_3[-1,0])+ " mJ/mm^2")

References
----------


.. [1] `Rice JR. A Path Independent Integral and the Approximate Analysis of Strain Concentration by Notches and Cracks. Journal of Applied Mechanics 1968;35(2):379–86. <https://doi.org/10.1115/1.3601206>`_
.. [2] `ASTM E1820–05 (2005) Standard test method for measurement of fracture toughness. In: Annual Book of ASTM Standards, vol 03.01. ASTM <https://www.astm.org/e1820-18.html>`_