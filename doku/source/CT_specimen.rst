.. _CT_specimen:


Example 2: CT-specimen
======================

This example uses Configurational Forces to determine the crack driving force. The sum of the CFs over a volume in crack direction represents the well-known J-integral [1]_.
The dimensions of the CT-specimen are shown in Fig.1. In the following two cases with differnet material behaviour are presented and compared with literature. 
In the linear elastic case (a) the classical J-integral is calculated from the configurational forces.
In the elastoplastic case (b) a incremental J-integral is evaluated. In contrast to the conventional J-integral, only the elastic energy density is considered.
The geometry and the mesh are in both cases the same.
The models are provided as Abaqus Input files and an Abaqus cae files and use fully-integrated bilinear quadrilateral plane strain elements (CPE4).


.. figure:: 101_ct_probe_002.svg 
   :alt: Geometry of the specimen
   
   Fig.1: Geometry of the specimen 

For the partitions F and G in Fig.1, a mapped mesh is used. In area F, the mesh size is 0.15 mm. In area G, the mesh size is 0.165 mm on the crack line. 
All other areas are meshed using a free meshing algorithm with a mesh size of 1 mm. 

a) Linear elastic case
----------------------

A linear elastic material with a Young's modulus of :math:`E_{1}=200\,\mathrm{GPa}` and a Poisson's ratio :math:`\nu` of 0.3 is used for the whole specimen.
For this example a displacement of :math:`u_{y}=0.5\,\mathrm{mm}` is used, which leads to a reaction force of :math:`P=56.32\,\mathrm{kN}`.
According to an ASTM standard [2]_, the J-integral for a CT-specimen can be calculated in the following way:

:math:`J=\cfrac{K^2 (1 - \nu^2)}{E}`

Where :math:`K` represents the fracture thoughness depending on the thickness :math:`B=B_N=25\,\mathrm{mm}` and the geometry factor :math:`f(a/W)`.


:math:`K=\cfrac{P}{(B B_N W)^{1/2}} f(a/W)`

The geometry factor :math:`f(a/W)` is a function of the crack length :math:`a` and the width :math:`W` of the specimen.

:math:`f(a/W)=\cfrac{(2+a/W)(0.886+4.64 a/W -13.32 (a/W)^2+14.72 (a/W)^3 -5.6 (a/W)^4)}{1-(a/W)^{3/2}}`

Summing up the configurational forces in crack direction over three countours around the crack tip yields the J-integral :math:`J_{3}=53.99\,\mathrm{kJ/m^2}`.
This can be considered as converged since 30 contours only leads to a minimal additional increase of 1% in the J-integral. Evaluating the analytic formula for the choosen geometry gives 
:math:`J=55.72\,\mathrm{kJ/m^2}`, which agrees well with ASTM.

b) Elastoplastic case
------------------

The presented algorithm is validated based on a similar example from Kolednik et al. [4]_.
The specimen is made from an annealed mild steel (S235) with a Young’s modulus of :math:`E_{1}=200\,\mathrm{GPa}`, Poisson’s ratio of 0.3 and a yield stress of :math:`\sigma_{y}=270\,\mathrm{MPa}`. 
The true stress versus strain curve can be found in the thesis of Schöngrundner [3]_.
Additionally, the data for this curve can be found in the Abaqus input file below the :code:`*Plastic` keyword.
The material is modelled using the incremental plasticity model provided by Abaqus. 
Therefore, the plastic behaviour is modelled using 10 measured points of the true stress-strain curve. 
Beyond the ultimate tensile strength of :math:`\sigma_{u}=426\,\mathrm{MPa}`, the material response is linearly extrapolated up to a strain of 300 % and stress of :math:`\sigma_{u}=1460\,\mathrm{MPa}`. 
To avoid large deformations at the load application points, the area around them is modelled with linear elastic material (same Young’s modulus and Poisson’s ratio as the elastoplastic material)
shown in Fig.1 as dark gray area.


On the :math:`25\,\mathrm{mm}` thick specimen, a displacement of :math:`u_{y}=0.5\,\mathrm{mm}` is applied in 10 time steps. 
The presented algorithm, calculates the incremental far-field J-integral :math:`J_{ep\,far}=29.58\,\mathrm{mJ/mm^2}` which is in good agreement to :math:`J_{ep\,far}^{ref}=30.4\,\mathrm{mJ/mm^2}` from Kolednik.
In contrast to the conventional J-integral, only the elastic energy density must be considered.
The results are consistent with this work, but it should be noted that there is some influence regarding the mesh size and the element type. 

Using the code
--------------

This example demonstrates the evaluation of Configurational Forces from an Abaqus output database. Alternatively, the results of the FE-calculation are provided within this documentation.
Therefore, the extraction of the FE results is optional.

Extraction of necessary input data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Abaqus Python script :code:`get_Data_from_abq.py` extracts all necessary data from the Abaqus output database. 
Additionally to the :ref:`Two_phase_bar`, a node set which defines the area where the J-integral is evaluated is exported. For the linear elastic specimen this is a 
node set defining 3 contours around the crack tip, the elastoplastic specimen this is a node set defining the far-field J-integral, see Fig.1.
The compiled functions for Configurational Force evaluation are element-dependent, therefore the results have to be extracted on a per-element basis.
The example uses only linear quadrilateral plane-strain elements.

The provided script outputs the data in a numpy :code:`.npz` file format. If some results are not available, e.g. plastic strain for a linear elastic calculation, 
an array filled with zeros with the shape described below has to be generated. The same applies for other inputs, e.g. the stress vector must always be of shape 6, regardless if the
model is 3d or 2d.

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


Calculation of Configurational Forces
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the following step, the script :code:`J_Integral.py` evaluates the configurational forces for all available time steps. The function call :code:`calc_Conf_Force_{Element type}` is element-dependent. 
If there are multiple element types in the model, the corresponding function for each element type has to be called seperately.
First, CFs are calculated at element nodal position. To combine the CFs on the nodes from all adjecent elements, nodal unique values can be calculated by calling the function :code:`calc_Nodal`. 
This function accepts both Numpy nd-arrays if only one element type is present in the model or a list of Numpy nd-arrays for multiple element types. Note that node labels within the model must be unique. 
A part/assembly structure is therefore only possible if a node label only occurs once in a model. 
In Abaqus CAE, the option "Do not use parts and assemblies in input files" is recommended to avoid this issue entirely.

**a) Linear elastic example**

    >>> Node_Labels=np.unique(Element_Connectivity)
    >>> # Calculate configurational forces
    >>> CF_Nodal=np.empty((Element_U.shape[0],Node_Labels.shape[0],3))
    >>> for i in range(Element_U.shape[0]):
    >>>     # Calculate on Element Nodal Position
    >>>     CF_Element_Nodal=cf.calc_Conf_Force_CPE4_static(Coords,Element_U[i],S_vec[i],PENER[i],SENER[i],method='dbf')
    >>>     # Get Nodal unique value
    >>>     Node_labels,CF_Nodal[i]=cf.calc_Nodal(Element_Connectivity,CF_Element_Nodal)

**b) Elastoplastic example**

To calculate the incremental J-integral afterwards the plastic energy density :code:`PENER` is set to zero. For comparison to the conventional J-integral this factor can be changed to 1.

    >>> Node_Labels=np.unique(Element_Connectivity)
    >>> # Calculate configurational forces
    >>> CF_Nodal=np.empty((Element_U.shape[0],Node_Labels.shape[0],3))
    >>> for i in range(Element_U.shape[0]):
    >>>     # Calculate on element nodal position
    >>>     CF_Element_Nodal=cf.calc_Conf_Force_CPE4_static(Coords,Element_U[i],S_vec[i],PENER[i]*0.,SENER[i],method='dbf')
    >>>     # Get nodal unique value
    >>>     Node_labels,CF_Nodal[i]=cf.calc_Nodal(Element_Connectivity,CF_Element_Nodal)

Evaluation of the J-integral
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this step the configurational forces in crack direction, which coincides in this case with the x-direction, are summed up to calculate the J-integral. 
The area, where this summation is performed is given by the node set `eval_Node_Labels`.

    >>> # Select the node labels which are in a given node set
    >>> idx=np.isin(Node_labels,eval_Node_Labels,assume_unique=True)
    >>> # Sum configurational forces of all selected nodes
    >>> J_dbf_3 = CF_Nodal[:,idx].sum(axis=1)
    >>>
    >>> # Output the result in x-direction (coincides with crack direction)
    >>> print("J-integral in evaluation region: "+ str(J_dbf_3[-1,0]*(-1.))+ " mJ/mm^2")


References
----------


.. [1] `Rice JR. A Path Independent Integral and the Approximate Analysis of Strain Concentration by Notches and Cracks. Journal of Applied Mechanics 1968;35(2):379–86. <https://doi.org/10.1115/1.3601206>`_
.. [2] `ASTM E1820–05 (2005) Standard test method for measurement of fracture toughness. In: Annual Book of ASTM Standards, vol 03.01. ASTM <https://www.astm.org/e1820-18.html>`_
.. [3] `Ronald Schöngrundner. Numerische Studien zur Ermittlung der risstreibenden Kraft in elastisch-plastischen Materialien bei unterschiedlichen Belastungsbedingungen [Dissertation]. Leoben: Montanuniverstät Leoben; 2010. <https://pure.unileoben.ac.at/portal/de/publications/numerische-studien-zur-ermittlung-der-risstreibenden-kraft-in-elastischplastischen-materialien-bei-unterschiedlichen-belastungsbedingungen(0c6ab65a-3702-4001-87a5-e53456916731).html?customType=theses>`_
.. [4] `Kolednik O, Schöngrundner R, Fischer FD. A new view on J-integrals in elastic–plastic materials. Int J Fract 2014;187(1):77–107. <https://doi.org/10.1007/s10704-013-9920-6>`_