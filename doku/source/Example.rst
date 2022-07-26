.. _Example:


.. py:currentmodule:: ConF3D.Auxiliary_functions

Symbolic implementation of configurational forces
=================================================

:math:`\newcommand{\mytensor}[1] {\boldsymbol{\mathrm{#1}}}`
:math:`\newcommand{\myjaci}[2]   {\displaystyle \sum^i\frac{\partial \mathrm{N}^{\,i}}{\partial #2} #1^{\,i}}`
:math:`\newcommand{\mynderiv}[2] {\displaystyle \frac{\partial\mathrm{N}^{#1}}{\partial\mathrm{r}_{#2}}}`

In the following, the main part of the derivation of configurational forces implementation for the elastic and elastic-plastic is shown.
Configurational forces :math:`\mytensor{g}` are defined depending on the shape functions :math:`\mytensor{N}`, the cartesian coordinates :math:`\mytensor{x}`, 
the configurational stress :math:`\mytensor{C}` and the Jacobi determinant :math:`\det\left(\mytensor{J}\right)` in the following way:

:math:`\mytensor{g}=\int_{V}{\frac{\partial\mytensor{N}}{\partial\mytensor{x}}{\mytensor{C}}^\mathrm{T}}\det\left(\mytensor{J}\right)dV`

In this formula the inner part of the integral must be derived symbolically, the integration is than performed using a automatically generated 
Gauss integration routine. This derivation is shown in detail in the following sections. At the end numeric code is generated for all supported elements.

For clarification each section contains the corresponding code snippets and formulas.


**Calculation of the shape function coefficients**

The elements are defined by their nodal coordinates and integration points in the natural coordinate system.
Shape functions are defined by the occuring polynomial exponents. Therefore, the coefficients of the shape functions have 
to be derived. Note that the node ordering of an element is Abaqus-specific needs to be adapted for different FEM-solvers. 

    >>> shapeFuncCoef_bild,shapeFuncCoef_int,shapeFuncCoef_bild_diff,shapeFuncCoef_int_diff=\
    >>> generate_shape_func_coeff(bild_points,int_points,poly_power)


**Definitions of neccesary symbols**

In this section, all necessary symbols are defined. Note that this definition also specifies the way in which the FE-results have to be passed to
the generated functions. For example, the order of the stress vector can be stated in this code segment.

    >>> num_nodes      = shapeFuncCoef_bild.shape[0]
    >>> num_int_points = int_points.shape[0]
    >>> r, s, t        = sy.symbols("r s t")
    >>> rst            = sy.Matrix((r, s, t))
    >>> coord          = sy.MatrixSymbol('coord', num_nodes, 3)
    >>> Element_U      = sy.MatrixSymbol('U', num_nodes, 3)
    >>> SENER,PENER    = sy.symbols("SENER PENER")
    >>> S_Ten          = sy.symbols("S11 S22 S33 S12 S13 S23")


**Generation of shape functions**

The shape functions :math:`\mytensor{N}` in symbolic represenation depending on the natural coordinates vector :math:`\mytensor{r}` with its components :math:`r \, s \, t`, 
are generated from their coefficients. Of course, the shape functions can be also directly 
provided at this point, but this would be more error-prone and less user user-friendly.

:math:`\mytensor{N}=\left[\begin{array}{ccc}\mathrm{N}^1 \\ \vdots \\ \mathrm{N}^i \end{array}\right]`

 >>> shapeFunc=get_shapeFunc(shapeFuncCoef_bild,r,s,t)


**Calculation of the Jacobi matrix**

The Jacobi matrix :math:`\mytensor{J}` is calculated as shown in the following formula. The :math:`\mytensor{r}` vector represents the natural coordinates of a node, 
the :math:`\mytensor{x}` vector represents the global coordinates.
Additionally, the inverse jacobian matrix and the jacobi determinant is calculated.


:math:`\mytensor{J}= \left[\begin{array}{cccc} \myjaci{x_1}{r_1} &  \myjaci{x_2}{r_1} &  \myjaci{x_3}{r_1} \\ \myjaci{x_1}{r_2} &  \myjaci{x_2}{r_2} &  \myjaci{x_3}{r_2} \\ \myjaci{x_1}{r_3} &  \myjaci{x_2}{r_3} &  \myjaci{x_3}{r_3} \\ \end{array} \right]`

    >>> jac = get_jac(shapeFunc,coord,r,s,t)
    >>> jac_inv,jacobi_det = get_det_and_inv(jac)


**Determining the derivivative of the displacements with respect to** :math:`\mytensor{x}` **and deformation gradient**

At first, the derivivative of the nodal displacements :math:`\mytensor{u}` with respect to the natural coordinates :math:`\mytensor{r}` is calculated simmilar to the Jacobi matrix above.

:math:`\displaystyle \frac{\partial\mytensor{u}}{\partial\mytensor{r}}= \left[\begin{array}{cccc} \myjaci{u_1}{r_1} &  \myjaci{u_2}{r_1} &  \myjaci{u_3}{r_1} \\ \myjaci{u_1}{r_2} &  \myjaci{u_2}{r_2} &  \myjaci{u_3}{r_2} \\ \myjaci{u_1}{r_3} &  \myjaci{u_2}{r_3} &  \myjaci{u_3}{r_3} \\ \end{array} \right]`

    >>> Jacobi_Element_U=get_jac(shapeFunc,Element_U,r,s,t)

Next, the derivivative of the displacement :math:`\mytensor{u}` with respect to :math:`\mytensor{x}` can be calculated.

:math:`\displaystyle \frac{\partial\mytensor{u}}{\partial\mytensor{x}}= \left ( \mytensor{J^{-1}} \displaystyle \frac{\partial\mytensor{u}}{\partial\mytensor{r}} \right )^\mathrm{T}`

    >>> dU_dx=(jac_inv*Jacobi_Element_U).T

The deformation gradient :math:`\mytensor{F}` is defined by the following equation, with :math:`\mytensor{I}` as the idenedity matrix.

:math:`\mytensor{F} = \mytensor{I}+\displaystyle \frac{\partial\mytensor{u}}{\partial\mytensor{x}}`

    >>> Def_grad=sy.Matrix(np.eye(3))+dU_dx


**Determining the derivivative of the shape functions** :math:`\mytensor{N}` **with respect to** :math:`\mytensor{x}` **from** :math:`\displaystyle \frac{\partial\mytensor{N}}{\partial\mytensor{r}}` **and the Jabobian**

:math:`\displaystyle \frac{\partial\mytensor{N}}{\partial\mytensor{r}} = \left[\begin{array}{cccc} \mynderiv{1}{1} & \mynderiv{1}{2} & \mynderiv{1}{3} \\ \vdots & \vdots & \vdots \\ \mynderiv{i}{1} & \mynderiv{i}{2} &  \mynderiv{i}{3} \end{array} \right]`

    >>> dN_rst=calc_dN_rst(shapeFunc,r,s,t)
    
:math:`\displaystyle \frac{\partial\mytensor{N}}{\partial\mytensor{x}} = \displaystyle \frac{\partial\mytensor{N}}{\partial\mytensor{r}} \mytensor{J}^\mathrm{-T}`

    >>> dN_dxyz=dN_rst*jac_inv.T


**Calculation of the first Piola Kirchhoff stress**

The first Piola Kirchhoff stress :math:`\mytensor{P}` is calculated from the the chauchy stress :math:`\mytensor{S}` and the deformation gradient :math:`\mytensor{F}` by the following equation.

:math:`\mytensor{P} = \mathrm{det}(\mytensor{F}) \, \mytensor{S} \,\mytensor{F}^{-T}`

    >>> S=gen_2D_Ten_from_vec(S_Ten)
    >>> Def_grad_inv,Def_grad_det=get_det_and_inv(Def_grad)
    >>> Piola_1= Def_grad_det*S*Def_grad_inv.T


**Calculation of the configurational stress**

The configurational stress :math:`\mytensor{C}` can be writen in a motion-based and a deformation-based formulation. Both are supported by this package.
The energy density :math:`\phi` represents the sum of the plastic energy density :math:`\phi_\mathrm{pl}` and the strain energy density :math:`\phi_\mathrm{el}`.
In the script which generates the numerical implementations for all supported element types, both the motion based and the deformation based formulations are generated.
In the interface they can be selected by the parameter :func:`method`.

:math:`\phi = \phi_\mathrm{pl} + \phi_\mathrm{el}`

    >>> ENER = SENER+PENER

The motion-based formulation defines the configurational stress :math:`\mytensor{C}` depending on the energy density :math:`\phi`, the idendity matrix :math:`\mytensor{I}`, 
the deformation gradient :math:`\mytensor{F}` and the first Piola Kirchhoff stress :math:`\mytensor{P}`.

:math:`\mytensor{C} = \mathrm{\phi} \, \mytensor{I} - \mytensor{F}^\mathrm{T} \, \boldsymbol{\mathrm{P}}`

    >>> C = ENER*sy.Matrix(np.eye(3))-Def_grad.T*Piola_1

The displacement-based formulation defines the configurational stress :math:`\mytensor{C}` depending on the energy density :math:`\phi`, the idendity matrix :math:`\mytensor{I}`, 
:math:`\displaystyle \frac{\partial\mytensor{u}}{\partial\mytensor{x}}` and the first Piola Kirchhoff stress :math:`\mytensor{P}`.

:math:`\mytensor{C} = \mathrm{\phi} \, \mytensor{I} - \displaystyle \left (\frac{\partial\mytensor{u}}{\partial\mytensor{x}} \right)^\mathrm{T} \, \boldsymbol{\mathrm{P}}`

    >>> C = ENER*sy.Matrix(np.eye(3))-dU_dx.T*Piola_1


**Determining inner part of the configurational force integral**

The inner part of the volume integral in a symbolic definition, shown in the beginning, represents the final result of the derivation. 
By calling the function :code:`lambdify_C`, a numerical implementation can be generated.

:math:`\mytensor{f} = \displaystyle \frac{\partial\mytensor{N}}{\partial\mytensor{x}} \, \mytensor{C}^\mathrm{T} \, \mathrm{det}(\mytensor{J})`

    >>> C_Force=dN_dxyz*C.T*jacobi_det



**Generate the Implementation**

The function :func:`gen_Integration_C_static` generates two C-functions. In the first step :func:`lambdify_C` is called to generate a numerical implementation of the provided expression.
In the second step a function to numerically integrate the expression is generated.

    >>> impl=gen_Integration_C_static(rst,coord,Element_U,S_Ten,PENER,SENER,int_points,int_weights,C_Force,expr_name)