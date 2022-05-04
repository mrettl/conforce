.. _Example:

Symbolic implementation of Configurational Forces for the elasto-platic case
============================================================================


**Calculation of the shape function coefficients**

    >>> shapeFuncCoef_bild,shapeFuncCoef_int,shapeFuncCoef_bild_diff,shapeFuncCoef_int_diff=\
    >>> generate_shape_func_coeff(bild_points,int_points,poly_power)



**Definitions of neccesary symbols**

    >>> num_nodes      = shapeFuncCoef_bild.shape[0]
    >>> num_int_points = int_points.shape[0]
    >>> r, s, t        = sy.symbols("r s t")
    >>> rst            = sy.Matrix((r, s, t))
    >>> coord          = sy.MatrixSymbol('coord', num_nodes, 3)
    >>> Element_U      = sy.MatrixSymbol('U', num_nodes, 3)
    >>> SENER,PENER    = sy.symbols("SENER PENER")
    >>> S_Ten          = sy.symbols("S11 S22 S33 S12 S13 S23")

:math:`\newcommand{\mytensor}[1] {\boldsymbol{\mathrm{#1}}}`
:math:`\newcommand{\myjaci}[2]   {\displaystyle \sum^i\frac{\partial \mathrm{N}^{\,i}}{\partial #2} #1^{\,i}}`
:math:`\newcommand{\mynderiv}[2] {\displaystyle \frac{\partial\mathrm{N}^{#1}}{\partial\mathrm{r}_{#2}}}`



**Generation of shape functions**

Definition of symbolic shape functions for the defined element.

:math:`\mytensor{N}=\left[\begin{array}{ccc}\mathrm{N}^1 \\ \vdots \\ \mathrm{N}^i \end{array}\right]`

 >>> shapeFunc=get_shapeFunc(shapeFuncCoef_bild,r,s,t)



**Calculation of the Jacobian matrix**

Additionally the inverse jacobian matrix and the jacobi determinant is calculated.


:math:`\mytensor{J}= \left[\begin{array}{cccc} \myjaci{x_1}{r_1} &  \myjaci{x_2}{r_1} &  \myjaci{x_3}{r_1} \\ \myjaci{x_1}{r_2} &  \myjaci{x_2}{r_2} &  \myjaci{x_3}{r_2} \\ \myjaci{x_1}{r_3} &  \myjaci{x_2}{r_3} &  \myjaci{x_3}{r_3} \\ \end{array} \right]`

    >>> jac = get_jac(shapeFunc,coord,r,s,t)
    >>> jac_inv,jacobi_det = get_det_and_inv(jac)



**Determining the derivivative of the displacements with respect to x and deformation gradient**

:math:`\displaystyle \frac{\partial\mytensor{u}}{\partial\mytensor{r}}= \left[\begin{array}{cccc} \myjaci{u_1}{r_1} &  \myjaci{u_2}{r_1} &  \myjaci{u_3}{r_1} \\ \myjaci{u_1}{r_2} &  \myjaci{u_2}{r_2} &  \myjaci{u_3}{r_2} \\ \myjaci{u_1}{r_3} &  \myjaci{u_2}{r_3} &  \myjaci{u_3}{r_3} \\ \end{array} \right]`

    >>> Jacobi_Element_U=get_jac(shapeFunc,Element_U,r,s,t)
:math:`\displaystyle \frac{\partial\mytensor{u}}{\partial\mytensor{x}}= \left ( \mytensor{J^{-1}} \displaystyle \frac{\partial\mytensor{u}}{\partial\mytensor{r}} \right )^\mathrm{T}`

    >>> dU_dx=(jac_inv*Jacobi_Element_U).T
:math:`\mytensor{F} = \mytensor{I}+\displaystyle \frac{\partial\mytensor{u}}{\partial\mytensor{x}}`

    >>> Def_grad=sy.Matrix(np.eye(3))+dU_dx



**Determining the derivivative of the shape functions with respect to x**

:math:`\displaystyle \frac{\partial\mytensor{N}}{\partial\mytensor{r}} = \left[\begin{array}{cccc} \mynderiv{1}{1} & \mynderiv{1}{2} & \mynderiv{1}{3} \\ \vdots & \vdots & \vdots \\ \mynderiv{i}{1} & \mynderiv{i}{2} &  \mynderiv{i}{3} \end{array} \right]`

    >>> dN_rst=calc_dN_rst(shapeFunc,r,s,t)
    
:math:`\displaystyle \frac{\partial\mytensor{N}}{\partial\mytensor{x}} = \displaystyle \frac{\partial\mytensor{N}}{\partial\mytensor{r}} \mytensor{J}^\mathrm{-T}`

    >>> dN_dxyz=dN_rst*jac_inv.T



**Calculation of the first Piola Kirchhoff stress**

:math:`\mytensor{P} = \mathrm{det}(\mytensor{F}) \, \mytensor{S} \,\mytensor{F}^{-T}`

    >>> S=gen_2D_Ten_from_vec(S_Ten)
    >>> Def_grad_inv,Def_grad_det=get_det_and_inv(Def_grad)
    >>> Piola_1= Def_grad_det*S*Def_grad_inv.T



**Calculation of the Configurational stress**

:math:`\phi = \phi_\mathrm{pl} + \phi_\mathrm{el}`

    >>> ENER = SENER+PENER

Motion based formulation:

:math:`\mytensor{C} = \mathrm{\phi} \, \mytensor{I} - \mytensor{F}^\mathrm{T} \, \boldsymbol{\mathrm{P}}`

    >>> C = ENER*sy.Matrix(np.eye(3))-Def_grad.T*Piola_1

Displacement based formulation:

:math:`\mytensor{C} = \mathrm{\phi} \, \mytensor{I} - \displaystyle \left (\frac{\partial\mytensor{u}}{\partial\mytensor{x}} \right)^\mathrm{T} \, \boldsymbol{\mathrm{P}}`

    >>> C = ENER*sy.Matrix(np.eye(3))-dU_dx.T*Piola_1



**Determining inner part of the integral**

:math:`\mytensor{f} = \displaystyle \frac{\partial\mytensor{N}}{\partial\mytensor{x}} \, \mytensor{C}^\mathrm{T} \, det(\mytensor{J})`

    >>> C_Force=dN_dxyz*C.T*jacobi_det