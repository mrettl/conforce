r"""
List of abbreviations:

- `d`: number of space dimensions (2D, 3D)
- `n`: number of nodes
- `ips`: number of integration points

- `R`: (d, 1)-Matrix of reference space coordinates
- `X`: (d, 1)-Matrix of real space coordinates

- `H`: (n, 1)-Matrix of shape functions
- `dH_dR`: (n, d)-Jacobian matrix of H with respect to R: :math:`dh\_dr_{ik} = \partial h_{i} / \partial r_{k}`
- `dX_dR`: (d, d)-Jacobian matrix of X with respect to R: :math:`dx\_dr_{ik}= \partial x_{i} / \partial r_{k}`
- `dH_dX`: (n, d)-Jacobian matrix of H with respect to X: :math:`dh\_dx_{ik} = \partial h_{i} / \partial x_{k}`
- `U`: (d,)-Matrix of displacements in the real space
- `U_at_nodes`: (n, d)-Matrix of displacements in the real space at the nodes
- `dU_dX`: (d, d)-Jacobian matrix of U with respect to X: :math:`du\_dx_{ik} = \partial u_{i} / \partial x_{k}`
- `S`: (d, d)-symmetric Cauchy stress tensor
- `F`: (d, d)-Deformation gradient
- `P`: (d, d)-First Piola-Kirchhoff stress tensor
- `e`: Helmholtz free energy density
- `CS`: (d, d)-Matrix of configurational stresses. Implemented formulations:

    - `mbf`: motion based formulation (Eshelby's formulation)
    - `dbf`: displacement based formulation

- `CF`: (d,)-Matrix of a configurational force in the real space
- `CF_at_nodes`: (n, d)-Matrix of configurational forces in the real space at the nodes
"""
import doctest

if __name__ == '__main__':
    doctest.testmod()
