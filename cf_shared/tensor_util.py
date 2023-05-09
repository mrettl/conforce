"""
This module contains methods to work with first-order tensors (vector) and second-order tensors.
"""

import numpy as np


def tensor_from_abaqus_notation(array):
    """
    Convert a second-order tensor from the abaqus notation into the matrix form.

    A tensor `T` in the abaqus notation is ordered:
        - with 6 components as `(T11, T22, T33, T12, T13, T23)`
        - with 4 components as `(T11, T22, T33, T12)`

    **Examples**

    The second order tensor can be converted to the abaqus notation
    with 6 and 4 components.
    By using 4 components, two components `T13` and `T23` are lost.

    >>> T = np.array([
    ...     [1.0, 0.2, 0.3],
    ...     [0.2, 2.0, 0.4],
    ...     [0.3, 0.4, 3.],
    ... ])
    >>> T_abaqus_6 = abaqus_notation_from_tensor(T, 6)
    >>> T_abaqus_6
    array([1. , 2. , 3. , 0.2, 0.3, 0.4])
    >>> T_abaqus_4 = abaqus_notation_from_tensor(T, 4)
    >>> T_abaqus_4
    array([1. , 2. , 3. , 0.2])

    For this reason, the tensor is reconstructed correctly,
    if 6 components are used, but not if 4 components are used.

    >>> tensor_from_abaqus_notation(T_abaqus_6)
    array([[1. , 0.2, 0.3],
           [0.2, 2. , 0.4],
           [0.3, 0.4, 3. ]])
    >>> tensor_from_abaqus_notation(T_abaqus_4)
    array([[1. , 0.2, 0. ],
           [0.2, 2. , 0. ],
           [0. , 0. , 3. ]])

    :param array: array of shape (..., 4) or (..., 6) in abaqus notation
    :return: array of shape (..., 3, 3) containing tensors in their matrix form.
    """
    array = np.asarray(array)
    dimensions = array.shape
    dim_array = dimensions[-1]
    dim = dimensions[:-1]

    tensor = np.zeros(list(dim) + [3, 3], dtype=float)

    if dim_array == 4:
        tensor[..., 0, 0] = array[..., 0]
        tensor[..., 1, 1] = array[..., 1]
        tensor[..., 2, 2] = array[..., 2]
        tensor[..., 0, 1] = tensor[..., 1, 0] = array[..., 3]

    elif dim_array == 6:
        tensor[..., 0, 0] = array[..., 0]
        tensor[..., 1, 1] = array[..., 1]
        tensor[..., 2, 2] = array[..., 2]
        tensor[..., 0, 1] = tensor[..., 1, 0] = array[..., 3]
        tensor[..., 0, 2] = tensor[..., 2, 0] = array[..., 4]
        tensor[..., 1, 2] = tensor[..., 2, 1] = array[..., 5]

    else:
        raise NotImplementedError("dim_array=" + str(dim_array))

    return tensor


def abaqus_notation_from_tensor(tensor, dim_vector):
    """
    Convert tensors into the abaqus notation.

    .. seealso:: :py:func:`tensor_from_abaqus_notation`

    :param tensor: array of shape (..., 3, 3) containing tensors in their matrix form
    :param dim_vector: int, number of components used in the abaqus notation (4 or 6).
    :return: array of shape (..., 4) or (..., 6) containing tensor in the abaqus notation
    """
    tensor = np.asarray(tensor)
    dimensions = tensor.shape
    dim = dimensions[:-2]

    vector = np.zeros(list(dim) + [dim_vector], dtype=float)

    if dim_vector == 4:
        vector[..., 0] = tensor[..., 0, 0]
        vector[..., 1] = tensor[..., 1, 1]
        vector[..., 2] = tensor[..., 2, 2]
        vector[..., 3] = tensor[..., 0, 1]
        vector[..., 3] = tensor[..., 1, 0]
    
    elif dim_vector == 6:
        vector[..., 0] = tensor[..., 0, 0]
        vector[..., 1] = tensor[..., 1, 1]
        vector[..., 2] = tensor[..., 2, 2]
        vector[..., 3] = tensor[..., 0, 1]
        vector[..., 3] = tensor[..., 1, 0]
        vector[..., 4] = tensor[..., 0, 2]
        vector[..., 4] = tensor[..., 2, 0]
        vector[..., 5] = tensor[..., 1, 2]
        vector[..., 5] = tensor[..., 2, 1]

    else:
        raise NotImplementedError("d_vec=" + str(dim_vector))

    return vector


def rotation_matrix_from_quaternion(Q):
    """
    Create rotation matrices out of quaternions.
    Quaternions have the form

    `Q[..., 0]*i + Q[..., 1]*j + Q[..., 2]*k + Q[..., 3]`.

    **Examples**

    Given the quaternion for a 90-degree rotation along the z-axis,
    the corresponding rotation matrix is computed.

    >>> Q90z = np.array([ 0, 0, 0.7071068, 0.7071068])
    >>> ROT90z = rotation_matrix_from_quaternion(Q90z)

    The rotation matrix is used to rotate the x-Vector forward and
    backward.

    >>> x_vector = np.array([1, 0, 0])
    >>> rotated_vector = do_rotation_vector(ROT90z, x_vector)
    >>> np.round(rotated_vector, 6)
    array([-0.,  1.,  0.])
    >>> np.round(undo_rotation_vector(ROT90z, rotated_vector), 6)
    array([1., 0., 0.])

    A tensor is rotated using :py:func:`do_rotation_tensor` and :py:func:`undo_rotation_tensor`.

    >>> tensor = np.array([
    ...     [1.0, 0.2, 0.3],
    ...     [0.2, 2.0, 0.4],
    ...     [0.3, 0.4, 3.0],
    ... ])
    >>> rotated_tensor = do_rotation_tensor(ROT90z, tensor)
    >>> np.round(rotated_tensor, 6)
    array([[ 2. , -0.2, -0.4],
           [-0.2,  1. ,  0.3],
           [-0.4,  0.3,  3. ]])
    >>> np.round(undo_rotation_tensor(ROT90z, rotated_tensor), 6)
    array([[1. , 0.2, 0.3],
           [0.2, 2. , 0.4],
           [0.3, 0.4, 3. ]])

    :param Q: array of shape (..., 4) containing the quaternions
    :return: array of shape (..., 3, 3) containing the rotation matrices.
    """
    Q = np.asarray(Q)
    dim = list(Q.shape[:-1])

    q0 = Q[..., 3].reshape(dim + [1, 1])
    q1 = Q[..., 0].reshape(dim + [1, 1])
    q2 = Q[..., 1].reshape(dim + [1, 1])
    q3 = Q[..., 2].reshape(dim + [1, 1])

    # First row of the rotation matrix
    r00 = 1 - 2*(q2*q2 + q3*q3)
    r01 = -2*q0*q3 + 2*q1*q2
    r02 = 2*q0*q2+2*q1*q3
    r0_ = np.concatenate([r00, r01, r02], axis=-1)

    # Second row of the rotation matrix
    r10 = 2*q0*q3 + 2*q1*q2
    r11 = 1 - 2*(q1*q1 + q3*q3)
    r12 = -2*q0*q1 + 2*q2*q3
    r1_ = np.concatenate([r10, r11, r12], axis=-1)

    # Third row of the rotation matrix
    r20 = -2*q0*q2+2*q1*q3
    r21 = 2*q0*q1 + 2*q2*q3
    r22 = 1 - 2*(q1*q1 + q2*q2)
    r2_ = np.concatenate([r20, r21, r22], axis=-1)

    # dim X 3 X 3 rotation matrix
    ROT = np.concatenate([
        r0_,
        r1_,
        r2_
    ], axis=-2)
    return ROT


def do_rotation_vector(ROT, vectors):
    """
    Rotate vectors using rotation matrices.

    :param ROT: array of shape (..., 3, 3) containing rotation matrices
    :param vectors: array of shape (..., 3) containing vectors
    :return: array of shape (..., 3) containing the rotated vectors
    """
    return np.einsum("...ij,...j", ROT, vectors)


def undo_rotation_vector(ROT, rotated_vectors):
    """
    Revert the rotation of vectors done by the given rotation matrices.

    :param ROT: array of shape (..., 3, 3) containing rotation matrices
    :param rotated_vectors: array of shape (..., 3) containing the rotated vectors
    :return: array of shape (..., 3) containing vectors
    """
    return np.einsum("...ji,...j", ROT, rotated_vectors)


def do_rotation_tensor(ROT, tensors):
    """
    Rotate second-order tensors using rotation matrices.

    :param ROT: array of shape (..., 3, 3) containing rotation matrices
    :param tensors: array of shape (..., 3) containing tensors
    :return: array of shape (..., 3, 3) containing the rotated tensors
    """
    return np.einsum("...ij,...jk,...lk", ROT, tensors, ROT)


def undo_rotation_tensor(ROT, rotated_tensors):
    """
    Revert the rotation of second-order tensors done by rotation matrices.

    :param ROT: array of shape (..., 3, 3) containing rotation matrices
    :param rotated_tensors: array of shape (..., 3) containing the rotated tensors
    :return: array of shape (..., 3, 3) containing tensors
    """
    return np.einsum("...ji,...jk,...kl", ROT, rotated_tensors, ROT)
