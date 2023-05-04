import numpy as np


def tensor_from_vector(vector):
    vector = np.asarray(vector)
    dimensions = vector.shape
    d_vec = dimensions[-1]
    dim = dimensions[:-1]

    if d_vec == 3:
        tensor = np.zeros(list(dim) + [2, 2], dtype=float)
        tensor[..., 0, 0] = vector[..., 0]
        tensor[..., 1, 1] = vector[..., 1]
        tensor[..., 0, 1] = vector[..., 2]
        tensor[..., 1, 0] = vector[..., 2]

    elif d_vec == 4:
        tensor = np.zeros(list(dim) + [3, 3], dtype=float)
        tensor[..., 0, 0] = vector[..., 0]
        tensor[..., 1, 1] = vector[..., 1]
        tensor[..., 2, 2] = vector[..., 2]
        tensor[..., 0, 1] = vector[..., 3]
        tensor[..., 1, 0] = vector[..., 3]

    elif d_vec == 6:
        tensor = np.zeros(list(dim) + [3, 3], dtype=float)
        tensor[..., 0, 0] = vector[..., 0]
        tensor[..., 1, 1] = vector[..., 1]
        tensor[..., 2, 2] = vector[..., 2]
        tensor[..., 0, 1] = vector[..., 3]
        tensor[..., 1, 0] = vector[..., 3]
        tensor[..., 0, 2] = vector[..., 4]
        tensor[..., 2, 0] = vector[..., 4]
        tensor[..., 1, 2] = vector[..., 5]
        tensor[..., 2, 1] = vector[..., 5]

    else:
        raise NotImplementedError("d_vec=" + str(d_vec))

    return tensor
