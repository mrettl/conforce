"""
Element Definitions
===================

This file contains all the necessary information about the supported Finite Elements in a dictionary.
The number of dimensions is 2 for 2D elements and 3 for 3D elements.

For the definition of an element the following arrays are necessary:

**Coordinates of the integration points in the reference coordinate system**
    (number of integration points, number of dimensions) nd-array

**Integration weights**
    (number of integration points) nd-array

**Coordinates of the nodes in the reference coordinate system**
    (number of nodes, number of dimensions) nd-array

**Polynomial powers of shape functions**
    (at least number of nodes, number of dimensions) nd-array

"""
import numpy as np

# Integration Points
int_points = {}
# Integration Weights
int_weights = {}
# Natural Coordinates of the Element Nodes
ref_nodes = {}
# Polynomial exponents of shape functions
poly_power = {}

# Calculation of gaussian weights
EdW3 = 1. / np.sqrt(3)
Ed3 = 1. / 3.
W3d5 = np.sqrt(3. / 5.)
tri41 = 1. / 4. + 3 * np.sqrt(5) / 20
tri42 = 1. / 4. - np.sqrt(5) / 20
tri31 = 1. / 6.
tri32 = 4. / 6.
B5d9 = 5. / 9.
B8d9 = 8. / 9.


###################################################################################
# 2D - elements ###################################################################
###################################################################################
# CPE4
int_points['CPE4'] = np.array((
    (-EdW3, -EdW3), (EdW3, -EdW3),
    (-EdW3, EdW3), (EdW3, EdW3)
))
int_weights['CPE4'] = np.array((1., 1., 1., 1.))
ref_nodes['CPE4'] = np.array((
    (-1., -1.), (1., -1.),
    (1., 1.), (-1., 1.)
))
poly_power['CPE4'] = np.array((
    (0, 0), (1, 0),
    (0, 1), (1, 1)
))

# CPE4R
int_points['CPE4R'] = np.array(((0., 0.),))
int_weights['CPE4R'] = np.array((4.,))
ref_nodes['CPE4R'] = ref_nodes['CPE4']
poly_power['CPE4R'] = poly_power['CPE4']

# CPE8
int_points['CPE8'] = np.array((
    (-W3d5, -W3d5), (0., -W3d5), (W3d5, -W3d5),
    (-W3d5, 0.), (0., 0.), (W3d5, 0.),
    (-W3d5, W3d5), (0., W3d5), (W3d5, W3d5)
))
int_weights['CPE8'] = np.array((
    B5d9 * B5d9, B8d9 * B5d9, B5d9 * B5d9,
    B8d9 * B5d9, B8d9 * B8d9, B8d9 * B5d9,
    B5d9 * B5d9, B8d9 * B5d9, B5d9 * B5d9
))
ref_nodes['CPE8'] = np.array((
    (-1., -1.), (1., -1.), (1., 1.), (-1., 1.),
    (0., -1.), (1., 0.), (0., 1.), (-1., 0.)
))
poly_power['CPE8'] = np.array((
    (0, 0), (1, 0), (0, 1), (1, 1),
    (2, 0), (0, 2), (2, 1), (1, 2)
))

# CPE8R
int_points['CPE8R'] = np.array((
    (-EdW3, -EdW3), (EdW3, -EdW3),
    (-EdW3, EdW3), (EdW3, EdW3)
))
int_weights['CPE8R'] = np.array((1., 1., 1., 1.))
ref_nodes['CPE8R'] = ref_nodes['CPE8']
poly_power['CPE8R'] = poly_power['CPE8']

# CPE3
int_points['CPE3'] = np.array(((Ed3, Ed3),))
int_weights['CPE3'] = np.array([0.5])
ref_nodes['CPE3'] = np.array((
    (0., 0.),
    (1., 0.),
    (0., 1.)
))
poly_power['CPE3'] = np.array((
    (0, 0),
    (1, 0),
    (0, 1)
))

# CPE6
int_points['CPE6'] = np.array((
    (tri31, tri31),
    (tri32, tri31),
    (tri31, tri32)
))
int_weights['CPE6'] = np.array((1. / 6, 1. / 6, 1. / 6))
ref_nodes['CPE6'] = np.array((
    (0., 0.), (1., 0.), (0., 1.),
    (0.5, 0.), (0.5, 0.5), (0., 0.5)
))
poly_power['CPE6'] = np.array((
    (0, 0),
    (1, 0), (0, 1), (1, 1),
    (0, 2), (2, 0)
))

###################################################################################
# 3D - elements ###################################################################
###################################################################################

# C3D8
int_points['C3D8'] = np.array((
    (-EdW3, -EdW3, -EdW3),
    (EdW3, -EdW3, -EdW3),
    (-EdW3, EdW3, -EdW3),
    (EdW3, EdW3, -EdW3),
    (-EdW3, -EdW3, EdW3),
    (EdW3, -EdW3, EdW3),
    (-EdW3, EdW3, EdW3),
    (EdW3, EdW3, EdW3)
))
int_weights['C3D8'] = np.array((1., 1., 1., 1., 1., 1., 1., 1.))
ref_nodes['C3D8'] = np.array((
    (-1., -1., -1.),
    (1., -1., -1.),
    (1., 1., -1.),
    (-1., 1., -1.),
    (-1., -1., 1.),
    (1., -1., 1.),
    (1., 1., 1.),
    (-1., 1., 1.)
))
poly_power['C3D8'] = np.array((
    (0, 0, 0),
    (1, 0, 0), (0, 1, 0), (0, 0, 1),
    (1, 1, 0), (1, 0, 1), (0, 1, 1),
    (1, 1, 1)
))

# C3D8R
int_points['C3D8R'] = np.array(((0., 0., 0.),))
int_weights['C3D8R'] = np.zeros(1) + 8

ref_nodes['C3D8R'] = ref_nodes['C3D8']
poly_power['C3D8R'] = poly_power['C3D8']

# C3D20
int_points['C3D20'] = np.array((
    (-W3d5, -W3d5, -W3d5),
    (0., -W3d5, -W3d5),
    (W3d5, -W3d5, -W3d5),
    (-W3d5, 0., -W3d5),
    (0., 0., -W3d5),
    (W3d5, 0., -W3d5),
    (-W3d5, W3d5, -W3d5),
    (0., W3d5, -W3d5),
    (W3d5, W3d5, -W3d5),
    (-W3d5, -W3d5, 0.),
    (0., -W3d5, 0.),
    (W3d5, -W3d5, 0.),
    (-W3d5, 0., 0.),
    (0., 0., 0.),
    (W3d5, 0., 0.),
    (-W3d5, W3d5, 0.),
    (0., W3d5, 0.),
    (W3d5, W3d5, 0.),
    (-W3d5, -W3d5, W3d5),
    (0., -W3d5, W3d5),
    (W3d5, -W3d5, W3d5),
    (-W3d5, 0., W3d5),
    (0., 0., W3d5),
    (W3d5, 0., W3d5),
    (-W3d5, W3d5, W3d5),
    (0., W3d5, W3d5),
    (W3d5, W3d5, W3d5)
))
int_weights['C3D20'] = np.array((
    B5d9 * B5d9 * B5d9,
    B8d9 * B5d9 * B5d9,
    B5d9 * B5d9 * B5d9,
    B5d9 * B8d9 * B5d9,
    B8d9 * B8d9 * B5d9,
    B5d9 * B8d9 * B5d9,
    B5d9 * B5d9 * B5d9,
    B8d9 * B5d9 * B5d9,
    B5d9 * B5d9 * B5d9,
    B5d9 * B5d9 * B8d9,
    B8d9 * B5d9 * B8d9,
    B5d9 * B5d9 * B8d9,
    B5d9 * B8d9 * B8d9,
    B8d9 * B8d9 * B8d9,
    B5d9 * B8d9 * B8d9,
    B5d9 * B5d9 * B8d9,
    B8d9 * B5d9 * B8d9,
    B5d9 * B5d9 * B8d9,
    B5d9 * B5d9 * B5d9,
    B8d9 * B5d9 * B5d9,
    B5d9 * B5d9 * B5d9,
    B5d9 * B8d9 * B5d9,
    B8d9 * B8d9 * B5d9,
    B5d9 * B8d9 * B5d9,
    B5d9 * B5d9 * B5d9,
    B8d9 * B5d9 * B5d9,
    B5d9 * B5d9 * B5d9
))
ref_nodes['C3D20'] = np.array((
    (-1., -1., -1.),
    (1., -1., -1.),
    (1., 1., -1.),
    (-1., 1., -1.),
    (-1., -1., 1.),
    (1., -1., 1.),
    (1., 1., 1.),
    (-1., 1., 1.),
    (0., -1., -1.),
    (1., 0., -1.),
    (0., 1., -1.),
    (-1., 0., -1.),
    (0., -1., 1.),
    (1., 0., 1.),
    (0., 1., 1.),
    (-1., 0., 1.),
    (-1., -1., 0.),
    (1., -1., 0.),
    (1., 1., 0.),
    (-1., 1., 0.)
))
poly_power['C3D20'] = np.array((
    (0, 0, 0),
    (1, 0, 0), (0, 1, 0), (0, 0, 1),
    (1, 1, 0), (1, 0, 1), (0, 1, 1),
    (1, 1, 1),
    (2, 0, 0), (0, 2, 0), (0, 0, 2),
    (2, 1, 0), (2, 0, 1), (1, 2, 0), (0, 2, 1), (1, 0, 2), (0, 1, 2),
    (2, 1, 1), (1, 2, 1), (1, 1, 2)
))

# C3D20R
int_points['C3D20R'] = np.array((
    (-EdW3, -EdW3, -EdW3),
    (EdW3, -EdW3, -EdW3),
    (-EdW3, EdW3, -EdW3),
    (EdW3, EdW3, -EdW3),
    (-EdW3, -EdW3, EdW3),
    (EdW3, -EdW3, EdW3),
    (-EdW3, EdW3, EdW3),
    (EdW3, EdW3, EdW3)
))
int_weights['C3D20R'] = np.array((1., 1., 1., 1., 1., 1., 1., 1.))
ref_nodes['C3D20R'] = ref_nodes['C3D20']
poly_power['C3D20R'] = poly_power['C3D20']

# C3D4
int_points['C3D4'] = np.array(((0.25, 0.25, 0.25),))
int_weights['C3D4'] = np.zeros(1) + 1. / 6.
ref_nodes['C3D4'] = np.array((
    (0., 0., 0.),
    (1., 0., 0.),
    (0., 1., 0.),
    (0., 0., 1.)
))
poly_power['C3D4'] = np.array(((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)))

# C3D10
int_points['C3D10'] = np.array((
    (tri42, tri42, tri42),
    (tri41, tri42, tri42),
    (tri42, tri41, tri42),
    (tri42, tri42, tri41)
))
int_weights['C3D10'] = np.array((1. / 24, 1. / 24, 1. / 24, 1. / 24))
ref_nodes['C3D10'] = np.array((
    (0., 0., 0.),
    (1., 0., 0.),
    (0., 1., 0.),
    (0., 0., 1.),
    (0.5, 0., 0.),
    (0.5, 0.5, 0.),
    (0., 0.5, 0.),
    (0., 0., 0.5),
    (0.5, 0., 0.5),
    (0., 0.5, 0.5)
))
poly_power['C3D10'] = np.array((
    (0, 0, 0),
    (1, 0, 0), (0, 1, 0), (0, 0, 1),
    (1, 1, 0), (1, 0, 1), (0, 1, 1),
    (2, 0, 0), (0, 2, 0), (0, 0, 2)
))

# C3D10R
int_points['C3D10R'] = np.array(((0.25, 0.25, 0.25),))
int_weights['C3D10R'] = np.zeros(1) + 1. / 6.
ref_nodes['C3D10R'] = ref_nodes['C3D10']
poly_power['C3D10R'] = poly_power['C3D10']
