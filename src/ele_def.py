"""
Element Definitions
===================

This file contains all the necessary information about the supported Finite Elements in a dictionary.

For the definiton of an element the following arrays are necessary:

**Coordinates of the integration points**
    (number of integration points,3) nd-array
    
    >>> int_points['CPE4'] = np.array(((-EdW3, -EdW3, 0),
    ...                               (EdW3, -EdW3, 0),
    ...                               (-EdW3, EdW3, 0),
    ...                               (EdW3, EdW3, 0)))

**Coordinates of the nodes in natural coordinates**
    (number of nodes,3) nd-array
    
    >>> bild_points['CPE4'] = np.array(((-1., -1., 0),
    ...                                    (1., -1., 0),
    ...                                    (1., 1., 0),
    ...                                    (-1., 1., 0)))

**Polynomial coefficents of shape functions**
    (at least number of nodes,3) nd-array
    
    >>> poly_power['CPE4'] = np.array(((0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0)))

**Integration weights**
    (number of integration points) nd-array
    
    >>> int_weights['CPE4'] = np.array((1., 1., 1., 1.))
"""
import numpy as np

# Integration Points
int_points = {}
# Integration Weights
int_weights = {}
# Natural Coordinates of the Element Nodes
bild_points = {}
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
int_points['CPE4'] = np.array(((-EdW3, -EdW3, 0),
                               (EdW3, -EdW3, 0),
                               (-EdW3, EdW3, 0),
                               (EdW3, EdW3, 0)))
bild_points['CPE4'] = np.array(((-1., -1., 0),
                                (1., -1., 0),
                                (1., 1., 0),
                                (-1., 1., 0)))
poly_power['CPE4'] = np.array(((0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0)))
int_weights['CPE4'] = np.array((1., 1., 1., 1.))

# CPE4R
int_points['CPE4R'] = np.array(((0., 0., 0),))
bild_points['CPE4R'] = np.array(((-1., -1., 0),
                                 (1., -1., 0),
                                 (1., 1., 0),
                                 (-1., 1., 0)))
poly_power['CPE4R'] = np.array(((0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0)))
int_weights['CPE4R'] = np.array((4.,))

# CPE8
int_points['CPE8'] = np.array(((-W3d5, -W3d5, 0),
                               (0., -W3d5, 0),
                               (W3d5, -W3d5, 0),
                               (-W3d5, 0., 0),
                               (0., 0., 0),
                               (W3d5, 0., 0),
                               (-W3d5, W3d5, 0),
                               (0., W3d5, 0),
                               (W3d5, W3d5, 0)))
int_weights['CPE8'] = np.array((B5d9 * B5d9,
                                B8d9 * B5d9,
                                B5d9 * B5d9,
                                B8d9 * B5d9,
                                B8d9 * B8d9,
                                B8d9 * B5d9,
                                B5d9 * B5d9,
                                B8d9 * B5d9,
                                B5d9 * B5d9))
bild_points['CPE8'] = np.array(((-1., -1., 0.),
                                (1., -1., 0.),
                                (1., 1., 0.),
                                (-1., 1., 0.),
                                (0., -1., 0.),
                                (1., 0., 0.),
                                (0., 1., 0.),
                                (-1., 0., 0.)))
poly_power['CPE8'] = np.array(
    ((0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0), (2, 0, 0), (0, 2, 0), (2, 1, 0), (1, 2, 0)))

# CPE8R
int_points['CPE8R'] = np.array(((-EdW3, -EdW3, 0),
                                (EdW3, -EdW3, 0),
                                (-EdW3, EdW3, 0),
                                (EdW3, EdW3, 0)))
bild_points['CPE8R'] = np.array(((-1., -1., 0),
                                 (1., -1., 0),
                                 (1., 1., 0),
                                 (-1., 1., 0),
                                 (0., -1., 0),
                                 (1., 0., 0),
                                 (0., 1., 0),
                                 (-1., 0., 0)))
poly_power['CPE8R'] = np.array(((0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0), (2, 0, 0), (0, 2, 0), (2, 1, 0), (1, 2, 0)))
int_weights['CPE8R'] = np.array((1., 1., 1., 1.))

# CPE3
int_points['CPE3'] = np.array(((Ed3, Ed3, 0),))
bild_points['CPE3'] = np.array(((0., 0., 0.),
                                (1., 0., 0.),
                                (0., 1., 0.)))
int_weights['CPE3'] = (np.zeros(1) + 0.5).reshape(1)
poly_power['CPE3'] = np.array(((0, 0, 0), (1, 0, 0), (0, 1, 0)))

# CPE6
int_points['CPE6'] = np.array(((tri31, tri31, 0),
                               (tri32, tri31, 0),
                               (tri31, tri32, 0)))
bild_points['CPE6'] = np.array(((0., 0., 0.),
                                (1., 0., 0.),
                                (0., 1., 0.),
                                (0.5, 0., 0.),
                                (0.5, 0.5, 0.),
                                (0., 0.5, 0.)))
poly_power['CPE6'] = np.array(((0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0), (0, 2, 0), (2, 0, 0)))
int_weights['CPE6'] = np.array((1 / 6, 1 / 6, 1 / 6))

###################################################################################
# 3D - elements ###################################################################
###################################################################################

# C3D8
int_points['C3D8'] = np.array(((-EdW3, -EdW3, -EdW3),
                               (EdW3, -EdW3, -EdW3),
                               (-EdW3, EdW3, -EdW3),
                               (EdW3, EdW3, -EdW3),
                               (-EdW3, -EdW3, EdW3),
                               (EdW3, -EdW3, EdW3),
                               (-EdW3, EdW3, EdW3),
                               (EdW3, EdW3, EdW3)))
bild_points['C3D8'] = np.array(((-1., -1., -1.),
                                (1., -1., -1.),
                                (1., 1., -1.),
                                (-1., 1., -1.),
                                (-1., -1., 1.),
                                (1., -1., 1.),
                                (1., 1., 1.),
                                (-1., 1., 1.),))
poly_power['C3D8'] = np.array(((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1)))
int_weights['C3D8'] = np.array((1., 1., 1., 1., 1., 1., 1., 1.))

# C3D8R
int_points['C3D8R'] = np.array(((0., 0., 0.),))
bild_points['C3D8R'] = np.array(((-1., -1., -1.),
                                 (1., -1., -1.),
                                 (1., 1., -1.),
                                 (-1., 1., -1.),
                                 (-1., -1., 1.),
                                 (1., -1., 1.),
                                 (1., 1., 1.),
                                 (-1., 1., 1.)))
poly_power['C3D8R'] = np.array(((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1)))
int_weights['C3D8R'] = np.zeros(1) + 8

# C3D20
int_points['C3D20'] = np.array(((-W3d5, -W3d5, -W3d5),
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
                                (W3d5, W3d5, W3d5)))
bild_points['C3D20'] = np.array(((-1., -1., -1.),
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
                                 (-1., 1., 0.)))
poly_power['C3D20'] = np.array(((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1),
                                (2, 0, 0), (0, 2, 0), (0, 0, 2), (2, 1, 0), (2, 0, 1), (1, 2, 0), (0, 2, 1), (1, 0, 2),
                                (0, 1, 2), (2, 1, 1), (1, 2, 1), (1, 1, 2)))
int_weights['C3D20'] = np.array((B5d9 * B5d9 * B5d9,
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
                                 B5d9 * B5d9 * B5d9))

# C3D20R
int_points['C3D20R'] = np.array(((-EdW3, -EdW3, -EdW3),
                                 (EdW3, -EdW3, -EdW3),
                                 (-EdW3, EdW3, -EdW3),
                                 (EdW3, EdW3, -EdW3),
                                 (-EdW3, -EdW3, EdW3),
                                 (EdW3, -EdW3, EdW3),
                                 (-EdW3, EdW3, EdW3),
                                 (EdW3, EdW3, EdW3)))
bild_points['C3D20R'] = np.array(((-1., -1., -1.),
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
                                  (-1., 1., 0.)))
poly_power['C3D20R'] = np.array(((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1),
                                 (2, 0, 0), (0, 2, 0), (0, 0, 2), (2, 1, 0), (2, 0, 1), (1, 2, 0), (0, 2, 1), (1, 0, 2),
                                 (0, 1, 2), (2, 1, 1), (1, 2, 1), (1, 1, 2)))
int_weights['C3D20R'] = np.array((1., 1., 1., 1., 1., 1., 1., 1.))

# C3D4
int_points['C3D4'] = np.array(((0.25, 0.25, 0.25),))
bild_points['C3D4'] = np.array(((0., 0., 0.),
                                (1., 0., 0.),
                                (0., 1., 0.),
                                (0., 0., 1.)))
poly_power['C3D4'] = np.array(((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)))
int_weights['C3D4'] = np.zeros(1) + 1. / 6.

# C3D10
int_points['C3D10'] = np.array(((tri42, tri42, tri42),
                                (tri41, tri42, tri42),
                                (tri42, tri41, tri42),
                                (tri42, tri42, tri41)))
bild_points['C3D10'] = np.array(((0., 0., 0.),
                                 (1., 0., 0.),
                                 (0., 1., 0.),
                                 (0., 0., 1.),
                                 (0.5, 0., 0.),
                                 (0.5, 0.5, 0.),
                                 (0., 0.5, 0.),
                                 (0., 0., 0.5),
                                 (0.5, 0., 0.5),
                                 (0., 0.5, 0.5)))
poly_power['C3D10'] = np.array(((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0),
                                (1, 0, 1), (0, 1, 1), (2, 0, 0), (0, 2, 0), (0, 0, 2)))
int_weights['C3D10'] = np.array((1. / 24, 1. / 24, 1. / 24, 1. / 24))

# C3D10R
int_points['C3D10R'] = np.array(((0.25, 0.25, 0.25),))
bild_points['C3D10R'] = np.array(((0., 0., 0.),
                                  (1., 0., 0.),
                                  (0., 1., 0.),
                                  (0., 0., 1.),
                                  (0.5, 0., 0.),
                                  (0.5, 0.5, 0.),
                                  (0., 0.5, 0.),
                                  (0., 0., 0.5),
                                  (0.5, 0., 0.5),
                                  (0., 0.5, 0.5)))
poly_power['C3D10R'] = np.array(((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0),
                                 (1, 0, 1), (0, 1, 1), (2, 0, 0), (0, 2, 0), (0, 0, 2)))
int_weights['C3D10R'] = np.zeros(1) + 1. / 6.
