"""
This module provides methods to read and create field outputs.
"""
import logging

import odbAccess
import abaqusConstants as abqConst

import numpy as np

from conforce import cf_c
from conforce.tensor_util import (
    tensor_from_abaqus_notation,
    abaqus_notation_from_tensor,
    rotation_matrix_from_quaternion
)
from conforce.element_type_mapping import map_abaqus_element_type_to_supported_element_type
import conforce_abq


LOGGER = conforce_abq.LOGGER.getChild(__name__)
"""default logger in this module"""


class FieldOutputReader(object):
    def __init__(self):
        """
        Reads field outputs and prepares
        data for the latter computation
        """
        # odb inst
        self._odb_inst = None
        self._element_labels_to_node_labels_for_type = None
        self._node_labels_to_coordinates = None
        self._X_el_labels_for_type = None
        self._X_at_nodes_for_type = None

        # element type
        self._element_type = None

        # fo_U
        self._fo_U = None
        self._U_el_labels_for_type = None
        self._U_at_nodes_for_type = None

        # fo_S
        self._fo_S = None
        self._S_el_labels = None
        self._S_at_int_points = None

        # fo_e
        self._fo_e = None
        self._e_el_labels = None
        self._e_at_int_points = None

        # masks
        self._el_labels = None
        self._X_mask = None
        self._U_mask = None
        self._e_mask = None
        self._S_mask = None

    def set_odb_inst(self, odb_inst):
        """
        Set the odb instance

        :param odb_inst: OdbInstance
        """
        self._odb_inst = odb_inst
        self._element_labels_to_node_labels_for_type = None
        self._node_labels_to_coordinates = None
        self._X_el_labels_for_type = None
        self._X_at_nodes_for_type = None

        self._U_el_labels_for_type = None
        self._U_at_nodes_for_type = None
        self._S_el_labels = None
        self._S_at_int_points = None
        self._e_el_labels = None
        self._e_at_int_points = None
        self._reset_masks()

    def set_element_type(self, element_type):
        """
        Define the element type

        :param element_type: str
        """
        self._element_type = element_type

        self._S_el_labels = None
        self._S_at_int_points = None
        self._e_el_labels = None
        self._e_at_int_points = None
        self._reset_masks()

    def set_fo_U(self, fo_U):
        """
        Set the displacement field output

        :param fo_U: FieldOutput
        """
        self._fo_U = fo_U
        self._U_el_labels_for_type = None
        self._U_at_nodes_for_type = None
        self._reset_masks()

    def set_fo_S(self, fo_S):
        """
        Set the stress field output

        :param fo_S: FieldOutput
        """
        self._fo_S = fo_S
        self._S_el_labels = None
        self._S_at_int_points = None
        self._reset_masks()

    def set_fo_e(self, fo_e):
        """
        Set the energy density field output

        :param fo_e: FieldOutput
        """
        self._fo_e = fo_e
        self._e_el_labels = None
        self._e_at_int_points = None
        self._reset_masks()

    def _reset_masks(self):
        self._el_labels = None
        self._X_mask = None
        self._U_mask = None
        self._e_mask = None
        self._S_mask = None

    @property
    def odb_inst(self):
        """
        Get the odb instance

        :return: OdbInstance
        """
        return self._odb_inst

    @property
    def element_type(self):
        """
        Get the defined element type

        :return: str
        """
        return self._element_type

    @property
    def fo_U(self):
        """
        Get the displacement field output at nodes of the defined instance

        :return: FieldOutput
        """
        return self._fo_U.getSubset(
            position=abqConst.NODAL,
            region=self._odb_inst
        )

    @property
    def fo_S(self):
        """
        Get the stress field output at integration points of
        the defined instance and element type

        :return: FieldOutput
        """
        return self._fo_S.getSubset(
            position=abqConst.INTEGRATION_POINT,
            region=self._odb_inst,
            elementType=self.element_type
        )

    @property
    def fo_e(self):
        """
        Get the energy density field output at integration points of
        the defined instance and element type

        :return: FieldOutput
        """
        return self._fo_e.getSubset(
            position=abqConst.INTEGRATION_POINT,
            region=self._odb_inst,
            elementType=self.element_type
        )

    @property
    def element_labels_to_node_labels_for_type(self):
        """
        Dictionary mapping the element type to a dictionary
        that maps element labels to the node labels of the element

        :return: dict
        """
        if self._element_labels_to_node_labels_for_type is None:
            odb_elements = self.odb_inst.elements

            el_to_n_label_for_type = dict()
            self._element_labels_to_node_labels_for_type = el_to_n_label_for_type
            for odb_element in odb_elements:
                el_type = odb_element.type
                label = odb_element.label
                connectivity = odb_element.connectivity

                if el_type not in el_to_n_label_for_type:
                    el_to_n_label_for_type[el_type] = dict()

                el_to_n_label_for_type[el_type][label] = (
                    connectivity)

        return self._element_labels_to_node_labels_for_type

    @property
    def element_labels_to_node_labels(self):
        """
        Dictionary that maps element labels to the node labels of the element.
        Works only for the defined element type.

        :return: dict
        """
        return self.element_labels_to_node_labels_for_type[self.element_type]

    @property
    def node_label_to_coordinates_mapping(self):
        """
        Dictionary mapping node labels to coordinates of the nodes

        :return: dict
        """
        if self._node_labels_to_coordinates is None:

            odb_nodes = self.odb_inst.nodes
            self._node_labels_to_coordinates = {
                odb_node.label: odb_node.coordinates
                for odb_node in odb_nodes
            }

        return self._node_labels_to_coordinates

    def _update_X_values(self):
        self._X_el_labels_for_type, self._X_at_nodes_for_type = self.map_nodes_to_element_nodal(
            self.node_label_to_coordinates_mapping
        )

    @property
    def X_el_labels_for_type(self):
        """
        Dictionary mapping the element type to an array of shape
        (num_elem,) containing the element labels.
        Element labels are in the same order as in :py:attr:`X_at_nodes_for_type`.

        :return: dict
        """
        if self._X_el_labels_for_type is None:
            self._update_X_values()

        return self._X_el_labels_for_type

    @property
    def X_el_labels(self):
        """
        Array of shape (num_elem,) containing the element labels
        of the defined element type
        Element labels are in the same order as in :py:attr:`X_at_nodes`.

        :return: np.ndarray
        """
        return self.X_el_labels_for_type[self.element_type]

    @property
    def X_at_nodes_for_type(self):
        """
        Dictionary mapping the element type to an array of shape
        (num_elem, n, 3/2) containing the coordinates.
        The number of nodes per element is n.

        :return: dict
        """
        if self._X_at_nodes_for_type is None:
            self._update_X_values()

        return self._X_at_nodes_for_type

    @property
    def X_at_nodes(self):
        """
        Array of shape (num_elem, n, 3/2) containing the coordinates
        of the defined element type.
        The number of nodes per element is n.

        :return: np.ndarray
        """
        return self.X_at_nodes_for_type[self.element_type]

    def _update_U_values(self):
        self._U_el_labels_for_type, self._U_at_nodes_for_type = self.map_nodes_to_element_nodal(
            self.create_node_label_to_bulk_data_mapping(self._fo_U.bulkDataBlocks)
        )

    @property
    def U_el_labels_for_type(self):
        """
        Array of shape (num_elem,) containing the element labels
        of the defined element type.
        Element labels are in the same order as in :py:attr:`U_at_nodes_for_type`.

        :return: np.ndarray
        """
        if self._U_el_labels_for_type is None:
            self._update_U_values()

        return self._U_el_labels_for_type

    @property
    def U_el_labels(self):
        """
        Array of shape (num_elem,) containing the element labels
        of the defined element type.
        Element labels are in the same order as in :py:attr:`U_at_nodes`.

        :return: np.ndarray
        """
        return self.U_el_labels_for_type[self.element_type]

    @property
    def U_at_nodes_for_type(self):
        """
        Dictionary mapping the element type to an array of shape
        (num_elem, n, 3/2) containing the displacements.
        The number of nodes per element is n.

        :return: dict
        """
        if self._U_at_nodes_for_type is None:
            self._update_U_values()

        return self._U_at_nodes_for_type

    @property
    def U_at_nodes(self):
        """
        Array of shape (num_elem, n, 3/2) containing the displacements
        of the defined element type.
        The number of nodes per element is n.

        :return: np.ndarray
        """
        return self.U_at_nodes_for_type[self.element_type]

    def _update_e_values(self):
        self._e_el_labels, self._e_at_int_points = self.extract_integration_points_values(
            self.fo_e.bulkDataBlocks
        )
        self._e_at_int_points = self._e_at_int_points.reshape((len(self._e_el_labels), -1))

    @property
    def e_el_labels(self):
        """
        Array of shape (num_elem,) containing the element labels
        of the defined element type.
        Element labels are in the same order as in :py:attr:`e_at_int_points`.

        :return: np.ndarray
        """
        if self._e_el_labels is None:
            self._update_e_values()

        return self._e_el_labels

    @property
    def e_at_int_points(self):
        """
        Array of shape (num_elem, ips) containing the energy densities
        at ips integration point for the defined element type.

        :return: np.ndarray
        """
        if self._e_at_int_points is None:
            self._update_e_values()

        return self._e_at_int_points

    def _update_S_values(self):
        self._S_el_labels, self._S_at_int_points = self.extract_integration_points_values(
            self.fo_S.bulkDataBlocks
        )
        self._S_at_int_points = tensor_from_abaqus_notation(self._S_at_int_points)

    @property
    def S_el_labels(self):
        """
        Array of shape (num_elem,) containing the element labels
        of the defined element type.
        Element labels are in the same order as in :py:attr:`S_at_int_points`.

        :return: np.ndarray
        """
        if self._S_el_labels is None:
            self._update_S_values()

        return self._S_el_labels

    @property
    def S_at_int_points(self):
        """
        Array of shape (num_elem, ips, 2/3, 2/3) containing the
        stress tensors at ips integration point for the defined element type.

        :return: np.ndarray
        """
        if self._S_at_int_points is None:
            self._update_S_values()

        return self._S_at_int_points

    def _update_masks(self):
        # filter elements not defined in all outputs
        self._el_labels = np.intersect1d(
            np.intersect1d(self.X_el_labels, self.U_el_labels),
            np.intersect1d(self.e_el_labels, self.S_el_labels)
        )
        self._X_mask = np.in1d(self.X_el_labels, self._el_labels)
        self._U_mask = np.in1d(self.U_el_labels, self._el_labels)
        self._e_mask = np.in1d(self.e_el_labels, self._el_labels)
        self._S_mask = np.in1d(self.S_el_labels, self._el_labels)

    @property
    def el_labels(self):
        """
        Labels of elements for which all data are present

        :return: np.ndarray
        """
        if self._el_labels is None:
            self._update_masks()

        return self._el_labels

    @property
    def X_mask(self):
        """
        Mask for :py:attr:`X_el_labels` and :py:attr:`X_at_nodes`
        such that the masked data correspond to :py:attr:`el_labels`

        :return: np.ndarray
        """
        if self._X_mask is None:
            self._update_masks()

        return self._X_mask

    @property
    def U_mask(self):
        """
        Mask for :py:attr:`U_el_labels` and :py:attr:`U_at_nodes`
        such that the masked data correspond to :py:attr:`el_labels`

        :return: np.ndarray
        """
        if self._U_mask is None:
            self._update_masks()

        return self._U_mask

    @property
    def e_mask(self):
        """
        Mask for :py:attr:`e_el_labels` and :py:attr:`e_at_int_points`
        such that the masked data correspond to :py:attr:`el_labels`

        :return: np.ndarray
        """
        if self._e_mask is None:
            self._update_masks()

        return self._e_mask

    @property
    def S_mask(self):
        """
        Mask for :py:attr:`S_el_labels` and :py:attr:`S_at_int_points`
        such that the masked data correspond to :py:attr:`el_labels`

        :return: np.ndarray
        """
        if self._S_mask is None:
            self._update_masks()

        return self._S_mask

    @staticmethod
    def create_node_label_to_bulk_data_mapping(bulk_data_blocks):
        """
        Create a dictionary mapping node labels to data

        :param bulk_data_blocks: sequence of bulk data objects
        :return: dict
        """
        data = np.concatenate([block.data for block in bulk_data_blocks])
        node_labels = np.concatenate([block.nodeLabels for block in bulk_data_blocks])

        return {
            label: data
            for label, data in zip(node_labels, data)
        }

    @staticmethod
    def extract_integration_points_values(bulk_data_blocks):
        """
        Read integration point values and put them into:

        - `el_labels_unique`: array of shape (num_elements,)
        - `reshaped_data`: array of shape (num_elements, ips, ?) in the same order as `el_labels_unique`

        :param bulk_data_blocks: sequence of bulk data objects
        :return: `el_labels_unique`, `reshaped_data`
        """
        el_labels = np.concatenate([block.elementLabels for block in bulk_data_blocks])
        int_points = np.concatenate([block.integrationPoints for block in bulk_data_blocks])
        data = np.concatenate([block.data for block in bulk_data_blocks])

        assert 1 <= len({block.baseElementType for block in bulk_data_blocks})

        # sort: ascending element labels and integration points
        idx = np.lexsort((int_points, el_labels))
        el_labels = el_labels[idx]
        int_points = int_points[idx]
        data = data[idx]

        # assign an index to integration points
        int_points_unique, int_points_indices = np.unique(int_points, return_inverse=True)
        num_ip = len(int_points_unique)

        # assign an index to element labels
        el_labels_unique, el_indices = np.unique(el_labels, return_inverse=True)
        num_el = len(el_labels_unique)

        # reshape data
        d = data.shape[1]
        reshaped_data = np.zeros(
            shape=(num_el, num_ip, d),
            dtype=float
        )
        reshaped_data[el_indices, int_points_indices] = data[:]
        return el_labels_unique, reshaped_data

    def map_nodes_to_element_nodal(self, node_label_to_data_mapping):
        """
        Map data given at nodes to element nodes.
        This implies that data is duplicated if more elements share one node.

        - `element_labels_for_type`: Dictionary mapping the element type to element labels
        - `element_nodal_data_for_type`: Dictionary mapping the element type to an array of shape
          (num_elements, n, ?) in the same order as in `element_labels_for_type`

        :param node_label_to_data_mapping: Dictionary mapping a node label to data
        :return: `element_labels_for_type`, `element_nodal_data_for_type`
        """
        el_to_n_label_for_type = self.element_labels_to_node_labels_for_type
        element_labels_for_type = dict()
        element_nodal_data_for_type = dict()

        for el_type, el_to_n_labels in el_to_n_label_for_type.items():
            element_nodal_data = list()
            element_labels = list()

            for element_label, node_labels in el_to_n_labels.items():
                data = list()
                for label in node_labels:
                    if label in node_label_to_data_mapping:
                        data.append(node_label_to_data_mapping[label])
                    else:
                        break
                else:
                    element_nodal_data.append(data)
                    element_labels.append(element_label)

            element_nodal_data = np.array(element_nodal_data, dtype=float)
            element_labels = np.array(element_labels, dtype=int)

            # sort
            idx = np.argsort(element_labels)
            element_nodal_data_for_type[el_type] = element_nodal_data[idx]
            element_labels_for_type[el_type] = element_labels[idx]

        return element_labels_for_type, element_nodal_data_for_type


class _FieldOutputWriter(object):
    def add(self, reader, supported_element_type):
        """
        Add data for an element type to the FieldOutput

        :param reader: FieldOutputReader
        :param supported_element_type: str, a supported element type
            as described in :py:attr:`conforce.cf_c.map_type_to_info`
        """

    def flush(self, odb_inst):
        """
        Call this after all data are added for all element types.
        This adds data to the FieldOutput that considers multiple element types.

        :param odb_inst: OdbInstance
        """


class FFieldOutputWriter(_FieldOutputWriter):
    def __init__(self, frame, d, name):
        """
        Create new FieldOutputs for each
        component of the deformation gradient.

        :param frame: OdbFrame to which FieldOutputs are added
        :param d: int, number of dimensions
        :param name: str, name of new FieldOutputs
        """
        _FieldOutputWriter.__init__(self)
        self.name = str(name)
        self._fo = [
            [
                frame.FieldOutput(
                    name=self.name + "_" + str(i + 1) + str(j + 1),
                    description="deformation gradient",
                    type=abqConst.SCALAR
                )
                for j in range(d)
            ]
            for i in range(d)
        ]
        self._d = d

    def add(self, reader, supported_element_type):
        f_data = cf_c.compute_F(
            X_at_nodes=reader.X_at_nodes[reader.X_mask, :, :self._d],
            U_at_nodes=reader.U_at_nodes[reader.U_mask, :, :self._d],
            element_type=supported_element_type,
        )

        for i in range(self._d):
            for j in range(self._d):
                self._fo[i][j].addData(
                    position=abqConst.INTEGRATION_POINT,
                    instance=reader.odb_inst,
                    labels=reader.el_labels,
                    data=np.ascontiguousarray(
                        f_data[:, :, i, j].reshape((-1, 1))
                    )
                )

    def flush(self, odb_inst):
        pass


class PFieldOutputWriter(_FieldOutputWriter):
    def __init__(self, frame, d, name):
        """
        Create new FieldOutputs for each
        component of the first Piola-Kirchhoff stress.

        :param frame: OdbFrame to which FieldOutputs are added
        :param d: int, number of dimensions
        :param name: str, name of new FieldOutputs
        """
        _FieldOutputWriter.__init__(self)
        self.name = str(name)
        self._fo = [
            [
                frame.FieldOutput(
                    name=self.name + "_" + str(i + 1) + str(j + 1),
                    description="First Piola-Kirchhoff stress tensor",
                    type=abqConst.SCALAR
                )
                for j in range(d)
            ]
            for i in range(d)
        ]
        self._d = d

    def add(self, reader, supported_element_type):
        p_data = cf_c.compute_P(
            X_at_nodes=reader.X_at_nodes[reader.X_mask, :, :self._d],
            U_at_nodes=reader.U_at_nodes[reader.U_mask, :, :self._d],
            S_at_int_points=reader.S_at_int_points[reader.S_mask, :, :self._d, :self._d],
            element_type=supported_element_type,
        )

        for i in range(self._d):
            for j in range(self._d):
                self._fo[i][j].addData(
                    position=abqConst.INTEGRATION_POINT,
                    instance=reader.odb_inst,
                    labels=reader.el_labels,
                    data=np.ascontiguousarray(
                        p_data[:, :, i, j].reshape((-1, 1))
                    )
                )


class CSFieldOutputWriter(_FieldOutputWriter):
    def __init__(self, frame, d, name, method):
        """
        Create new FieldOutputs for each
        component of the first Piola-Kirchhoff stress.

        :param frame: OdbFrame to which FieldOutputs are added
        :param d: int, number of dimensions
        :param name: str, name of new FieldOutputs
        :param method: Method as described in :py:func:`conforce.cf_c.compute_CS`
        """
        _FieldOutputWriter.__init__(self)
        self.name = str(name)
        self._fo = [
            [
                frame.FieldOutput(
                    name=self.name + "_" + str(i + 1) + str(j + 1),
                    description="configurational stresses",
                    type=abqConst.SCALAR
                )
                for j in range(d)
            ]
            for i in range(d)
        ]

        self._d = d
        self._method = method

    def add(self, reader, supported_element_type):
        # compute configurational stresses
        cs_data = cf_c.compute_CS(
            e_at_int_points=reader.e_at_int_points[reader.e_mask],
            X_at_nodes=reader.X_at_nodes[reader.X_mask, :, :self._d],
            U_at_nodes=reader.U_at_nodes[reader.U_mask, :, :self._d],
            S_at_int_points=reader.S_at_int_points[reader.S_mask, :, :self._d, :self._d],
            element_type=supported_element_type,
            method=self._method
        )

        for i in range(self._d):
            for j in range(self._d):
                self._fo[i][j].addData(
                    position=abqConst.INTEGRATION_POINT,
                    instance=reader.odb_inst,
                    labels=reader.el_labels,
                    data=np.ascontiguousarray(
                        cs_data[:, :, i, j].reshape((-1, 1))
                    )
                )


class CFFieldOutputWriter(_FieldOutputWriter):
    def __init__(self, frame, d, name, method):
        """
        Create one FieldOutput for the configurational forces.

        :param frame: OdbFrame to which the FieldOutput is added
        :param d: int, number of dimensions
        :param name: str, name of the new FieldOutput
        :param method: Method as described in :py:func:`conforce.cf_c.compute_CS`
        """
        _FieldOutputWriter.__init__(self)
        self.name = str(name)
        self._fo = frame.FieldOutput(
            name=self.name,
            description="configurational forces",
            type=abqConst.VECTOR,
            validInvariants=[
                abqConst.MAGNITUDE
            ]
        )
        self._d = d
        self._method = method

        self._CF_at_nodes = dict()

    def add(self, reader, supported_element_type):
        # compute configurational forces
        cf_data = cf_c.compute_CF(
            e_at_int_points=reader.e_at_int_points[reader.e_mask],
            X_at_nodes=reader.X_at_nodes[reader.X_mask, :, :self._d],
            U_at_nodes=reader.U_at_nodes[reader.U_mask, :, :self._d],
            S_at_int_points=reader.S_at_int_points[reader.S_mask, :, :self._d, :self._d],
            element_type=supported_element_type,
            method=self._method
        )

        # create datastructure abaqus understands
        el_to_n_label = reader.element_labels_to_node_labels_for_type[reader.element_type]
        CF_at_nodes = self._CF_at_nodes
        for CF_el_data, el_label in zip(cf_data, reader.el_labels):
            for CF_el_node_data, node_label in zip(CF_el_data, el_to_n_label[el_label]):
                if node_label not in CF_at_nodes:
                    CF_at_node = CF_at_nodes[node_label] = np.zeros(self._d)

                else:
                    CF_at_node = CF_at_nodes[node_label]

                CF_at_node += CF_el_node_data

    def flush(self, odb_inst):
        # add data to field output
        self._fo.addData(
            position=abqConst.NODAL,
            instance=odb_inst,
            labels=self._CF_at_nodes.keys(),
            data=np.array(self._CF_at_nodes.values(), dtype=float)
        )

        # reset CF_at_nodes
        self._CF_at_nodes = dict()


def get_present_element_types_in(bulk_data_blocks):
    """
    Return all element types that occur in at least one bulk data block.
    A bulk data block contains data and is associated to one element type.

    :param bulk_data_blocks: A sequence of bulk data blocks
    :return: a set of strings defining the element types
    """
    return {
        block.baseElementType
        for block in bulk_data_blocks
    }


def eval_field_output_expression(field_outputs, expression):
    """
    Create a new temporary anonymous field output that
    is computed out of the existing field outputs as described by the expression.

    The expression may contain:

    - values:

        - floating point or integer values
        - names of existing field outputs

    - operators: +,-,*,/

    .. note:

        A requirement for the expression is, that all used field outputs have:

            - the same type (Scalar, Vector, Tensor)
            - are evaluated on the same position (Integration point, Nodes, Elements)

    :param field_outputs: a dictionary like object mapping field output names to field outputs
    :param expression: str, an expression how to compute the new field output
    :return: new temporary anonymous field output

    :raises KeyError: If a value in the expressions is neither a scalar value
        nor a name of an existing field output
    """
    fo = None
    for addition in expression.split("+"):
        fo_subtract = None
        for subtraction in addition.split("-"):
            fo_multiplication = None
            for multiplication in subtraction.split("*"):
                fo_division = None
                for division in multiplication.split("/"):
                    division = division.strip()
                    try:
                        new_term = float(division)
                    except ValueError:
                        if division not in field_outputs.keys():
                            raise KeyError(division, expression)

                        new_term = field_outputs[division]

                    if fo_division is None:
                        fo_division = new_term
                    else:
                        fo_division = fo_division / new_term

                if fo_multiplication is None:
                    fo_multiplication = fo_division
                else:
                    fo_multiplication = fo_multiplication * fo_division

            if fo_subtract is None:
                fo_subtract = fo_multiplication
            else:
                fo_subtract = fo_subtract - fo_multiplication

        if fo is None:
            fo = fo_subtract
        else:
            fo = fo + fo_subtract

    return fo


def rotate_field_output_to_global_coordinate_system(frame, field_output, name, description, logger=None):
    """
    Create a new field output with the given name and description
    and saves the new field output into `frame.fieldOutputs`.
    The new field output contains data of the given `field_output`
    that are rotated to the global coordinate system.

    Supported data types are:

    - Scalar values,
    - Vectors with 2 and 3 components
    - Tensors supported by :py:func:`conforce.tensor_util.tensor_from_abaqus_notation`

    :param frame: The new field output is saved to this frame's fieldOutputs
    :param field_output: The existing field output whose values should be rotated
    :param name: The name of the new field output
    :param description: The description of the new field output
    :param logger: Logging instance, default uses :py:attr:`LOGGER`
    :return: A new field output with data rotated to the global coordinate system
    """
    if logger is None:
        logger = LOGGER

    if field_output.type == abqConst.SCALAR:
        return field_output

    new_field_output = frame.FieldOutput(
        name=name,
        description=description,
        type=field_output.type,
        validInvariants=field_output.validInvariants
    )

    for block in field_output.bulkDataBlocks:
        if block.position == abqConst.NODAL:
            labels = block.nodeLabels
        elif block.position == abqConst.INTEGRATION_POINT:
            labels = block.elementLabels
            labels, indices = np.unique(labels, return_index=True)
            labels = labels[np.argsort(indices)]
        else:
            raise NotImplementedError("not supported position " + str(block.position))

        if block.localCoordSystem is None:
            data = block.data

        elif field_output.type == abqConst.VECTOR:
            ROT = rotation_matrix_from_quaternion(block.localCoordSystem)

            d = block.data.shape[1]
            assert d in (3, 2)
            if d == 3:
                local_vectors = block.data

            else:
                # create 3D vector out of 2D vector
                local_vectors = np.zeros((block.data.shape[0], 3), dtype=float)
                local_vectors[:, :2] = block.data

            data = np.einsum("...ji,...j", ROT, local_vectors)

        else:  # TENSOR
            d_quaternion = block.localCoordSystem.shape[1]
            assert d_quaternion in (4, 2)
            if d_quaternion == 4:
                Q = block.localCoordSystem
            else:
                # Create full quaternions from the given `Q[..., 2]*k + Q[..., 3]`
                Q = np.zeros((block.localCoordSystem.shape[0], 4), dtype=float)
                Q[:, 2:] = block.localCoordSystem

            ROT = rotation_matrix_from_quaternion(Q)

            local_vectors = block.data
            local_tensors = tensor_from_abaqus_notation(local_vectors)
            global_tensors = np.einsum("...ji,...jk,...kl", ROT, local_tensors, ROT)
            data = abaqus_notation_from_tensor(global_tensors, local_vectors.shape[-1])

        # add data
        if block.instance is not None:
            new_field_output.addData(
                position=block.position,
                instance=block.instance,
                labels=np.ascontiguousarray(labels),
                data=np.ascontiguousarray(data)
            )

        else:
            logger.warning("skip in field output %s labels %s (not associated with an instance)", name, labels[:10])

    return new_field_output


def add_field_outputs(
        odb, 
        method="mbf",
        e_expression="SENER+PENER",
        name_U_global_csys="U_GLOBAL_CSYS",
        name_S_global_csys="S_GLOBAL_CSYS",
        request_F=True,
        name_F="DEF_GRAD",
        request_P=True,
        name_P="FIRST_PIOLA_STRESS",
        request_CS=True,
        name_CS="CONF_STRESS",
        request_CF=True,
        name_CF="CONF_FORCE",
        odb_instances=None,
        odb_set=None,
        el_type_mapping=None,
        logger=None
):
    """
    Add field outputs to each OdbInstance, OdbStep and OdbFrame in the given `odb`.
    The newly created field outputs are named according to the parameters
    `name_F`, `name_P`, `name_CS`, and `name_CF`.
    The field outputs are only computed and saved to the `odb`,
    if they are requested by setting the corresponding parameters
    `request_F`, `request_P`, `request_CS`, and `request_CF` to True.

    For the configurational stresses and forces, an energy density is required.
    Use the parameter `e_expression` to define how the energy density is computed.
    `e_expression` can refer to any scalar field output.
    Assure, that the field output exists in the odb by requesting e.g. `ENER` for all steps.
    Furthermore, the configurational stresses and forces can be computed with
    various `methods`.

    Two field outputs named as defined in `name_U_global_csys` and `name_S_global_csys` are created.
    These two field outputs contain the displacements and stresses in the global coordinate system,
    regardless whether a transformation has been applied or not.

    Not supported element types can be replaced by supported element types.
    See :py:attr:`conforce.cf_c.map_type_to_info` for a list of supported element types.
    The parameter `el_type_mapping` maps element types to supported element types.
    The default mapping of `el_type_mapping` is
    :py:attr:`conforce.element_type_mapping.map_abaqus_element_type_to_supported_element_type`.

    .. note::

        The odb is closed and reopened.
        Consequently, the input `odb` object is not the same at the resulting `odb` object.
        After this function call the input `odb` is invalidated and must not be used anymore.

    :param odb: Odb into which the FieldOutput objects are written.
    :type odb: Odb
    :param method: see :py:func:`conforce.cf_c.compute_CS`
    :type method: str
    :param e_expression: see :py:func:`eval_field_output_expression`
    :type e_expression: str
    :param name_U_global_csys: name of the newly generated field output containing displacements
        in the global coordinate system
    :type name_U_global_csys: str
    :param name_S_global_csys: name of the newly generated field output containing stresses
        in the global coordinate system
    :type name_S_global_csys: str
    :param request_F: True to create field outputs for the deformation gradients.
    :type request_F: bool
    :param name_F: name of the field outputs that contain the deformation gradient.
        These field outputs are only generated if "F" is in `fields`.
    :type name_F: str
    :param request_P: True to create field outputs for the first Piola-Kirchhoff stresses.
    :type request_P: bool
    :param name_P: name of the field outputs that contain the first Piola-Kirchhoff stresses.
        These field outputs are only generated if "P" is in `fields`.
    :type name_P: str
    :param request_CS: True to create field outputs for the configurational stresses.
    :type request_CS: bool
    :param name_CS: name of the field outputs that contain the configurational stresses.
        These field outputs are only generated if "CS" is in `fields`.
    :type name_CS: str
    :param request_CF: True to create field outputs for the configurational forces.
    :type request_CF: bool
    :param name_CF: name of the field output that contains the configurational forces.
        This field outputs is only generated if "CF" is in `fields`.
    :type name_CF: str
    :param odb_instances: Compute quantities only for these instances
    :type odb_instances: Sequence of OdbInstance
    :param odb_set: Compute quantities only for nodes and element defined in this set.
    :type odb_set: OdbSet
    :param el_type_mapping: Maps element types to supported element types.
    :type el_type_mapping: Dict[str, str]
    :param logger: Print logging messages, default uses :py:attr:`LOGGER`
    :type logger: logging.Logger
    :return: modified Odb or None if Odb was read-only
    :rtype: Odb
    """

    if el_type_mapping is None:
        el_type_mapping = map_abaqus_element_type_to_supported_element_type

    if logger is None:
        logger = LOGGER

    path = odb.path
    is_read_only = odb.isReadOnly
    if is_read_only:
        logger.error("Odb is read-only. Reopen the odb and assure the Read-only checkbox is not checked.")
        return None

    fo_reader = FieldOutputReader()

    if odb_instances is None:
        odb_instances = odb.rootAssembly.instances.values()

    for odb_inst in odb_instances:
        fo_reader.set_odb_inst(odb_inst)

        if len(fo_reader.element_labels_to_node_labels_for_type) == 0:
            logger.info("skip instance %s (contains no elements)", odb_inst.name)
            continue

        for step in odb.steps.values():
            for frame in step.frames:
                msg = (
                        "instance=" + str(odb_inst.name)
                        + "; step=" + str(step.name)
                        + "; frame=" + str(frame.frameId)
                )

                # get field outputs and rotate them to the global coordinate system
                fo = frame.fieldOutputs

                # energy density
                try:
                    fo_e = eval_field_output_expression(fo, e_expression)
                except KeyError as error:
                    logger.error("invalid field output %s in expression %s (%s)", error.args[0], error.args[1], msg)
                    break

                # displacements in global coordinate system
                if name_U_global_csys in fo.keys():
                    logger.info("found field output %s (%s)", name_U_global_csys, msg)
                    fo_U = fo[name_U_global_csys]

                else:
                    logger.info("create field output %s (%s)", name_U_global_csys, msg)
                    fo_U = fo["U"]
                    fo_U = rotate_field_output_to_global_coordinate_system(
                        frame,
                        fo_U,
                        name_U_global_csys,
                        "Displacement in global coordinate system",
                        logger=logger
                    )

                # stresses in global coordinate system
                if name_S_global_csys in fo.keys():
                    logger.info("found field output %s (%s)", name_S_global_csys, msg)
                    fo_S = fo[name_S_global_csys]

                else:
                    logger.info("create field output %s (%s)", name_S_global_csys, msg)

                    fo_S = fo["S"]
                    fo_S = rotate_field_output_to_global_coordinate_system(
                        frame,
                        fo_S,
                        name_S_global_csys,
                        "Stresses in global coordinate system",
                        logger=logger
                    )

                # compute values only for the given odb set
                if odb_set is not None:
                    if odb_set.nodes is not None:
                        fo_U = fo_U.getSubset(region=odb_set)

                    elif odb_set.elements is not None:
                        fo_e = fo_e.getSubset(region=odb_set)
                        fo_S = fo_S.getSubset(region=odb_set)

                    else:
                        logger.warning(
                            "odb_set %s does not contain nodes or elements and is ignored.",
                            odb_set.name
                        )

                # read odb output
                fo_reader.set_fo_U(fo_U)
                fo_reader.set_fo_e(fo_e)
                fo_reader.set_fo_S(fo_S)

                # write computed fields to odb
                d = fo_U.bulkDataBlocks[0].data.shape[1]
                fo_writers = list()
                fo_keys = set(fo.keys())

                if request_F and (name_F + "_11") not in fo_keys:
                    logger.info("create field output %s_ij (%s)", name_F, msg)
                    fo_writers.append(FFieldOutputWriter(frame, d, name_F))
                elif request_F:
                    logger.warning("skip field output %s_ij (%s)", name_F, msg)

                if request_P and (name_P + "_11") not in fo_keys:
                    logger.info("create field output %s_ij (%s)", name_P, msg)
                    fo_writers.append(PFieldOutputWriter(frame, d, name_P))
                elif request_P:
                    logger.warning("skip field output %s_ij (%s)", name_P, msg)

                if request_CS and (name_CS + "_11") not in fo_keys:
                    logger.info("create field output %s_ij (%s)", name_CS, msg)
                    fo_writers.append(CSFieldOutputWriter(frame, d, name=name_CS, method=method))
                elif request_CS:
                    logger.warning("skip field output %s_ij (%s)", name_CS, msg)

                if request_CF and name_CF not in fo_keys:
                    logger.info("create field output %s (%s)", name_CF, msg)
                    fo_writers.append(CFFieldOutputWriter(frame, d, name=name_CF, method=method))
                elif request_CF:
                    logger.warning("skip field output %s (%s)", name_CF, msg)

                # add data for all element types
                for element_type in get_present_element_types_in(fo["S"].bulkDataBlocks):
                    fo_reader.set_element_type(element_type)

                    # get supported element type (CF, ... is implemented for the element type)
                    if element_type in el_type_mapping:
                        supported_el_type = el_type_mapping[element_type]

                        if supported_el_type != element_type:
                            logger.warning(
                                "compute element type %s as if it were %s. (%s is not supported)",
                                element_type, supported_el_type, element_type)

                    else:
                        logger.warning(
                            "skip not supported element type %s (can not be replaced by a similar element type)",
                            element_type)
                        continue

                    # add data
                    for fo_writer in fo_writers:
                        fo_writer.add(fo_reader, supported_el_type)

                # consolidate data of all element types
                for fo_writer in fo_writers:
                    fo_writer.flush(odb_inst)

    odb.save()
    odb.close()
    return odbAccess.openOdb(path, readOnly=False)
