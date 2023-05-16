"""
This module requires Abaqus.
"""
import odbAccess
import abaqusConstants as abqConst

import numpy as np

from cf_shared import cf_c
from cf_shared.tensor_util import (
    tensor_from_abaqus_notation,
    abaqus_notation_from_tensor,
    rotation_matrix_from_quaternion
)
from cf_shared.element_type_mapping import map_abaqus_element_type_to_supported_element_type
import cf_abq


LOGGER = cf_abq.LOGGER.getChild(__name__)


class FieldOutputReader(object):
    map_el_info_to_el_type = {
        info: type
        for type, info in cf_c.map_type_to_info.items()
    }

    def __init__(self):
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
        self._element_type = element_type

        self._S_el_labels = None
        self._S_at_int_points = None
        self._e_el_labels = None
        self._e_at_int_points = None
        self._reset_masks()

    def set_fo_U(self, fo_U):
        self._fo_U = fo_U
        self._U_el_labels_for_type = None
        self._U_at_nodes_for_type = None
        self._reset_masks()

    def set_fo_S(self, fo_S):
        self._fo_S = fo_S
        self._S_el_labels = None
        self._S_at_int_points = None
        self._reset_masks()

    def set_fo_e(self, fo_e):
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
        return self._odb_inst

    @property
    def element_type(self):
        return self._element_type

    @property
    def fo_U(self):
        return self._fo_U.getSubset(
            position=abqConst.NODAL,
            region=self._odb_inst
        )

    @property
    def fo_S(self):
        return self._fo_S.getSubset(
            position=abqConst.INTEGRATION_POINT,
            region=self._odb_inst,
            elementType=self.element_type
        )

    @property
    def fo_e(self):
        return self._fo_e.getSubset(
            position=abqConst.INTEGRATION_POINT,
            region=self._odb_inst,
            elementType=self.element_type
        )

    @property
    def element_labels_to_node_labels_for_type(self):
        if self._element_labels_to_node_labels_for_type is None:
            odb_elements = self._odb_inst.elements

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
        return self.element_labels_to_node_labels_for_type[self.element_type]

    @property
    def node_label_to_coordinates_mapping(self):
        if self._node_labels_to_coordinates is None:

            odb_nodes = self._odb_inst.nodes
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
        if self._X_el_labels_for_type is None:
            self._update_X_values()

        return self._X_el_labels_for_type

    @property
    def X_el_labels(self):
        return self.X_el_labels_for_type[self.element_type]

    @property
    def X_at_nodes_for_type(self):
        if self._X_at_nodes_for_type is None:
            self._update_X_values()

        return self._X_at_nodes_for_type

    @property
    def X_at_nodes(self):
        return self.X_at_nodes_for_type[self.element_type]

    def _update_U_values(self):
        self._U_el_labels_for_type, self._U_at_nodes_for_type = self.map_nodes_to_element_nodal(
            self.create_node_label_to_bulk_data_mapping(self._fo_U.bulkDataBlocks)
        )

    @property
    def U_el_labels_for_type(self):
        if self._U_el_labels_for_type is None:
            self._update_U_values()

        return self._U_el_labels_for_type

    @property
    def U_el_labels(self):
        return self.U_el_labels_for_type[self.element_type]

    @property
    def U_at_nodes_for_type(self):
        if self._U_at_nodes_for_type is None:
            self._update_U_values()

        return self._U_at_nodes_for_type

    @property
    def U_at_nodes(self):
        return self.U_at_nodes_for_type[self.element_type]

    def _update_e_values(self):
        self._e_el_labels, self._e_at_int_points = self.extract_integration_points_values(
            self.fo_e.bulkDataBlocks
        )
        self._e_at_int_points = self._e_at_int_points.reshape((len(self._e_el_labels), -1))

    @property
    def e_el_labels(self):
        if self._e_el_labels is None:
            self._update_e_values()

        return self._e_el_labels

    @property
    def e_at_int_points(self):
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
        if self._S_el_labels is None:
            self._update_S_values()

        return self._S_el_labels

    @property
    def S_at_int_points(self):
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
        if self._el_labels is None:
            self._update_masks()

        return self._el_labels

    @property
    def X_mask(self):
        if self._X_mask is None:
            self._update_masks()

        return self._X_mask

    @property
    def U_mask(self):
        if self._U_mask is None:
            self._update_masks()

        return self._U_mask

    @property
    def e_mask(self):
        if self._e_mask is None:
            self._update_masks()

        return self._e_mask

    @property
    def S_mask(self):
        if self._S_mask is None:
            self._update_masks()

        return self._S_mask

    @staticmethod
    def create_node_label_to_bulk_data_mapping(bulk_data_blocks):
        data = np.concatenate([block.data for block in bulk_data_blocks])
        node_labels = np.concatenate([block.nodeLabels for block in bulk_data_blocks])

        return {
            label: data
            for label, data in zip(node_labels, data)
        }

    @staticmethod
    def extract_integration_points_values(bulk_data_blocks):
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


class FFieldOutputWriter(object):
    def __init__(self, frame, d, name):
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
        # compute deformation gradient
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


class PFieldOutputWriter(object):
    def __init__(self, frame, d, name):
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
        # compute first Piola-Kirchhoff stress tensor
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

    def flush(self, odb_inst):
        pass


class CSFieldOutputWriter(object):
    def __init__(self, frame, d, name, method):
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

    def flush(self, odb_inst):
        pass


class CFFieldOutputWriter(object):
    def __init__(self, frame, d, name, method):
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


def element_types(bulk_data_blocks):
    return {
        block.baseElementType
        for block in bulk_data_blocks
    }


def field_output_expression(field_outputs, expression):
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


def rotate_field_output_to_global_coordinate_system(frame, field_output, name, description, logger=LOGGER):
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
            repeat_labels = np.unique(labels, return_counts=True)[1][0]
            labels = labels[::repeat_labels]
        else:
            raise NotImplementedError("not supported position " + str(block.position))

        if block.localCoordSystem is None:
            data = block.data

        elif field_output.type == abqConst.VECTOR:
            ROT = rotation_matrix_from_quaternion(block.localCoordSystem)

            d = block.data.shape[1]
            assert d == 3 or d == 2
            if d == 3:
                local_vectors = block.data

            else:
                # create 3D vector out of 2D vector
                local_vectors = np.zeros((block.data.shape[0], 3), dtype=float)
                local_vectors[:, :2] = block.data

            data = np.einsum("...ji,...j", ROT, local_vectors)

        else:  # TENSOR
            d_quaternion = block.localCoordSystem.shape[1]
            assert d_quaternion == 4 or d_quaternion == 2
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
        fields=("F", "P", "CS", "CF"), 
        method="mbf", 
        e_expression="SENER+PENER",
        name_U_global_csys="U_GLOBAL_CSYS",
        name_S_global_csys="S_GLOBAL_CSYS",
        name_F="DEF_GRAD",
        name_P="FIRST_PIOLA_STRESS",
        name_CS="CONF_STRESS",
        name_CF="CONF_FORCE",
        el_type_mapping=map_abaqus_element_type_to_supported_element_type,
        logger=LOGGER
):
    path = odb.path
    is_read_only = odb.isReadOnly
    if is_read_only:
        logger.error("Odb is read-only. Reopen the odb and assure the Read-only checkbox is not checked.")
        return None

    fo_reader = FieldOutputReader()

    for odb_inst in odb.rootAssembly.instances.values():
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
                    fo_e = field_output_expression(fo, e_expression)
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

                # read odb output
                fo_reader.set_fo_U(fo_U)
                fo_reader.set_fo_e(fo_e)
                fo_reader.set_fo_S(fo_S)

                # write computed fields to odb
                d = fo_U.bulkDataBlocks[0].data.shape[1]
                fo_writers = list()
                fo_keys = set(fo.keys())

                if "F" in fields and (name_F + "_11") not in fo_keys:
                    logger.info("create field output %s_ij (%s)", name_F, msg)
                    fo_writers.append(FFieldOutputWriter(frame, d, name_F))
                elif "F" in fields:
                    logger.warning("skip field output %s_ij (%s)", name_F, msg)

                if "P" in fields and (name_P + "_11") not in fo_keys:
                    logger.info("create field output %s_ij (%s)", name_P, msg)
                    fo_writers.append(PFieldOutputWriter(frame, d, name_P))
                elif "P" in fields:
                    logger.warning("skip field output %s_ij (%s)", name_P, msg)

                if "CS" in fields and (name_CS + "_11") not in fo_keys:
                    logger.info("create field output %s_ij (%s)", name_CS, msg)
                    fo_writers.append(CSFieldOutputWriter(frame, d, name=name_CS, method=method))
                elif "CS" in fields:
                    logger.warning("skip field output %s_ij (%s)", name_CS, msg)

                if "CF" in fields and name_CF not in fo_keys:
                    logger.info("create field output %s_ij (%s)", name_CF, msg)
                    fo_writers.append(CFFieldOutputWriter(frame, d, name=name_CF, method=method))
                elif "CF" in fields:
                    logger.warning("skip field output %s_ij (%s)", name_CF, msg)

                # add data for all element types
                for element_type in element_types(fo["S"].bulkDataBlocks):
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
