import abaqus as abq
import abaqusConstants as abqConst

import numpy as np

from cf import cf_c
from cf.abaqus_util import tensor_from_vector


class FieldOutputReader(object):
    def __init__(self):
        # odb inst
        self._odb_inst = None
        self._element_labels_to_node_labels = None
        self._node_labels_to_coordinates = None
        self._X_el_labels = None
        self._X_at_nodes = None

        # element type
        self._element_type = None

        # fo_U
        self._fo_U = None
        self._U_el_labels = None
        self._U_at_nodes = None

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
        self._element_labels_to_node_labels = None
        self._node_labels_to_coordinates = None
        self._X_el_labels = None
        self._X_at_nodes = None

        self._U_el_labels = None
        self._U_at_nodes = None
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
        self._U_el_labels = None
        self._U_at_nodes = None
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
    def element_labels_to_node_labels(self):
        if self._element_labels_to_node_labels is None:
            odb_elements = self._odb_inst.elements

            self._element_labels_to_node_labels = {
                odb_element.label: odb_element.connectivity
                for odb_element in odb_elements
            }

        return self._element_labels_to_node_labels

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
        self._X_el_labels, self._X_at_nodes = self.map_nodes_to_element_nodal(
            self.node_label_to_coordinates_mapping
        )

    @property
    def X_el_labels(self):
        if self._X_el_labels is None:
            self._update_X_values()

        return self._X_el_labels

    @property
    def X_at_nodes(self):
        if self._X_at_nodes is None:
            self._update_X_values()

        return self._X_at_nodes

    def _update_U_values(self):
        self._U_el_labels, self._U_at_nodes = self.map_nodes_to_element_nodal(
            self.create_node_label_to_bulk_data_mapping(self._fo_U.bulkDataBlocks)
        )

    @property
    def U_el_labels(self):
        if self._U_el_labels is None:
            self._update_U_values()

        return self._U_el_labels

    @property
    def U_at_nodes(self):
        if self._U_at_nodes is None:
            self._update_U_values()

        return self._U_at_nodes

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
        self._S_at_int_points = tensor_from_vector(self._S_at_int_points)

    @property
    def S_el_labels(self):
        if self._S_el_labels is None:
            self._update_S_values()

        return self._S_el_labels

    @property
    def S_at_int_points(self):
        if self._S_at_int_points is None:
            self._update_e_values()

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
        element_nodal_data = list()
        element_labels = list()
        for element_label, node_labels in self.element_labels_to_node_labels.items():
            data = list()
            for label in node_labels:
                if label in node_label_to_data_mapping:
                    data.append(node_label_to_data_mapping[label])
                else:
                    break
            else:
                element_nodal_data.append(data)
                element_labels.append(element_label)

        element_nodal_data = np.array(element_nodal_data)
        element_labels = np.array(element_labels)

        # sort
        idx = np.argsort(element_labels)
        element_labels = element_labels[idx]
        element_nodal_data = element_nodal_data[idx]

        return element_labels, element_nodal_data


class FFieldOutputWriter(object):
    def __init__(self, frame, d):
        self._fo = [
            [
                frame.FieldOutput(
                    name="DEF_GRAD_" + str(i + 1) + str(j + 1),
                    description="deformation gradient",
                    type=abqConst.SCALAR
                )
                for j in range(d)
            ]
            for i in range(d)
        ]
        self._d = d

    def add(self, reader):
        # compute deformation gradient
        f_data = cf_c.compute_F(
            X_at_nodes=reader.X_at_nodes[reader.X_mask, :, :self._d],
            U_at_nodes=reader.U_at_nodes[reader.U_mask, :, :self._d],
            element_type=reader.element_type,
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


class PFieldOutputWriter(object):
    def __init__(self, frame, d):
        self._fo = [
            [
                frame.FieldOutput(
                    name="FIRST_PIOLA_STRESS_" + str(i + 1) + str(j + 1),
                    description="First Piola-Kirchhoff stress tensor",
                    type=abqConst.SCALAR
                )
                for j in range(d)
            ]
            for i in range(d)
        ]
        self._d = d

    def add(self, reader):
        # compute first Piola-Kirchhoff stress tensor
        p_data = cf_c.compute_P(
            X_at_nodes=reader.X_at_nodes[reader.X_mask, :, :self._d],
            U_at_nodes=reader.U_at_nodes[reader.U_mask, :, :self._d],
            S_at_int_points=reader.S_at_int_points[reader.S_mask, :, :self._d, :self._d],
            element_type=reader.element_type,
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


class CSFieldOutputWriter(object):
    def __init__(self, frame, d, method):
        self._fo = [
            [
                frame.FieldOutput(
                    name="CONF_STRESS_" + str(i + 1) + str(j + 1),
                    description="configurational stresses",
                    type=abqConst.SCALAR
                )
                for j in range(d)
            ]
            for i in range(d)
        ]

        self._d = d
        self._method = method

    def add(self, reader):
        # compute configurational stresses
        cs_data = cf_c.compute_CS(
            e_at_int_points=reader.e_at_int_points[reader.e_mask],
            X_at_nodes=reader.X_at_nodes[reader.X_mask, :, :self._d],
            U_at_nodes=reader.U_at_nodes[reader.U_mask, :, :self._d],
            S_at_int_points=reader.S_at_int_points[reader.S_mask, :, :self._d, :self._d],
            element_type=reader.element_type,
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


class CFFieldOutputWriter(object):
    def __init__(self, frame, d, method):
        self._fo = frame.FieldOutput(
            name="CONF_FORCE",
            description="configurational forces",
            type=abqConst.VECTOR,
            validInvariants=[
                abqConst.MAGNITUDE
            ]
        )
        self._d = d
        self._method = method

    def add(self, reader):
        # compute configurational forces
        cf_data = cf_c.compute_CF(
            e_at_int_points=reader.e_at_int_points[reader.e_mask],
            X_at_nodes=reader.X_at_nodes[reader.X_mask, :, :self._d],
            U_at_nodes=reader.U_at_nodes[reader.U_mask, :, :self._d],
            S_at_int_points=reader.S_at_int_points[reader.S_mask, :, :self._d, :self._d],
            element_type=reader.element_type,
            method=self._method
        )

        # create datastructure abaqus understands
        cf_nodes = dict()
        for cf_el_data, el_label in zip(cf_data, reader.el_labels):
            for cf_el_node_data, node_label in zip(cf_el_data, reader.element_labels_to_node_labels[el_label]):
                cf_node = cf_nodes.setdefault(node_label, np.zeros(self._d))
                cf_node += cf_el_node_data

        # add data to field output
        self._fo.addData(
            position=abqConst.NODAL,
            instance=reader.odb_inst,
            labels=cf_nodes.keys(),
            data=np.array(cf_nodes.values())
        )


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
                    try:
                        new_term = float(division)
                    except ValueError:
                        new_term = field_outputs[division]

                    if fo_division is None:
                        fo_division = new_term
                    else:
                        fo_division = fo_division - new_term

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


def add_field_outputs(odb, fields=("F", "P", "CS", "CF"), method="mbf", e_name="SENER+PENER"):
    path = odb.path
    is_read_only = odb.isReadOnly
    if is_read_only:
        odb.save()
        odb.close()
        odb = abq.session.openOdb(path, readOnly=False)

    fo_reader = FieldOutputReader()
    csys = odb.rootAssembly.DatumCsysByThreePoints(
        name="global",
        coordSysType=abqConst.CARTESIAN,
        origin=(0, 0, 0),
        point1=(1, 0, 0),
        point2=(0, 1, 0))

    for odb_inst in odb.rootAssembly.instances.values():
        fo_reader.set_odb_inst(odb_inst)

        for step in odb.steps.values():
            for frame in step.frames:
                fo_e = field_output_expression(frame.fieldOutputs, e_name)
                fo_U = frame.fieldOutputs["U"]
                fo_S = frame.fieldOutputs["S"]

                d = fo_U.bulkDataBlocks[0].data.shape[1]

                if d == 3:
                    print("local coordinate systems are not supported")  # TODO: logger
                    """
                    fo_U = fo_U.getTransformedField(csys)
                    fo_S = fo_S.getTransformedField(csys)
                    """
                elif len(odb.rootAssembly.datumCsyses) > 0:
                    print("local coordinate systems are not supported for 2d")  # TODO: logger

                fo_reader.set_fo_U(fo_U.getTransformedField(csys))
                fo_reader.set_fo_e(fo_e)
                fo_reader.set_fo_S(fo_S.getTransformedField(csys))

                fo_writers = list()
                if "F" in fields:
                    fo_writers.append(FFieldOutputWriter(frame, d))
                if "P" in fields:
                    fo_writers.append(PFieldOutputWriter(frame, d))
                if "CS" in fields:
                    fo_writers.append(CSFieldOutputWriter(frame, d, method=method))
                if "CF" in fields:
                    fo_writers.append(CFFieldOutputWriter(frame, d, method=method))

                for element_type in element_types(frame.fieldOutputs["S"].bulkDataBlocks):
                    fo_reader.set_element_type(element_type)

                    for fo_writer in fo_writers:
                        fo_writer.add(fo_reader)

    odb.save()
    odb.close()
    return abq.session.openOdb(path, readOnly=False)