"""
This module provides methods to export computed results to Python objects.
"""
import numpy as np

from conforce import cf_c
from conforce_abq.field_output_util import (
    FieldOutputReader,
    eval_field_output_expression,
    get_present_element_types_in
)
from conforce.element_type_mapping import map_abaqus_element_type_to_supported_element_type
import conforce_abq


LOGGER = conforce_abq.LOGGER.getChild(__name__)
"""default logger in this module"""


class ExportCF(object):
    def __init__(
            self,
            odb,
            odb_instance,
            odb_set,
            method="mbf",
            e_expression="SENER+PENER",
            name_U_global_csys="U_GLOBAL_CSYS",
            name_S_global_csys="S_GLOBAL_CSYS",
            el_type_mapping=None,
            logger=None
    ):
        """
        Compute Configurational forces without writing them to the ODB as FieldOutput.
        Set a frame before calling the export functions.

        :param odb: Odb into which the FieldOutput objects are written.
        :type odb: Odb
        :param odb_instance: Compute quantities only for this instance.
        :type odb_instance: OdbInstance
        :param odb_set: Compute quantities only for nodes and elements defined in this set.
        :type odb_set: OdbSet
        :param method: see :py:func:`conforce.cf_c.compute_CS`
        :type method: str
        :param e_expression: see :py:func:`conforce_abq.field_output_util.eval_field_output_expression`
        :type e_expression: str
        :param name_U_global_csys: name of an existing field output containing displacements
            in the global coordinate system. Use "U" if the local coordinate system is the global coordinate system.
        :type name_U_global_csys: str
        :param name_S_global_csys: name of an existing field output containing stresses
            in the global coordinate system. Use "S" if the local coordinate system is the global coordinate system.
        :param el_type_mapping: Maps element types to supported element types.
        :type el_type_mapping: Dict[str, str]
        :param logger: Print logging messages, default uses :py:attr:`LOGGER`
        :type logger: logging.Logger
        """

        self._odb = odb
        self._reader = FieldOutputReader(logger=logger)
        self._reader.set_odb_inst(odb_instance)
        self._odb_set = odb_set

        self._method = method
        self._e_expression = e_expression
        self._name_U_global_csys = name_U_global_csys
        self._name_S_global_csys = name_S_global_csys

        if el_type_mapping is None:
            el_type_mapping = map_abaqus_element_type_to_supported_element_type

        self._el_type_mapping = el_type_mapping

        if logger is None:
            logger = LOGGER

        self._logger = logger

        self._step = None
        self._frame = None
        self._dimensions = None
        self._element_types = None

    def set_frame(self, step, frame):
        """
        The Configurational forces are computed for the given frame in the given step.

        :param step: Step to compute forces for.
        :type step: OdbStep
        :param frame: Frame to compute forces for.
        :type frame: OdbFrame
        """
        self._step = step
        self._frame = frame

        fo = self._frame.fieldOutputs

        fo_U = fo[self._name_U_global_csys]
        fo_S = fo[self._name_S_global_csys]
        fo_e = eval_field_output_expression(fo, self._e_expression)

        if self._odb_set is not None and self._odb_set.nodes is not None:
            fo_U = fo_U.getSubset(region=self._odb_set)

        if self._odb_set is not None and self._odb_set.elements is not None:
            fo_S = fo_S.getSubset(region=self._odb_set)
            fo_e = fo_e.getSubset(region=self._odb_set)

        self._dimensions = fo_U.bulkDataBlocks[0].data.shape[1]
        self._element_types = get_present_element_types_in(fo_S.bulkDataBlocks)

        self._reader.set_fo_e(fo_e)
        self._reader.set_fo_U(fo_U)
        self._reader.set_fo_S(fo_S)

    def export_CF_contributions(self):
        """
        Compute element contributions to the nodal configurational forces for the frame defined in :py:meth:`set_frame`.
        The resulting nodal configurational forces are the sum of all element contributions.
        To get the contribution of element :code:`el_label` to the nodal configurational forces at node :code:`node_label`,
        use the following code snippet:

        .. code-block:: python

            CF_contributions, el_labels, node_labels, step_name, frame_id = exportCF.export_CF_contributions()
            el_label = el_labels[i]
            node_label = node_labels[i]
            [CF_x, CF_y, CF_z] = CF_contributions[i]

        :return: (CF_contributions, el_labels, node_labels, step_name, frame_id)
        :rtype: (List[np.ndarray], List[int], list[int], str, int)
        """
        if self._frame is None:
            raise ValueError("frame not defined. Set a frame first before extracting data")

        el_labels = list()
        node_labels = list()
        CF_contributions = list()
        for element_type in self._element_types:
            self._reader.set_element_type(element_type)
            supported_element_type = self._el_type_mapping[element_type]

            if len(self._reader.el_labels) == 0:
                # no elements of this type are contained in the OdbSet
                continue

            # compute configurational forces
            cf_data = cf_c.compute_CF(
                e_at_int_points=self._reader.e_at_int_points[self._reader.e_mask],
                X_at_nodes=self._reader.X_at_nodes[self._reader.X_mask, :, :self._dimensions],
                U_at_nodes=self._reader.U_at_nodes[self._reader.U_mask, :, :self._dimensions],
                S_at_int_points=self._reader.S_at_int_points[self._reader.S_mask, :, :self._dimensions, :self._dimensions],
                element_type=supported_element_type,
                method=self._method
            )

            # organize CF contributions
            el_to_n_label = self._reader.element_labels_to_node_labels
            for CF_el_data, el_label in zip(cf_data, self._reader.el_labels):
                for CF_el_node_data, node_label in zip(CF_el_data, el_to_n_label[el_label]):
                    el_labels.append(el_label)
                    node_labels.append(node_label)
                    CF_contributions.append(CF_el_node_data)

        return CF_contributions, el_labels, node_labels, self._step.name, self._frame.frameId

    def export_resulting_CF(self):
        """
        Compute nodal configurational forces for the frame defined in :py:meth:`set_frame`.
        To get the nodal configurational forces at node :code:`node_label`,
        use the following code snippet:

        .. code-block:: python

            CF, node_labels, step_name, frame_id = exportCF.export_resulting_CF()
            el_label = el_labels[i]
            node_label = node_labels[i]
            [CF_x, CF_y, CF_z] = CF[i]

        :return: (CF, node_labels, step_name, frame_id)
        :rtype: (List[np.ndarray], list[int], str, int)
        """
        CF_contributions, el_labels, node_labels, step_name, frame_id = self.export_CF_contributions()
        return np.sum(CF_contributions, axis=0), node_labels, step_name, frame_id
