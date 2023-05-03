import abaqus as abq
import abaqusConstants as abqConst
import caeModules as cae
import regionToolset

import numpy as np

from cf import cf_c


def build():
    # reset model
    abq.Mdb()

    # close odbs
    for odb in abq.session.odbs.values():
        odb.close()

    model = abq.mdb.Model(name="test")
    sketch = model.ConstrainedSketch(
        name='__profile__',
        sheetSize=200.0
    )
    sketch.Line(point1=(-10.0, 15.0), point2=(10.0, 15.0))
    sketch.Line(point1=(10.0, 15.0), point2=(10.0, 5.0))
    sketch.Line(point1=(-10.0, 15.0), point2=(-10.0, 5.0))
    sketch.ArcByCenterEnds(
        center=(-10.0, 0.0), point1=(-10.0, 5.0),
        point2=(-5.0, 0.0), direction=abqConst.CLOCKWISE)
    sketch.ArcByCenterEnds(
        center=(10.0, 0.0), point1=(10.0, 5.0),
        point2=(5.0, 0.0), direction=abqConst.COUNTERCLOCKWISE)
    sketch.Line(point1=(-5.0, 0.0), point2=(-5.0, -25.0))
    sketch.Line(point1=(-5.0, -25.0), point2=(5.0, -25.0))
    sketch.Line(point1=(5.0, -25.0), point2=(5.0, 0.0))

    part = model.Part(
        name='Part-1',
        dimensionality=abqConst.THREE_D,
        type=abqConst.DEFORMABLE_BODY
    )

    part.BaseSolidExtrude(sketch=sketch, depth=5.0)

    material = model.Material(name="Material-1")
    material.Elastic(table=((1000.0, 0.3), ))
    section = model.HomogeneousSolidSection(
        name='Section-1',
        material=material.name,
        thickness=None)

    part.SectionAssignment(
        region=(part.cells, ),
        sectionName=section.name
    )

    part.seedPart(size=2.5, deviationFactor=0.1, minSizeFactor=0.1)
    elemType1 = cae.mesh.ElemType(elemCode=abqConst.C3D8, elemLibrary=abqConst.STANDARD)
    elemType2 = cae.mesh.ElemType(elemCode=abqConst.C3D6, elemLibrary=abqConst.STANDARD)
    elemType3 = cae.mesh.ElemType(elemCode=abqConst.C3D4, elemLibrary=abqConst.STANDARD)
    part.setElementType(
        regions=(part.cells, ),
        elemTypes=(elemType1, elemType2, elemType3))
    part.generateMesh()

    assembly = model.rootAssembly
    instance = assembly.Instance(name='Part-1-1', part=part, dependent=abqConst.ON)

    step = model.StaticStep(name='Step-1', previous='Initial', nlgeom=abqConst.ON)
    model.fieldOutputRequests['F-Output-1'].setValues(variables=(
        'U', 'S', 'COORD', 'ENER', 'ELEN', 'ELEDEN'))

    model.YsymmBC(
        name="BC-1",
        createStepName="Initial",
        region=(instance.faces.findAt(((0., 15., 2.5), )),))
    model.ZsymmBC(
        name="BC-2",
        createStepName="Initial",
        region=(instance.edges.findAt(((0., 15., 0.),)),))
    model.PinnedBC(
        name="BC-3",
        createStepName="Initial",
        region=(instance.vertices.findAt(((-10., 15., 0.),)),))
    pressure_region = regionToolset.Region(
        side1Faces=instance.faces.findAt(((0., -25., 2.5),)))
    model.Pressure(
        name="Load",
        createStepName=step.name,
        magnitude=-5.0,
        region=pressure_region)

    job = abq.mdb.Job(
        name='Job-1', model=model.name, type=abqConst.ANALYSIS,
        explicitPrecision=abqConst.DOUBLE, nodalOutputPrecision=abqConst.FULL
    )
    job.submit()
    job.waitForCompletion()

    return abq.session.openOdb(job.name + ".odb", readOnly=False)


class FieldOutputAdder(object):

    def __init__(self, odb_inst):
        self._odb_inst = odb_inst

        self._element_labels_to_node_labels = None

    @property
    def element_labels_to_node_labels(self):
        if self._element_labels_to_node_labels is None:
            odb_elements = self._odb_inst.elements

            self._element_labels_to_node_labels = {
                odb_element.label: odb_element.connectivity
                for odb_element in odb_elements
            }

        return self._element_labels_to_node_labels

    def create_node_label_to_coordinates_mapping(self):
        odb_nodes = self._odb_inst.nodes

        return {
            odb_node.label: odb_node.coordinates
            for odb_node in odb_nodes
        }

    @staticmethod
    def create_node_label_to_bulk_data_mapping(bulk_data_blocks):
        data = np.concatenate([block.data for block in bulk_data_blocks])
        node_labels = np.concatenate([block.nodeLabels for block in bulk_data_blocks])

        return {
            label: data
            for label, data in zip(node_labels, data)
        }

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

    def extract_element_nodal_coordinates(self):
        return self.map_nodes_to_element_nodal(
            self.create_node_label_to_coordinates_mapping())

    def extract_element_nodal_values(self, bulk_data_blocks):
        return self.map_nodes_to_element_nodal(
            self.create_node_label_to_bulk_data_mapping(bulk_data_blocks))

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

    def add_CF_for_element_type(self, element_type, fo_U, fo_e, fo_S, fo_CF):
        #
        X_el_labels, X_at_nodes = self.extract_element_nodal_coordinates()

        #
        U_blocks = fo_U.getSubset(
            position=abqConst.NODAL,
            region=self._odb_inst

        ).bulkDataBlocks
        U_el_labels, U_at_nodes = self.extract_element_nodal_values(U_blocks)

        #
        e_blocks = fo_e.getSubset(
            position=abqConst.INTEGRATION_POINT,
            region=self._odb_inst,
            elementType=element_type
        ).bulkDataBlocks
        e_el_labels, e_at_int_points = self.extract_integration_points_values(e_blocks)
        e_at_int_points = e_at_int_points.reshape((len(e_el_labels), -1))

        #
        S_blocks = fo_S.getSubset(
            position=abqConst.INTEGRATION_POINT,
            region=self._odb_inst,
            elementType=element_type
        ).bulkDataBlocks
        S_el_labels, S_at_int_points = self.extract_integration_points_values(S_blocks)
        S_at_int_points = tensor_from_vector(S_at_int_points)

        # filter elements not defined in all outputs
        el_labels = np.intersect1d(
            np.intersect1d(X_el_labels, U_el_labels),
            np.intersect1d(e_el_labels, S_el_labels)
        )
        X_mask = np.in1d(X_el_labels, el_labels)
        U_mask = np.in1d(U_el_labels, el_labels)
        e_mask = np.in1d(e_el_labels, el_labels)
        S_mask = np.in1d(S_el_labels, el_labels)

        # compute configurational forces
        cf_data = cf_c.compute_CF(
            e_at_int_points=e_at_int_points[e_mask],
            X_at_nodes=X_at_nodes[X_mask],
            U_at_nodes=U_at_nodes[U_mask],
            S_at_int_points=S_at_int_points[S_mask],
            element_type=element_type,
            method="mbf"
        )

        #
        cf_nodes = dict()
        for cf_el_data, el_label in zip(cf_data, el_labels):
            for cf_el_node_data, node_label in zip(cf_el_data, self.element_labels_to_node_labels[el_label]):
                cf_node = cf_nodes.setdefault(node_label, np.zeros(3))
                cf_node += cf_el_node_data
        fo_CF.addData(
            position=abqConst.NODAL,
            instance=self._odb_inst,
            labels=cf_nodes.keys(),
            data=np.array(cf_nodes.values())
        )

    def add_CF(self, fo_U, fo_e, fo_S, fo_CF):
        blocks = fo_e.getSubset(
            position=abqConst.INTEGRATION_POINT,
            region=self._odb_inst
        ).bulkDataBlocks

        element_types = {
            block.baseElementType
            for block in blocks
        }

        for element_type in element_types:
            self.add_CF_for_element_type(element_type, fo_U, fo_e, fo_S, fo_CF)

    def add_CS_for_element_type(self, element_type, fo_U, fo_e, fo_S, fo_CS):
        #
        X_el_labels, X_at_nodes = self.extract_element_nodal_coordinates()

        #
        U_blocks = fo_U.getSubset(
            position=abqConst.NODAL,
            region=self._odb_inst

        ).bulkDataBlocks
        U_el_labels, U_at_nodes = self.extract_element_nodal_values(U_blocks)

        #
        e_blocks = fo_e.getSubset(
            position=abqConst.INTEGRATION_POINT,
            region=self._odb_inst,
            elementType=element_type
        ).bulkDataBlocks
        e_el_labels, e_at_int_points = self.extract_integration_points_values(e_blocks)
        e_at_int_points = e_at_int_points.reshape((len(e_el_labels), -1))

        #
        S_blocks = fo_S.getSubset(
            position=abqConst.INTEGRATION_POINT,
            region=self._odb_inst,
            elementType=element_type
        ).bulkDataBlocks
        S_el_labels, S_at_int_points = self.extract_integration_points_values(S_blocks)
        S_at_int_points = tensor_from_vector(S_at_int_points)

        # filter elements not defined in all outputs
        el_labels = np.intersect1d(
            np.intersect1d(X_el_labels, U_el_labels),
            np.intersect1d(e_el_labels, S_el_labels)
        )
        X_mask = np.in1d(X_el_labels, el_labels)
        U_mask = np.in1d(U_el_labels, el_labels)
        e_mask = np.in1d(e_el_labels, el_labels)
        S_mask = np.in1d(S_el_labels, el_labels)

        # compute configurational forces
        cs_data = cf_c.compute_CS(
            e_at_int_points=e_at_int_points[e_mask],
            X_at_nodes=X_at_nodes[X_mask],
            U_at_nodes=U_at_nodes[U_mask],
            S_at_int_points=S_at_int_points[S_mask],
            element_type=element_type,
            method="mbf"
        )

        #
        cs_vectors = list()
        for cs_at_integration_points, el_label in zip(cs_data, el_labels):
            for cs_at_integration_point in cs_at_integration_points:
                cs_vectors.append([
                    cs_at_integration_point[0, 0],
                    cs_at_integration_point[1, 1],
                    cs_at_integration_point[2, 2],
                    cs_at_integration_point[0, 1],
                    cs_at_integration_point[0, 2],
                    cs_at_integration_point[1, 2],
                ])
                # TODO: cs is not symmetric

        fo_CS.addData(
            position=abqConst.INTEGRATION_POINT,
            instance=self._odb_inst,
            labels=el_labels,
            data=cs_vectors
        )

    def add_CS(self, fo_U, fo_e, fo_S, fo_CS):
        blocks = fo_e.getSubset(
            position=abqConst.INTEGRATION_POINT,
            region=self._odb_inst
        ).bulkDataBlocks

        element_types = {
            block.baseElementType
            for block in blocks
        }

        for element_type in element_types:
            self.add_CS_for_element_type(element_type, fo_U, fo_e, fo_S, fo_CS)


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


if __name__ == '__main__':
    """
    odb = build()
    odb_inst = odb.rootAssembly.instances['PART-1-1']
    frame = odb.steps['Step-1'].frames[1]

    fo_adder = FieldOutputAdder(
        odb_inst
    )

    frame.FieldOutput(
        name="CONF_FORCE",
        description="configurational forces",
        type=abqConst.VECTOR,
        validInvariants=[
            abqConst.MAGNITUDE
        ]
    )
    fo_adder.add_CF(
        frame.fieldOutputs["U"],
        frame.fieldOutputs["SENER"],
        frame.fieldOutputs["S"],
        frame.fieldOutputs["CONF_FORCE"]
    )

    frame.FieldOutput(
        name="CONF_STRESS",
        description="configurational stresses",
        type=abqConst.TENSOR_3D_FULL,
        validInvariants=[
            abqConst.MAX_PRINCIPAL,
            abqConst.MID_PRINCIPAL,
            abqConst.MIN_PRINCIPAL,
        ]
    )"""
    fo_adder.add_CS(
        frame.fieldOutputs["U"],
        frame.fieldOutputs["SENER"],
        frame.fieldOutputs["S"],
        frame.fieldOutputs["CONF_STRESS"]
    )
