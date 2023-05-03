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


if __name__ == '__main__':
    # odb = build()

    odb_inst = odb.rootAssembly.instances['PART-1-1']
    odb_nodes = odb_inst.nodes
    coordinates = {
        odb_node.label: odb_node.coordinates
        for odb_node in odb_nodes
    }

    odb_elements = odb_inst.elements
    n = len(odb_elements)
    element_node_labels = {
        odb_element.label: odb_element.connectivity
        for odb_element in odb_elements
    }

    frame = odb.steps['Step-1'].frames[1]
    fo = frame.fieldOutputs
    subset = dict(
        region=odb_inst,
        elementType="C3D8"
    )

    #
    e = fo["SENER"].getSubset(position=abqConst.INTEGRATION_POINT, **subset)
    e_bulk_data_blocks = e.bulkDataBlocks
    e_element_type = e_bulk_data_blocks[0].baseElementType
    e_data = e_bulk_data_blocks[0].data
    e_el_labels = e_bulk_data_blocks[0].elementLabels  # TODO: label to id
    unique_e_el_labels, idx_e_el_labels = np.unique(e_el_labels, return_inverse=True)
    e_int_point_ids = e_bulk_data_blocks[0].integrationPoints - 1
    e_el = np.zeros(
        shape=(len(unique_e_el_labels), 8),
        dtype=float
    )
    e_el[idx_e_el_labels, e_int_point_ids] = e_data[:, 0]

    #
    s = fo["S"].getSubset(position=abqConst.INTEGRATION_POINT, **subset)
    s_bulk_data_blocks = s.bulkDataBlocks
    s_element_type = s_bulk_data_blocks[0].baseElementType
    s_data = s_bulk_data_blocks[0].data
    s_el_labels = s_bulk_data_blocks[0].elementLabels  # TODO: label to id
    unique_s_el_labels, idx_s_el_labels = np.unique(e_el_labels, return_inverse=True)
    s_int_point_ids = s_bulk_data_blocks[0].integrationPoints - 1
    s_symmetric = np.zeros(
        shape=(len(unique_s_el_labels), 8, 3, 3),
        dtype=float
    )
    s_symmetric[idx_s_el_labels, s_int_point_ids, 0, 0] = s_data[:, 0]
    s_symmetric[idx_s_el_labels, s_int_point_ids, 1, 1] = s_data[:, 1]
    s_symmetric[idx_s_el_labels, s_int_point_ids, 2, 2] = s_data[:, 2]
    s_symmetric[idx_s_el_labels, s_int_point_ids, 0, 1] \
        = s_symmetric[idx_s_el_labels, s_int_point_ids, 0, 1] \
        = s_data[:, 3]
    s_symmetric[idx_s_el_labels, s_int_point_ids, 0, 2] = \
        s_symmetric[idx_s_el_labels, s_int_point_ids, 0, 2] \
        = s_data[:, 4]
    s_symmetric[idx_s_el_labels, s_int_point_ids, 1, 2] = \
        s_symmetric[idx_s_el_labels, s_int_point_ids, 1, 2] \
        = s_data[:, 5]

    #
    u = fo["U"].getSubset(position=abqConst.NODAL, region=odb_inst)
    u_bulk_data_blocks = u.bulkDataBlocks
    u_element_type = u_bulk_data_blocks[0].baseElementType
    u_data = u_bulk_data_blocks[0].data
    u_node_labels = u_bulk_data_blocks[0].nodeLabels
    u_map = {
        label: data
        for label, data in zip(u_node_labels, u_data)
    }
    u_data_el = list()
    u_el_labels = list()
    for element_label, node_labels in element_node_labels.items():
        data = list()
        for label in node_labels:
            if label in u_map:
                data.append(u_map[label])
            else:
                break
        else:
            u_data_el.append(data)
            u_el_labels.append(element_label)

    u_data_el = np.array(u_data_el)
    u_el_labels = np.array(u_el_labels)

    #
    coord_data_el = list()
    coord_el_labels = list()
    for element_label, node_labels in element_node_labels.items():
        data = list()
        for label in node_labels:
            if label in coordinates:
                data.append(coordinates[label])
            else:
                break
        else:
            coord_data_el.append(data)
            coord_el_labels.append(element_label)

    coord_data_el = np.array(coord_data_el)
    coord_el_labels = np.array(coord_el_labels)

    #
    el_labels = np.intersect1d(
        np.intersect1d(e_el_labels, s_el_labels),
        np.intersect1d(u_el_labels, coord_el_labels)
    )
    e_mask = np.in1d(unique_e_el_labels, el_labels)
    s_mask = np.in1d(unique_s_el_labels, el_labels)
    u_mask = np.in1d(u_el_labels, el_labels)
    coord_mask = np.in1d(coord_el_labels, el_labels)

    #
    cf_data = cf_c.compute_CF(
        e_at_int_points=e_el[e_mask],
        X_at_nodes=coord_data_el[coord_mask],
        U_at_nodes=u_data_el[u_mask],
        S_at_int_points=s_symmetric[s_mask],
        element_type=subset["elementType"],
        method="mbf"
    )

    cf_nodes = dict()
    for cf_el_data, el_label in zip(cf_data, el_labels):
        for cf_el_node_data, node_label in zip(cf_el_data, element_node_labels[el_label]):
            cf_node = cf_nodes.setdefault(node_label, np.zeros(3))
            cf_node += cf_el_node_data

    #
    cf_fo = frame.FieldOutput(
        name="CF",
        description="configurational forces",
        type=abqConst.VECTOR
    )
    cf_fo.addData(
        position=abqConst.NODAL,
        instance=odb_inst,
        labels=cf_nodes.keys(),
        data=np.array(cf_nodes.values())
    )
