import abaqus as abq
import abaqusConstants as abqConst
import caeModules as cae
import regionToolset

import numpy as np

from cf.abaqus_field_output_utils import add_field_outputs


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
    material.Elastic(table=((1000.0, 0.3),))
    section = model.HomogeneousSolidSection(
        name='Section-1',
        material=material.name,
        thickness=None)

    part.SectionAssignment(
        region=(part.cells,),
        sectionName=section.name
    )
    part.MaterialOrientation(
        region=(part.cells,),
        localCsys=None,
        axis=abqConst.AXIS_3,
        angle=30.0,
        stackDirection=abqConst.STACK_3,
        orientationType=abqConst.SYSTEM,
        additionalRotationType=abqConst.ROTATION_ANGLE
    )

    part.seedPart(size=2.5, deviationFactor=0.1, minSizeFactor=0.1)
    elemType1 = cae.mesh.ElemType(elemCode=abqConst.C3D8, elemLibrary=abqConst.STANDARD)
    elemType2 = cae.mesh.ElemType(elemCode=abqConst.C3D6, elemLibrary=abqConst.STANDARD)
    elemType3 = cae.mesh.ElemType(elemCode=abqConst.C3D4, elemLibrary=abqConst.STANDARD)
    part.setElementType(
        regions=(part.cells,),
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
        region=(instance.faces.findAt(((0., 15., 2.5),)),))
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
        nodalOutputPrecision=abqConst.FULL
    )
    job.submit()
    job.waitForCompletion()

    return abq.session.openOdb(job.name + ".odb", readOnly=False)


if __name__ == '__main__':
    odb = build()
    odb = add_field_outputs(odb)
