from __future__ import print_function

import sys
import os
import json

import numpy as np

import abaqus as abq
import abaqusConstants as abqConst
import caeModules as cae

# append path to the folder where conforce_abq is located and import conforce_abq afterward
sys.path.append(os.path.abspath("../.."))
from conforce_abq.main import apply

# parameters
C10 = 86.
D1 = 1.2e-3


def main(
        name="bending_beam_model",
        lx=100., ly=50., lz=5.,
        cx=50., cy=40., r=5.,
        r1=6., r2=8.,
        uy=10.
):
    """
    creates a model with a cylindrical hole that is fixed on the left side.
    A displacement is applied on the right side.

    :param name: name of the model, the job and odb files
    :param lx: length of the model
    :param ly: height of the model
    :param lz: thickness of the model
    :param cx: x-position of the cylindrical hole
    :param cy: y-position of the cylindrical hole
    :param r: radius of the cylindrical hole
    :param r1: radius (r1 > r) of region "1" around the cylindrical hole for the evaluation of computational forces
    :param r2: radius (r2 > r1) of region "2" around the cylindrical hole for the evaluation of computational forces
    :param uy: applied displacement on the right side of the model
    :return: dictionary containing the strain energy (ALLSE) and the configurational forces of region "1" and "2" (CF_1 and CF_2)
    """

    # put all results in this dictionary
    results = dict(
        name=name,
        lx=lx,
        ly=ly,
        lz=lz,
        cx=cx,
        cy=cy,
        r=r,
        r1=r1,
        r2=r2,
        uy=uy,
        C10=C10,
        D1=D1
    )

    # Clear model database
    abq.Mdb()

    # Create model
    model = abq.mdb.Model(name=name)
    del abq.mdb.models["Model-1"]

    # create part
    part = model.Part(
        name='Part-1',
        dimensionality=abqConst.THREE_D,
        type=abqConst.DEFORMABLE_BODY
    )

    sketch = model.ConstrainedSketch(name='sketch', sheetSize=lx)
    sketch.rectangle(point1=(0.0, 0.0), point2=(lx, ly))
    sketch.CircleByCenterPerimeter(center=(cx, cy), point1=(cx+r, cy))

    part.BaseSolidExtrude(sketch=sketch, depth=lz)

    # partition part
    partition_sketch = model.ConstrainedSketch(name='partition', sheetSize=lx)
    partition_sketch.CircleByCenterPerimeter(center=(cx, cy), point1=(cx+r1, cy))
    partition_sketch.CircleByCenterPerimeter(center=(cx, cy), point1=(cx+r2, cy))

    part.PartitionFaceBySketch(
        faces=part.faces.findAt(((cx, cy*0.2, 0), )),
        sketch=partition_sketch
    )
    part.PartitionCellByExtrudeEdge(
        cells=part.cells[:],
        edges=part.edges.findAt(
            ((cx, cy+r1, 0), ),
        ),
        line=part.edges.findAt((0., 0., lz*0.2), ),
        sense=abqConst.FORWARD
    )
    part.PartitionCellByExtrudeEdge(
        cells=part.cells[:],
        edges=part.edges.findAt(
            ((cx, cy+r2, 0), ),
        ),
        line=part.edges.findAt((0., 0., lz*0.2), ),
        sense=abqConst.FORWARD
    )

    # Set
    part_all_set = part.Set(name="all", cells=part.cells)

    # create material, section, section assignments
    material = model.Material(name="material")
    material.Hyperelastic(
        type=abqConst.NEO_HOOKE,
        materialType=abqConst.ISOTROPIC,
        testData=abqConst.OFF,
        volumetricResponse=abqConst.VOLUMETRIC_DATA,
        table=((C10, D1), ))

    section = model.HomogeneousSolidSection(
        name="section",
        material="material",
        thickness=None
    )

    part.SectionAssignment(
        region=part_all_set,
        sectionName=section.name,
    )

    # mesh part
    part.seedPart(size=lx/50, deviationFactor=0.1, minSizeFactor=0.1)
    num_nodes_circle = int((2*r2*np.pi) / (lx/100))
    part.seedEdgeByNumber(
        edges=part.edges.findAt(((cx, cy+r, 0), ),),
        number=num_nodes_circle
    )
    part.seedEdgeByNumber(
        edges=part.edges.findAt(((cx, cy+r1, 0), ),),
        number=num_nodes_circle
    )
    part.seedEdgeByNumber(
        edges=part.edges.findAt(((cx, cy+r2, 0), ),),
        number=num_nodes_circle
    )
    part.setMeshControls(
        regions=part.cells.findAt(((cx, cy + (r + r1) / 2, lz*0.2),)),
        algorithm=abqConst.MEDIAL_AXIS
    )

    part.setElementType(
        regions=part_all_set,
        elemTypes=(
            cae.mesh.ElemType(elemCode=abqConst.C3D20R, elemLibrary=abqConst.STANDARD),
            cae.mesh.ElemType(elemCode=abqConst.C3D15, elemLibrary=abqConst.STANDARD),
            cae.mesh.ElemType(elemCode=abqConst.C3D10, elemLibrary=abqConst.STANDARD)
        )
    )
    part.generateMesh()

    # create instance
    assembly = model.rootAssembly
    instance = assembly.Instance(name='Part-1-1', part=part, dependent=abqConst.ON)

    rp_feature = assembly.ReferencePoint(point=(lx, ly/2, lz/2))
    rp = assembly.referencePoints[rp_feature.id]

    # create sets
    assembly.Set(
        name="REGION_A",
        cells=instance.cells.findAt(((cx, cy + (r + r1) / 2, lz*0.2),))
    )

    assembly.Set(
        name="REGION_B",
        cells=instance.cells.findAt(
            ((cx, cy + (r + r1) / 2, lz*0.2),),
            ((cx, cy + (r1 + r2) / 2, lz*0.2),),
        )
    )

    node_rp_set = assembly.Set(
        name="node_rp",
        referencePoints=(rp, )
    )

    node_000_set = assembly.Set(
        name="node_000",
        vertices=instance.vertices.findAt(((0., 0., 0.),))
    )

    face_0_set = assembly.Set(
        name="face_0",
        faces=instance.faces.findAt(((0., ly*0.2, lz*0.2),))
    )
    face_1_set = assembly.Set(
        name="face_1",
        faces=instance.faces.findAt(((lx, ly*0.2, lz*0.2),))
    )

    surface_1 = assembly.Surface(
        name='face_1',
        side1Faces=face_1_set.faces
    )

    # coupling
    model.Coupling(
        name='Constraint-1',
        surface=surface_1,
        controlPoint=node_rp_set,
        influenceRadius=abqConst.WHOLE_SURFACE,
        couplingType=abqConst.KINEMATIC,
        u1=abqConst.ON, u2=abqConst.ON, u3=abqConst.ON,
        ur1=abqConst.ON, ur2=abqConst.ON, ur3=abqConst.ON
    )

    # create step
    step = model.StaticStep(name='Step-1', previous='Initial', nlgeom=abqConst.ON)
    model.fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'LE', 'U', 'ENER'))

    # create boundary conditions
    model.PinnedBC(name='pinned', createStepName='Initial', region=face_0_set, localCsys=None)
    model.DisplacementBC(name='BC', createStepName=step.name, region=node_rp_set, u2=-uy)

    # save cae
    abq.mdb.saveAs(pathName=name + ".cae")

    #################################
    # SIMULATION AND POSTPROCESSING #
    #################################

    # create a job, start the simulation and wait until the simulation completed
    job = abq.mdb.Job(
        name=name,
        model=model.name,
        nodalOutputPrecision=abqConst.FULL
    )
    job.submit()
    job.waitForCompletion()

    # open odb with **readOnly=False**
    odb = abq.session.openOdb(job.name + ".odb", readOnly=False)

    # apply conforce plugin and request configurational forces as field output
    odb = apply(
        odb,
        request_CF=True,
        CF_name="CONF_FORCE",
        method="mbf",
        e_expression="SENER"
    )

    # compute results

    # configurational forces
    fo_CF = odb.steps['Step-1'].frames[-1].fieldOutputs["CONF_FORCE"]
    results["CF_A"] = np.sum([
        value.data
        for value in fo_CF.getSubset(region=odb.rootAssembly.nodeSets["REGION_A"]).values
    ], axis=0).tolist()
    results["CF_B"] = np.sum([
        value.data
        for value in fo_CF.getSubset(region=odb.rootAssembly.nodeSets["REGION_B"]).values
    ], axis=0).tolist()

    # strain energy
    results["ALLSE"] = odb.steps['Step-1'].historyRegions['Assembly ASSEMBLY'].historyOutputs['ALLSE'].data[-1][-1]

    return results


if __name__ == '__main__':
    name = str(sys.argv[-3])
    cx = float(sys.argv[-2])
    cy = float(sys.argv[-1])
    result = main(
        name=name,
        cx=cx,
        cy=cy
    )

    with open(name + "_results.json", "w") as fh:
        json.dump(result, fh, indent=1)
