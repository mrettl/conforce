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
E_MPa = 210000
nu = 0.3
KI_MPa_m_05 = 20
KII_MPa_m_05 = 10

ascending_radii_mm, ascending_el_size_mm = np.array([
    [2, 0.05],
    [4, 0.1],
    [7, 0.5],
    [13, 2.],
    [25, 4.],
    [50, 8.],
]).T


def u_mm_KI(x, y):
    radius = (x**2 + y**2)**0.5
    theta = np.arctan2(y, x)

    G_MPa = E_MPa/(2*(1+nu))
    kappa = 3 - 4*nu  # plane strain
    # kappa = (3-nu)/(1+nu)  # plane stress

    f = (
            KI_MPa_m_05 * 1000 ** 0.5
            / (2 * G_MPa)
            * (radius / (2 * np.pi)) ** 0.5
    )

    ux = (
        f
        * np.cos(theta / 2)
        * (kappa - 1 + 2 * np.sin(theta / 2)**2)
    )
    uy = (
        f
        * np.sin(theta / 2)
        * (kappa + 1 - 2 * np.cos(theta / 2)**2)
    )
    return ux, uy


def u_mm_KII(x, y):
    radius = (x**2 + y**2)**0.5
    theta = np.arctan2(y, x)

    G_MPa = E_MPa/(2*(1+nu))
    kappa = 3 - 4*nu  # plane strain
    # kappa = (3-nu)/(1+nu)  # plane stress

    f = (
            KII_MPa_m_05 * 1000 ** 0.5
            / (2 * G_MPa)
            * (radius / (2 * np.pi)) ** 0.5
    )

    ux = (
        f
        * np.sin(theta / 2)
        * (kappa + 1 + 2 * np.cos(theta / 2)**2)
    )
    uy = (
        -f
        * np.cos(theta / 2)
        * (kappa - 1 - 2 * np.sin(theta / 2)**2)
    )
    return ux, uy


def ux_mm(x, y):
    return u_mm_KI(x, y)[0] + u_mm_KII(x, y)[0]


def uy_mm(x, y):
    return u_mm_KI(x, y)[1] + u_mm_KII(x, y)[1]


def main():
    # Clear model database
    abq.Mdb()

    # Create model
    model = abq.mdb.Model(name="K_field_model")
    del abq.mdb.models["Model-1"]

    # create part
    part = model.Part(
        name='Part-1',
        dimensionality=abqConst.TWO_D_PLANAR,
        type=abqConst.DEFORMABLE_BODY
    )

    sketch = model.ConstrainedSketch(name='sketch', sheetSize=200)
    sketch.CircleByCenterPerimeter(
        center=(0.0, 0.0),
        point1=(-ascending_radii_mm[-1], 0.0)
    )
    part.BaseShell(sketch=sketch)

    # partition part
    partition_sketch = model.ConstrainedSketch(name='partition', sheetSize=200)
    partition_sketch.Line(point1=(-ascending_radii_mm[-1], 0.0), point2=(0.0, 0.0))
    partition_sketch.Line(point1=(0.0, 0.0), point2=(ascending_radii_mm[-1], 0.0))

    for radius in ascending_radii_mm[1:-1]:
        partition_sketch.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(-radius, 0.0))

    # smallest radii is as rectangle
    inner_region_length_mm = ascending_radii_mm[0] * np.sqrt(2)
    partition_sketch.rectangle(
        point1=(-inner_region_length_mm/2, -inner_region_length_mm/2),
        point2=(inner_region_length_mm/2, inner_region_length_mm/2)
    )

    part.PartitionFaceBySketch(
        faces=part.faces,
        sketch=partition_sketch
    )

    # partition edge
    part.PartitionEdgeByParam(
        edges=part.edges.findAt(((0, ascending_radii_mm[-1], 0),), ),
        parameter=0.95
    )
    part.PartitionEdgeByParam(
        edges=part.edges.findAt(((0, -ascending_radii_mm[-1], 0),), ),
        parameter=0.05
    )

    # Set
    part_all_set = part.Set(name="all", faces=part.faces)

    # create material, section, section assignments
    material = model.Material(name="material")
    material.Elastic(table=((E_MPa, nu),))

    section = model.HomogeneousSolidSection(
        name="section",
        material="material",
        thickness=None
    )

    part.SectionAssignment(
        region=part_all_set,
        sectionName=section.name,
    )

    # create instance
    assembly = model.rootAssembly
    instance = assembly.Instance(name='Part-1-1', part=part, dependent=abqConst.OFF)

    # create sets
    crack_tip = assembly.Set(
        name="crack_tip",
        vertices=instance.vertices.findAt(((0, 0, 0),))
    )

    outer_contour = assembly.Set(
        name="outer_contour",
        edges=instance.edges.findAt(
            ((0, ascending_radii_mm[-1], 0),),
            ((0, -ascending_radii_mm[-1], 0),)
        )
    )

    middle_radii = (np.append(0, ascending_radii_mm[:-1]) + ascending_radii_mm) / 2
    crack_edges = assembly.Set(
        name="crack",
        edges=instance.edges.findAt(*[
            ((-radius, 0, 0), )
            for radius in middle_radii
        ])
    )

    region_sets = [
        assembly.Set(
            name="region_" + str(i + 1),
            faces=instance.faces.findAt(*[
                ((0, sign * radius, 0),)
                for radius in middle_radii[:i + 1]
                for sign in [1, -1]
            ])
        )
        for i in range(len(middle_radii))
    ]

    # add crack
    assembly.engineeringFeatures.assignSeam(
        regions=crack_edges
    )
    crack_contour_integral_1 = assembly.engineeringFeatures.ContourIntegral(
        name='Crack-1',
        crackFront=crack_tip,
        crackTip=crack_tip,
        symmetric=abqConst.OFF,
        extensionDirectionMethod=abqConst.Q_VECTORS,
        qVectors=((
            (0., 0., 0.),
            (1., 0., 0.)
        ),)
    )
    crack_contour_integral_2 = assembly.engineeringFeatures.ContourIntegral(
        name='Crack-2',
        crackFront=crack_tip,
        crackTip=crack_tip,
        symmetric=abqConst.OFF,
        extensionDirectionMethod=abqConst.Q_VECTORS,
        qVectors=((
                      (0., 0., 0.),
                      (0., 1., 0.)
                  ),)
    )

    # mesh
    assembly.setElementType(
        regions=(instance.faces,),
        elemTypes=(
            cae.mesh.ElemType(elemCode=abqConst.CPE4, elemLibrary=abqConst.STANDARD),
            cae.mesh.ElemType(elemCode=abqConst.CPE3, elemLibrary=abqConst.STANDARD)
        )
    )
    assembly.seedPartInstance(
        regions=(instance,),
        size=ascending_el_size_mm[-1]
    )
    for radius, el_size in zip(ascending_radii_mm[1:], ascending_el_size_mm[1:]):
        assembly.seedEdgeBySize(
            edges=instance.edges.findAt(
                ((0, radius, 0),),
                ((0, -radius, 0),)
            ),
            size=el_size
        )

    assembly.seedEdgeBySize(
        edges=instance.edges.findAt(
            ((0, inner_region_length_mm/2, 0),),
            ((0, -inner_region_length_mm/2, 0),),
            ((inner_region_length_mm/2, -inner_region_length_mm/4, 0),),
            ((inner_region_length_mm/2, inner_region_length_mm/4, 0),),
            ((-inner_region_length_mm/2, -inner_region_length_mm/4, 0),),
            ((-inner_region_length_mm/2, inner_region_length_mm/4, 0),),
        ),
        size=ascending_el_size_mm[0]
    )

    for mid_radius, min_el_size, max_el_size in zip(
        middle_radii,
        np.append(ascending_el_size_mm[0], ascending_el_size_mm[:-1]),
        ascending_el_size_mm
    ):
        assembly.seedEdgeByBias(
            end2Edges=instance.edges.findAt(((-mid_radius, 0, 0),),),
            minSize=min_el_size,
            maxSize=max_el_size,
            biasMethod=abqConst.SINGLE
        )
        assembly.seedEdgeByBias(
            end1Edges=instance.edges.findAt(((mid_radius, 0, 0),),),
            minSize=min_el_size,
            maxSize=max_el_size,
            biasMethod=abqConst.SINGLE
        )

    assembly.setMeshControls(
        regions=region_sets[0].faces,
        elemShape=abqConst.QUAD
    )
    assembly.generateMesh(regions=(instance,))

    # create step
    step_1 = model.StaticStep(
        name='Step-1',
        previous='Initial',
        nlgeom=abqConst.ON
    )

    # request output
    model.fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'E', 'U', 'ENER', 'COORD'))
    num_contours_region_1 = int(np.floor((2 / np.sqrt(2)) / 0.05))
    model.HistoryOutputRequest(
        name='H-Output-2',
        createStepName=step_1.name,
        contourIntegral=crack_contour_integral_1.name,
        sectionPoints=abqConst.DEFAULT,
        rebar=abqConst.EXCLUDE,
        numberOfContours=num_contours_region_1
    )
    model.HistoryOutputRequest(
        name='H-Output-3',
        createStepName=step_1.name,
        contourIntegral=crack_contour_integral_2.name,
        sectionPoints=abqConst.DEFAULT,
        rebar=abqConst.EXCLUDE,
        numberOfContours=num_contours_region_1
    )

    #
    nodes = instance.nodes
    contours = list()
    for i in range(num_contours_region_1):
        distance = ascending_el_size_mm[0]*(i + 0.5)
        nodes_within = nodes.getByBoundingBox(
            xMin=-distance,
            yMin=-distance,
            xMax=distance,
            yMax=distance
        )
        contour_name = "CONTOUR_" + str(i)
        contours.append(dict(
            r=-nodes_within.getBoundingBox()["low"][0],
            name=contour_name
        ))

        assembly.Set(
            name=contour_name,
            nodes=nodes_within
        )

    # field
    ux_field = model.ExpressionField(name='UX', expression='ux_mm(X, Y)')
    uy_field = model.ExpressionField(name='UY', expression='uy_mm(X, Y)')

    # boundary conditions
    model.DisplacementBC(
        name='BC_UX',
        distributionType=abqConst.FIELD,
        fieldName=ux_field.name,
        u1=1.0,
        region=outer_contour,
        createStepName=step_1.name,
    )
    model.DisplacementBC(
        name='BC_UY',
        distributionType=abqConst.FIELD,
        fieldName=uy_field.name,
        u2=1.0,
        region=outer_contour,
        createStepName=step_1.name,
    )

    #################################
    # SIMULATION AND POSTPROCESSING #
    #################################

    # create a job, start the simulation and wait until the simulation completed
    job = abq.mdb.Job(
        name='K_field_model',
        model='K_field_model',
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
        request_CS=True,
        CF_name="CONF_FORCE",
        CS_name="CONF_STRESS",
        method="mbf",
        e_expression="SENER"
    )

    # put all results in this dictionary
    results = dict()
    fo_CF_a0 = odb.steps['Step-1'].frames[-1].fieldOutputs["CONF_FORCE"]

    # regions
    results["regions"] = regions = list()
    for region_id in range(len(ascending_radii_mm)):
        region = "REGION_"+str(1+region_id)
        region_result = dict()
        regions.append(region_result)

        # extract CF field output
        region_result["CF_mbf_mJ_mm2"] = np.sum([
            value.data
            for value in fo_CF_a0.getSubset(region=odb.rootAssembly.nodeSets[region]).values
        ], axis=0).tolist()

        # radius
        if region_id == 0:
            region_result["R_mm"] = inner_region_length_mm / 2
        else:
            region_result["R_mm"] = ascending_radii_mm[region_id]

    # contours
    results["contours"] = contours_results = list()
    for contour in contours:
        contours_results.append({
            "R_mm": contour["r"],
            "CF_mbf_mJ_mm2": np.sum([
                value.data
                for value in fo_CF_a0.getSubset(region=odb.rootAssembly.nodeSets[contour["name"]]).values
            ], axis=0).tolist()
        })

    # save J-Integral
    ho = odb.steps['Step-1'].historyRegions.values()[1].historyOutputs
    results["J1_mJ_mm2"] = [
        J1_integral.data[-1][1]
        for key, J1_integral in ho.items()
        if "H-OUTPUT-2" in key
    ]

    results["J2_mJ_mm2"] = [
        J2_integral.data[-1][1]
        for key, J2_integral in ho.items()
        if "H-OUTPUT-3" in key
    ]

    # model parameters
    results["el_size_mm"] = ascending_el_size_mm.tolist()
    results["E_MPa"] = E_MPa
    results["nu"] = nu
    results["KI_MPa_m_05"] = KI_MPa_m_05
    results["KII_MPa_m_05"] = KII_MPa_m_05

    # save results
    with open("results.json", "w") as fh:
        json.dump(results, fh, indent=4)

    print(results)

    return odb


if __name__ == '__main__':
    main()
