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
KII_MPa_m_05 = 0

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
    uy = - (
        f
        * np.cos(theta / 2)
        * (kappa - 1 - 2 * np.cos(theta / 2)**2)
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
    crack_contour_integral = assembly.engineeringFeatures.ContourIntegral(
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

    # find nodes before crack tip
    def _nodes_before_crack_tip():
        nodes_at_crack = list(assembly.sets["crack"].nodes)
        idx_node_before_crack_tip = np.argsort([node.coordinates[0] for node in nodes_at_crack])[-2]
        node_before_crack_tip = nodes_at_crack[idx_node_before_crack_tip]
        _da = -node_before_crack_tip.coordinates[0]

        # nodes that should be closed
        crack_closure_nodes = instance.nodes.getByBoundingSphere(
            center=node_before_crack_tip.coordinates,
            radius=1e-6
        )
        side_1 = assembly.Set(
            name="nodes_before_crack_tip_1",
            nodes=instance.nodes.sequenceFromLabels((crack_closure_nodes[0].label,))
        )
        side_2 = assembly.Set(
            name="nodes_before_crack_tip_2",
            nodes=instance.nodes.sequenceFromLabels((crack_closure_nodes[1].label,))
        )

        # reference point at arbitrary position
        rp_feature = assembly.ReferencePoint(
            point=(-ascending_radii_mm[-1], -ascending_radii_mm[-1], 0.0)
        )
        rp = assembly.referencePoints[rp_feature.id]
        _crack_closure_displacement = assembly.Set(
            name='crack_closure_displacement',
            referencePoints=(rp,)
        )

        # constrain crack closure nodes
        model.Equation(
            name='Crack_closure',
            terms=(
                (1.0, 'nodes_before_crack_tip_1', 2),  # first side of crack face
                (-1.0, 'nodes_before_crack_tip_2', 2),  # second side of crack face
                (1.0, 'crack_closure_displacement', 2),  # slack variable
            )
        )

        return _da, _crack_closure_displacement

    da, crack_closure_displacement = _nodes_before_crack_tip()

    # create step
    step_1 = model.StaticStep(
        name='Step-1',
        previous='Initial',
        nlgeom=abqConst.ON
    )
    step_2 = model.StaticStep(
        name='Step-2',
        previous=step_1.name,
        nlgeom=abqConst.ON
    )

    # request output
    model.fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'E', 'U', 'ENER', 'COORD'))
    model.HistoryOutputRequest(
        name='H-Output-2',
        createStepName=step_2.name,
        contourIntegral=crack_contour_integral.name,
        sectionPoints=abqConst.DEFAULT,
        rebar=abqConst.EXCLUDE,
        numberOfContours=int(np.floor((2 / np.sqrt(2)) / 0.05))
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
    bc_crack_closure = model.DisplacementBC(
        name='crack_closure',
        createStepName=step_1.name,
        region=crack_closure_displacement,
        u1=0.0,
        u2=0.0,
        u3=0.0
    )
    bc_crack_closure.deactivate(step_2.name)

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
        CF_name="CONF_FORCE",
        method="mbf",
        e_expression="SENER"
    )
    odb = apply(
        odb,
        request_CF=True,
        CF_name="CONF_FORCE_DBF",
        method="dbf",
        e_expression="SENER"
    )

    # put all results in this dictionary
    results = dict()

    # extract CF field output
    fo_CF_a0 = odb.steps['Step-2'].frames[-1].fieldOutputs["CONF_FORCE"]
    results["CF_mbf_mJ_mm2"] = [
        np.sum([
            value.data
            for value in fo_CF_a0.getSubset(region=odb.rootAssembly.nodeSets["REGION_"+str(1+i)]).values
        ], axis=0).tolist()
        for i in range(len(ascending_radii_mm))
    ]

    fo_CF_a0 = odb.steps['Step-2'].frames[-1].fieldOutputs["CONF_FORCE_DBF"]
    results["CF_dbf_mJ_mm2"] = [
        np.sum([
            value.data
            for value in fo_CF_a0.getSubset(region=odb.rootAssembly.nodeSets["REGION_"+str(1+i)]).values
        ], axis=0).tolist()
        for i in range(len(ascending_radii_mm))
    ]

    # save J-Integral
    J_integrals = odb.steps['Step-2'].historyRegions.values()[1].historyOutputs
    results["J_mJ_mm2"] = [
        J_integral.data[-1][1]
        for J_integral in J_integrals.values()
    ]

    # compute energy release rate from crack closure
    SE_000 = odb.steps['Step-1'].historyRegions['Assembly ASSEMBLY'].historyOutputs['ALLSE'].data[-1][1]
    SE_001 = odb.steps['Step-2'].historyRegions['Assembly ASSEMBLY'].historyOutputs['ALLSE'].data[-1][1]
    results["G_mJ_mm2"] = (SE_000 - SE_001) / da

    # save model parameters
    results["da_mm"] = da
    results["R_mm"] = ascending_radii_mm.tolist()

    results["inner_region_length_mm"] = inner_region_length_mm
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
