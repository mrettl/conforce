# Model with crack tip, plastic and CC
# ------------------------------------------------------------------
# Martin Pletz, 2023-05

from abaqus import *
from abaqusConstants import *
from caeModules import *
import os,shutil
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(8, r'c:/SIMULIA/CAE/plugins/2020/cf')
from cf_abq import main as cf_abq_main

DIR0 = os.path.abspath('')
TOL = 1e-4
session.journalOptions.setValues(replayGeometry=COORDINATE)

# general model functions
# ------------------------------------------------------------------

def remove_files(dir0,type_list=('com','sim','prt','msg','mdl','reg','7','simlog',
                                 'exception','023','simdir','SMABulk','stt')):
    # get all files and folders in dir0
    file_list = os.listdir(dir0)
    # go through the files and delete some of them
    for file in file_list:
        if file.split('.')[-1] in type_list:
            try:
                os.remove(dir0+'/'+file)                    
            except:
                shutil.rmtree(dir0+'/'+file)
                # print('file '+file+' could not be deleted!')
    return

def make_dir(dir_name, if_change=0, if_clear=0):
    # change to subfolder and create that folder if it does not already exist
    #dir_abs = os.path.abspath('')
    if os.path.exists(dir_name) == 0:
        os.mkdir(dir_name)
    else:
        if if_clear:
            shutil.rmtree(dir_name)
            os.mkdir(dir_name)
    #dir1 = dir_abs + "//" + dir_name
    if if_change:
        os.chdir(dir_name)
    return dir_name

def get_rtheta(x, y):
    """convert polar coordinates to carthesian coordinates
    """
    r, theta = (x**2 + y**2)**0.5, np.arctan2(y,x)
    return r,theta

def get_u(K_I, E, nu, pos, is_pe=1):
    """Calculate (u_x,u_y) [mm] around crack tip according to the formulas from Andersson.

    Args:
        K_I (float >= 0): Mode I K-factor [MPa m**0.5]
        E (float): Young's modulus [MPa]
        nu (float): Poisson's ratio [1]
        pos (tuple): (r, theta), position for getting u [mm,rad]
        is_pe (boolean): if plane strain (1) or plane stress (0)
    """
    r, theta = pos

    G = E/(2*(1+nu))
    if is_pe:
        kappa = 3 - 4*nu
    else:
        kappa = (3-nu)/(1+nu)

    fac0 = K_I*1000**0.5/(2*G)*(r/(2*np.pi))**0.5
    ux = fac0 * np.cos(theta/2)*(kappa-1+2*np.sin(theta/2)**2)
    uy = fac0 * np.sin(theta/2)*(kappa+1-2*np.cos(theta/2)**2)
    return ux, uy

def make_sym_nsets(p,a0=0):
    # create node sets in the part p for opening them

    p.Set(name='nodes_open', nodes=p.nodes.getByBoundingBox(yMax=TOL, xMin=a0-TOL))
    nodes_open = sorted([(node.coordinates[0],node.label) for node in p.sets['nodes_open'].nodes])

    # open_list: (y0,dy,set_name)
    crack_dict = {'sets':[], 'c_vals':[]}
    #
    # make sets of 2 nodes
    for i_nodes in range(len(nodes_open)):
        nkey_1 = nodes_open[i_nodes][1]
        c_temp = nodes_open[i_nodes][0]
        #
        crack_dict['sets'] += ['nopen-set-'+str(i_nodes).zfill(3)]
        crack_dict['c_vals'] += [c_temp]
        p.Set(name='nopen-set-'+str(i_nodes).zfill(3), nodes=p.nodes[nkey_1-1:nkey_1])
    #
    crack_dict['c_vals'] = np.array(crack_dict['c_vals'])
    crack_dict['sets'] = np.array(crack_dict['sets'])
    #
    return crack_dict

# crack model functions
# ------------------------------------------------------------------

def make_mesh_mat(model,is_pe,is_lin,(E,nu,sig_y0,B,n)):
    
    # set the element types
    # --------------------------
    p = model.parts.values()[0]
    
    if is_pe:
        if is_lin:
            p.setElementType(regions=(p.faces,), elemTypes=(mesh.ElemType(elemCode=CPE4R, 
                    elemLibrary=STANDARD), mesh.ElemType(elemCode=CPE3, elemLibrary=STANDARD)))
        else:
            p.setElementType(regions=(p.faces,), elemTypes=(mesh.ElemType(elemCode=CPE8R, 
                    elemLibrary=STANDARD), mesh.ElemType(elemCode=CPE6, elemLibrary=STANDARD)))
    else:
        if is_lin:
            p.setElementType(regions=(p.faces,), elemTypes=(mesh.ElemType(elemCode=CPS4R, 
                    elemLibrary=STANDARD), mesh.ElemType(elemCode=CPS3, elemLibrary=STANDARD)))
        else:
            p.setElementType(regions=(p.faces,), elemTypes=(mesh.ElemType(elemCode=CPS8R,# CPS8R
                    elemLibrary=STANDARD), mesh.ElemType(elemCode=CPS6, elemLibrary=STANDARD)))

    # create instance and sections
    # --------------------------
    ass = model.rootAssembly
    inst = ass.Instance(name='crack', part=p, dependent=ON)

    # create material
    mat = model.Material(name='steel')
    mat.Elastic(table=((E,nu), ))

    # some conditions for plastic part of material model
    if sig_y0 != -1:
        if B == 0:
            mat.Plastic(table=((sig_y0, 0.0),))
        else:
            mat.Plastic(hardening=JOHNSON_COOK,
                table=((sig_y0, B, n, 0.0, 1000, 500), ))

    # Section erstellen und zuweisen
    model.HomogeneousSolidSection(material='steel', name='steel', thickness=None)
    p.SectionAssignment(region=p.sets['all'], sectionName='steel',
                        thicknessAssignment=FROM_SECTION)
    return inst


def make_geom_circ(model,r_list,mesh_sizes):
    """create the geometry of the crack model

    Args:
        model (abaqus model): the model
        R (float): Radius of the crack model
        mesh_size,coarse_fac (float,float): Mesh size in fine region and the factor that the coarse region should be more coarse
        E,nu,sig_y0,B,n (all float): Young's modulus, Poisson's ratio, initial yield stress (if -1, linear-el. material), hardening factor B (if 0, ideal plastic) and hardening exponent n. (sig_y = sig_y0 + B*eps_pl**n)
        is_pe (int, optional): If the 2D model should be plane strain (1) or plane stress (0). Defaults to 1.

    Returns:
        abaqus instance: the crack instance.
    """
    R = r_list[-1]

    # create crack part
    s = model.ConstrainedSketch(name='crack', sheetSize=200.0)

    s.ArcByCenterEnds(center=(0,0), point1=(-R,0), point2=(R,0), 
                      direction=CLOCKWISE)

    s.Line(point1=(-R,0), point2=(0,0))
    s.Line(point1=(0,0), point2=(R,0))

    p = model.Part(name='crack', dimensionality=TWO_D_PLANAR, 
                   type=DEFORMABLE_BODY)

    p.BaseShell(sketch=s)

    # Sets and Surfaces
    p.Set(name='all', faces=p.faces)
    p.Set(name='load_u', edges=p.edges.findAt(coordinates=((0,R,0),)))
    #p.Set(name='nodes_open', edges=p.edges.getByBoundingBox(yMax=TOL, xMin=-TOL))
    p.Set(name='crack_tip', vertices=p.vertices.findAt(coordinates=((0,0,0),)))
    
    s_part = model.ConstrainedSketch(name='partition', sheetSize=200.0)

    r_fine = r_list[0]
    # rectangle partition
    s_part.Line(point1=(-r_fine/4., 0), point2=(-r_fine/4., r_fine/4.))
    s_part.Line(point1=(r_fine/2., 0), point2=(r_fine/2., r_fine/4.))
    s_part.Line(point1=(-r_fine/4., r_fine/4.), point2=(r_fine/2., r_fine/4.))


    for r_fine in r_list[:-1]:
        s_part.CircleByCenterPerimeter(center=(0, 0), point1=(r_fine, 0))

    p.PartitionFaceBySketch(faces=p.faces, sketch=s_part)
    
    
    p.Set(name='cf-eval', faces=p.faces.getByBoundingBox(yMax=r_list[0]/4.+TOL))
    
    #mdb.saveAs('test')
    #raise ValueError()
    
    p.seedPart(size=mesh_sizes[1])

    # seed the rectangle
    p.seedEdgeBySize(edges=p.edges.getByBoundingBox(yMax=r_fine/2.+TOL,xMin=-r_fine/2.-TOL,xMax=r_fine+TOL),
                     size=mesh_sizes[0], constraint=FINER)    

    for i_part in range(len(r_list)):
        
        r_temp = r_list[i_part]
        m_size = mesh_sizes[i_part+1]

        p.seedEdgeBySize(edges=p.edges.findAt(((0,r_temp,0),)),
                        size=m_size, constraint=FINER)

        if i_part > 0:
            mean_size = (m_size+mesh_sizes[i_part-1])/2.

            p.seedEdgeBySize(edges=p.edges.findAt(((r_temp-TOL,0,0),)),
                        size=mean_size, constraint=FINER)
            p.seedEdgeBySize(edges=p.edges.findAt(((-r_temp+TOL,0,0),)),
                        size=mean_size, constraint=FINER)
            
    p.setMeshControls(regions=p.faces, elemShape=QUAD)
    p.generateMesh()

    return 

def make_geom_rect(model,b_list,h_list,a0,mesh_sizes):
    """create the geometry of the crack model

    Args:
        model (abaqus model): the model
        R (float): Radius of the crack model
        mesh_size,coarse_fac (float,float): Mesh size in fine region and the factor that the coarse region should be more coarse
        E,nu,sig_y0,B,n (all float): Young's modulus, Poisson's ratio, initial yield stress (if -1, linear-el. material), hardening factor B (if 0, ideal plastic) and hardening exponent n. (sig_y = sig_y0 + B*eps_pl**n)
        is_pe (int, optional): If the 2D model should be plane strain (1) or plane stress (0). Defaults to 1.

    Returns:
        abaqus instance: the crack instance.
    """

    B,H = b_list[-1], h_list[-1]

    # create crack part
    s = model.ConstrainedSketch(name='plate', sheetSize=200.0)

    s.rectangle(point1=(0,0), point2=(B,H))

    p = model.Part(name='plate', dimensionality=TWO_D_PLANAR, 
                   type=DEFORMABLE_BODY)

    p.BaseShell(sketch=s)
    
    rp = p.ReferencePoint(point=(0,H,0))
    p.Set(name='RP', referencePoints=(p.referencePoints[rp.id],)) 

    # Sets and Surfaces
    p.Set(name='all', faces=p.faces)
    p.Set(name='top-edge', edges=p.edges.getByBoundingBox(yMin=H-TOL))
    p.Set(name='left-edge', edges=p.edges.getByBoundingBox(xMax=TOL))

    s_part = model.ConstrainedSketch(name='partition', sheetSize=200.0)

    s_part.Line(point1=(0,h_list[0]), point2=(b_list[0],h_list[0]))
    s_part.Line(point1=(b_list[0],0), point2=(b_list[0],h_list[0]))

    for b_i,h_i in zip(b_list,h_list)[1:-1]:
        s_part.Line(point1=(0,h_i), point2=(b_i,h_i))
        s_part.ArcByCenterEnds(center=(b_i,0), point1=(b_i,h_i), point2=(b_i+h_i,0), 
                               direction=CLOCKWISE)

    p.PartitionFaceBySketch(faces=p.faces, sketch=s_part)

    p.seedPart(size=mesh_sizes[-1])

    # seed the part regions (outside to inside)
    for i_part in range(len(h_list)-1)[::-1]:
        
        b_temp = b_list[i_part]
        h_temp = h_list[i_part]
        m_size = mesh_sizes[i_part]

        p.seedEdgeBySize(edges=p.edges.getByBoundingBox(xMax=b_temp+h_temp+TOL,yMax=h_temp+TOL),
                        size=m_size, constraint=FINER)
    p.generateMesh()

    return


def make_loads(model,inst,is_k,K_I,a0,n_open,E,nu,is_pe=1,is_lin=0,uy=0,n_incs=10,appl_u=1):
    
    p = model.parts.values()[0]
    # make the node sets for opening the crack
    crack_dict = make_sym_nsets(p,a0=a0)
    
    if not is_k:
        # couple top nodes to RP
        model.Coupling(name='Couple-Y-top', controlPoint=inst.sets['RP'], surface=inst.sets['top-edge'],
                    influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, u1=OFF, u2=ON, ur3=OFF)

    # Step und Loads
    model.StaticStep(name='load-000', previous='Initial', maxNumInc=1000, initialInc=0.05, 
                    minInc=1e-08, maxInc=1, nlgeom=ON)
    
    if is_k:
        # apply the ux,uy from the analytical K_I field

        nodes_K = sorted([(node.coordinates[0],node.coordinates[1],node.label)
                        for node in inst.sets['load_u'].nodes])
        
        for x,y,node_label in nodes_K:
            r, theta = get_rtheta(x,y)
            ux, uy =  get_u(K_I, E, nu, (r, theta), is_pe)
            # create set and load this node
            p.Set(name='appl-K-n'+str(node_label), nodes=p.nodes[node_label-1:node_label])

            model.DisplacementBC(createStepName='load-000', name='appl-K-n'+str(node_label),
                                region=inst.sets['appl-K-n'+str(node_label)], u1=ux, u2=uy)
    else:
        # x-syymetry on the left edge
        model.DisplacementBC(createStepName='Initial', name='xsym',
                                region=inst.sets['left-edge'], u1=0)
        
        model.HistoryOutputRequest(name='H-Output-3', createStepName='load-000',
                                   variables=('U2', 'RF2'), region=inst.sets['RP'])
        
        if appl_u:
            # apply the vertical displacement on the top edge
            model.DisplacementBC(createStepName='load-000', name='appl-uyy',
                                    region=inst.sets['RP'], u1=0, ur3=0, u2=uy)
        else:
            None

    # open the crack node by node until n_open nodes are opened
    # --------------------------------------------------------------

    if is_k:
        # leave the last so that the u from K_I can be applied there
        for set_name in crack_dict['sets'][:-1]:
            model.DisplacementBC(createStepName='Initial', name=set_name,
                                region=inst.sets[set_name], u2=0)
    else:
        # all sets for the rect
        for set_name in crack_dict['sets']:
            model.DisplacementBC(createStepName='Initial', name=set_name,
                                region=inst.sets[set_name], u2=0)        

    c_list = []
    
    for i_open in range(n_open):
        step_name = 'load-'+str(i_open+1).zfill(3)
        
        model.StaticStep(name=step_name, nlgeom=ON, previous='load-'+str(i_open).zfill(3),
                         initialInc=1./n_incs, maxInc=1./n_incs)

        if is_lin:
            model.boundaryConditions[crack_dict['sets'][i_open]].deactivate(step_name)
            c_list += [crack_dict['c_vals'][i_open]]
        else:
            model.boundaryConditions[crack_dict['sets'][i_open*2]].deactivate(step_name)
            model.boundaryConditions[crack_dict['sets'][i_open*2+1]].deactivate(step_name)
            c_list += [crack_dict['c_vals'][i_open*2]]

    np.savetxt(model.name+'_a_list.dat',np.array(c_list))

    # output the energy densities
    model.FieldOutputRequest(name='F-Output-2', createStepName='load-000', variables=('ENER', 'COORD', ))
    
    # run the job
    # --------------------------

    mdb.saveAs(model.name)
    job = mdb.Job(name=model.name, model=model.name, type=ANALYSIS)

    job.submit()
    job.waitForCompletion()
    return

def eval_sig_energies_circ(job_name,R,K_I,E,nu,is_pe=1,sig_y=300.,mesh_size=0.03,if_cf=1):
    # History Output: Diagram and dat file
    if if_cf:
        odb = session.openOdb(name=job_name+'.odb', readOnly=False)
        session.viewports['Viewport: 1'].setValues(displayedObject=odb)
    
        cf_abq_main.main(odb_name=odb.name, method='motion based formulation',
                        e_expression='SENER+PENER', request_F=False, request_P=False,
                        request_CS=False, CS_name='CON_S', request_CF=True, CF_name='CON_F_SP')
        
        odb = session.openOdb(name=job_name+'.odb', readOnly=False)
        session.viewports['Viewport: 1'].setValues(displayedObject=odb)
        
        cf_abq_main.main(odb_name=odb.name, method='motion based formulation',
                        e_expression='SENER', request_F=False, request_P=False,
                        request_CS=False, CS_name='CON_S', request_CF=True, CF_name='CON_F_S')
        odb = session.openOdb(name=job_name+'.odb')
    else:
        odb = session.openOdb(name=job_name+'.odb')

    out_arr = []
    cf_s_list,cf_sp_list = [],[]

    # define CF evaluation region and coords for CF output
    eval_region = odb.rootAssembly.instances['CRACK'].nodeSets['CF-EVAL']
    coord_vals = odb.steps.values()[0].frames[0].fieldOutputs['COORD'].getSubset(region=eval_region).values
    
    coord_res = np.array([[coord.nodeLabel]+list(coord.data) for coord in coord_vals])
    np.savetxt(job_name+'_coord.dat',coord_res,header='label, x, y')

    for i_step,step in enumerate(odb.steps.values()):
        hr_ass = step.historyRegions.values()[0]

        t = i_step
        allse = hr_ass.historyOutputs['ALLSE'].data[-1][-1] 
        allpd = hr_ass.historyOutputs['ALLPD'].data[-1][-1]
        #
        out_arr += [[t,allse,allpd]]
        
        # die CF in ein file schreiben
        if if_cf:
            # read coord, cf and write to dat file
        
            cf_vals_sp = step.frames[-1].fieldOutputs['CON_F_SP'].getSubset(region=eval_region).values
            cf_vals_s = step.frames[-1].fieldOutputs['CON_F_S'].getSubset(region=eval_region).values

            cf_res = np.array([[cf_sp.nodeLabel]+list(cf_sp.data)+list(cf_s.data) for cf_sp, cf_s in zip(cf_vals_sp,cf_vals_s)])
            
            cf_s_list += [-2*np.sum(cf_res[:,-2])]
            cf_sp_list += [-2*np.sum(cf_res[:,1])]
            
            np.savetxt(job_name+'_cf_'+step.name+'.dat',cf_res,header='label, CFx-SP, CFy-SP, CFx-S, CFy-S')
    
    out_arr = np.array(out_arr)
    
    np.savetxt(job_name+'_res.dat',out_arr,delimiter=',',
            header='t (s), ALLSE (Nmm), ALLPD (Nmm)')

    # Pfad fuer die Spannung
    vp = session.viewports['Viewport: 1']
    vp.setValues(displayedObject=odb)

    session.Path(name='pfad', type=POINT_LIST, expression=((0,0,0), (R,0,0)))
    # field output to be evaluated along path
    vp.odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF, ))
    vp.odbDisplay.setPrimaryVariable(variableLabel='S', 
            outputPosition=INTEGRATION_POINT,refinement=(COMPONENT, 'S22'))
    vp.odbDisplay.setFrame(step=1, frame=0)
    
    # output for path
    xy = xyPlot.XYDataFromPath(path=session.paths['pfad'], includeIntersections=True, 
                projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                projectionTolerance=0, shape=UNDEFORMED, labelType=TRUE_DISTANCE_X, 
                removeDuplicateXYPairs=True, includeAllElements=False)
    
    xy = [[float(i),float(j)] for i,j in xy]
    #
    np.savetxt(job_name+'_s_res.dat',xy,delimiter=',',
            header='x (mm), sYY (MPa)')
    
    try:
        session.Path(name='pfad-pe', type=POINT_LIST, expression=((0,mesh_size/2.,0), (R,mesh_size/2,0)))
        # PE22 ausgeben
        vp.odbDisplay.setPrimaryVariable(variableLabel='PE', outputPosition=INTEGRATION_POINT,
                                        refinement=(COMPONENT, 'PE22'), ) 
        vp.odbDisplay.setPrimaryVariable(variableLabel='PEEQ', outputPosition=INTEGRATION_POINT) 

        xy = xyPlot.XYDataFromPath(path=session.paths['pfad-pe'], includeIntersections=True, 
                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                    projectionTolerance=0, shape=UNDEFORMED, labelType=TRUE_DISTANCE_X, 
                    removeDuplicateXYPairs=True, includeAllElements=False)
        
        xy = [[float(i),float(j)] for i,j in xy]
        #
        np.savetxt(job_name+'_peeq_res.dat',xy,delimiter=',',
                header='x (mm), PEEQ (1)')
    except:
        print('anscheinend lin.el. Rechnung')

    # TODO: Rissverlaengerung plastisch mitnehmen?!?, also alle Knoten auf der Linie, die zu oeffnen ist :-)

    # draw the diagram of the stress and the energy release rate
    # --------------------------------------------
    res_e = np.loadtxt(job_name+'_res.dat', delimiter=',', skiprows=1)
    res_a = np.loadtxt(job_name+'_a_list.dat', delimiter=',')
    res_s = np.loadtxt(job_name+'_s_res.dat', delimiter=',', skiprows=1)

    # --------- this part only works for Abaqus 2020 or later -----------
    if is_pe:
        r_pl = 1000/(6*np.pi)*K_I**2/sig_y**2
    else:
        r_pl = 1000/(2*np.pi)*K_I**2/sig_y**2

    plt.figure(figsize=(4,4))
    # compute G_inc
    G_tot_inc = []
    G_wo_pl_inc = []

    for i in range(1,len(res_a)):
        a_temp = res_a[i]
        dE_el = res_e[i,1] - res_e[0,1]
        dE_pl = res_e[i,2] - res_e[0,2]

        G_tot_inc += [-dE_el/a_temp*2]
        G_wo_pl_inc += [-(dE_el+dE_pl)*2/a_temp]

    # plot, plot, plot
    if is_pe:
        G_pe = K_I**2*(1-nu**2)/E*1000
    else:
        G_pe = K_I**2/E*1000

    plt.plot(res_a[1::],G_tot_inc, '.-',label='incremental $G_{tot}$')
    plt.plot(res_a[1::],G_wo_pl_inc, '.-', label='incremental $G_{nopl}$')

    plt.plot(res_a,cf_sp_list[:-1],'.-',label='CF (SE+PE)')
    plt.plot(res_a,cf_s_list[:-1],'.-',label='CF (SE)')

    plt.plot([0,max(res_a)],[G_pe,G_pe],'--',c='Grey', label='global $G$')

    plt.plot((r_pl,r_pl),(0,max(G_tot_inc)),'--',c='c')

    plt.plot(res_s[:,0],res_s[:,1]/1000, '.-', label='stress $\sigma_{yy}$ [GPa]')

    try:
        res_eps = np.loadtxt(job_name+'_peeq_res.dat', delimiter=',', skiprows=1)
        # PE: fast keine pl. Dehnungen!!
        # ersten Punkt nciht plotten
        take_none0 = res_eps[:,0]>TOL

        plt.plot(res_eps[take_none0,0],res_eps[take_none0,1]*100, '.-', label=r'PEEQ [%]')
    except:
        print('noch immer keine plastischen Dehnungen!')
    
    plt.legend()

    plt.xlabel('crack extension $\Delta a$ (mm)')
    plt.ylabel('Energy release rate (kJ/m$^2$)')

    plt.xlim(xmin=0, xmax=max(res_a)) #0.6)
    plt.ylim(ymin=0)

    plt.title('$K$ = '+str(int(K_I))+' MPa m$^{0.5}$')
    plt.tight_layout()

    plt.savefig(job_name+'_sig_G.pdf')
    plt.savefig(job_name+'_sig_G.png',dpi=400)
    plt.clf()
    #plt.show()
    #
    odb.close()
    return

def eval_sig_energies_rect(job_name,B,a0,eps_yy,mesh_size=0.03):
    # History Output: Diagram and dat file
    odb = session.openOdb(name=job_name+'.odb')

    out_arr = []

    for i_step,step in enumerate(odb.steps.values()):
        hr_ass = step.historyRegions.values()[0]

        t = i_step
        allse = hr_ass.historyOutputs['ALLSE'].data[-1][-1] 
        allpd = hr_ass.historyOutputs['ALLPD'].data[-1][-1]
        #
        out_arr += [[t,allse,allpd]]
    
    out_arr = np.array(out_arr)
    
    np.savetxt(job_name+'_res.dat',out_arr,delimiter=',',
            header='t (s), ALLSE (Nmm), ALLPD (Nmm)')

    # Pfad fuer die Spannung
    vp = session.viewports['Viewport: 1']
    vp.setValues(displayedObject=odb)

    session.Path(name='pfad', type=POINT_LIST, expression=((a0,0,0), (B,0,0)))
    # field output to be evaluated along path
    vp.odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF, ))
    vp.odbDisplay.setPrimaryVariable(variableLabel='S', 
            outputPosition=INTEGRATION_POINT,refinement=(COMPONENT, 'S22'))
    vp.odbDisplay.setFrame(step=1, frame=0)
    
    # output for path
    xy = xyPlot.XYDataFromPath(path=session.paths['pfad'], includeIntersections=True, 
                projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                projectionTolerance=0, shape=UNDEFORMED, labelType=TRUE_DISTANCE_X, 
                removeDuplicateXYPairs=True, includeAllElements=False)
    
    xy = [[float(i),float(j)] for i,j in xy]
    #
    np.savetxt(job_name+'_s_res.dat',xy,delimiter=',',
            header='x (mm), sYY (MPa)')
    
    try:
        session.Path(name='pfad-pe', type=POINT_LIST, expression=((a0,mesh_size/2.,0), (B,mesh_size/2,0)))
        # PE22 ausgeben
        vp.odbDisplay.setPrimaryVariable(variableLabel='PE', outputPosition=INTEGRATION_POINT,
                                        refinement=(COMPONENT, 'PE22'), ) 
        vp.odbDisplay.setPrimaryVariable(variableLabel='PEEQ', outputPosition=INTEGRATION_POINT) 

        xy = xyPlot.XYDataFromPath(path=session.paths['pfad-pe'], includeIntersections=True, 
                    projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
                    projectionTolerance=0, shape=UNDEFORMED, labelType=TRUE_DISTANCE_X, 
                    removeDuplicateXYPairs=True, includeAllElements=False)
        
        xy = [[float(i),float(j)] for i,j in xy]
        #
        np.savetxt(job_name+'_peeq_res.dat',xy,delimiter=',',
                header='x (mm), PEEQ (1)')
    except:
        print('anscheinend lin.el. Rechnung')

    # TODO: Rissverlaengerung plastisch mitnehmen?!?, also alle Knoten auf der Linie, die zu oeffnen ist :-)

    # draw the diagram of the stress and the energy release rate
    # --------------------------------------------
    res_e = np.loadtxt(job_name+'_res.dat', delimiter=',', skiprows=1)
    res_a = np.loadtxt(job_name+'_a_list.dat', delimiter=',')
    res_s = np.loadtxt(job_name+'_s_res.dat', delimiter=',', skiprows=1)

    # --------- this part only works for Abaqus 2020 or later -----------

    plt.figure(figsize=(4,4))
    # compute G_inc
    G_tot_inc = []
    G_wo_pl_inc = []

    for i in range(1,len(res_a)):
        a_temp = res_a[i]
        dE_el = res_e[i,1] - res_e[0,1]
        dE_pl = res_e[i,2] - res_e[0,2]

        G_tot_inc += [-dE_el/(a_temp-a0)*2]
        G_wo_pl_inc += [-(dE_el+dE_pl)*2/(a_temp-a0)]

    # plot, plot, plot


    #plt.plot(res_a[1::],G_tot_inc, '.-',label='incremental $G_{tot}$')
    plt.plot(res_a[1::],G_wo_pl_inc, '.-', label='incremental $G_*$')

    plt.plot(res_s[:,0],res_s[:,1]/1000, '.-', label='stress $\sigma_{yy}$ [GPa]')

    try:
        res_eps = np.loadtxt(job_name+'_peeq_res.dat', delimiter=',', skiprows=1)
        # PE: fast keine pl. Dehnungen!!
        # ersten Punkt nciht plotten
        take_none0 = res_eps[:,0]>TOL

        plt.plot(res_eps[take_none0,0],res_eps[take_none0,1]*100, '.-', label=r'PEEQ [%]')
    except:
        print('noch immer keine plastischen Dehnungen!')
    
    plt.legend()

    plt.xlabel('crack extension $\Delta a$ (mm)')
    plt.ylabel('Energy release rate (kJ/m$^2$)')

    plt.xlim(xmin=a0, xmax=max(res_a)) #0.6)
    plt.ylim(ymin=0)

    plt.title('a_0 = '+str(a0)+ ' mm') #'$K$ = '+str(int(eps_yy))+' MPa m$^{0.5}$')
    plt.tight_layout()

    plt.savefig(job_name+'_sig_G.pdf')
    plt.savefig(job_name+'_sig_G.png',dpi=400)
    plt.clf()
    #plt.show()
    #
    odb.close()
    return

# call the models

def make_crack_model_circ(dir_name,model_name,r_list,(E,nu,sig_y0,B,n),mesh_sizes,
                     K_I,n_open,is_pe,is_lin):
    
    # reset model and give it the model_name
    Mdb()
    mdb.models.changeKey(fromName='Model-1', toName=model_name)
    model = mdb.models[model_name]

    # change into the subfolder for running model
    make_dir(dir_name,if_change=1)

    # create geometry and mesh it
    make_geom_circ(model,r_list,mesh_sizes)
    inst = make_mesh_mat(model,is_pe,is_lin,(E,nu,sig_y0,B,n))

    # create load steps with opening nodes and run the model
    make_loads(model,inst,1,K_I,0,n_open,E,nu,is_pe,is_lin)

    # evaluate the stress path and incremental G, plot it already
    # (working for Abq 2020 or later)
    eval_sig_energies_circ(model_name,r_list[1],K_I,E,nu,is_pe,sig_y0,mesh_sizes[0])
    
    remove_files(DIR0+'/'+dir_name)
    os.chdir(DIR0)
    return

def make_crack_model_rect(dir_name,model_name,b_list,h_list,a0,(E,nu,sig_y0,B,n),mesh_sizes,
                     eps_yy,n_open,is_pe,is_lin):
    
    # reset model and give it the model_name
    Mdb()
    mdb.models.changeKey(fromName='Model-1', toName=model_name)
    model = mdb.models[model_name]

    # change into the subfolder for running model
    make_dir(dir_name,if_change=1)

    # create geometry and mesh it
    make_geom_rect(model,b_list,h_list,a0,mesh_sizes)
    inst = make_mesh_mat(model,is_pe,is_lin,(E,nu,sig_y0,B,n))

    # create load steps with opening nodes and run the model
    make_loads(model,inst,0,0,a0,n_open,E,nu,is_pe,is_lin,uy=eps_yy*h_list[-1])

    # evaluate the stress path and incremental G, plot it already
    # (working for Abq 2020 or later)
    eval_sig_energies_rect(model_name,b_list[-1],a0,eps_yy,mesh_sizes[0])
    
    remove_files(DIR0+'/'+dir_name)
    os.chdir(DIR0)
    return

# model parameters (N-mm-s)
# ------------------------------------------------------------------

r_list = (5,8,25,70,200,500,1200,3500)
mesh_sizes = (0.02,0.1,0.5,3,8,25,50,100,200) # (0.03,0.1,0.5,3,8,25,50,100,200)

E, nu = 210000., 0.3
k_list = (40.,) # (10.,20.,40.)

# list of (sig_y0, B, n), lin.el. for sig_y = -1
load_cases = ((-1, E/10., 1),
              (300., E/10., 1),
              (300., E/20., 3),
              (300., E/100., 0.7),
              (300., 0, 1.))

load_cases = ((300., E/10., 1.),)

# run PE and PS, all material cases, and all applied Ks
# ------------------------------------------------------------------

sig_y0, B, n = load_cases[0]

#make_crack_model_rect('test-rect-eval-PE-04mm-mesh0015','crack-rect',(3,6,25,100),(1.4,8,20,100),0.4,(E,nu,sig_y0,B,n),
#                      (0.015,0.15,0.8,4,10),0.002,10,1,0)

for is_pe in (1,): # (1,0):
    for sig_y0, B, n in load_cases:
        
        if sig_y0 == -1:
            dir_name = 'closeEel-'+{0:'PS',1:'PE'}[is_pe]
        else:
            dir_name = ('cf4-closeEpl'+str(int(B/E*100)).zfill(3)+'-n'+
                        str(int(n*10)).zfill(2)+'-'+{0:'PS',1:'PE'}[is_pe])

        for K_I in k_list:
            if K_I == 10.:
                n_open = 30
            elif K_I == 20.:
                n_open = 5 # 60
            else:
                n_open = 3 #100
            #
            make_crack_model_circ(dir_name,'crack-k'+str(int(K_I)),r_list,(E,nu,sig_y0,B,n),
                                  mesh_sizes,K_I,n_open,is_pe,is_lin=1)

remove_files(DIR0)