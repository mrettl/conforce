from odbAccess import *
import numpy as np
import sys

def get_section_density(odb,instance):
    """
    Reads the section density from odb
    Output: 
          dictionary[section_name:density]
          dictionary[section_name:region_object]
    """
    mat={}
    for matName in odb.materials.keys():
        mat[matName]=odb.materials[matName].density.table[0][0]
    
    sec_density={}
    sec_region ={}
    for secName in odb.sections.keys():
        sec_density[secName]=mat[odb.sections[secName].material]
        for sec in instance.sectionAssignments:
            if sec.sectionName == secName:
                sec_region[secName]=sec.region
    return sec_density,sec_region

def get_Mesh_Data(instance):
    """
    Get Data from Mesh
    Input:
        instance object
    Output:
        NodeLables            array(numNodes)
        NodeCoordinates       array(numNodes,3)
        Element_Labels        dict['elemtype'][array(num_elements)]
        Element_Connectivity  dict['elemtype'][array(num_elements,numNodesPerElement)]
    """
    NodeObj = instance.nodes
    NodeLables=np.array([i.label for i in NodeObj],np.int64)
    NodeCoordinates=np.array([i.coordinates for i in NodeObj],np.float64)
    idx=np.argsort(NodeLables)
    NodeLables=NodeLables[idx]
    NodeCoordinates=NodeCoordinates[idx,:]
    
    Element_Labels={}
    Element_Connectivity={}
    Element_types=set()
    ElementObj=instance.elements
    for elem in ElementObj:
        type = elem.type
        if type not in Element_types:
            Element_Labels[type]=[elem.label,]
            Element_Connectivity[type]=[elem.connectivity,]
            Element_types.add(type)
        else:
            Element_Labels[type].append(elem.label)
            Element_Connectivity[type].append(elem.connectivity)
    for key in Element_Labels.keys():
        Labels=np.array(Element_Labels[key],dtype=np.int64)
        Connectivity=np.array(Element_Connectivity[key],dtype=np.int64)
        idx=np.argsort(Labels)
        Element_Labels[key]=Labels[idx]
        Element_Connectivity[key]=Connectivity[idx,:]
    return NodeLables,NodeCoordinates,Element_Labels,Element_Connectivity

def read_frame_time(odb):
    """
    Get the step time from ODB
    Output:
        Time array(numOfSteps*numFrames)
    """
    time=[]
    Steps=odb.steps.keys()
    for i in range(len(Steps)):
        Step=odb.steps[Steps[i]]
        for j in range(len(Step.frames)):
            time.append(Step.frames[j].frameValue)
    S=np.array(time)
    return time


def Nodal_to_ElementNodal(U_label,U_data,Element_Connectivity_inv,U):
    for elemtype in Element_Connectivity_inv.keys():
        valid_Labels=Element_Connectivity_inv[elemtype]['NodeLabels']
        inverse=Element_Connectivity_inv[elemtype]['inverse']
        idx=np.isin(U_label,valid_Labels,assume_unique=True)
        valid_data=U_data[idx]
        U[elemtype].append(valid_data[inverse])
    return U


def get_Nodal_field(U_Field,Element_Connectivity_inv,U):
    #Sorted Unique Nodal values
    U_data=np.copy(U_Field.bulkDataBlocks[0].data)
    U_label=np.copy(U_Field.bulkDataBlocks[0].nodeLabels)
    for ii in range(1,len(U_Field.bulkDataBlocks)):
        U_data=np.concatenate((U_data,U_Field.bulkDataBlocks[ii].data))
        U_label=np.concatenate((U_label,U_Field.bulkDataBlocks[ii].nodeLabels))
    idx=np.argsort(U_label)
    U_data=U_data[idx,:]
    U_label=U_label[idx]
    
    U=Nodal_to_ElementNodal(U_label,U_data,Element_Connectivity_inv,U)
    
    return U


def get_Int_Point_field(S_Field,Element_Connectivity,S):
    for elemtype in Element_Connectivity.keys():
        num_Elements=Element_Connectivity[elemtype].shape[0]
        
        S_Field_Elem=S_Field.getSubset(elementType=elemtype)
        S_data=np.copy(S_Field_Elem.bulkDataBlocks[0].data)
        S_label=np.copy(S_Field_Elem.bulkDataBlocks[0].elementLabels)
        S_int_points=np.copy(S_Field_Elem.bulkDataBlocks[0].integrationPoints)
        for ii in range(1,len(S_Field_Elem.bulkDataBlocks)):
            S_data=np.concatenate((S_data,S_Field_Elem.bulkDataBlocks[ii].data))
            S_label=np.concatenate((S_label,S_Field_Elem.bulkDataBlocks[ii].elementLabels))
            S_int_points=np.concatenate((S_int_points,S_Field_Elem.bulkDataBlocks[ii].integrationPoints))
        
        S_label_sort=S_label*100+S_int_points
        idx=np.argsort(S_label_sort)
        S_data=S_data[idx]
        
        if S_Field_Elem.bulkDataBlocks[0].position == INTEGRATION_POINT:
            S_data=S_data.reshape(num_Elements,-1,S_data.shape[-1])
        
        S[elemtype].append(S_data)
    return S


def read_static_data(odb,instance,Element_Connectivity,NodeLables,NodeCoordinates):
    """
    Get field Output at nodes
    Output:
        array of shape (numSteps*num_Frames,num_Elements,numNodes,numValues)
    """
    #Inititiation
    U={}
    S={}
    SENER={}
    PENER={}
    Coord={}
    for elemtype in Element_Connectivity.keys():
        S[elemtype]=[]
        SENER[elemtype]=[]
        PENER[elemtype]=[]
        U[elemtype]=[]
        Coord[elemtype]=[]
    
    #inverse connectivity
    Element_Connectivity_inv={}
    for elemtype in Element_Connectivity.keys():
        NodeLabels_elem,unique_inverse=np.unique(Element_Connectivity[elemtype],return_inverse=True)
        unique_inverse=unique_inverse.reshape(Element_Connectivity[elemtype].shape)
        Element_Connectivity_inv[elemtype]={}
        Element_Connectivity_inv[elemtype]['NodeLabels']=NodeLabels_elem
        Element_Connectivity_inv[elemtype]['inverse']=unique_inverse
    
    Coord=Nodal_to_ElementNodal(NodeLables,NodeCoordinates,Element_Connectivity_inv,Coord)
    
    #Read field Outputs
    Steps=odb.steps.keys()
    for i in range(len(Steps)):
        Step=odb.steps[Steps[i]]
        for j in range(len(Step.frames)):
            U_Field=Step.frames[j].fieldOutputs['U'].getSubset(region=instance)
            U=get_Nodal_field(U_Field,Element_Connectivity_inv,U)
            
            S_Field=Step.frames[j].fieldOutputs['S'].getSubset(region=instance)
            S=get_Int_Point_field(S_Field,Element_Connectivity,S)
            
            SENER_Field=Step.frames[j].fieldOutputs['SENER'].getSubset(region=instance)
            SENER=get_Int_Point_field(SENER_Field,Element_Connectivity,SENER)
            
            PENER_Field=Step.frames[j].fieldOutputs['PENER'].getSubset(region=instance)
            PENER=get_Int_Point_field(PENER_Field,Element_Connectivity,PENER)
    
    for elemtype in Element_Connectivity.keys():
        U[elemtype]=np.array(U[elemtype],dtype=np.float64)
        S[elemtype]=np.array(S[elemtype],dtype=np.float64)
        SENER[elemtype]=np.array(SENER[elemtype],dtype=np.float64)
        SENER[elemtype]=SENER[elemtype].reshape(SENER[elemtype].shape[:-1])
        PENER[elemtype]=np.array(PENER[elemtype],dtype=np.float64)
        PENER[elemtype]=PENER[elemtype].reshape(PENER[elemtype].shape[:-1])
        Coord[elemtype]=np.squeeze(np.array(Coord[elemtype],dtype=np.float64))
    return U,S,SENER,PENER,Coord

###########################################
###########################################
odbPath = "Two_phase_bar.odb"
odb=openOdb(odbPath,readOnly=False)
instance=odb.rootAssembly.instances['BAR']

NodeLables,NodeCoordinates,Element_Labels,Element_Connectivity=get_Mesh_Data(instance)
U,S,SENER,PENER,Coords=read_static_data(odb,instance,Element_Connectivity,NodeLables,NodeCoordinates)

##########################################
#Extract Data from ODB
##########################################
U_data=U['CPE4R']
S_data=S['CPE4R']
SENER=SENER['CPE4R']
PENER=PENER['CPE4R']
Coords=Coords['CPE4R']
Element_Connectivity=Element_Connectivity['CPE4R']

#To 3D
Element_U=np.zeros(U_data.shape[:-1]+(3,))
Element_U[:,:,:,0:2]=U_data
S_vec=np.zeros(S_data.shape[:-1]+(6,)) 
S_vec[:,:,:,0:4]=S_data

odb.close()
#########################################
#Save Data
#########################################
np.savez("Data.npz",Coords=Coords,Element_Connectivity=Element_Connectivity,Element_U=Element_U,S_vec=S_vec,PENER=PENER,SENER=SENER)


