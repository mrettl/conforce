import numpy as np
import Conf_Forces_py as cf

#Load Data
data                  = np.load("Data.npz")
Coords                = data["Coords"]
Element_U             = data["Element_U"]
S_vec                 = data["S_vec"]
PENER                 = data["PENER"]
SENER                 = data["SENER"]
Element_Connectivity  = data["Element_Connectivity"]
eval_Node_Labels      = data["eval_Node_Labels"]
t                     = data['t']

#Calculate the dispalement
v_ll=np.round(0.5*t,3)
Node_Labels=np.unique(Element_Connectivity)

########################################################
#Calculate configurational forces without plastic energy
########################################################
CF_Nodal=np.empty((Element_U.shape[0],Node_Labels.shape[0],3))
for i in range(Element_U.shape[0]):
    # Calculate on element nodal position
    CF_Element_Nodal=cf.calc_Conf_Force_CPE4_static(Coords,Element_U[i],S_vec[i],PENER[i]*0.,SENER[i],method='dbf')
    # Get nodal unique value
    Node_labels,CF_Nodal[i]=cf.calc_Nodal(Element_Connectivity,CF_Element_Nodal)

#Select the node labels which are in a given node set
idx=np.isin(Node_labels,eval_Node_Labels,assume_unique=True)
#Sum configurational forces of all selected nodes
J_dbf_far = CF_Nodal[:,idx].sum(axis=1)

#Output the result in x-direction
print("J_far_ep: "+ str(J_dbf_far[-1,0])+ " mJ/mm^2")