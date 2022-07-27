import numpy as np
import Conf_Forces_py as cf

# Load Data
data = np.load("Data.npz")
Coords = data["Coords"]
Element_U = data["Element_U"]
S_vec = data["S_vec"]
PENER = data["PENER"]
SENER = data["SENER"]
Element_Connectivity = data["Element_Connectivity"]
eval_Node_Labels = data["eval_Node_Labels"]

Node_Labels = np.unique(Element_Connectivity)

# Calculate Configurational Forces
CF_Nodal = np.empty((Element_U.shape[0], Node_Labels.shape[0], 3))
for i in range(Element_U.shape[0]):
    # Calculate on Element Nodal Position
    CF_Element_Nodal = cf.calc_Conf_Force_CPE4_static(Coords, Element_U[i], S_vec[i], PENER[i], SENER[i], method='dbf')
    # Get Nodal unique value
    Node_labels, CF_Nodal[i] = cf.calc_Nodal(Element_Connectivity, CF_Element_Nodal)

np.savez('Configurational_Forces.npz', Node_labels=Node_labels, CF_Nodal=CF_Nodal)

# Select the node labels which are in a given node set
idx = np.isin(Node_labels, eval_Node_Labels, assume_unique=True)
# Sum configurational forces of all selected nodes
J_dbf_3 = CF_Nodal[:, idx].sum(axis=1)

# Output the result in x-direction
print("J-integral over 3 Contours: " + str(J_dbf_3[-1, 0] * -1) + " mJ/mm^2")
