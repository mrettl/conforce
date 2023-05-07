import io
import os
import subprocess
import json

import numpy as np


def generate_abaqus_input(element_type, R_at_nodes, U_at_nodes):
    new_line = "\n"

    # node block
    buffer = io.StringIO()
    R_at_nodes = np.asarray(R_at_nodes, dtype=float)
    rows, columns = R_at_nodes.shape

    node_ids = np.arange(1, rows+1).reshape((-1, 1))
    node_set_names = [f"NODE_SET_{node_id}" for node_id in node_ids.flat]

    node_data = np.column_stack([node_ids.astype(float), R_at_nodes])
    np.savetxt(buffer, node_data, fmt=["%5d"] + ["%10.16g"]*columns, delimiter=",")
    buffer.seek(0)

    return f"""*Heading
** Generated abaqus input file
*Preprint, echo=NO, model=NO, history=NO, contact=NO
**
** PART
**
*Part, name=PART-1
*Node
{buffer.read()}**
*Element, type={element_type}
1, {', '.join([str(node_id) for node_id in node_ids.flat])}
*Elset, elset=ELEMENT_SET
1,
*Solid Section, elset=ELEMENT_SET, material=MATERIAL-1
,
*End Part
**
**  ASSEMBLY
**
*Assembly, name=Assembly
**
*Instance, name=PART-1-1, part=PART-1
*End Instance
**
{new_line.join([
        f'*Nset, nset={set_name}, instance=PART-1-1' + new_line + 
        f'{node_id},'
        for set_name, node_id in zip(node_set_names, node_ids.flat)
])}
*End Assembly
**
** MATERIAL
**
*Material, name=MATERIAL-1
*Elastic
1000., 0.3
** 
** STEP
**
*Step, name=STEP-1, nlgeom=NO
*Static
1., 1., 1e-05, 1.
**
** BOUNDARY CONDITIONS
**
{
    new_line.join([
        new_line.join([
            f"*Boundary{new_line}"
            f"{set_name}, {axis}, {axis}, {U}"
            for U, axis in zip(U_at_node, [1, 2, 3])
        ])
        for set_name, U_at_node in zip(node_set_names, U_at_nodes)
    ])}
** 
** OUTPUT REQUESTS
** 
** FIELD OUTPUT
** 
*Output, field
*Node Output
COORD, RF, U
*Element Output, directions=YES
COORD, S, E, SENER, ESEDEN, EVOL, IVOL
** 
** HISTORY OUTPUT
** 
*Output, history, variable=PRESELECT
*End Step
"""


def simulate_one_element(R_at_nodes, U_at_nodes, element_type: str, load_name: str, folder: str):
    job_name = f"{element_type}_{load_name}"
    folder = os.path.abspath(folder)
    file = f"{folder}/{job_name}.inp"
    if not os.path.exists(folder):
        os.makedirs(folder)

    result_file_path = f"{folder}/{job_name}_result.json"
    if not os.path.exists(result_file_path):
        with open(file, "w") as f:
            f.write(generate_abaqus_input(
                element_type=element_type,
                R_at_nodes=R_at_nodes,
                U_at_nodes=U_at_nodes,
            ))

        subprocess.call([
            "abaqus", "cae", r"noGUI=cf\abaqus_one_element_script.py", "--", file
        ], shell=True)

        for file in {"abaqus.rpy"} | {
            f"{folder}/{job_name}.{extension}"
            for extension
            in ["com", "dat", "inp", "log", "msg", "odb", "prt", "sim", "sta"]
        }:
            if os.path.exists(file):
                os.remove(file)

    with open(result_file_path, "r") as fh:
        return json.load(fh)
