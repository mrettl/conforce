"""
This module contains methods to simulate one element models.
"""
import io
import os
import subprocess
import json

import numpy as np


def generate_abaqus_input(element_type, X_at_nodes, U_at_nodes, E: float = 1000., nu: float = 0.3):
    r"""
    Generate a model containing one element of the given element type
    with the given nodal coordinates and nodal displacements.
    The element of the given element type has d=2 or d=3 dimensions and n nodes.

    :param element_type: str, name of the element type
    :param X_at_nodes: array of shape (n, d) containing nodal coordinates
    :param U_at_nodes: array of shape (n, d) containing nodal displacements
    :param E: Young's modulus
    :param nu: Poisson's ratio
    :return: str, valid abaqus input format
    """

    new_line = "\n"

    # node block
    buffer = io.StringIO()
    X_at_nodes = np.asarray(X_at_nodes, dtype=float)
    rows, columns = X_at_nodes.shape

    node_ids = np.arange(1, rows+1).reshape((-1, 1))
    node_set_names = [f"NODE_SET_{node_id}" for node_id in node_ids.flat]

    node_data = np.column_stack([node_ids.astype(float), X_at_nodes])
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
{E}, {nu}
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


def simulate_one_element(X_at_nodes, U_at_nodes, element_type: str, load_name: str, folder: str,
                         E: float = 1000., nu: float = 0.3) -> dict:
    """
    Generate and simulate a one element model.
    The input is generated using the method :py:func`generate_abaqus_input`.
    For the simulation the script `one_element_script.py` is executed in Abaqus.
    The model is simulated in the given folder.
    After the simulation has completed successfully,
    a result file in the json format is writen into the folder.

    This result file contains:

        - model data like the strain energy (ALLSE),
        - element data like the element volume (EVOL),
        - nodal data like the nodal coordinates (COORD),
        - integration point data like the integration point coordinates (COORD),

    .. note::

        If the result file already exists, no simulation is performed.
        Instead, the content of the result file is returned.

    **Examples**

    The result contains the strain energy (ALLSE),
    reaction forces (RF), strain energy densities (SENER) and
    the element volume (EVOL)

    >>> result = simulate_one_element(
    ...     X_at_nodes=np.array([[0, 0], [1, 0], [0.5, 1]]),
    ...     U_at_nodes=np.array([[0, 0], [0, 0], [0, 1]]),
    ...     element_type="CPE3",
    ...     load_name="triangle_tension",
    ...     folder="res/tests/triangle_tension"
    ... )
    >>> np.around(result["model"]["ALLSE"], 3)
    336.538
    >>> np.around(result["nodes"]["RF"], 3)
    array([[-288.462, -336.538],
           [ 288.462, -336.538],
           [  -0.   ,  673.077]])
    >>> np.around(result["integration_points"]["SENER"], 3)
    array([[673.077]])
    >>> np.around(result["element"]["EVOL"], 3)
    0.5

    :param X_at_nodes: array of shape (n, d) containing nodal coordinates
    :param U_at_nodes: array of shape (n, d) containing nodal displacements
    :param element_type: str, name of the element type
    :param load_name: str, a valid load name in abaqus
    :param folder: str, simulation folder
    :param E: Young's modulus
    :param nu: Poisson's ratio
    :return: dict, containing the data of the result file
    """
    job_name = f"{element_type}_{load_name}"
    folder = os.path.abspath(folder)
    file = f"{folder}/{job_name}.inp"
    if not os.path.exists(folder):
        os.makedirs(folder)

    result_file_path = f"{folder}/{job_name}_result.json"
    if not os.path.exists(result_file_path):
        with open(file, "w", encoding="utf-8") as f:
            f.write(generate_abaqus_input(
                element_type=element_type,
                X_at_nodes=X_at_nodes,
                U_at_nodes=U_at_nodes,
                E=E,
                nu=nu
            ))

        subprocess.call([
            "abaqus", "cae", r"noGUI=conforce_abq\one_element_script.py", "--", file
        ], shell=True)

        for file in {"abaqus.rpy"} | {
            f"{folder}/{job_name}.{extension}"
            for extension
            in ["com", "dat", "inp", "log", "msg", "odb", "prt", "sim", "sta"]
        }:
            if os.path.exists(file):
                os.remove(file)

    with open(result_file_path, "r", encoding="utf-8") as fh:
        return json.load(fh)
