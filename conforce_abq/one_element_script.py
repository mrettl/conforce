"""
This Abaqus script simulates an input file and writes results to a json file.

Call this script using :py:func:`conforce_3.one_element_runner.simulate_one_element`.
Or call this script directly using a shell with:

    ``abaqus cae noGUI="{path to this script file}" -- {path to the abaqus input file}``

After the successful simulation a json file named "{inp_name}_result.json" is placed
in the same folder as the Abaqus input file.
"""
import sys
import os
import json

import abaqus as abq
import abaqusConstants as abqConst


def simulate(inp_file_path):
    """
    Simulate the given Abaqus input
    and write a json file named "{inp_name}_result.json" containing the results of the simulation.
    Ths json file is placed into the current working directory.

    :param inp_file_path: path to an Abaqus input file
    """
    job_name = os.path.basename(inp_file_path).split(".")[0]
    job = abq.mdb.JobFromInputFile(
        name=job_name,
        inputFileName=inp_file_path,
        nodalOutputPrecision=abqConst.FULL
    )
    job.submit()
    job.waitForCompletion()
    odb = abq.session.openOdb(job_name + ".odb", readOnly=False)

    step = odb.steps.values()[0]
    ho = step.historyRegions.values()[0]
    fo = step.frames[-1].fieldOutputs

    result = {
        "model": {
            "ALLSE": float(ho.historyOutputs['ALLSE'].data[-1][-1])
        },
        "element": {
            "ESEDEN": float(fo["ESEDEN"].bulkDataBlocks[0].data[0][0]),
            "EVOL": float(fo["EVOL"].bulkDataBlocks[0].data[0][0])
        },
        "nodes": {
            key: fo[key].getSubset(position=abqConst.NODAL).bulkDataBlocks[0].data.tolist()
            for key in ["COORD", "RF", "U"]
        },
        "integration_points": {
            key: fo[key].getSubset(position=abqConst.INTEGRATION_POINT).bulkDataBlocks[0].data.tolist()
            for key in ["COORD", "S", "E", "SENER", "IVOL"]
        },
        "faces": {
            str(face.faces[0][0]): [
                int(node.label)
                for node in face.nodes[0]
            ]
            for face in odb.rootAssembly.surfaces.values()
        }
    }

    with open(job_name + "_result.json", "w") as fh:
        json.dump(result, fh, indent=4)


if __name__ == '__main__':
    # global constants
    INP_FILE_PATH = os.path.abspath(sys.argv[-1])
    WORKING_DIRECTORY = os.path.abspath(os.path.join(INP_FILE_PATH, os.pardir))
    HOME_DIRECTORY = os.path.abspath(".")

    # reset model database
    abq.Mdb()

    os.chdir(WORKING_DIRECTORY)
    try:
        simulate(INP_FILE_PATH)

    finally:
        os.chdir(HOME_DIRECTORY)

