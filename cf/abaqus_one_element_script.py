"""
This module requires Abaqus (including Abaqus cae).
"""
from __future__ import print_function

import sys
import os
import json

import abaqus as abq
import abaqusConstants as abqConst


# global constants
INP_FILE_PATH = os.path.abspath(sys.argv[-1])
WORKING_DIRECTORY = os.path.abspath(os.path.join(INP_FILE_PATH, os.pardir))
HOME_DIRECTORY = os.path.abspath(".")

# reset model database
abq.Mdb()


def main():
    job_name = os.path.basename(INP_FILE_PATH).split(".")[0]
    job = abq.mdb.JobFromInputFile(
        name=job_name,
        inputFileName=INP_FILE_PATH,
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
        }
    }

    with open(job_name + "_result.json", "w") as fh:
        json.dump(result, fh, indent=4)


if __name__ == '__main__':
    os.chdir(WORKING_DIRECTORY)
    try:
        main()

    finally:
        os.chdir(HOME_DIRECTORY)
