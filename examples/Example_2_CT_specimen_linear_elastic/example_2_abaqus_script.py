"""
This is an Abaqus script for the example in :py:mod:`example_2_CT_specimen_linear_elastic`.

To simulated and write a result.json file open a shell and type:

.. code-block:: console

    abaqus cae noGUI="example_2_abaqus_script.py"

"""
from __future__ import print_function

import sys
import os
import json

import numpy as np

import abaqus as abq

# append path to the folder where conforce_abq is located and import conforce_abq afterward
sys.path.append(os.path.abspath("../.."))
from conforce_abq.main import apply


def simulate_and_save_results(inp_file_path="CT_specimen_CPE4.inp"):
    """
    This function:
    - simulate the input file
    - computes the configurational forces in the odb
    - writes results to a file called "results.json"

    :param inp_file_path: str, path to input file to simulate
    :return: updated Odb with configurational forces as field output
    """

    # create a job from the input file, start the simulation and wait until the simulation completed
    job_name = os.path.basename(inp_file_path).split(".")[0]
    job = abq.mdb.JobFromInputFile(
        name=job_name,
        inputFileName=inp_file_path
    )
    job.submit()
    job.waitForCompletion()

    # open odb with **readOnly=False**
    odb = abq.session.openOdb(job_name + ".odb", readOnly=False)

    # apply conforce plugin and request configurational forces as field output
    odb = apply(
        odb,
        request_CF=True,
        CF_name="CONF_FORCE_LIN_EL",
        method="mbf",
        e_expression="SENER"
    )

    # put all results in this dictionary
    results = dict()

    # extract CF field output
    fo_CF = odb.steps['Loading'].frames[-1].fieldOutputs["CONF_FORCE_LIN_EL"]

    results["CF"] = dict()
    for set_name, region in odb.rootAssembly.nodeSets.items():
        fo_CF_in_region = fo_CF.getSubset(
            region=region
        )

        results["CF"][set_name] = np.sum([
            value.data
            for value in fo_CF_in_region.values
        ], axis=0).tolist()

    # extract history outputs
    (ho_assembly, ho_J_integral, ho_bc_lower, ho_bc_upper) = [
        history_region.historyOutputs
        for history_region in odb.steps['Loading'].historyRegions.values()
    ]

    # extract J-Integrals at the last frame
    results["J"] = {
        key: value.data[-1][-1]
        for key, value in ho_J_integral.items()
    }

    # extract reaction forces
    results["reaction_force"] = ho_bc_upper['RF2'].data[-1][-1]

    # save results
    with open("results.json", "w") as fh:
        json.dump(results, fh, indent=4)

    return odb


if __name__ == '__main__':
    # delegate output to console
    stdout_default, stderr_default = sys.stdout, sys.stderr
    sys.stdout, sys.stderr = sys.__stdout__, sys.__stderr__
    try:
        # simulate
        simulate_and_save_results()

    finally:
        # restore default output streams
        sys.stdout, sys.stderr = stdout_default, stderr_default

