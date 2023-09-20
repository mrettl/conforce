import sys
import os
import json

import numpy as np

import abaqus as abq

# append path to the folder where conforce_abq is located and import conforce_abq afterward
sys.path.append(os.path.abspath("../../.."))

from conforce_abq import field_output_util as fou


def compute_CF(odb):
    # CF for whole model
    odb = fou.add_field_outputs(
        odb,
        request_F=False,
        request_P=False,
        request_CS=False,
        name_CF="CF_ALL"
    )

    # CF for odb_inst
    odb_inst = odb.rootAssembly.instances['PART-1-2']
    odb = fou.add_field_outputs(
        odb,
        request_F=False,
        request_P=False,
        request_CS=False,
        name_CF="CF_INSTANCE",
        odb_instances=[odb_inst]
    )

    # CF for element set
    odb_element_set = odb.rootAssembly.elementSets['SET-1']
    odb = fou.add_field_outputs(
        odb,
        request_F=False,
        request_P=False,
        request_CS=False,
        name_CF="CF_ELSET",
        odb_set=odb_element_set
    )

    # CF for node set
    odb_node_set = odb.rootAssembly.nodeSets['SET-1']
    odb = fou.add_field_outputs(
        odb,
        request_F=False,
        request_P=False,
        request_CS=False,
        name_CF="CF_NSET",
        odb_set=odb_node_set
    )

    return odb


def simulate():
    inp_file_path = "TestModel.inp"

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

    return odb


def export_results(odb):
    results = {
        CF_keys: np.sum(
            np.concatenate([
                block.data
                for block
                in odb.steps['Step-1'].frames[1].fieldOutputs[CF_keys].bulkDataBlocks
            ]),
            axis=0
        ).tolist()
        for CF_keys in [
            "CF_ALL",
            "CF_INSTANCE",
            "CF_ELSET",
            "CF_NSET"
        ]
    }

    with open("results.json", "w") as f:
        json.dump(results, f)


if __name__ == '__main__':
    export_results(compute_CF(simulate()))
