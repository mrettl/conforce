import sys
import os

import abaqus as abq

# append path to the folder where conforce_abq is located and import conforce_abq afterward
sys.path.append(os.path.abspath("../.."))

from conforce_abq import field_output_util as fou


def compute_CF(odb):
    # CF for whole model
    odb = fou.add_field_outputs(
        odb,
        request_F=False,
        request_P=False,
        request_CS=False,
        name_CF="CF_a"
    )

    # CF for odb_inst
    odb_inst = odb.rootAssembly.instances['PART-1-2']
    odb = fou.add_field_outputs(
        odb,
        request_F=False,
        request_P=False,
        request_CS=False,
        name_CF="CF_b",
        odb_instances=[odb_inst]
    )

    # CF for element set
    odb_element_set = odb.rootAssembly.elementSets['SET-1']
    odb = fou.add_field_outputs(
        odb,
        request_F=False,
        request_P=False,
        request_CS=False,
        name_CF="CF_c",
        odb_set=odb_element_set
    )

    # CF for node set
    odb_node_set = odb.rootAssembly.nodeSets['SET-1']
    odb = fou.add_field_outputs(
        odb,
        request_F=False,
        request_P=False,
        request_CS=False,
        name_CF="CF_d",
        odb_set=odb_node_set
    )

    return odb


def simulate():
    inp_file_path = "Job-1.inp"

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


if __name__ == '__main__':
    compute_CF(simulate())
