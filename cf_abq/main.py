import logging
import sys
from io import StringIO

import abaqus as abq

from cf_abq.field_output_util import add_field_outputs


def main(
        odb_name=None,
        request_F=False,
        request_P=False,
        request_CS=False,
        request_CF=False,
        CS_name="CONF_STRESS",
        CF_name="CONF_FORCE",
        method="mbf",
        e_expression="SENER+PENER"
):
    log = StringIO()
    handler_log = logging.StreamHandler(log)
    handler_log.setFormatter(logging.Formatter(fmt=u'%(levelname)s: %(message)s'))

    handler_stdout = logging.StreamHandler(sys.stdout)
    handler_stdout.setFormatter(logging.Formatter(fmt=u'%(levelname)s: %(message)s'))

    logger = logging.Logger("cf_plugin")
    logger.addHandler(handler_log)
    logger.addHandler(handler_stdout)

    session = abq.session
    if odb_name is None:
        odb_name = session.currentViewportName
        vp = session.viewports[odb_name]
        odb = vp.displayedObject

    else:
        odb = session.odbs[odb_name]

    fields = [
        field
        for field, request in zip(
            ["F", "P", "CS", "CF"],
            [request_F, request_P, request_CS, request_CF]
        )
        if bool(request)
    ]

    method_map = {
        "mbf": "mbf",
        "motion based formulation": "mbf",
        "dbf": "dbf",
        "deformation based formulation": "dbf"
    }
    method = method_map[method]

    if odb is not None:
        add_field_outputs(
            odb=odb,
            fields=fields,
            method=method,
            e_expression=e_expression,
            name_CS=CS_name,
            name_CF=CF_name,
            logger=logger
        )

    else:
        logger.error("No odb found. Open a odb in the current viewport.")

    # Summary
    print(" ")
    print("CF - log summary")
    print("----------------")

    log.seek(0)
    messages = log.read()
    level_msgs = dict()
    for line in messages.splitlines():
        level = line.split(":")[0]
        level_msgs.setdefault(level, list()).append(line)

    for level, msgs in level_msgs.items():
        print("{0:3d} x {1:s}".format(len(msgs), level))

    print("----------------\n ")
