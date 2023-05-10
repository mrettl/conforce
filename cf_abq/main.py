import logging
import sys
from io import StringIO

import abaqus as abq

from field_output_util import add_field_outputs


def main(
        odb_name=None,
        request_F=False,
        request_P=False,
        request_CS=False,
        request_CF=False,
        method="mbf",
        e_expression="SENER+PENER"
):
    log = StringIO()
    handler = logging.StreamHandler(log)
    handler.setFormatter(logging.Formatter(fmt=u'%(levelname)s: %(message)s'))

    logger = logging.Logger("cf_plugin")
    logger.addHandler(handler)
    logger.addHandler(logging.StreamHandler(sys.stdout))

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
            logger=logger
        )

    else:
        logger.error("No odb found. Open a odb in the current viewport.")

    # Summary
    print("\n\n")
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

    print("\n")
