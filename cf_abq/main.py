import logging
import sys
from io import StringIO

import abaqus as abq

from cf_abq.field_output_util import add_field_outputs


def main(
        odb_or_odb_name=None,
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
    if odb_or_odb_name is None:
        # use current odb
        odb_or_odb_name = session.currentViewportName
        vp = session.viewports[odb_or_odb_name]
        odb = vp.displayedObject

    elif odb_or_odb_name in session.odbs.keys():
        # lookup odb_name
        odb = session.odbs[odb_or_odb_name]

    else:
        # this has to be an odb
        odb = odb_or_odb_name

    if odb is not None and odb.__class__.__name__ == "Odb":
        method_map = {
            "mbf": "mbf",
            "motion based formulation": "mbf",
            "dbf": "dbf",
            "deformation based formulation": "dbf"
        }
        method = method_map[method]

        odb = add_field_outputs(
            odb=odb,
            request_F=request_F,
            request_P=request_P,
            request_CS=request_CS,
            name_CS=CS_name,
            request_CF=request_CF,
            name_CF=CF_name,
            e_expression=e_expression,
            method=method,
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

    return odb
