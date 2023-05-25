"""
This module contains a function the configurational force computation to an odb.
"""
import logging
import sys
from io import StringIO

import abaqus as abq

from conforce_abq.field_output_util import add_field_outputs


def apply(
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
    """
    Apply the configurational force computation to an odb.
    Additionally, this function prints a log summary after the computation.

    .. seealso::
        The computation is done by :py:func:`conforce_abq.field_output_util.add_field_outputs`


    :param odb_or_odb_name: If None, use current odb.
                            If a string, lookup the name of the odb in the list of opened odb and use this odb.
                            If an Odb, use the given Odb object
    :type odb_or_odb_name: Union[None, str, Odb]
    :param request_F: True, to create a FieldOutput for the deformation gradient
    :type request_F: bool
    :param request_P: True, to create a FieldOutput for the First Piola-Kirchhoff stresses
    :type request_P: bool
    :param request_CS: True, to create a FieldOutput for the configurational stresses
    :type request_CS: bool
    :param request_CF: True, to create a FieldOutput for the configurational forces
    :type request_CF: bool
    :param CS_name: name-prefix of the configurational stress FieldOutput (only considered if `request_CS == True`)
    :type CS_name: str
    :param CF_name: name of the configurational force FieldOutput (only considered if `request_CF == True`)
    :type CF_name: str
    :param method: see :py:func:`conforce_shared.cf_c.compute_CS`
    :type method: str
    :param e_expression: see :py:func:`conforce_abq.field_output_util.eval_field_output_expression`
    :type e_expression: str
    :return: Odb with new requested FieldOutputs
    :rtype: Odb
    """
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
