import logging
import sys
from io import StringIO

import abaqus as abq

from field_output_util import add_field_outputs


def main():
    log = StringIO()
    handler = logging.StreamHandler(log)
    handler.setFormatter(logging.Formatter(fmt=u'%(message)s'))

    logger = logging.Logger("cf_plugin")
    logger.addHandler(handler)
    logger.addHandler(logging.StreamHandler(sys.stdout))

    session = abq.session
    vp = session.viewports[abq.session.currentViewportName]
    odb = vp.displayedObject

    if odb is not None:
        add_field_outputs(
            odb=odb,
            logger=logger
        )

    else:
        logger.error("No odb found. Open a odb in the current viewport.")

    log.seek(0)
    print(log.read())
