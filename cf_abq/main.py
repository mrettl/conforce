import abaqus as abq

from field_output_util import add_field_outputs


def main():
    session = abq.session
    vp = session.viewports[abq.session.currentViewportName]
    odb = vp.displayedObject
    add_field_outputs(odb=odb)

