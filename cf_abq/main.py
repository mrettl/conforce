import abaqus as abq

from field_output_util import add_field_outputs


def main():
    for odb in abq.session.odbs.values():
        add_field_outputs(odb=odb)
