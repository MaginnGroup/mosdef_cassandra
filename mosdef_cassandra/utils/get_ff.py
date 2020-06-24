import os
from pkg_resources import resource_filename


def get_ff_path(ff_name):
    """Get the path to a force field xml file in the examples
    directory"""
    ff_path = resource_filename(
        "mosdef_cassandra", os.path.join("examples/ff_files", ff_name + ".xml")
    )

    return ff_path
