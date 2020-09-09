import os
from pkg_resources import resource_filename


def get_example_ff_path(ff_name):
    """Get the path to a force field xml file in the examples
    directory"""
    ff_path = resource_filename(
        "mosdef_cassandra", os.path.join("examples/ff_files", ff_name + ".xml")
    )

    return ff_path


def get_example_mol2_path(molecule_name):
    """Get the path to a mol2 file in the examples
    directory"""
    mol2_path = resource_filename(
        "mosdef_cassandra",
        os.path.join("examples/mol2_files", molecule_name + ".mol2"),
    )

    return mol2_path
