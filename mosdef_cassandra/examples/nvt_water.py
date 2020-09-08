import mbuild
import foyer
import numpy as np
import mosdef_cassandra as mc
import unyt as u
from scipy.spatial.transform import Rotation
from mosdef_cassandra.utils.get_files import get_ff_path, get_mol2_path


def run_water_nvt(custom_args={}):
    molecule = mbuild.load(get_mol2_path("spce"))

    # Create an empty mbuild.Box
    box = mbuild.Box(lengths=[3.0, 3.0, 3.0])

    # Load forcefields
    spce = foyer.Forcefield(get_ff_path("spce"))

    # Use foyer to apply forcefields
    molecule_ff = spce.apply(molecule)

    # Create box and species list
    box_list = [box]
    species_list = [molecule_ff]

    # Use Cassandra to insert some initial number of species
    mols_to_add = [[50]]

    # Define the system object
    system = mc.System(box_list, species_list, mols_to_add=mols_to_add)
    # Get the move probabilities
    moveset = mc.MoveSet("nvt", species_list)

    # Run a simulation with at 300 K with 10000 MC moveset
    mc.run(
        system=system,
        moveset=moveset,
        run_type="equilibration",
        run_length=10000,
        temperature=300.0 * u.K,
        **custom_args,
    )


if __name__ == "__main__":
    run_water_nvt()
