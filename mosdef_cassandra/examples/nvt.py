import mbuild
import foyer
import mosdef_cassandra as mc


def run_nvt(custom_args={}):

    # Use mbuild to create molecules
    methane = mbuild.load("C", smiles=True)

    # Create an empty mbuild.Box
    box = mbuild.Box(lengths=[3.0, 3.0, 3.0])

    # Load forcefields
    oplsaa = foyer.forcefields.load_OPLSAA()

    # Use foyer to apply forcefields
    methane_ff = oplsaa.apply(methane)

    # Create box and species list
    box_list = [box]
    species_list = [methane_ff]

    # Use Cassandra to insert some initial number of species
    mols_to_add = [[50]]

    # Define the system object
    system = mc.System(box_list, species_list, mols_to_add=mols_to_add)
    # Get the move probabilities
    moves = mc.Moves("nvt", species_list)

    # Run a simulation with at 300 K with 10000 MC moves
    mc.run(
        system=system,
        moves=moves,
        run_type="equilibration",
        run_length=10000,
        temperature=300.0,
        **custom_args,
    )


if __name__ == "__main__":
    run_nvt()
