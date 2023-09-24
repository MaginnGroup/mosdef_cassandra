import mbuild
import foyer
import mosdef_cassandra as mc
import unyt as u


def run_nvt(**custom_args):

    # Use mBuild to create a methane molecule
    methane = mbuild.load("C", smiles=True)

    # Create an empty mbuild.Box
    box = mbuild.Box(lengths=[3.0, 3.0, 3.0])

    # Load force field
    oplsaa = foyer.forcefields.load_OPLSAA()

    # Use foyer to apply force field to methane
    methane_ff = oplsaa.apply(methane)

    # Create box and species list
    box_list = [box]
    species_list = [methane_ff]

    # Use Cassandra to insert some initial number of methane molecules
    mols_to_add = [[50]]

    # Define the System
    system = mc.System(box_list, species_list, mols_to_add=mols_to_add)
    # Define the MoveSet
    moveset = mc.MoveSet("nvt", species_list)

    # Run a simulation at 300 K for 10000 MC moves
    mc.run(
        system=system,
        moveset=moveset,
        run_type="equilibration",
        run_name="nvt_parmed",
        run_length=10000,
        temperature=300.0 * u.K,
        **custom_args,
    )


if __name__ == "__main__":
    run_nvt()
