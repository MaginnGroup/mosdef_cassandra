import mbuild
import foyer
import mosdef_cassandra as mc
import unyt as u


def run_nvt_mixture(custom_args=None):
    # If no custom args are passed, assign empty dictionary
    if custom_args is None:
        custom_args = {}

    # Use mbuild to create molecules
    methane = mbuild.load("C", smiles=True)
    propane = mbuild.load("CCC", smiles=True)

    # Create an empty mbuild.Box
    box = mbuild.Box(lengths=[3.0, 3.0, 3.0])

    # Load forcefields
    oplsaa = foyer.forcefields.load_OPLSAA()

    # Use foyer to apply forcefields
    typed_methane = oplsaa.apply(methane)
    typed_propane = oplsaa.apply(propane)

    # Create box and species list
    box_list = [box]
    species_list = [typed_methane, typed_propane]

    # Use Cassandra to insert some initial number of species
    mols_to_add = [[100, 50]]

    system = mc.System(box_list, species_list, mols_to_add=mols_to_add)
    moveset = mc.MoveSet("nvt", species_list)

    mc.run(
        system=system,
        moveset=moveset,
        run_type="equilibration",
        run_length=10000,
        temperature=200.0 * u.K,
        **custom_args,
    )


if __name__ == "__main__":
    run_nvt_mixture()
