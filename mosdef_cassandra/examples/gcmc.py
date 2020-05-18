import mbuild
import foyer
import mosdef_cassandra as mc


def run_gcmc(custom_args={}):

    # Use mbuild to create molecules
    methane = mbuild.load("C", smiles=True)

    # Create an empty mbuild.Box
    box = mbuild.Box(lengths=[10.0, 10.0, 10.0])

    # Load forcefields
    oplsaa = foyer.forcefields.load_OPLSAA()

    # Use foyer to apply forcefields
    methane_ff = oplsaa.apply(methane)

    # Create box and species list
    box_list = [box]
    species_list = [methane_ff]

    mols_to_add = [[100]]

    system = mc.System(box_list, species_list, mols_to_add=mols_to_add)
    moveset = mc.MoveSet("gcmc", species_list)

    default_args = {
        "chemical_potentials": [-35.0],
        "prop_freq": 100,
    }

    # Combine default/custom args and override default
    custom_args = {**default_args, **custom_args}

    mc.run(
        system=system,
        moveset=moveset,
        run_type="equilibration",
        run_length=1000,
        temperature=300.0,
        **custom_args,
    )


if __name__ == "__main__":
    run_gcmc()
