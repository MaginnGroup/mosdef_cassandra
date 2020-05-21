import mbuild
import foyer
import mosdef_cassandra as mc
import unyt as u


def run_gcmc_restricted(custom_args={}):

    # Use mbuild to create molecules
    methane = mbuild.load("C", smiles=True)

    # Create an empty mbuild.Box
    box = mbuild.Box(lengths=[5.0, 5.0, 5.0])

    # Load forcefields
    oplsaa = foyer.forcefields.load_OPLSAA()

    # Use foyer to apply forcefields
    methane_ff = oplsaa.apply(methane)

    # Create box and species list
    box_list = [box]
    species_list = [methane_ff]

    mols_to_add = [[10]]

    system = mc.System(box_list, species_list, mols_to_add=mols_to_add)
    moveset = mc.MoveSet("gcmc", species_list)

    # Specify restricted insertions
    moveset.add_restricted_insertions(
        species_list, [["sphere"]], [[20 * u.angstrom]]
    )

    default_args = {
        "chemical_potentials": [-35.0 * (u.kJ / u.mol)],
        "prop_freq": 10,
    }

    # Combine default/custom args and override default
    custom_args = {**default_args, **custom_args}

    mc.run(
        system=system,
        moveset=moveset,
        run_type="equilibration",
        run_length=100,
        temperature=300.0 * u.K,
        **custom_args,
    )


if __name__ == "__main__":
    run_gcmc_restricted()
