import mbuild
import foyer
import mosdef_cassandra as mc
import unyt as u


def run_gemc(**custom_args):

    # Use mbuild to create molecules
    methane = mbuild.Compound(name="_CH4")

    # Create two empty mbuild.Box
    # (vapor = larger, liquid = smaller)
    liquid_box = mbuild.Box(lengths=[3.0, 3.0, 3.0])
    vapor_box = mbuild.Box(lengths=[4.0, 4.0, 4.0])

    # Load forcefields
    trappe = foyer.forcefields.load_TRAPPE_UA()

    # Use foyer to apply forcefields
    typed_methane = trappe.apply(methane)

    # Create box and species list
    box_list = [liquid_box, vapor_box]
    species_list = [typed_methane]

    mols_to_add = [[350], [100]]

    system = mc.System(box_list, species_list, mols_to_add=mols_to_add)
    moveset = mc.MoveSet("gemc", species_list)

    moveset.prob_volume = 0.010
    moveset.prob_swap = 0.11

    thermo_props = [
        "energy_total",
        "energy_intervdw",
        "pressure",
        "volume",
        "nmols",
        "mass_density",
    ]

    default_args = {
        "run_name": "equil",
        "charge_style": "none",
        "rcut_min": 2.0 * u.angstrom,
        "vdw_cutoff": 14.0 * u.angstrom,
        "units": "sweeps",
        "steps_per_sweep": 450,
        "coord_freq": 50,
        "prop_freq": 10,
        "properties": thermo_props,
    }

    # Combine default/custom args and override default
    custom_args = {**default_args, **custom_args}

    mc.run(
        system=system,
        moveset=moveset,
        run_type="equilibration",
        run_length=250,
        temperature=151.0 * u.K,
        **custom_args,
    )

    mc.restart(
        restart_from="equil",
        run_name="prod",
        run_type="production",
        total_run_length=750,
    )


if __name__ == "__main__":
    run_gemc()
