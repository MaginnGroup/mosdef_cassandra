import mbuild
import foyer
import mosdef_cassandra as mc


def run_gemc():

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
    moves = mc.Moves("gemc", species_list)

    moves.prob_volume = 0.010
    moves.prob_swap = 0.11

    thermo_props = [
        "energy_total",
        "energy_intervdw",
        "pressure",
        "volume",
        "nmols",
        "mass_density",
    ]

    custom_args = {
        "run_name": "equil",
        "charge_style": "none",
        "rcut_min": 2.0,
        "vdw_cutoff": 14.0,
        "units": "sweeps",
        "steps_per_sweep": 450,
        "coord_freq": 50,
        "prop_freq": 10,
        "properties": thermo_props,
    }

    mc.run(
        system=system,
        moves=moves,
        run_type="equilibration",
        run_length=250,
        temperature=151.0,
        **custom_args,
    )

    # Set max translate and volume for production
    moves.max_translate = [[0.5], [14.0]]
    moves.max_volume = [700.0]

    # Update run_name and restart_name
    custom_args["run_name"] = "prod"
    custom_args["restart_name"] = "equil"

    mc.restart(
        system=system,
        moves=moves,
        run_type="production",
        run_length=750,
        temperature=151.0,
        **custom_args,
    )


if __name__ == "__main__":
    run_gemc()
