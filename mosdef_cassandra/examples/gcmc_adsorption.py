import mbuild
import foyer
import mosdef_cassandra as mc
import unyt as u


from mosdef_cassandra.utils.get_files import (
    get_example_ff_path,
    get_example_cif_path,
)


def run_gcmc_adsorption(**custom_args):

    # Use mbuild to create a zeolite supercell from CIF
    lattice = mbuild.lattice.load_cif(get_example_cif_path("TON"))
    compound_dict = {
        "Si": mbuild.Compound(name="Si"),
        "O": mbuild.Compound(name="O"),
    }
    ton = lattice.populate(compound_dict, 3, 3, 6)

    # Create a coarse-grained methane
    methane = mbuild.Compound(name="_CH4")

    # Load forcefields
    trappe_zeo = foyer.Forcefield(get_example_ff_path("trappe_zeo"))
    trappe = foyer.forcefields.load_TRAPPE_UA()

    # Use foyer to apply forcefields
    ton_ff = trappe_zeo.apply(ton)
    methane_ff = trappe.apply(methane)

    # Create box and species list
    box_list = [ton]
    species_list = [ton_ff, methane_ff]

    # Since we have an occupied box we need to specify
    # the number of each species present in the intial config
    mols_in_boxes = [[1, 0]]

    system = mc.System(box_list, species_list, mols_in_boxes=mols_in_boxes)
    moveset = mc.MoveSet("gcmc", species_list)

    default_args = {
        "chemical_potentials": ["none", -30.0 * (u.kJ / u.mol)],
        "rcut_min": 0.5 * u.angstrom,
        "vdw_cutoff": 14.0 * u.angstrom,
        "charge_cutoff": 14.0 * u.angstrom,
        "coord_freq": 100,
        "prop_freq": 10,
    }

    # Combine default/custom args and override default
    custom_args = {**default_args, **custom_args}

    mc.run(
        system=system,
        moveset=moveset,
        run_type="equilibration",
        run_length=10000,
        temperature=300.0 * u.K,
        **custom_args,
    )


if __name__ == "__main__":
    run_gcmc_adsorption()
