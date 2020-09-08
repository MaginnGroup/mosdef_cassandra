import mbuild
import foyer
import mosdef_cassandra as mc
import unyt as u
from mosdef_cassandra.examples.structures import carbon_lattice


def run_gcmc_adsorption(custom_args=None):
    # If no custom args are passed, assign empty dictionary
    if custom_args is None:
        custom_args = {}

    # Use mbuild to create molecules
    lattice = carbon_lattice()
    methane = mbuild.load("C", smiles=True)

    # Load forcefields
    trappe = foyer.forcefields.load_TRAPPE_UA()
    oplsaa = foyer.forcefields.load_OPLSAA()

    # Use foyer to apply forcefields
    typed_lattice = trappe.apply(lattice)
    methane_ff = oplsaa.apply(methane)

    # Create box and species list
    box_list = [lattice]
    species_list = [typed_lattice, methane_ff]

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
