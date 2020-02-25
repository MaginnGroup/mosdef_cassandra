import mbuild
import foyer
import mosdef_cassandra as mc
from mosdef_cassandra.examples.structures import carbon_lattice


def run_gcmc_adsorption():

    # Use mbuild to create molecules
    lattice = carbon_lattice()
    methane = mbuild.load("C", smiles=True)

    # Load forcefields
    trappe = foyer.forcefields.load_TRAPPE_UA()
    oplsaa = foyer.forcefields.load_OPLSAA()

    # Use foyer to apply forcefields
    typed_lattice = trappe.apply(lattice)
    typed_methane = oplsaa.apply(methane)

    # Create box and species list
    box_list = [lattice]
    species_list = [typed_lattice, typed_methane]

    # Since we have an occupied box we need to specify
    # the number of each species present in the intial config
    mols_in_boxes = [[1, 0]]

    system = mc.System(box_list, species_list, mols_in_boxes=mols_in_boxes)
    moves = mc.Moves("gcmc", species_list)

    custom_args = {
        "chemical_potentials": ["none", -30.0],
        "rcut_min": 0.5,
        "vdw_cutoff": 14.0,
        "charge_cutoff": 14.0,
        "coord_freq": 100,
        "prop_freq": 10,
    }

    mc.run(system, moves, 300.0, "equilibration", 10000, **custom_args)


if __name__ == "__main__":
    run_gcmc_adsorption()
