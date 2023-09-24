import mbuild
import foyer
import mosdef_cassandra as mc
import unyt as u


def run_nvt_mbuild(fix_bonds, **custom_args):

    dme = mbuild.load("COC", smiles=True)
    dee = mbuild.load("CCOCC", smiles=True)

    # Create an empty mbuild.Box
    box = mbuild.Box(lengths=[3.0, 3.0, 3.0])

    # fill box
    box = mbuild.fill_box([dme, dee], n_compounds=[10, 10], box=box)

    # Load forcefields
    ff = foyer.forcefields.load_OPLSAA()

    # Use foyer to apply forcefields
    dme_ff = ff.apply(dme)
    dee_ff = ff.apply(dee)

    # Create box and species list
    box_list = [box]
    species_list = [dme_ff, dee_ff]

    # Use Cassandra to insert some initial number of species
    mols_in_boxes = [[10, 10]]

    # Define the system object
    system = mc.System(
        box_list,
        species_list,
        mols_in_boxes=mols_in_boxes,
        fix_bonds=fix_bonds,
    )
    # Get the move probabilities
    moveset = mc.MoveSet("nvt", species_list)

    # Run a simulation with at 300 K with 10000 MC moves
    mc.run(
        system=system,
        moveset=moveset,
        run_type="equilibration",
        run_length=100,
        temperature=300.0 * u.K,
        **custom_args,
    )


if __name__ == "__main__":
    run_nvt_mbuild(fix_bonds=True)
