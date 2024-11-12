import mbuild
import mosdef_cassandra as mc
import unyt as u
from gmso.core.forcefield import ForceField
from gmso.external import from_mbuild
import forcefield_utilities as ffutils
from gmso.parameterization import apply


def run_nvt_gmso(**custom_args):

    # Use mBuild to create a methane molecule
    methane = mbuild.load("C", smiles=True)
    methane_top = from_mbuild(methane)
    methane_top.identify_connections()
    ff = ffutils.FoyerFFs().load("oplsaa").to_gmso_ff()
    methane_top = apply(methane_top, ff, remove_untyped=True)
    methane_top.save("methane.mcf", overwrite=True)

    # Create an empty mbuild.Box
    box = mbuild.Box(lengths=[3.0, 3.0, 3.0])
    packed_system = mbuild.fill_box(compound=methane, n_compounds=10, box=box)

    # Create box and species list
    box_list = [packed_system]
    species_list = [methane_top]

    # Define the System
    system = mc.System(box_list, species_list, mols_in_boxes=[[10]])
    # Define the MoveSet
    moveset = mc.MoveSet("nvt", species_list)

    # Run a simulation at 300 K for 10000 MC moves
    mc.run(
        system=system,
        moveset=moveset,
        run_type="equilibration",
        run_length=10000,
        temperature=300.0 * u.K,
        seeds=[12345, 12345],
        run_name="nvt_gmso",
        **custom_args,
    )


#

if __name__ == "__main__":
    run_nvt()
