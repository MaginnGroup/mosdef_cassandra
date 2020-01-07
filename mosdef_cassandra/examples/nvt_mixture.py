
import mbuild
import foyer
import mosdef_cassandra as mc

def run_nvt_mixture():
    # Use mbuild to create molecules
    methane = mbuild.load('C',smiles=True)
    propane = mbuild.load('CCC',smiles=True)

    # Create an empty mbuild.Box
    box = mbuild.Box(lengths=[3.,3.,3.])

    # Load forcefields
    oplsaa = foyer.forcefields.load_OPLSAA()

    # Use foyer to apply forcefields
    typed_methane = oplsaa.apply(methane)
    typed_propane = oplsaa.apply(propane)

    # Create box and species list
    box_list = [box]
    species_list = [typed_methane, typed_propane]

    # Use Cassandra to insert some initial number of species
    species_to_add = [[100,50]]

    system = mc.System(box_list,species_list,
                       species_to_add=species_to_add)
    moves = mc.Moves('nvt', species_list)

    mc.run(system,moves,200.0,'equilibration',10000)

if __name__ == "__main__":
    run_nvt_mixture()
