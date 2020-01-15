
import mbuild
import foyer
import mosdef_cassandra as mc

def run_gcmc():

    # Use mbuild to create molecules
    methane = mbuild.load('C',smiles=True)

    # Create an empty mbuild.Box
    box = mbuild.Box(lengths=[10.,10.,10.])

    # Load forcefields
    oplsaa = foyer.forcefields.load_OPLSAA()

    # Use foyer to apply forcefields
    typed_methane = oplsaa.apply(methane)

    # Create box and species list
    box_list = [box]
    species_list = [typed_methane]

    mols_to_add = [[100]]

    system = mc.System(box_list,species_list,
                       mols_to_add=mols_to_add)
    moves = mc.Moves('gcmc', species_list)

    mc.run(system,moves,300.0,'equilibration',1000,
            chemical_potentials=[-35.],
            prop_freq=100)

if __name__ == "__main__":
    run_gcmc()
