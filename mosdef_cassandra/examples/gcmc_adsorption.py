
import mbuild
import foyer
import mosdef_cassandra as mc
from mosdef_cassandra.examples.structures import carbon_lattice

# Use mbuild to create molecules
lattice = carbon_lattice()
methane = mbuild.load('C',smiles=True)

# Use foyer to apply forcefield
trappe = foyer.forcefields.load_TRAPPE_UA()
oplsaa = foyer.forcefields.load_OPLSAA()

typed_lattice = trappe.apply(lattice)
typed_methane = oplsaa.apply(methane)

box_list = [lattice]
species_list = [typed_lattice,typed_methane]

system = mc.System(box_list,species_list,species_in_boxes=[[1,0]])
moves = mc.Moves('gcmc', species_list)

mc.run(system,moves,300.0,'equilibration',10000,chemical_potentials=[0,100],coord_freq=100)

