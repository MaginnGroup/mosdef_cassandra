
import mbuild

def carbon_lattice():
    carbon = mbuild.Compound(name='_CH4')
    angles = [90.,90.,90.]
    carbon_locations = [[0.,0.,0.],
                        [1./2.,1./2.,0.],
                        [0.,1./2.,1./2.],
                        [1./2.,0.,0.],
                        [0.,0.,1./2.]]
    basis = {'C' : carbon_locations }
    carbon_dict = {'C' : carbon}
    lattice = mbuild.Lattice(lattice_spacing = [0.746,0.746,0.746],
                             angles = angles,
                             lattice_points = basis)

    system = lattice.populate(
                        compound_dict=carbon_dict,x=4,y=4,z=4)

    return system

