
import pytest
import mbuild
import foyer
import mosdef_cassandra as mc

class BaseTest:

    @pytest.fixture
    def methane_oplsaa(self):
        oplsaa = foyer.forcefields.load_OPLSAA()
        methane = mbuild.load('C',smiles=True)
        methane = oplsaa.apply(methane)
        return methane

    @pytest.fixture
    def methane_trappe(self):
        trappe = foyer.forcefields.load_TRAPPE_UA()
        methane = mbuild.Compound(pos=[0.,0.,0.],name='_CH4')
        methane = trappe.apply(methane)
        return methane

    @pytest.fixture
    def ethane_oplsaa(self):
        oplsaa = foyer.forcefields.load_OPLSAA()
        ethane = mbuild.load('CC',smiles=True)
        ethane = oplsaa.apply(ethane)
        return ethane

    @pytest.fixture
    def fixed_lattice(self):
        trappe = foyer.forcefields.load_TRAPPE_UA()
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
        
        typed_system = trappe.apply(system)

        return typed_system

    @pytest.fixture
    def box(self):
        box = mbuild.Box(lengths=[5.,5.,5.])
        return box

    @pytest.fixture
    def wrongbox(self):
        box = mbuild.Compound()
        box.boundingbox.lengths = [5.,5.,5.]
        return box

    @pytest.fixture
    def methane_single(self):
        methane = mbuild.load('C',smiles=True)
        methane.boundingbox.lengths = [5.,5.,5.]
        return methane


