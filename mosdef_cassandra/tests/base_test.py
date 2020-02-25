import pytest
import mbuild
import foyer
import mosdef_cassandra as mc


class BaseTest:
    @pytest.fixture
    def methane_oplsaa(self):
        oplsaa = foyer.forcefields.load_OPLSAA()
        methane = mbuild.load("C", smiles=True)
        methane = oplsaa.apply(methane)
        return methane

    @pytest.fixture
    def butane_oplsaa(self):
        oplsaa = foyer.forcefields.load_OPLSAA()
        butane = mbuild.load("CCCC", smiles=True)
        butane = oplsaa.apply(butane)
        return butane

    @pytest.fixture
    def methane_trappe(self):
        trappe = foyer.forcefields.load_TRAPPE_UA()
        methane = mbuild.Compound(pos=[0.0, 0.0, 0.0], name="_CH4")
        methane = trappe.apply(methane)
        return methane

    @pytest.fixture
    def fixed_lattice_compound(self):
        carbon = mbuild.Compound(name="_CH4")
        angles = [90.0, 90.0, 90.0]
        carbon_locations = [
            [0.0, 0.0, 0.0],
            [1.0 / 2.0, 1.0 / 2.0, 0.0],
            [0.0, 1.0 / 2.0, 1.0 / 2.0],
            [1.0 / 2.0, 0.0, 0.0],
            [0.0, 0.0, 1.0 / 2.0],
        ]
        basis = {"C": carbon_locations}
        carbon_dict = {"C": carbon}
        lattice_spacing = [0.746, 0.746, 0.746]
        lattice = mbuild.Lattice(
            lattice_spacing=lattice_spacing,
            angles=angles,
            lattice_points=basis,
        )
        system = lattice.populate(compound_dict=carbon_dict, x=4, y=4, z=4)

        return system

    @pytest.fixture
    def fixed_lattice_trappe(self, fixed_lattice_compound):
        trappe = foyer.forcefields.load_TRAPPE_UA()
        lattice = fixed_lattice_compound
        typed_system = trappe.apply(lattice)

        return typed_system

    @pytest.fixture
    def box(self):
        box = mbuild.Box(lengths=[5.0, 5.0, 5.0])
        return box

    @pytest.fixture
    def wrongbox(self):
        box = mbuild.Compound()
        box.boundingbox.lengths = [5.0, 5.0, 5.0]
        return box

    @pytest.fixture
    def methane_single(self):
        methane = mbuild.load("C", smiles=True)
        methane.boundingbox.lengths = [5.0, 5.0, 5.0]
        return methane
