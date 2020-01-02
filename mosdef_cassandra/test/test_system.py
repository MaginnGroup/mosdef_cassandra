
import pytest
import mosdef_cassandra as mc

from mosdef_cassandra.test.base_test import BaseTest

class TestSystem(BaseTest):

    def test_invalid_box_list(self,methane_oplsaa,box):
        with pytest.raises(TypeError,match=r'"boxes" should be a list'):
            system = mc.System(box,[methane_oplsaa],
                               species_to_add=[[10]])

    def test_invalid_box_type(self,methane_oplsaa,box):
        with pytest.raises(TypeError,match=r'Each box should be'):
            system = mc.System([methane_oplsaa],[methane_oplsaa],
                               species_to_add=[[10]])

    def test_invalid_species_list(self,methane_oplsaa,box):
        with pytest.raises(TypeError,match=r'"species_topologies" should '
                                            'be a list'):
            system = mc.System([box],methane_oplsaa,
                               species_to_add=[[10]])

    def test_invalid_species_type(self,box):
        with pytest.raises(TypeError,match=r'Each species should be '
                                            'a parmed.Structure'):
            system = mc.System([box],[1.0],
                               species_to_add=[[10]])

    def test_invalid_add_species_list(self,methane_oplsaa,box):
        with pytest.raises(TypeError,match=r'"species_to_add" should be'):
            system = mc.System([box],[methane_oplsaa],
                               species_to_add=10)
        with pytest.raises(TypeError,match=r'with one list for each box'):
            system = mc.System([box],[methane_oplsaa],
                               species_to_add=[10])
        with pytest.raises(TypeError,match=r'added to each box must'):
            system = mc.System([box],[methane_oplsaa],
                               species_to_add=[[10.0]])

    def test_invalid_current_species_list(self,methane_oplsaa,box):
        with pytest.raises(TypeError,match=r'"species_in_boxes" should be'):
            system = mc.System([box],[methane_oplsaa],
                               species_in_boxes=10)
        with pytest.raises(TypeError,match=r'with one list for each box'):
            system = mc.System([box],[methane_oplsaa],
                               species_in_boxes=[10])
        with pytest.raises(TypeError,match=r'in each box must be'):
            system = mc.System([box],[methane_oplsaa],
                               species_in_boxes=[[10.0]])

    def test_nboxes_mismatch(self,methane_oplsaa,box):
        with pytest.raises(ValueError,match=r'The number of boxes inferred '
                                'from the length of "species_in_boxes"'):
            system = mc.System([box],[methane_oplsaa],
                               species_in_boxes=[[0],[0]])
        with pytest.raises(ValueError,match=r'The number of boxes inferred '
                                'from the length of "species_to_add"'):
            system = mc.System([box],[methane_oplsaa],
                               species_to_add=[[10],[10]])

    def test_natoms_mismatch(self,methane_single,methane_oplsaa,box,
                             wrongbox):
        with pytest.raises(ValueError,match=r'is an mbuild.Box'):
            system = mc.System([box],[methane_oplsaa],
                               species_in_boxes=[[10]])
        with pytest.raises(ValueError,match=r'number of atoms'):
            system = mc.System([methane_single],[methane_oplsaa],
                               species_in_boxes=[[10]])
        with pytest.raises(ValueError,match=r'NOTE'):
            system = mc.System([wrongbox],[methane_oplsaa],
                               species_in_boxes=[[10]])


    def test_edit_boxes(self,methane_oplsaa,box):
        with pytest.raises(AttributeError,match=r'cannot be modified'):
            system = mc.System([box],[methane_oplsaa],
                    species_to_add=[[10]])
            system.boxes = [box]

    def test_edit_topologies(self,methane_oplsaa,box):
        with pytest.raises(AttributeError,match=r'cannot be modified'):
            system = mc.System([box],[methane_oplsaa],
                    species_to_add=[[10]])
            system.species_topologies = [methane_oplsaa]

    def test_edit_nspecies(self,methane_oplsaa,box):
        with pytest.raises(AttributeError,match=r'cannot be modified'):
            system = mc.System([box],[methane_oplsaa],
                    species_to_add=[[10]])
            system.species_in_boxes = [[0]]

    def test_edit_add_nspecies(self,methane_oplsaa,box):
        system = mc.System([box],[methane_oplsaa],
                    species_to_add=[[10]])
        system.species_to_add = [[100]]

        assert system.species_to_add == [[100]]
