
import pytest
import mosdef_cassandra as mc

from mosdef_cassandra.test.base_test import BaseTest

class TestMoves(BaseTest):

    def test_invalid_ensemble(self,methane_oplsaa):
        with pytest.raises(ValueError,match=r'Invalid ensemble'):
            moves = mc.Moves('nvtn',[methane_oplsaa])

    def test_invalid_species_list(self,methane_oplsaa):
        with pytest.raises(TypeError,match=r'species_topologies should '
                                            'be a list of species'):
            moves = mc.Moves('nvt',methane_oplsaa)

    def test_invalid_species_type(self):
        with pytest.raises(TypeError,match=r'each species should be '
                                            'a parmed.Structure'):
            moves = mc.Moves('nvt',[1])

    def test_ensemble_nvt(self,methane_oplsaa):
        moves = mc.Moves('nvt',[methane_oplsaa])
        assert moves.ensemble == 'nvt'
        assert moves.prob_translate == 0.35
        assert moves.prob_rotate == 0.35
        assert moves.prob_regrow == 0.30
        assert moves.prob_vol == 0.0
        assert moves.prob_angle == 0.0
        assert moves.prob_dihedral == 0.0
        assert moves.prob_insert == 0.0
        assert moves.prob_swap == 0.0

        # Per box attributes
        assert len(moves.max_translate) == 1
        assert len(moves.max_rotate) == 1
        assert len(moves.max_volume) == 1
        assert len(moves.prob_swap_from_box) == 1
        assert moves.prob_swap_from_box[0] == 1.0
        assert moves.max_volume[0] == 0.
        # Per species-per-box
        assert len(moves.max_translate[0]) == 1
        assert len(moves.max_rotate[0]) == 1
        # Per species attributes
        assert len(moves.sp_insertable) == 1
        assert len(moves.sp_prob_swap) == 1
        assert len(moves.sp_prob_regrow) == 1

        # Should be insertable and regrow-able
        assert moves.sp_insertable[0] == True
        assert moves.sp_prob_regrow[0] == 1.0

    def test_ensemble_npt(self,methane_oplsaa):
        moves = mc.Moves('npt',[methane_oplsaa])
        assert moves.ensemble == 'npt'
        assert moves.prob_translate == 0.34
        assert moves.prob_rotate == 0.34
        assert moves.prob_regrow == 0.30
        assert moves.prob_vol == 0.02
        assert moves.prob_angle == 0.0
        assert moves.prob_dihedral == 0.0
        assert moves.prob_insert == 0.0
        assert moves.prob_swap == 0.0

        # Per box attributes
        assert len(moves.max_translate) == 1
        assert len(moves.max_rotate) == 1
        assert len(moves.max_volume) == 1
        assert len(moves.prob_swap_from_box) == 1
        assert moves.prob_swap_from_box[0] == 1.0
        assert moves.max_volume[0] == 500.
        # Per species-per-box
        assert len(moves.max_translate[0]) == 1
        assert len(moves.max_rotate[0]) == 1
        # Per species attributes
        assert len(moves.sp_insertable) == 1
        assert len(moves.sp_prob_swap) == 1
        assert len(moves.sp_prob_regrow) == 1

        # Should be insertable and regrow-able
        assert moves.sp_insertable[0] == True
        assert moves.sp_prob_regrow[0] == 1.0

    def test_ensemble_gcmc(self,methane_oplsaa):
        moves = mc.Moves('gcmc',[methane_oplsaa])
        assert moves.ensemble == 'gcmc'
        assert moves.prob_translate == 0.30
        assert moves.prob_rotate == 0.30
        assert moves.prob_regrow == 0.30
        assert moves.prob_vol == 0.0
        assert moves.prob_angle == 0.0
        assert moves.prob_dihedral == 0.0
        assert moves.prob_insert == 0.1
        assert moves.prob_swap == 0.0

        # Per box attributes
        assert len(moves.max_translate) == 1
        assert len(moves.max_rotate) == 1
        assert len(moves.max_volume) == 1
        assert len(moves.prob_swap_from_box) == 1
        assert moves.prob_swap_from_box[0] == 1.0
        assert moves.max_volume[0] == 0.
        # Per species-per-box
        assert len(moves.max_translate[0]) == 1
        assert len(moves.max_rotate[0]) == 1
        # Per species attributes
        assert len(moves.sp_insertable) == 1
        assert len(moves.sp_prob_swap) == 1
        assert len(moves.sp_prob_regrow) == 1

        # Should be insertable and regrow-able
        assert moves.sp_insertable[0] == True
        assert moves.sp_prob_regrow[0] == 1.0

    def test_ensemble_gemc(self,methane_oplsaa):
        moves = mc.Moves('gemc',[methane_oplsaa])
        assert moves.ensemble == 'gemc'
        assert moves.prob_translate == 0.29
        assert moves.prob_rotate == 0.29
        assert moves.prob_regrow == 0.30
        assert moves.prob_vol == 0.02
        assert moves.prob_angle == 0.0
        assert moves.prob_dihedral == 0.0
        assert moves.prob_insert == 0.0
        assert moves.prob_swap == 0.1

        # Per box attributes
        assert len(moves.max_translate) == 2
        assert len(moves.max_rotate) == 2
        assert len(moves.prob_swap_from_box) == 2
        assert moves.prob_swap_from_box[0] == 0.5
        assert moves.prob_swap_from_box[1] == 0.5
        assert len(moves.max_volume) == 1
        assert moves.max_volume[0] == 500.
        # Per species-per-box
        assert len(moves.max_translate[0]) == 1
        assert len(moves.max_rotate[0]) == 1
        assert len(moves.max_translate[1]) == 1
        assert len(moves.max_rotate[1]) == 1
        # Per species attributes
        assert len(moves.sp_insertable) == 1
        assert len(moves.sp_prob_swap) == 1
        assert len(moves.sp_prob_regrow) == 1

        # Should be insertable and regrow-able
        assert moves.sp_insertable[0] == True
        assert moves.sp_prob_regrow[0] == 1.0

    def test_ensemble_gemcnpt(self,methane_oplsaa):
        moves = mc.Moves('gemc_npt',[methane_oplsaa])
        assert moves.ensemble == 'gemc_npt'
        assert moves.prob_translate == 0.29
        assert moves.prob_rotate == 0.29
        assert moves.prob_regrow == 0.30
        assert moves.prob_vol == 0.02
        assert moves.prob_angle == 0.0
        assert moves.prob_dihedral == 0.0
        assert moves.prob_insert == 0.0
        assert moves.prob_swap == 0.1

        # Per box attributes
        assert len(moves.max_translate) == 2
        assert len(moves.max_rotate) == 2
        assert len(moves.prob_swap_from_box) == 2
        assert moves.prob_swap_from_box[0] == 0.5
        assert moves.prob_swap_from_box[1] == 0.5
        assert len(moves.max_volume) == 2
        assert moves.max_volume[0] == 500.
        assert moves.max_volume[1] == 5000.
        # Per species-per-box
        assert len(moves.max_translate[0]) == 1
        assert len(moves.max_rotate[0]) == 1
        assert len(moves.max_translate[1]) == 1
        assert len(moves.max_rotate[1]) == 1
        # Per species attributes
        assert len(moves.sp_insertable) == 1
        assert len(moves.sp_prob_swap) == 1
        assert len(moves.sp_prob_regrow) == 1

        # Should be insertable and regrow-able
        assert moves.sp_insertable[0] == True
        assert moves.sp_prob_regrow[0] == 1.0

    def test_single_site_nvt(self,methane_trappe):

        moves = mc.Moves('nvt',[methane_trappe])
        assert moves.ensemble == 'nvt'
        assert moves.prob_rotate == 0.0
        assert moves.prob_regrow == 0.0
        assert moves.prob_translate == 1.0
        assert moves.prob_vol == 0.0
        assert moves.prob_angle == 0.0
        assert moves.prob_dihedral == 0.0
        assert moves.prob_insert == 0.0
        assert moves.prob_swap == 0.0

        # Per box attributes
        assert len(moves.max_translate) == 1
        assert len(moves.max_rotate) == 1
        assert len(moves.max_volume) == 1
        assert len(moves.prob_swap_from_box) == 1
        assert moves.prob_swap_from_box[0] == 1.0
        assert moves.max_volume[0] == 0.
        # Per species-per-box
        assert len(moves.max_translate[0]) == 1
        assert len(moves.max_rotate[0]) == 1
        # Per species attributes
        assert len(moves.sp_insertable) == 1
        assert len(moves.sp_prob_swap) == 1
        assert len(moves.sp_prob_regrow) == 1

        # Should be insertable and NOT regrow-able
        assert moves.sp_insertable[0] == True
        assert moves.sp_prob_regrow[0] == 0.0


    def test_gcmc_lattice(self,fixed_lattice_trappe,methane_trappe):

        moves = mc.Moves('gcmc',[fixed_lattice_trappe,methane_trappe])
        assert moves.ensemble == 'gcmc'
        assert moves.prob_translate == pytest.approx(0.9)
        assert moves.prob_rotate == 0.0
        assert moves.prob_regrow == 0.0
        assert moves.prob_vol == 0.0
        assert moves.prob_angle == 0.0
        assert moves.prob_dihedral == 0.0
        assert moves.prob_insert == 0.1
        assert moves.prob_swap == pytest.approx(0.0)

        # Per box attributes
        assert len(moves.max_translate) == 1
        assert len(moves.max_rotate) == 1
        assert len(moves.max_volume) == 1
        assert len(moves.prob_swap_from_box) == 1
        assert moves.prob_swap_from_box[0] == 1.0
        assert moves.max_volume[0] == 0.
        # Per species-per-box
        assert len(moves.max_translate[0]) == 2
        assert len(moves.max_rotate[0]) == 2
        # Per species attributes
        assert len(moves.sp_insertable) == 2
        assert len(moves.sp_prob_swap) == 2
        assert len(moves.sp_prob_regrow) == 2

        # Lattice should not be insertable or regrow-able
        assert moves.sp_insertable[0] == False
        assert moves.sp_prob_regrow[0] == 0.0
        # Methane should be insertable and not regrow-able
        assert moves.sp_insertable[1] == True
        assert moves.sp_prob_regrow[1] == 0.0



