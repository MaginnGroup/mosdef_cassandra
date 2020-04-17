import pytest

import mosdef_cassandra as mc
import warnings
from mosdef_cassandra.tests.base_test import BaseTest


class TestMoves(BaseTest):
    def test_invalid_ensemble(self, methane_oplsaa):
        with pytest.raises(ValueError, match=r"Invalid ensemble"):
            moves = mc.Moves("nvtn", [methane_oplsaa])

    def test_invalid_species_list(self, methane_oplsaa):
        with pytest.raises(
            TypeError,
            match=r"species_topologies should " "be a list of species",
        ):
            moves = mc.Moves("nvt", methane_oplsaa)

    def test_invalid_species_type(self):
        with pytest.raises(
            TypeError, match=r"each species should be " "a parmed.Structure"
        ):
            moves = mc.Moves("nvt", [1])

    def test_ensemble_nvt(self, methane_oplsaa):
        moves = mc.Moves("nvt", [methane_oplsaa])
        assert moves.ensemble == "nvt"
        assert moves.prob_translate == 0.35
        assert moves.prob_rotate == 0.35
        assert moves.prob_regrow == 0.30
        assert moves.prob_volume == 0.0
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
        assert moves.max_volume[0] == 0.0
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

    def test_ensemble_npt(self, methane_oplsaa):
        moves = mc.Moves("npt", [methane_oplsaa])
        assert moves.ensemble == "npt"
        assert moves.prob_translate == 0.34
        assert moves.prob_rotate == 0.34
        assert moves.prob_regrow == 0.30
        assert moves.prob_volume == 0.02
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
        assert moves.max_volume[0] == 500.0
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

    def test_ensemble_gcmc(self, methane_oplsaa):
        moves = mc.Moves("gcmc", [methane_oplsaa])
        assert moves.ensemble == "gcmc"
        assert moves.prob_translate == 0.25
        assert moves.prob_rotate == 0.25
        assert moves.prob_regrow == 0.30
        assert moves.prob_volume == 0.0
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
        assert moves.max_volume[0] == 0.0
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

    @pytest.mark.parametrize(
        "typ,value",
        [
            ("slitpore", 1),
            ("cylinder", 1),
            ("sphere", 1),
            ("interface", [1, 2]),
        ],
    )
    def test_restricted_gcmc(self, methane_oplsaa, typ, value):
        moves = mc.Moves("gcmc", [methane_oplsaa])
        moves.add_restricted_insertions([methane_oplsaa], [[typ]], [[value]])

        assert moves._restricted_type == [[typ]]
        assert moves._restricted_value == [[value]]

    def test_ensemble_gemc(self, methane_oplsaa):
        moves = mc.Moves("gemc", [methane_oplsaa])
        assert moves.ensemble == "gemc"
        assert moves.prob_translate == 0.29
        assert moves.prob_rotate == 0.29
        assert moves.prob_regrow == 0.30
        assert moves.prob_volume == 0.02
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
        assert moves.max_volume[0] == 500.0
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

    @pytest.mark.parametrize(
        "typ,value",
        [
            ("slitpore", 1),
            ("cylinder", 1),
            ("sphere", 1),
            ("interface", [1, 2]),
        ],
    )
    def test_restricted_gemc(self, methane_oplsaa, typ, value):
        moves = mc.Moves("gemc", [methane_oplsaa])
        moves.add_restricted_insertions(
            [methane_oplsaa], [[None], [typ]], [[None], [value]]
        )

        assert moves._restricted_type == [[None], [typ]]
        assert moves._restricted_value == [[None], [value]]

    def test_ensemble_gemcnpt(self, methane_oplsaa):
        moves = mc.Moves("gemc_npt", [methane_oplsaa])
        assert moves.ensemble == "gemc_npt"
        assert moves.prob_translate == 0.29
        assert moves.prob_rotate == 0.29
        assert moves.prob_regrow == 0.30
        assert moves.prob_volume == 0.02
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
        assert moves.max_volume[0] == 500.0
        assert moves.max_volume[1] == 5000.0
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

    def test_restricted_gemc_npt(self, methane_oplsaa):
        moves = mc.Moves("gemc_npt", [methane_oplsaa])
        moves.add_restricted_insertions(
            [methane_oplsaa], [[None], ["slitpore"]], [[None], [3]]
        )

        assert moves._restricted_type == [[None], ["slitpore"]]
        assert moves._restricted_value == [[None], [3]]

    def test_single_site_nvt(self, methane_trappe):

        moves = mc.Moves("nvt", [methane_trappe])
        assert moves.ensemble == "nvt"
        assert moves.prob_rotate == 0.0
        assert moves.prob_regrow == 0.0
        assert moves.prob_translate == 1.0
        assert moves.prob_volume == 0.0
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
        assert moves.max_volume[0] == 0.0
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

    def test_single_site_gemc(self, methane_trappe):

        moves = mc.Moves("gemc", [methane_trappe])
        assert moves.ensemble == "gemc"
        assert moves.prob_rotate == 0.0
        assert moves.prob_regrow == 0.0
        assert moves.prob_translate == pytest.approx(0.88)
        assert moves.prob_volume == 0.02
        assert moves.prob_angle == 0.0
        assert moves.prob_dihedral == 0.0
        assert moves.prob_insert == 0.0
        assert moves.prob_swap == 0.1

        # Per box attributes
        assert len(moves.max_translate) == 2
        assert len(moves.max_rotate) == 2
        assert len(moves.max_volume) == 1
        assert len(moves.prob_swap_from_box) == 2
        assert moves.prob_swap_from_box[0] == 0.5
        assert moves.prob_swap_from_box[1] == 0.5
        assert moves.max_volume[0] == 500.0
        # Per species-per-box
        assert len(moves.max_translate[0]) == 1
        assert len(moves.max_translate[1]) == 1
        assert len(moves.max_rotate[0]) == 1
        assert len(moves.max_rotate[1]) == 1
        # Per species attributes
        assert len(moves.sp_insertable) == 1
        assert len(moves.sp_prob_swap) == 1
        assert len(moves.sp_prob_regrow) == 1

        # Should be insertable and NOT regrow-able
        assert moves.sp_insertable[0] == True
        assert moves.sp_prob_regrow[0] == 0.0

    def test_gcmc_lattice(self, fixed_lattice_trappe, methane_trappe):

        moves = mc.Moves("gcmc", [fixed_lattice_trappe, methane_trappe])
        assert moves.ensemble == "gcmc"
        assert moves.prob_translate == pytest.approx(0.8)
        assert moves.prob_rotate == 0.0
        assert moves.prob_regrow == 0.0
        assert moves.prob_volume == 0.0
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
        assert moves.max_volume[0] == 0.0
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

    def test_prob_setters(self, methane_oplsaa):
        moves = mc.Moves("nvt", [methane_oplsaa])
        with pytest.raises(TypeError, match=r"prob_translate must"):
            moves.prob_translate = []
        with pytest.raises(TypeError, match=r"prob_translate must"):
            moves.prob_translate = True
        with pytest.raises(ValueError, match=r"Probability must"):
            moves.prob_translate = -0.2

        with pytest.raises(TypeError, match=r"prob_rotate must"):
            moves.prob_rotate = []
        with pytest.raises(TypeError, match=r"prob_rotate must"):
            moves.prob_rotate = True
        with pytest.raises(ValueError, match=r"Probability must"):
            moves.prob_rotate = -0.2

        with pytest.raises(TypeError, match=r"prob_angle must"):
            moves.prob_angle = []
        with pytest.raises(TypeError, match=r"prob_angle must"):
            moves.prob_angle = True
        with pytest.raises(ValueError, match=r"Probability must"):
            moves.prob_angle = -0.2

        with pytest.raises(TypeError, match=r"prob_dihedral must"):
            moves.prob_dihedral = []
        with pytest.raises(TypeError, match=r"prob_dihedral must"):
            moves.prob_dihedral = True
        with pytest.raises(ValueError, match=r"Probability must"):
            moves.prob_dihedral = -0.2

        with pytest.raises(TypeError, match=r"prob_regrow must"):
            moves.prob_regrow = []
        with pytest.raises(TypeError, match=r"prob_regrow must"):
            moves.prob_regrow = True
        with pytest.raises(ValueError, match=r"Probability must"):
            moves.prob_regrow = -0.2

        with pytest.raises(TypeError, match=r"prob_volume must"):
            moves.prob_volume = []
        with pytest.raises(TypeError, match=r"prob_volume must"):
            moves.prob_volume = True
        with pytest.raises(ValueError, match=r"Probability must"):
            moves.prob_volume = -0.2
        with pytest.raises(ValueError, match=r"Ensemble is nvt."):
            moves.prob_volume = 0.02

        moves = mc.Moves("gemc", [methane_oplsaa])
        with pytest.raises(ValueError, match=r"Ensemble is gemc."):
            moves.prob_volume = 0.0

        moves = mc.Moves("nvt", [methane_oplsaa])
        with pytest.raises(TypeError, match=r"prob_insert must"):
            moves.prob_insert = []
        with pytest.raises(TypeError, match=r"prob_insert must"):
            moves.prob_insert = True
        with pytest.raises(ValueError, match=r"Probability must"):
            moves.prob_insert = -0.2
        with pytest.raises(ValueError, match=r"Ensemble is nvt."):
            moves.prob_insert = 0.2
        moves = mc.Moves("gcmc", [methane_oplsaa])
        with pytest.raises(ValueError, match=r"Ensemble is gcmc."):
            moves.prob_insert = 0.0

        moves = mc.Moves("nvt", [methane_oplsaa])
        with pytest.raises(TypeError, match=r"prob_swap must"):
            moves.prob_swap = []
        with pytest.raises(TypeError, match=r"prob_swap must"):
            moves.prob_swap = True
        with pytest.raises(ValueError, match=r"Probability must"):
            moves.prob_swap = -0.2
        with pytest.raises(ValueError, match=r"Ensemble is nvt."):
            moves.prob_swap = 0.2
        moves = mc.Moves("gemc", [methane_oplsaa])
        with pytest.raises(ValueError, match=r"Ensemble is gemc."):
            moves.prob_swap = 0.0

    def test_maxval_setters(self, methane_oplsaa):
        moves = mc.Moves("gemc", [methane_oplsaa])
        with pytest.raises(ValueError, match=r"must be a list"):
            moves.max_translate = 1.0
        with pytest.raises(ValueError, match=r"must be a list"):
            moves.max_translate = [[1.0, 1.0], [1.0]]
        with pytest.raises(TypeError, match=r"of type float"):
            moves.max_translate = [[1.0], [True]]
        with pytest.raises(ValueError, match=r"cannot be less than zero"):
            moves.max_translate = [[1.0], [-1.0]]
        moves.max_translate = [[1.0], [1]]
        assert moves.max_translate[1][0] == 1.0

        with pytest.raises(ValueError, match=r"must be a list"):
            moves.max_rotate = 1.0
        with pytest.raises(ValueError, match=r"must be a list"):
            moves.max_rotate = [[1.0, 1.0], [1.0]]
        with pytest.raises(TypeError, match=r"of type float"):
            moves.max_rotate = [[1.0], [True]]
        with pytest.raises(ValueError, match=r"cannot be less than zero"):
            moves.max_rotate = [[1.0], [-1.0]]
        moves.max_rotate = [[1.0], [1]]
        assert moves.max_rotate[1][0] == 1.0

        with pytest.raises(ValueError, match=r"must be a list"):
            moves.max_dihedral = 1.0
        with pytest.raises(ValueError, match=r"must be a list"):
            moves.max_dihedral = [1.0, 1.0, 1.0]
        with pytest.raises(TypeError, match=r"of type float"):
            moves.max_dihedral = [True]
        with pytest.raises(ValueError, match=r"cannot be less than zero"):
            moves.max_dihedral = [-1.0]
        moves.max_dihedral = [1.0]
        assert moves.max_dihedral[0] == 1.0

        with pytest.raises(ValueError, match=r"must be a list"):
            moves.prob_swap_from_box = 1.0
        with pytest.raises(ValueError, match=r"must be a list"):
            moves.prob_swap_from_box = [1.0, 1.0, 1.0]
        with pytest.raises(TypeError, match=r"of type float"):
            moves.prob_swap_from_box = [0.5, True]
        with pytest.raises(ValueError, match=r"cannot be less than zero"):
            moves.prob_swap_from_box = [0.5, -0.5]
        moves.prob_swap_from_box = [0.5, 0.5]
        assert moves.prob_swap_from_box[0] == 0.5

        moves = mc.Moves("gemc", [methane_oplsaa])
        with pytest.raises(ValueError, match=r"must be a list"):
            moves.max_volume = 1.0
        with pytest.raises(ValueError, match=r"must be a list"):
            moves.max_volume = [1.0, 1.0, 1.0]
        with pytest.raises(TypeError, match=r"of type float"):
            moves.max_volume = [True]
        with pytest.raises(ValueError, match=r"cannot be less than zero"):
            moves.max_volume = [-100]
        moves.max_volume = [5000.0]
        assert moves.max_volume[0] == 5000.0
        moves = mc.Moves("gemc_npt", [methane_oplsaa])
        with pytest.raises(ValueError, match=r"must be a list"):
            moves.max_volume = 1.0
        with pytest.raises(ValueError, match=r"must be a list"):
            moves.max_volume = [1.0, 1.0, 1.0]
        with pytest.raises(TypeError, match=r"of type float"):
            moves.max_volume = [True, 50000.0]
        with pytest.raises(ValueError, match=r"cannot be less than zero"):
            moves.max_volume = [-100, 100.0]
        moves.max_volume = [5000.0, 50000.0]
        assert moves.max_volume[1] == 50000.0

        moves = mc.Moves("gemc", [methane_oplsaa])
        with pytest.raises(ValueError, match=r"must be a list"):
            moves.sp_insertable = 1.0
        with pytest.raises(ValueError, match=r"must be a list"):
            moves.sp_insertable = [1.0, 1.0, 1.0]
        with pytest.raises(TypeError, match=r"as a boolean type"):
            moves.sp_insertable = [2.0]
        moves.sp_insertable = [True]
        assert moves.sp_insertable[0] == True

        with pytest.raises(ValueError, match=r"must be a list"):
            moves.sp_prob_swap = 1.0
        with pytest.raises(ValueError, match=r"must be a list"):
            moves.sp_prob_swap = [1.0, 1.0, 1.0]
        with pytest.raises(TypeError, match=r"of type float"):
            moves.sp_prob_swap = [True]
        with pytest.raises(ValueError, match=r"cannot be less than zero"):
            moves.sp_prob_swap = [-1.0]
        moves.sp_prob_swap = [1.0]
        assert moves.sp_prob_swap[0] == 1.0

        with pytest.raises(ValueError, match=r"must be a list"):
            moves.sp_prob_regrow = 1.0
        with pytest.raises(ValueError, match=r"must be a list"):
            moves.sp_prob_regrow = [1.0, 1.0, 1.0]
        with pytest.raises(TypeError, match=r"of type float"):
            moves.sp_prob_regrow = [True]
        with pytest.raises(ValueError, match=r"cannot be less than zero"):
            moves.sp_prob_regrow = [-1.0]
        moves.sp_prob_regrow = [1.0]
        assert moves.sp_prob_regrow[0] == 1.0

    def test_print_moves(self, methane_oplsaa):
        """Simple test to make sure moves object is printed"""

        moves = mc.Moves("gemc", [methane_oplsaa])
        moves.print()

    def test_invalid_ensemble_and_restriction(self, methane_oplsaa):
        moves = mc.Moves("nvt", [methane_oplsaa])
        with pytest.raises(ValueError, match=r"only valid"):
            moves.add_restricted_insertions(
                [methane_oplsaa], [["slitpore"]], [[1]]
            )

    @pytest.mark.parametrize(
        "typ,value",
        [("slitpore", [[1], [2]]), ("cylinder", [[None]]), (None, [[1]]),],
    )
    def test_value_error_restricted_type_and_value(
        self, methane_oplsaa, typ, value
    ):
        moves = mc.Moves("gcmc", [methane_oplsaa])
        with pytest.raises(ValueError):
            moves.add_restricted_insertions([methane_oplsaa], [[typ]], value)

    def test_type_error_restricted_type_and_value(self, methane_oplsaa):
        moves = mc.Moves("gcmc", [methane_oplsaa])
        with pytest.raises(TypeError):
            moves.add_restricted_insertions(
                [methane_oplsaa], ["cylinder"], [3]
            )

    def test_invalid_restricted_type_and_species(self, methane_oplsaa):
        moves = mc.Moves("gcmc", [methane_oplsaa])
        with pytest.raises(ValueError, match=r"Length of 'species'"):
            moves.add_restricted_insertions(
                [methane_oplsaa], [["slitpore", None]], [[1, None]]
            )

    def test_add_multiple_restricted_insertions(self, methane_oplsaa):
        moves = mc.Moves("gcmc", [methane_oplsaa])
        moves.add_restricted_insertions(
            [methane_oplsaa], [["slitpore"]], [[3]]
        )
        with pytest.warns(None) as record:
            moves.add_restricted_insertions(
                [methane_oplsaa], [["cylinder"]], [[3]]
            )
