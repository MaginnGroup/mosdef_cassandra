import pytest
import warnings
import unyt as u
import mosdef_cassandra as mc

from mosdef_cassandra.tests.base_test import BaseTest
from unyt.exceptions import IterableUnitCoercionError


class TestMoveSet(BaseTest):
    def test_invalid_ensemble(self, methane_oplsaa):
        with pytest.raises(ValueError, match=r"Invalid ensemble"):
            moveset = mc.MoveSet("nvtn", [methane_oplsaa])

    def test_invalid_species_list(self, methane_oplsaa):
        with pytest.raises(
            TypeError,
            match=r"species_topologies should " "be a list of species",
        ):
            moveset = mc.MoveSet("nvt", methane_oplsaa)

    def test_invalid_species_type(self):
        with pytest.raises(
            TypeError,
            match=r"Each species should be a "
            "parmed.Structure or gmso.Topology"
            "and must be of the same type",
        ):
            moveset = mc.MoveSet("nvt", [1])

    def test_ensemble_nvt(self, methane_oplsaa):
        moveset = mc.MoveSet("nvt", [methane_oplsaa])
        assert moveset.ensemble == "nvt"
        assert moveset.prob_translate == 0.33
        assert moveset.prob_rotate == 0.33
        assert moveset.prob_regrow == 0.34
        assert moveset.prob_volume == 0.0
        assert moveset.prob_angle == 0.0
        assert moveset.prob_dihedral == 0.0
        assert moveset.prob_insert == 0.0
        assert moveset.prob_swap == 0.0

        # Per box attributes
        assert len(moveset.max_translate) == 1
        assert len(moveset.max_rotate) == 1
        assert len(moveset.max_volume) == 1
        assert len(moveset.prob_swap_from_box) == 1
        assert moveset.prob_swap_from_box[0] == 1.0
        assert moveset.max_volume[0] == 0.0 * (u.angstrom**3)
        # Per species-per-box
        assert len(moveset.max_translate[0]) == 1
        assert len(moveset.max_rotate[0]) == 1
        # Per species attributes
        assert len(moveset.insertable) == 1
        assert len(moveset.prob_swap_species) == 1
        assert len(moveset.prob_regrow_species) == 1
        # Should be regrow-able but not insertable
        assert moveset.prob_regrow_species[0] == 1.0
        assert moveset.insertable[0] == False

    def test_ensemble_npt(self, methane_oplsaa):
        moveset = mc.MoveSet("npt", [methane_oplsaa])
        assert moveset.ensemble == "npt"
        assert moveset.prob_translate == 0.33
        assert moveset.prob_rotate == 0.33
        assert moveset.prob_regrow == 0.335
        assert moveset.prob_volume == 0.005
        assert moveset.prob_angle == 0.0
        assert moveset.prob_dihedral == 0.0
        assert moveset.prob_insert == 0.0
        assert moveset.prob_swap == 0.0

        # Per box attributes
        assert len(moveset.max_translate) == 1
        assert len(moveset.max_rotate) == 1
        assert len(moveset.max_volume) == 1
        assert len(moveset.prob_swap_from_box) == 1
        assert moveset.prob_swap_from_box[0] == 1.0
        assert moveset.max_volume[0] == 500.0 * (u.angstrom**3)
        # Per species-per-box
        assert len(moveset.max_translate[0]) == 1 * u.angstrom
        assert len(moveset.max_rotate[0]) == 1 * u.degree
        # Per species attributes
        assert len(moveset.insertable) == 1
        assert len(moveset.prob_swap_species) == 1
        assert len(moveset.prob_regrow_species) == 1
        # Should be regrow-able but not insertable
        assert moveset.prob_regrow_species[0] == 1.0
        assert moveset.insertable[0] == False

    def test_ensemble_gcmc(self, methane_oplsaa):
        moveset = mc.MoveSet("gcmc", [methane_oplsaa])
        assert moveset.ensemble == "gcmc"
        assert moveset.prob_translate == 0.25
        assert moveset.prob_rotate == 0.25
        assert moveset.prob_regrow == 0.30
        assert moveset.prob_volume == 0.0
        assert moveset.prob_angle == 0.0
        assert moveset.prob_dihedral == 0.0
        assert moveset.prob_insert == 0.1
        assert moveset.prob_swap == 0.0

        # Per box attributes
        assert len(moveset.max_translate) == 1
        assert len(moveset.max_rotate) == 1
        assert len(moveset.max_volume) == 1
        assert len(moveset.prob_swap_from_box) == 1
        assert moveset.prob_swap_from_box[0] == 1.0
        assert moveset.max_volume[0] == 0.0 * (u.angstrom**3)
        # Per species-per-box
        assert len(moveset.max_translate[0]) == 1
        assert len(moveset.max_rotate[0]) == 1
        # Per species attributes
        assert len(moveset.insertable) == 1
        assert len(moveset.prob_swap_species) == 1
        assert len(moveset.prob_regrow_species) == 1
        # Should be insertable and regrow-able
        assert moveset.insertable[0] == True
        assert moveset.prob_regrow_species[0] == 1.0

    @pytest.mark.parametrize(
        "typ,value",
        [
            ("slitpore", 1 * u.angstrom),
            ("cylinder", 1 * u.angstrom),
            ("sphere", 1 * u.angstrom),
            ("interface", [1 * u.angstrom, 2 * u.angstrom]),
        ],
    )
    def test_restricted_gcmc(self, methane_oplsaa, typ, value):
        moveset = mc.MoveSet("gcmc", [methane_oplsaa])
        moveset.add_restricted_insertions([methane_oplsaa], [[typ]], [[value]])

        assert moveset._restricted_type == [[typ]]
        assert moveset._restricted_value == [[value]]

    def test_ensemble_gemc(self, methane_oplsaa):
        moveset = mc.MoveSet("gemc", [methane_oplsaa])
        assert moveset.ensemble == "gemc"
        assert moveset.prob_translate == 0.30
        assert moveset.prob_rotate == 0.30
        assert moveset.prob_regrow == 0.295
        assert moveset.prob_volume == 0.005
        assert moveset.prob_angle == 0.0
        assert moveset.prob_dihedral == 0.0
        assert moveset.prob_insert == 0.0
        assert moveset.prob_swap == 0.1

        # Per box attributes
        assert len(moveset.max_translate) == 2
        assert len(moveset.max_rotate) == 2
        assert len(moveset.prob_swap_from_box) == 2
        assert moveset.prob_swap_from_box[0] == 0.5
        assert moveset.prob_swap_from_box[1] == 0.5
        assert len(moveset.max_volume) == 1
        assert moveset.max_volume[0] == 500.0 * (u.angstrom**3)
        # Per species-per-box
        assert len(moveset.max_translate[0]) == 1
        assert len(moveset.max_rotate[0]) == 1
        assert len(moveset.max_translate[1]) == 1
        assert len(moveset.max_rotate[1]) == 1
        # Per species attributes
        assert len(moveset.insertable) == 1
        assert len(moveset.prob_swap_species) == 1
        assert len(moveset.prob_regrow_species) == 1
        # Should be insertable and regrow-able
        assert moveset.insertable[0] == True
        assert moveset.prob_regrow_species[0] == 1.0

    @pytest.mark.parametrize(
        "typ,value",
        [
            ("slitpore", 1 * u.angstrom),
            ("cylinder", 1 * u.angstrom),
            ("sphere", 1 * u.angstrom),
            ("interface", [1 * u.angstrom, 2 * u.angstrom]),
        ],
    )
    def test_restricted_gemc(self, methane_oplsaa, typ, value):
        moveset = mc.MoveSet("gemc", [methane_oplsaa])
        moveset.add_restricted_insertions(
            [methane_oplsaa], [[None], [typ]], [[None], [value]]
        )

        assert moveset._restricted_type == [[None], [typ]]
        assert moveset._restricted_value == [[None], [value]]

    def test_ensemble_gemcnpt(self, methane_oplsaa):
        moveset = mc.MoveSet("gemc_npt", [methane_oplsaa])
        assert moveset.ensemble == "gemc_npt"
        assert moveset.prob_translate == 0.30
        assert moveset.prob_rotate == 0.30
        assert moveset.prob_regrow == 0.295
        assert moveset.prob_volume == 0.005
        assert moveset.prob_angle == 0.0
        assert moveset.prob_dihedral == 0.0
        assert moveset.prob_insert == 0.0
        assert moveset.prob_swap == 0.1

        # Per box attributes
        assert len(moveset.max_translate) == 2
        assert len(moveset.max_rotate) == 2
        assert len(moveset.prob_swap_from_box) == 2
        assert moveset.prob_swap_from_box[0] == 0.5
        assert moveset.prob_swap_from_box[1] == 0.5
        assert len(moveset.max_volume) == 2
        assert moveset.max_volume[0] == 500.0 * (u.angstrom**3)
        assert moveset.max_volume[1] == 5000.0 * (u.angstrom**3)
        # Per species-per-box
        assert len(moveset.max_translate[0]) == 1
        assert len(moveset.max_rotate[0]) == 1
        assert len(moveset.max_translate[1]) == 1
        assert len(moveset.max_rotate[1]) == 1
        # Per species attributes
        assert len(moveset.insertable) == 1
        assert len(moveset.prob_swap_species) == 1
        assert len(moveset.prob_regrow_species) == 1
        # Should be insertable and regrow-able
        assert moveset.insertable[0] == True
        assert moveset.prob_regrow_species[0] == 1.0

    def test_restricted_gemc_npt(self, methane_oplsaa):
        moveset = mc.MoveSet("gemc_npt", [methane_oplsaa])
        moveset.add_restricted_insertions(
            [methane_oplsaa],
            [[None], ["slitpore"]],
            [[None], [3 * u.angstrom]],
        )

        assert moveset._restricted_type == [[None], ["slitpore"]]
        assert moveset._restricted_value == [[None], [3 * u.angstrom]]

    def test_single_site_nvt(self, methane_trappe):

        moveset = mc.MoveSet("nvt", [methane_trappe])
        assert moveset.ensemble == "nvt"
        assert moveset.prob_rotate == 0.0
        assert moveset.prob_regrow == 0.0
        assert moveset.prob_translate == 1.0
        assert moveset.prob_volume == 0.0
        assert moveset.prob_angle == 0.0
        assert moveset.prob_dihedral == 0.0
        assert moveset.prob_insert == 0.0
        assert moveset.prob_swap == 0.0

        # Per box attributes
        assert len(moveset.max_translate) == 1
        assert len(moveset.max_rotate) == 1
        assert len(moveset.max_volume) == 1
        assert len(moveset.prob_swap_from_box) == 1
        assert moveset.prob_swap_from_box[0] == 1.0
        assert moveset.max_volume[0] == 0.0
        # Per species-per-box
        assert len(moveset.max_translate[0]) == 1
        assert len(moveset.max_rotate[0]) == 1
        # Per species attributes
        assert len(moveset.insertable) == 1
        assert len(moveset.prob_swap_species) == 1
        assert len(moveset.prob_regrow_species) == 1
        # Should NOT be insertable or regrow-able
        assert moveset.insertable[0] == False
        assert moveset.prob_regrow_species[0] == 0.0

    def test_single_site_gemc(self, methane_trappe):

        moveset = mc.MoveSet("gemc", [methane_trappe])
        assert moveset.ensemble == "gemc"
        assert moveset.prob_rotate == 0.0
        assert moveset.prob_regrow == 0.0
        assert moveset.prob_translate == pytest.approx(0.895)
        assert moveset.prob_volume == 0.005
        assert moveset.prob_angle == 0.0
        assert moveset.prob_dihedral == 0.0
        assert moveset.prob_insert == 0.0
        assert moveset.prob_swap == 0.1

        # Per box attributes
        assert len(moveset.max_translate) == 2
        assert len(moveset.max_rotate) == 2
        assert len(moveset.max_volume) == 1
        assert len(moveset.prob_swap_from_box) == 2
        assert moveset.prob_swap_from_box[0] == 0.5
        assert moveset.prob_swap_from_box[1] == 0.5
        assert moveset.max_volume[0] == 500.0
        # Per species-per-box
        assert len(moveset.max_translate[0]) == 1
        assert len(moveset.max_translate[1]) == 1
        assert len(moveset.max_rotate[0]) == 1
        assert len(moveset.max_rotate[1]) == 1
        # Per species attributes
        assert len(moveset.insertable) == 1
        assert len(moveset.prob_swap_species) == 1
        assert len(moveset.prob_regrow_species) == 1
        # Should be insertable and NOT regrow-able
        assert moveset.insertable[0] == True
        assert moveset.prob_regrow_species[0] == 0.0

    def test_gcmc_lattice(self, fixed_lattice_trappe, methane_trappe):

        moveset = mc.MoveSet("gcmc", [fixed_lattice_trappe, methane_trappe])
        assert moveset.ensemble == "gcmc"
        assert moveset.prob_translate == pytest.approx(0.8)
        assert moveset.prob_rotate == 0.0
        assert moveset.prob_regrow == 0.0
        assert moveset.prob_volume == 0.0
        assert moveset.prob_angle == 0.0
        assert moveset.prob_dihedral == 0.0
        assert moveset.prob_insert == 0.1
        assert moveset.prob_swap == pytest.approx(0.0)

        # Per box attributes
        assert len(moveset.max_translate) == 1
        assert len(moveset.max_rotate) == 1
        assert len(moveset.max_volume) == 1
        assert len(moveset.prob_swap_from_box) == 1
        assert moveset.prob_swap_from_box[0] == 1.0
        assert moveset.max_volume[0] == 0.0
        # Per species-per-box
        assert len(moveset.max_translate[0]) == 2
        assert len(moveset.max_rotate[0]) == 2
        # Per species attributes
        assert len(moveset.insertable) == 2
        assert len(moveset.prob_swap_species) == 2
        assert len(moveset.prob_regrow_species) == 2
        # Lattice should not be insertable or regrow-able
        assert moveset.insertable[0] == False
        assert moveset.prob_regrow_species[0] == 0.0
        # Methane should be insertable and not regrow-able
        assert moveset.insertable[1] == True
        assert moveset.prob_regrow_species[1] == 0.0

    def test_prob_setters(self, methane_oplsaa):
        moveset = mc.MoveSet("nvt", [methane_oplsaa])
        with pytest.raises(TypeError, match=r"prob_translate must"):
            moveset.prob_translate = []
        with pytest.raises(TypeError, match=r"prob_translate must"):
            moveset.prob_translate = True
        with pytest.raises(ValueError, match=r"be between 0.0 and 1.0"):
            moveset.prob_translate = -0.2

        with pytest.raises(TypeError, match=r"prob_rotate must"):
            moveset.prob_rotate = []
        with pytest.raises(TypeError, match=r"prob_rotate must"):
            moveset.prob_rotate = True
        with pytest.raises(ValueError, match=r"be between 0.0 and 1.0"):
            moveset.prob_rotate = -0.2

        with pytest.raises(TypeError, match=r"prob_angle must"):
            moveset.prob_angle = []
        with pytest.raises(TypeError, match=r"prob_angle must"):
            moveset.prob_angle = True
        with pytest.raises(ValueError, match=r"be between 0.0 and 1.0"):
            moveset.prob_angle = -0.2

        with pytest.raises(TypeError, match=r"prob_dihedral must"):
            moveset.prob_dihedral = []
        with pytest.raises(TypeError, match=r"prob_dihedral must"):
            moveset.prob_dihedral = True
        with pytest.raises(ValueError, match=r"be between 0.0 and 1.0"):
            moveset.prob_dihedral = -0.2

        with pytest.raises(TypeError, match=r"prob_regrow must"):
            moveset.prob_regrow = []
        with pytest.raises(TypeError, match=r"prob_regrow must"):
            moveset.prob_regrow = True
        with pytest.raises(ValueError, match=r"be between 0.0 and 1.0"):
            moveset.prob_regrow = -0.2

        with pytest.raises(TypeError, match=r"prob_volume must"):
            moveset.prob_volume = []
        with pytest.raises(TypeError, match=r"prob_volume must"):
            moveset.prob_volume = True
        with pytest.raises(ValueError, match=r"be between 0.0 and 1.0"):
            moveset.prob_volume = -0.2
        with pytest.raises(ValueError, match=r"Ensemble is nvt."):
            moveset.prob_volume = 0.02

        moveset = mc.MoveSet("gemc", [methane_oplsaa])
        with pytest.raises(ValueError, match=r"Ensemble is gemc."):
            moveset.prob_volume = 0.0

        moveset = mc.MoveSet("nvt", [methane_oplsaa])
        with pytest.raises(TypeError, match=r"prob_insert must"):
            moveset.prob_insert = []
        with pytest.raises(TypeError, match=r"prob_insert must"):
            moveset.prob_insert = True
        with pytest.raises(ValueError, match=r"be between 0.0 and 1.0"):
            moveset.prob_insert = -0.2
        with pytest.raises(ValueError, match=r"Ensemble is nvt."):
            moveset.prob_insert = 0.2
        moveset = mc.MoveSet("gcmc", [methane_oplsaa])
        with pytest.raises(ValueError, match=r"Ensemble is gcmc."):
            moveset.prob_insert = 0.0

        moveset = mc.MoveSet("nvt", [methane_oplsaa])
        with pytest.raises(TypeError, match=r"prob_swap must"):
            moveset.prob_swap = []
        with pytest.raises(TypeError, match=r"prob_swap must"):
            moveset.prob_swap = True
        with pytest.raises(ValueError, match=r"be between 0.0 and 1.0"):
            moveset.prob_swap = -0.2
        with pytest.raises(ValueError, match=r"Ensemble is nvt."):
            moveset.prob_swap = 0.2
        moveset = mc.MoveSet("gemc", [methane_oplsaa])
        with pytest.raises(ValueError, match=r"Ensemble is gemc."):
            moveset.prob_swap = 0.0

    def test_maxval_setters(self, methane_oplsaa):
        moveset = mc.MoveSet("gemc", [methane_oplsaa])
        with pytest.raises(TypeError, match=r"must be a list"):
            moveset.max_translate = 1.0 * u.angstrom
        with pytest.raises(TypeError, match=r"must be a list"):
            moveset.max_translate = [
                [1.0 * u.angstrom, 1.0 * u.angstrom],
                [1.0 * u.angstrom],
            ]
        with pytest.raises(TypeError, match=r"unyt_array"):
            moveset.max_translate = [[1.0], [True]]
        with pytest.raises(ValueError, match=r"cannot be less than zero"):
            moveset.max_translate = [[1.0 * u.angstrom], [-1.0 * u.angstrom]]
        moveset.max_translate = [[1.0 * u.angstrom], [1 * u.angstrom]]
        assert moveset.max_translate[1][0] == 1.0 * u.angstrom

        with pytest.raises(TypeError, match=r"must be a list"):
            moveset.max_rotate = 1.0 * u.degree
        with pytest.raises(TypeError, match=r"must be a list"):
            moveset.max_rotate = [
                [1.0 * u.degree, 1.0 * u.degree],
                [1.0 * u.degree],
            ]
        with pytest.raises(TypeError, match=r"unyt_array"):
            moveset.max_rotate = [[1.0], [True]]
        with pytest.raises(ValueError, match=r"must be between"):
            moveset.max_rotate = [[1.0 * u.degree], [-1.0 * u.degree]]
        moveset.max_rotate = [[1.0 * u.degree], [1 * u.degree]]
        assert moveset.max_rotate[1][0] == 1.0 * u.degree

        with pytest.raises(TypeError, match=r"must be a list"):
            moveset.max_dihedral = 1.0 * u.degree
        with pytest.raises(TypeError, match=r"must be a list"):
            moveset.max_dihedral = [
                1.0 * u.degree,
                1.0 * u.degree,
                1.0 * u.degree,
            ]
        with pytest.raises(TypeError, match=r"unyt_array"):
            moveset.max_dihedral = [True]
        with pytest.raises(ValueError, match=r"must be between"):
            moveset.max_dihedral = [-1.0 * u.degree]
        with pytest.raises(ValueError, match=r"must be between"):
            moveset.max_dihedral = [370.0 * u.degree]
        moveset.max_dihedral = [1.0 * u.degree]
        assert moveset.max_dihedral[0] == 1.0 * u.degree

        with pytest.raises(TypeError, match=r"must be a list"):
            moveset.prob_swap_from_box = 1.0
        with pytest.raises(TypeError, match=r"must be a list"):
            moveset.prob_swap_from_box = [1.0, 1.0, 1.0]
        with pytest.raises(TypeError, match=r"of type float"):
            moveset.prob_swap_from_box = [0.5, True]
        with pytest.raises(ValueError, match=r"must be between"):
            moveset.prob_swap_from_box = [0.5, -0.5]
        moveset.prob_swap_from_box = [0.5, 0.5]
        assert moveset.prob_swap_from_box[0] == 0.5

        moveset = mc.MoveSet("gemc", [methane_oplsaa])
        with pytest.raises(TypeError, match=r"must be a list"):
            moveset.max_volume = 1.0
        with pytest.raises(TypeError, match=r"must be a list"):
            moveset.max_volume = [
                1.0 * (u.angstrom**3),
                1.0 * (u.angstrom**3),
                1.0 * (u.angstrom**3),
            ]
        with pytest.raises(TypeError, match=r"unyt_array"):
            moveset.max_volume = [True]
        with pytest.raises(ValueError, match=r"cannot be less than zero"):
            moveset.max_volume = [-100 * (u.angstrom**3)]
        moveset.max_volume = [5000.0 * (u.angstrom**3)]
        assert moveset.max_volume[0] == 5000.0 * (u.angstrom**3)
        moveset = mc.MoveSet("gemc_npt", [methane_oplsaa])
        moveset.max_volume = 1.0 * (u.angstrom**3)
        assert len(moveset.max_volume) == 2
        assert moveset.max_volume[0] == 1.0 * u.angstrom**3
        assert moveset.max_volume[1] == 1.0 * u.angstrom**3
        with pytest.raises(TypeError, match=r"must be a list"):
            moveset.max_volume = [
                1.0 * (u.angstrom**3),
                1.0 * (u.angstrom**3),
                1.0 * (u.angstrom**3),
            ]
        with pytest.raises(TypeError, match=r"unyt_array"):
            moveset.max_volume = [True, 50000.0]
        with pytest.raises(ValueError, match=r"cannot be less than zero"):
            moveset.max_volume = [
                -100 * (u.angstrom**3),
                100.0 * (u.angstrom**3),
            ]
        moveset.max_volume = [
            5000.0 * (u.angstrom**3),
            50000.0 * (u.angstrom**3),
        ]
        assert moveset.max_volume[1] == 50000.0 * (u.angstrom**3)

        moveset = mc.MoveSet("gemc", [methane_oplsaa])
        with pytest.raises(TypeError, match=r"must be a list"):
            moveset.insertable = 1.0
        with pytest.raises(TypeError, match=r"must be a list"):
            moveset.insertable = [1.0, 1.0, 1.0]
        with pytest.raises(TypeError, match=r"as a boolean type"):
            moveset.insertable = [2.0]
        moveset.insertable = [True]
        assert moveset.insertable[0] == True

        with pytest.raises(TypeError, match=r"must be a list"):
            moveset.prob_swap_species = 1.0
        with pytest.raises(TypeError, match=r"must be a list"):
            moveset.prob_swap_species = [1.0, 1.0, 1.0]
        with pytest.raises(TypeError, match=r"of type float"):
            moveset.prob_swap_species = [True]
        with pytest.raises(ValueError, match=r"must be between"):
            moveset.prob_swap_species = [-1.0]
        moveset.prob_swap_species = [1.0]
        assert moveset.prob_swap_species[0] == 1.0

        with pytest.raises(TypeError, match=r"must be a list"):
            moveset.prob_regrow_species = 1.0
        with pytest.raises(TypeError, match=r"must be a list"):
            moveset.prob_regrow_species = [1.0, 1.0, 1.0]
        with pytest.raises(TypeError, match=r"of type float"):
            moveset.prob_regrow_species = [True]
        with pytest.raises(ValueError, match=r"must be between"):
            moveset.prob_regrow_species = [-1.0]
        moveset.prob_regrow_species = [1.0]
        assert moveset.prob_regrow_species[0] == 1.0

    def test_cbmc_setters(self, methane_oplsaa):
        moveset = mc.MoveSet("gemc", [methane_oplsaa])
        with pytest.raises(TypeError, match=r"must be of type int"):
            moveset.cbmc_n_insert = 12.2
        with pytest.raises(ValueError, match=r"must be greater than zero"):
            moveset.cbmc_n_insert = -2
        moveset.cbmc_n_insert = 20
        assert moveset.cbmc_n_insert == 20
        with pytest.raises(TypeError, match=r"must be of type int"):
            moveset.cbmc_n_dihed = 12.2
        with pytest.raises(ValueError, match=r"must be greater than zero"):
            moveset.cbmc_n_dihed = -2
        moveset.cbmc_n_dihed = 20
        assert moveset.cbmc_n_dihed == 20
        with pytest.raises(TypeError, match=r"unyt_array"):
            moveset.cbmc_rcut = [3.0 * u.angstrom]
        with pytest.raises(TypeError, match=r"unyt_array"):
            moveset.cbmc_rcut = 3.0
        with pytest.raises(ValueError, match=r"less than zero"):
            moveset.cbmc_rcut = [3.0 * u.angstrom, -3.0 * u.angstrom]
        with pytest.raises(TypeError, match=r"cbmc_rcut must be a list"):
            moveset.cbmc_rcut = [0.4 * u.nm, 8.0 * u.gram]
        moveset.cbmc_rcut = 5.0 * u.angstrom
        assert len(moveset.cbmc_rcut) == 2
        assert moveset.cbmc_rcut[0].to_value("angstrom") == 5.0
        assert moveset.cbmc_rcut[1].to_value("angstrom") == 5.0
        moveset = mc.MoveSet("nvt", [methane_oplsaa])
        moveset.cbmc_rcut = 7.0 * u.angstrom
        assert len(moveset.cbmc_rcut) == 1
        assert moveset.cbmc_rcut[0].to_value("angstrom") == 7.0

    def test_print_moveset(self, methane_oplsaa):
        """Simple test to make sure moveset object is printed"""

        moveset = mc.MoveSet("gemc", [methane_oplsaa])
        moveset.print()

    def test_invalid_ensemble_and_restriction(self, methane_oplsaa):
        moveset = mc.MoveSet("nvt", [methane_oplsaa])
        with pytest.raises(ValueError, match=r"only valid"):
            moveset.add_restricted_insertions(
                [methane_oplsaa], [["slitpore"]], [[1]]
            )

    @pytest.mark.parametrize(
        "typ,value",
        [
            ("slitpore", [[1], [2]]),
            ("cylinder", [[None]]),
            (None, [[1]]),
        ],
    )
    def test_value_error_restricted_type_and_value(
        self, methane_oplsaa, typ, value
    ):
        moveset = mc.MoveSet("gcmc", [methane_oplsaa])
        with pytest.raises(ValueError):
            moveset.add_restricted_insertions([methane_oplsaa], [[typ]], value)

    def test_type_error_restricted_type_and_value(self, methane_oplsaa):
        moveset = mc.MoveSet("gcmc", [methane_oplsaa])
        with pytest.raises(TypeError):
            moveset.add_restricted_insertions(
                [methane_oplsaa], ["cylinder"], [3]
            )

    def test_invalid_restricted_type_and_species(self, methane_oplsaa):
        moveset = mc.MoveSet("gcmc", [methane_oplsaa])
        with pytest.raises(ValueError, match=r"Length of 'species'"):
            moveset.add_restricted_insertions(
                [methane_oplsaa], [["slitpore", None]], [[1, None]]
            )

    @pytest.mark.skip(reason="The purpose of this test is not clear.")
    def test_add_multiple_restricted_insertions(self, methane_oplsaa):
        moveset = mc.MoveSet("gcmc", [methane_oplsaa])
        moveset.add_restricted_insertions(
            [methane_oplsaa], [["slitpore"]], [[3 * u.angstrom]]
        )
        with pytest.warns(None) as record:
            moveset.add_restricted_insertions(
                [methane_oplsaa], [["cylinder"]], [[3 * u.angstrom]]
            )

    def test_change_ensemble(self, methane_oplsaa):
        moveset = mc.MoveSet("gcmc", [methane_oplsaa])
        with pytest.raises(
            AttributeError, match=r"Ensemble cannot be changed"
        ):
            moveset.ensemble = "nvt"
