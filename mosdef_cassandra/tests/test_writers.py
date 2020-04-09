import pytest
from copy import deepcopy

import mosdef_cassandra as mc
from mosdef_cassandra.tests.base_test import BaseTest
from mosdef_cassandra.writers.inp_functions import generate_input


class TestInpFunctions(BaseTest):
    @pytest.fixture
    def onecomp_system(self, methane_oplsaa, box):
        system = mc.System([box], [methane_oplsaa], mols_to_add=[[10]])
        moves = mc.Moves("nvt", [methane_oplsaa])
        return system, moves

    @pytest.fixture
    def twocomp_system(self, methane_oplsaa, butane_oplsaa, box):
        system = mc.System(
            [box], [methane_oplsaa, butane_oplsaa], mols_to_add=[[10, 100]]
        )
        moves = mc.Moves("nvt", [methane_oplsaa, butane_oplsaa])
        return system, moves

    @pytest.fixture
    def twobox_system(self, methane_oplsaa, box):
        system = mc.System(
            [box, box], [methane_oplsaa], mols_to_add=[[10], [5]]
        )
        moves = mc.Moves("gemc", [methane_oplsaa])
        return system, moves

    @pytest.fixture
    def twocomptwobox_system(self, methane_oplsaa, butane_oplsaa, box):
        system = mc.System(
            [box, box],
            [methane_oplsaa, butane_oplsaa],
            mols_to_add=[[10, 100], [1, 5]],
        )
        moves = mc.Moves("gemc_npt", [methane_oplsaa, butane_oplsaa])
        return system, moves

    @pytest.fixture
    def gcmc_system(
        self, methane_oplsaa, fixed_lattice_compound, fixed_lattice_trappe
    ):
        box_list = [fixed_lattice_compound]
        species_list = [fixed_lattice_trappe, methane_oplsaa]
        system = mc.System(
            box_list,
            species_list,
            mols_in_boxes=[[1, 0]],
            mols_to_add=[[0, 10]],
        )
        moves = mc.Moves("gcmc", species_list)
        return system, moves

    def test_invalid_kwargs(self, onecomp_system):
        (system, moves) = onecomp_system
        with pytest.raises(ValueError, match=r"Invalid input argument"):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=500,
                temperature=300.0,
                random_arg=1,
            )

    def test_run_name(self, onecomp_system):
        (system, moves) = onecomp_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            run_name="test name",
        )

        assert "# Run_Name\ntest-name.out" in inp_data

        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            run_name="test_name",
        )

        assert "# Run_Name\ntest_name.out" in inp_data

        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
        )

        assert "# Run_Name\nnvt.out" in inp_data

        with pytest.raises(TypeError, match=r"must be a string"):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=500,
                temperature=300.0,
                run_name=1,
            )

    def test_sim_type(self, onecomp_system):
        (system, moves) = onecomp_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
        )

        assert "# Sim_Type\nnvt" in inp_data

        with pytest.raises(ValueError, match=r"Unsupported sim_type"):
            inp_data = mc.writers.inp_functions.get_sim_type("gccmc")

    def test_nbr_species(self, onecomp_system, twocomp_system):
        (system, moves) = onecomp_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
        )
        assert "# Nbr_Species\n1" in inp_data
        (system, moves) = twocomp_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
        )
        assert "# Nbr_Species\n2" in inp_data

    def test_vdw_style(self, twocomp_system, twobox_system):
        (system, moves) = twocomp_system

        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
        )
        assert "# VDW_Style\nlj cut_tail 12.0" in inp_data

        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            vdw_style="none",
        )
        assert "# VDW_Style\nnone\n" in inp_data

        with pytest.raises(ValueError, match=r"Unsupported vdw_style"):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=500,
                temperature=300.0,
                vdw_style="cutoff",
                vdw_cutoff=12.0,
            )

        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            cutoff_style="cut",
            vdw_cutoff=15.0,
        )
        assert "# VDW_Style\nlj cut 15.0" in inp_data

        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            cutoff_style="cut_shift",
            vdw_cutoff=15.0,
        )
        assert "# VDW_Style\nlj cut_shift 15.0" in inp_data

        with pytest.raises(ValueError, match=r"Only one box"):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=500,
                temperature=300.0,
                vdw_cutoff_box2=10.0,
            )

        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            cutoff_style="cut_switch",
            vdw_cutoff=[12.0, 15.0],
        )
        assert "# VDW_Style\nlj cut_switch 12.0 15.0" in inp_data

        with pytest.raises(ValueError, match=r"requires an inner"):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=500,
                temperature=300.0,
                cutoff_style="cut_switch",
                vdw_cutoff=12.0,
            )

        (system, moves) = twobox_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
        )
        assert "# VDW_Style\nlj cut_tail 12.0\nlj cut_tail 12.0" in inp_data

        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            cutoff_style="cut_switch",
            vdw_cutoff_box1=[12.0, 15.0],
            vdw_cutoff_box2=[11.0, 13.0],
        )
        assert (
            "# VDW_Style\nlj cut_switch 12.0 15.0\nlj cut_switch 11.0 13.0"
            in inp_data
        )

        with pytest.raises(ValueError, match=r"Unsupported cutoff style"):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=500,
                temperature=300.0,
                cutoff_style="cutoff",
                vdw_cutoff=12.0,
            )

    def test_charge_style(self, twocomp_system, twobox_system):
        (system, moves) = twocomp_system

        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
        )
        assert "# Charge_Style\ncoul ewald 12.0 1e-05\n" in inp_data

        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            charge_style="cut",
        )
        assert "# Charge_Style\ncoul cut 12.0\n" in inp_data

        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            charge_style="dsf",
        )
        assert "# Charge_Style\ncoul dsf 12.0\n" in inp_data

        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            charge_style="dsf",
            dsf_damping=0.2,
        )
        assert "# Charge_Style\ncoul dsf 12.0 0.2\n" in inp_data

        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            charge_style="none",
        )
        assert "# Charge_Style\nnone\n" in inp_data

        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            charge_cutoff=15.0,
            ewald_accuracy=5e-6,
        )
        assert "# Charge_Style\ncoul ewald 15.0 5e-06\n" in inp_data

        with pytest.raises(ValueError, match=r"Only one box"):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=500,
                temperature=300.0,
                charge_cutoff_box2=1.0,
            )

        (system, moves) = twobox_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            charge_cutoff_box2=30.0,
            ewald_accuracy=5e-6,
        )
        assert (
            "# Charge_Style\ncoul ewald 12.0 5e-06\ncoul ewald 30.0 5e-06\n"
            in inp_data
        )

    def test_mixing_rule(self, onecomp_system):
        (system, moves) = onecomp_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
        )
        assert "# Mixing_Rule\nlb\n" in inp_data

        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            mixing_rule="geometric",
        )
        assert "# Mixing_Rule\ngeometric\n" in inp_data

        mixing_dict = {"ls_138_s1 ls_140_s1": "1.0 1.0"}
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            mixing_rule="custom",
            custom_mixing_dict=mixing_dict,
        )
        assert (
            "# Mixing_Rule\ncustom\nls_138_s1 ls_140_s1 1.0 1.0\n" in inp_data
        )

        with pytest.raises(
            ValueError, match=r"Custom mixing rule requested but"
        ):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=500,
                temperature=300.0,
                mixing_rule="custom",
            )

        with pytest.raises(ValueError, match=r"Unsupported mixing rule"):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=500,
                temperature=300.0,
                mixing_rule="other",
            )

    def test_seeds(self, onecomp_system):
        (system, moves) = onecomp_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
        )

        assert "# Seed_Info\n" in inp_data

        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            seeds=[1, 2],
        )

        assert "# Seed_Info\n1 2\n" in inp_data

        with pytest.raises(TypeError, match=r"argument should be a list"):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=500,
                temperature=300.0,
                seeds=100,
            )

        with pytest.raises(ValueError, match=r"must be integers"):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=500,
                temperature=300.0,
                seeds=[100, -1],
            )

    def test_rcut_min(self, onecomp_system):
        (system, moves) = onecomp_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
        )

        assert "# Rcutoff_Low\n1.0\n" in inp_data

        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            rcut_min=10.0,
        )

        assert "# Rcutoff_Low\n10.0\n" in inp_data

        with pytest.raises(TypeError, match=r"be of type float"):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=500,
                temperature=300.0,
                rcut_min="hello",
            )

    def test_pair_energy(self, onecomp_system):
        (system, moves) = onecomp_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            rcut_min=10.0,
        )

        assert "# Pair_Energy\ntrue\n" in inp_data

        with pytest.raises(TypeError, match=r"be of type boolean"):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=500,
                temperature=300.0,
                pair_energy=1,
            )

    def test_max_molecules(self, twocomp_system, gcmc_system):
        (system, moves) = twocomp_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
        )

        assert (
            "# Molecule_Files\nspecies1.mcf 10\nspecies2.mcf 100" in inp_data
        )

        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            max_molecules=[100, 1000],
        )

        assert (
            "# Molecule_Files\nspecies1.mcf 100\nspecies2.mcf 1000" in inp_data
        )

        (system, moves) = gcmc_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            chemical_potentials=["none", 10.0],
        )

        assert (
            "# Molecule_Files\nspecies1.mcf 1\nspecies2.mcf 510\n" in inp_data
        )

        (system, moves) = twocomp_system
        with pytest.raises(TypeError, match=r"should be a list"):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=500,
                temperature=300.0,
                max_molecules=100,
            )

        with pytest.raises(ValueError, match=r"Length of list specified"):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=500,
                temperature=300.0,
                max_molecules=[100],
            )

    def test_boxes(self, onecomp_system, twobox_system, gcmc_system):
        (system, moves) = onecomp_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
        )

        assert "# Box_Info\n1\ncubic\n50.0\n" in inp_data

        (system, moves) = twobox_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
        )

        assert "# Box_Info\n2\ncubic\n50.0\n\ncubic\n50.0\n" in inp_data

        (system, moves) = gcmc_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            chemical_potentials=["none", 10.0],
        )

        assert "# Box_Info\n1\ncubic\n29.84\n" in inp_data

    def test_temperature(self, onecomp_system, twobox_system):
        (system, moves) = onecomp_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=200.0,
        )

        assert "# Temperature_Info\n200.0\n" in inp_data

        (system, moves) = twobox_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=200.0,
        )

        assert "# Temperature_Info\n200.0\n200.0\n" in inp_data

        with pytest.raises(ValueError, match=r"less than zero"):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=500,
                temperature=-300.0,
            )

        with pytest.raises(TypeError, match=r"of type float"):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=500,
                temperature="hi",
            )

    def test_pressure(self, twocomptwobox_system):
        (system, moves) = twocomptwobox_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            pressure=2.0,
        )

        assert "# Pressure_Info\n2.0\n2.0\n" in inp_data

        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            pressure=2.0,
            pressure_box2=10.0,
        )

        assert "# Pressure_Info\n2.0\n10.0\n" in inp_data

        with pytest.raises(ValueError, match=r"Pressure must be specified"):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=500,
                temperature=300.0,
            )

        with pytest.raises(TypeError, match=r"of type float"):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=500,
                temperature=300.0,
                pressure="string",
            )

    def test_chempot(self, gcmc_system):
        (system, moves) = gcmc_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            chemical_potentials=["none", 10.0],
        )

        assert "# Chemical_Potential_Info\nnone 10.0 \n" in inp_data

        with pytest.raises(
            ValueError, match=r"Chemical potential information"
        ):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=500,
                temperature=300.0,
            )

        with pytest.raises(TypeError, match=r"of type float"):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=500,
                temperature=300.0,
                chemical_potentials=["none", "string"],
            )

    def test_moves_formatting(self, onecomp_system):
        # Invalid keyword
        with pytest.raises(
            ValueError, match="Invalid probability info section"
        ):
            fake_prob_dict = {"trans": "test"}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        # Translate
        with pytest.raises(TypeError, match="not formatted properly"):
            fake_prob_dict = {"translate": "test"}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        with pytest.raises(TypeError, match="not formatted properly"):
            fake_prob_dict = {"translate": [0.1, 1.0]}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        with pytest.raises(TypeError, match="not formatted properly"):
            fake_prob_dict = {"translate": [0.1, ["test"]]}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        fake_prob_dict = {"translate": [0.1, [5.0]]}
        inp_data = mc.writers.inp_functions.get_move_probability_info(
            **fake_prob_dict
        )
        # Rotate
        with pytest.raises(TypeError, match="not formatted properly"):
            fake_prob_dict = {"rotate": "test"}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        with pytest.raises(TypeError, match="not formatted properly"):
            fake_prob_dict = {"rotate": [0.1, 1.0]}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        with pytest.raises(TypeError, match="not formatted properly"):
            fake_prob_dict = {"rotate": [0.1, ["test"]]}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        fake_prob_dict = {"rotate": [0.1, [5.0]]}
        inp_data = mc.writers.inp_functions.get_move_probability_info(
            **fake_prob_dict
        )
        # Angle
        with pytest.raises(TypeError, match="not formatted properly"):
            fake_prob_dict = {"angle": [14.0]}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        fake_prob_dict = {"angle": 14.0}
        inp_data = mc.writers.inp_functions.get_move_probability_info(
            **fake_prob_dict
        )
        # Dihedral
        with pytest.raises(TypeError, match="not formatted properly"):
            fake_prob_dict = {"dihed": "test"}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        with pytest.raises(TypeError, match="not formatted properly"):
            fake_prob_dict = {"dihed": [0.1, 1.0]}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        with pytest.raises(TypeError, match="not formatted properly"):
            fake_prob_dict = {"dihed": [0.1, ["test"]]}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        fake_prob_dict = {"dihed": [0.1, [5.0]]}
        inp_data = mc.writers.inp_functions.get_move_probability_info(
            **fake_prob_dict
        )
        # Regrow
        with pytest.raises(TypeError, match="not formatted properly"):
            fake_prob_dict = {"regrow": "test"}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        with pytest.raises(TypeError, match="not formatted properly"):
            fake_prob_dict = {"regrow": ["test", 0.1, 0.2]}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        with pytest.raises(TypeError, match="must be a floating"):
            fake_prob_dict = {"regrow": ["test", [1.0]]}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        with pytest.raises(TypeError, match="not formatted properly"):
            fake_prob_dict = {"regrow": [0.3, 1.0]}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        with pytest.raises(TypeError, match="must be a floating"):
            fake_prob_dict = {"regrow": [0.3, ["string"]]}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        fake_prob_dict = {"regrow": [0.3, [1.0]]}
        # Vol
        inp_data = mc.writers.inp_functions.get_move_probability_info(
            **fake_prob_dict
        )
        with pytest.raises(TypeError, match="not formatted properly"):
            fake_prob_dict = {"volume": "test"}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        with pytest.raises(TypeError, match="not formatted properly"):
            fake_prob_dict = {"volume": [0.1, 100.0, 0.2]}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        with pytest.raises(TypeError, match="must be a floating point"):
            fake_prob_dict = {"volume": ["test", [100.0]]}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        with pytest.raises(TypeError, match="not formatted properly"):
            fake_prob_dict = {"volume": [0.1, 100.0]}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        with pytest.raises(TypeError, match="must be a floating point"):
            fake_prob_dict = {"volume": [0.1, ["test"]]}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        fake_prob_dict = {"volume": [0.1, [100.0]]}
        inp_data = mc.writers.inp_functions.get_move_probability_info(
            **fake_prob_dict
        )
        # Insertable
        with pytest.raises(TypeError, match="not formatted properly"):
            fake_prob_dict = {"insert": "test"}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        with pytest.raises(TypeError, match="not formatted properly"):
            fake_prob_dict = {"insert": [0.1, True, True]}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        with pytest.raises(TypeError, match="must be a floating point"):
            fake_prob_dict = {"insert": ["test", [True]]}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        with pytest.raises(TypeError, match="not formatted properly"):
            fake_prob_dict = {"insert": [0.1, True]}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        with pytest.raises(TypeError, match="must be a boolean value"):
            fake_prob_dict = {"insert": [0.1, [1.0]]}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        fake_prob_dict = {"insert": [0.1, [True]]}
        inp_data = mc.writers.inp_functions.get_move_probability_info(
            **fake_prob_dict
        )
        # Swap
        with pytest.raises(TypeError, match="not formatted properly"):
            fake_prob_dict = {"swap": "test"}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        with pytest.raises(TypeError, match="not formatted properly"):
            fake_prob_dict = {"swap": [0.1, [True], [0.5]]}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        with pytest.raises(TypeError, match="must be a floating point"):
            fake_prob_dict = {"swap": ["test", [True], [0.5], [1.0]]}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        with pytest.raises(TypeError, match="not formatted properly"):
            fake_prob_dict = {"swap": [0.1, True, [0.5], [1.0]]}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        with pytest.raises(TypeError, match="must be a boolean value"):
            fake_prob_dict = {"swap": [0.1, [1.0], [0.5], [1.0]]}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        with pytest.raises(TypeError, match="not formatted properly"):
            fake_prob_dict = {"swap": [0.1, [True], 0.5, [1.0]]}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        with pytest.raises(TypeError, match="must be a floating point"):
            fake_prob_dict = {"swap": [0.1, [True], ["test"], [1.0]]}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        with pytest.raises(TypeError, match="not formatted properly"):
            fake_prob_dict = {"swap": [0.1, [True], [0.5], 1.0]}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )
        with pytest.raises(TypeError, match="must be a floating point"):
            fake_prob_dict = {"swap": [0.1, [True], [0.5], ["test"]]}
            inp_data = mc.writers.inp_functions.get_move_probability_info(
                **fake_prob_dict
            )

        fake_prob_dict = {"swap": [0.1, [True], [0.5], [1.0]]}
        inp_data = mc.writers.inp_functions.get_move_probability_info(
            **fake_prob_dict
        )
        fake_prob_dict = {"swap": [0.1, [True], [0.5], None]}
        inp_data = mc.writers.inp_functions.get_move_probability_info(
            **fake_prob_dict
        )
        fake_prob_dict = {"swap": [0.1, [True], None, None]}
        inp_data = mc.writers.inp_functions.get_move_probability_info(
            **fake_prob_dict
        )

    def test_moves_onecomp(self, onecomp_system):

        (system, moves) = onecomp_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
        )

        assert "# Move_Probability_Info" in inp_data
        assert "# Done_Probability_Info" in inp_data
        assert "# Prob_Translation\n0.35\n2.0 \n" in inp_data
        assert "# Prob_Rotation\n0.35\n30.0 \n" in inp_data
        assert "# Prob_Angle" not in inp_data
        assert "# Prob_Dihedral" not in inp_data
        assert "# Prob_Regrowth\n0.3\n1.0 \n" in inp_data
        assert "# Prob_Volume" not in inp_data
        assert "# Prob_Insertion" not in inp_data
        assert "# Prob_Deletion" not in inp_data
        assert "# Prob_Swap" not in inp_data
        assert "# Prob_Ring" not in inp_data

        moves.prob_angle = 0.1
        moves.prob_translate = 0.3
        moves.prob_rotate = 0.3
        moves.max_translate[0][0] = 10.0
        moves.max_rotate[0][0] = 10.0

        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
        )

        assert "# Move_Probability_Info" in inp_data
        assert "# Done_Probability_Info" in inp_data
        assert "# Prob_Translation\n0.3\n10.0 \n" in inp_data
        assert "# Prob_Rotation\n0.3\n10.0 \n" in inp_data
        assert "# Prob_Angle\n0.1\n" in inp_data
        assert "# Prob_Dihedral" not in inp_data
        assert "# Prob_Regrowth\n0.3\n1.0 \n" in inp_data
        assert "# Prob_Volume" not in inp_data
        assert "# Prob_Insertion" not in inp_data
        assert "# Prob_Deletion" not in inp_data
        assert "# Prob_Swap" not in inp_data
        assert "# Prob_Ring" not in inp_data

    def test_moves_twocomp(self, twocomp_system):

        (system, moves) = twocomp_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
        )

        assert "# Move_Probability_Info" in inp_data
        assert "# Done_Probability_Info" in inp_data
        assert "# Prob_Translation\n0.35\n2.0 2.0 \n" in inp_data
        assert "# Prob_Rotation\n0.35\n30.0 30.0 \n" in inp_data
        assert "# Prob_Angle" not in inp_data
        assert "# Prob_Dihedral" not in inp_data
        assert "# Prob_Regrowth\n0.3\n0.5 0.5 \n" in inp_data
        assert "# Prob_Volume" not in inp_data
        assert "# Prob_Insertion" not in inp_data
        assert "# Prob_Deletion" not in inp_data
        assert "# Prob_Swap" not in inp_data
        assert "# Prob_Ring" not in inp_data

        moves.prob_angle = 0.1
        moves.prob_translate = 0.3
        moves.prob_rotate = 0.3
        moves.max_translate[0][0] = 10.0
        moves.max_rotate[0][0] = 10.0

        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
        )

        assert "# Move_Probability_Info" in inp_data
        assert "# Done_Probability_Info" in inp_data
        assert "# Prob_Translation\n0.3\n10.0 2.0 \n" in inp_data
        assert "# Prob_Rotation\n0.3\n10.0 30.0 \n" in inp_data
        assert "# Prob_Angle\n0.1\n" in inp_data
        assert "# Prob_Dihedral" not in inp_data
        assert "# Prob_Regrowth\n0.3\n0.5 0.5 \n" in inp_data
        assert "# Prob_Volume" not in inp_data
        assert "# Prob_Insertion" not in inp_data
        assert "# Prob_Deletion" not in inp_data
        assert "# Prob_Swap" not in inp_data
        assert "# Prob_Ring" not in inp_data

    def test_moves_twobox(self, twobox_system):

        (system, moves) = twobox_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
        )

        assert "# Move_Probability_Info" in inp_data
        assert "# Done_Probability_Info" in inp_data
        assert "# Prob_Translation\n0.29\n2.0 \n2.0 \n" in inp_data
        assert "# Prob_Rotation\n0.29\n30.0 \n30.0 \n" in inp_data
        assert "# Prob_Angle" not in inp_data
        assert "# Prob_Dihedral" not in inp_data
        assert "# Prob_Regrowth\n0.3\n1.0 \n" in inp_data
        assert "# Prob_Volume\n0.02\n500.0\n" in inp_data
        assert "# Prob_Insertion" not in inp_data
        assert "# Prob_Deletion" not in inp_data
        assert (
            "# Prob_Swap\n0.1\ncbmc \nprob_swap_species 1.0 \nprob_swap_from_box 0.5 0.5 \n"
            in inp_data
        )
        assert "# Prob_Ring" not in inp_data

    def test_moves_twocomptwobox(self, twocomptwobox_system):

        (system, moves) = twocomptwobox_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            pressure=1.0,
        )
        assert "# Move_Probability_Info" in inp_data
        assert "# Done_Probability_Info" in inp_data
        assert "# Prob_Translation\n0.29\n2.0 2.0 \n2.0 2.0 \n" in inp_data
        assert "# Prob_Rotation\n0.29\n30.0 30.0 \n30.0 30.0 \n" in inp_data
        assert "# Prob_Angle" not in inp_data
        assert "# Prob_Dihedral" not in inp_data
        assert "# Prob_Regrowth\n0.3\n0.5 0.5 \n" in inp_data
        assert "# Prob_Volume\n0.02\n500.0\n5000.0\n" in inp_data
        assert "# Prob_Insertion" not in inp_data
        assert "# Prob_Deletion" not in inp_data
        assert (
            "# Prob_Swap\n0.1\ncbmc cbmc \nprob_swap_species 0.5 0.5 \nprob_swap_from_box 0.5 0.5 \n"
            in inp_data
        )
        assert "# Prob_Ring" not in inp_data

    def test_moves_gcmc(self, gcmc_system):

        (system, moves) = gcmc_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            chemical_potentials=["none", 1.0],
        )

        assert "# Move_Probability_Info" in inp_data
        assert "# Done_Probability_Info" in inp_data
        assert "# Prob_Translation\n0.25\n0.0 2.0 \n" in inp_data
        assert "# Prob_Rotation\n0.25\n0.0 30.0 \n" in inp_data
        assert "# Prob_Angle" not in inp_data
        assert "# Prob_Dihedral" not in inp_data
        assert "# Prob_Regrowth\n0.3\n0.0 1.0 \n" in inp_data
        assert "# Prob_Volume" not in inp_data
        assert "# Prob_Insertion\n0.1\nnone cbmc" in inp_data
        assert "# Prob_Deletion\n0.1\n" in inp_data
        assert "# Prob_Swap" not in inp_data
        assert "# Prob_Ring" not in inp_data

    def test_start_type(
        self,
        onecomp_system,
        twocomp_system,
        twobox_system,
        twocomptwobox_system,
        gcmc_system,
    ):

        (system, moves) = onecomp_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
        )

        assert "# Start_Type\nmake_config 10\n" in inp_data

        (system, moves) = twocomp_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
        )

        assert "# Start_Type\nmake_config 10 100\n" in inp_data

        (system, moves) = twobox_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
        )

        assert "# Start_Type\nmake_config 10\nmake_config 5\n" in inp_data

        (system, moves) = twocomptwobox_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            pressure=1.0,
        )
        assert (
            "# Start_Type\nmake_config 10 100\nmake_config 1 5\n" in inp_data
        )

        (system, moves) = gcmc_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            chemical_potentials=["none", 1.0],
        )
        assert "# Start_Type\nadd_to_config 1 0 box1.in.xyz 0 10\n" in inp_data

        # HACK to test read config
        system_copy = deepcopy(system)
        system_copy._mols_to_add = [[0, 0], [0, 0]]
        inp_data = generate_input(
            system=system_copy,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            chemical_potentials=["none", 1.0],
        )

        assert "# Start_Type\nread_config 1 0 box1.in.xyz\n" in inp_data

    def test_run_type(self, onecomp_system, twobox_system):
        (system, moves) = onecomp_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
        )
        assert "# Run_Type\nequilibration 1000 \n" in inp_data

        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="production",
            run_length=500,
            temperature=300.0,
        )
        assert "# Run_Type\nproduction 1000 \n" in inp_data
        with pytest.raises(ValueError, match=r"Invalid run type"):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="pro",
                run_length=500,
                temperature=300.0,
            )

        (system, moves) = twobox_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
        )

        assert "# Run_Type\nequilibration 1000 100\n" in inp_data

        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            thermal_stat_freq=100,
            vol_stat_freq=50,
        )

        assert "# Run_Type\nequilibration 100 50\n" in inp_data

        with pytest.raises(ValueError, match=r"must be an integer"):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=500,
                temperature=300.0,
                thermal_stat_freq=10.2,
                vol_stat_freq=50,
            )

        with pytest.raises(ValueError, match=r"must be an integer"):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=500,
                temperature=300.0,
                thermal_stat_freq=10,
                vol_stat_freq=1.2,
            )

    def test_length_info(self, onecomp_system):
        (system, moves) = onecomp_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
        )

        assert (
            "# Simulation_Length_Info\nunits steps\nprop_freq 500\ncoord_freq 5000\nrun 500"
            in inp_data
        )

        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            steps_per_sweep=10,
            units="sweeps",
        )

        assert (
            "# Simulation_Length_Info\nunits sweeps\nprop_freq 500\ncoord_freq 5000\nrun 500\nsteps_per_sweep 10\n"
            in inp_data
        )
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            block_avg_freq=10,
        )
        assert (
            "# Simulation_Length_Info\nunits steps\nprop_freq 500\ncoord_freq 5000\nrun 500\nblock_averages 10\n"
            in inp_data
        )

        with pytest.raises(ValueError, match=r"Invalid units"):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=500,
                temperature=300.0,
                units="stweeps",
            )
        with pytest.raises(ValueError, match=r"must be an integer"):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=500,
                temperature=300.0,
                prop_freq=1.2,
            )
        with pytest.raises(ValueError, match=r"must be an integer"):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=500,
                temperature=300.0,
                coord_freq=1.2,
            )
        with pytest.raises(ValueError, match=r"must be an integer"):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=5.2,
                temperature=300.0,
            )
        with pytest.raises(ValueError, match=r"must be an integer"):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=500,
                temperature=300.0,
                block_avg_freq=10.2,
            )
        with pytest.raises(ValueError, match=r"must be an integer"):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=500,
                temperature=300.0,
                steps_per_sweep=10.2,
            )

    def test_property_info(self, onecomp_system, twobox_system):
        (system, moves) = onecomp_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
        )

        assert (
            "# Property_Info 1\nenergy_total\nenergy_intra\nenergy_inter\nenthalpy\npressure\nvolume\nnmols\nmass_density\n"
            in inp_data
        )

        (system, moves) = twobox_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
        )

        assert (
            "# Property_Info 1\nenergy_total\nenergy_intra\nenergy_inter\nenthalpy\npressure\nvolume\nnmols\nmass_density\n\n# Property_Info 2\nenergy_total\nenergy_intra\nenergy_inter\nenthalpy\npressure\nvolume\nnmols\nmass_density\n"
            in inp_data
        )

        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            properties=["energy_total", "enthalpy", "density"],
        )

        assert (
            "# Property_Info 1\nenergy_total\nenthalpy\ndensity\n\n# Property_Info 2\nenergy_total\nenthalpy\ndensity\n"
            in inp_data
        )

        with pytest.raises(ValueError, match=r"Invalid property"):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=500,
                temperature=300.0,
                properties=["temperature"],
            )

    def test_fragment_files(self, onecomp_system):
        (system, moves) = onecomp_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
        )

        assert "# Fragment_Files\n" in inp_data

    def test_verbose_log(self, onecomp_system):
        (system, moves) = onecomp_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            verbose_log=True,
        )

        assert "# Verbose_Logfile\ntrue\n" in inp_data

        with pytest.raises(TypeError, match=r"Verbosity must be"):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=500,
                temperature=300.0,
                verbose_log="true",
            )

    def test_cbmc_info(self, onecomp_system, twobox_system):
        (system, moves) = onecomp_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
        )

        assert (
            "# CBMC_Info\nkappa_ins 10\nkappa_dih 10\nrcut_cbmc 6.0\n"
            in inp_data
        )

        (system, moves) = twobox_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
        )

        assert (
            "# CBMC_Info\nkappa_ins 10\nkappa_dih 10\nrcut_cbmc 6.0 6.0\n"
            in inp_data
        )

        (system, moves) = onecomp_system
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            cbmc_kappa_ins=2,
            cbmc_kappa_dih=5,
            cbmc_rcut=4.5,
        )

        assert (
            "# CBMC_Info\nkappa_ins 2\nkappa_dih 5\nrcut_cbmc 4.5\n"
            in inp_data
        )

        with pytest.raises(TypeError, match=r"must be an integer"):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=500,
                temperature=300.0,
                cbmc_kappa_ins=2.5,
                cbmc_kappa_dih=5,
                cbmc_rcut=4.5,
            )

        with pytest.raises(TypeError, match=r"must be an integer"):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=500,
                temperature=300.0,
                cbmc_kappa_ins=2,
                cbmc_kappa_dih=5.5,
                cbmc_rcut=4.5,
            )

        with pytest.raises(TypeError, match=r"must be a float"):
            inp_data = generate_input(
                system=system,
                moves=moves,
                run_type="equilibration",
                run_length=500,
                temperature=300.0,
                cbmc_kappa_ins=2,
                cbmc_kappa_dih=5,
                cbmc_rcut=[],
            )

    def test_write_restricted_gcmc(self, gcmc_system):
        (system, moves) = gcmc_system
        moves.restricted_type = [[None, 'sphere']]
        moves.restricted_value = [[None, 3]]
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            chemical_potentials=["none", 10.0],
        )

        assert "\nrestricted_insertion sphere 3\n" in inp_data

    def test_write_restricted_gemc_npt(self, twocomptwobox_system):
        (system, moves) = twocomptwobox_system
        moves.restricted_type = [[None, None], [None, 'sphere']]
        moves.restricted_value = [[None, None], [None, 3]]
        inp_data = generate_input(
            system=system,
            moves=moves,
            run_type="equilibration",
            run_length=500,
            temperature=300.0,
            pressure=1,
            chemical_potentials=["none", 10.0],
        )

        assert "\nrestricted_insertion sphere 3\n" in inp_data
