import pytest
from copy import deepcopy
from pathlib import Path

import mosdef_cassandra as mc
import unyt as u
from mosdef_cassandra.tests.base_test import BaseTest
from mosdef_cassandra.writers.inp_functions import generate_input
from mosdef_cassandra.writers.writers import _generate_restart_inp
from mosdef_cassandra.writers.writers import write_input
from mosdef_cassandra.writers.writers import write_restart_input
from mosdef_cassandra.writers.writers import write_mcfs
from mosdef_cassandra.utils.tempdir import *


class TestInpFunctions(BaseTest):

    @staticmethod
    def check_start_type_header(inp_contents):
        """Helper function to find the # Start_Type header index."""
        for idx, line in enumerate(inp_contents):
            if "# Start_Type" in line:
                return idx
        raise AssertionError("Missing '# Start_Type' header")
    
    @staticmethod
    def check_checkpoint_line(inp_contents, start_idx):
        """Helper function to check the checkpoint line format."""
        checkpoint_line = inp_contents[start_idx + 1].strip()
        parts = checkpoint_line.split()
        assert len(parts) == 2, "The line following '# Start_Type' should have exactly two entries"
        assert parts[0] == "checkpoint", "The first entry should be 'checkpoint'"

    @staticmethod
    def check_only_comments_or_whitespace(inp_contents, start_idx):
        """Helper function to check for only comments or whitespace until the next header."""
        for line in inp_contents[start_idx + 2:]:
            line = line.strip()
            if line.startswith("#"):
                break
            assert line == "" or line.startswith("!"), \
                "Only spaces or comments are allowed between the checkpoint line and the next header"

    @pytest.fixture
    def onecomp_system(self, methane_oplsaa, box):
        system = mc.System([box], [methane_oplsaa], mols_to_add=[[10]])
        moveset = mc.MoveSet("nvt", [methane_oplsaa])
        return system, moveset

    @pytest.fixture
    def twocomp_system(self, methane_oplsaa, butane_oplsaa, box):
        system = mc.System(
            [box], [methane_oplsaa, butane_oplsaa], mols_to_add=[[10, 100]]
        )
        moveset = mc.MoveSet("nvt", [methane_oplsaa, butane_oplsaa])
        return system, moveset

    @pytest.fixture
    def twobox_system(self, methane_oplsaa, box):
        system = mc.System(
            [box, box], [methane_oplsaa], mols_to_add=[[10], [5]]
        )
        moveset = mc.MoveSet("gemc", [methane_oplsaa])
        return system, moveset

    @pytest.fixture
    def twocomptwobox_system(self, methane_oplsaa, butane_oplsaa, box):
        system = mc.System(
            [box, box],
            [methane_oplsaa, butane_oplsaa],
            mols_to_add=[[10, 100], [1, 5]],
        )
        moveset = mc.MoveSet("gemc_npt", [methane_oplsaa, butane_oplsaa])
        return system, moveset

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
        moveset = mc.MoveSet("gcmc", species_list)
        return system, moveset

    def test_invalid_kwargs(self, onecomp_system):
        (system, moveset) = onecomp_system
        with pytest.raises(ValueError, match=r"Invalid input argument"):
            inp_data = generate_input(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=500,
                temperature=300.0 * u.K,
                random_arg=1,
            )

    def test_run_name(self, onecomp_system):
        (system, moveset) = onecomp_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            run_name="test name",
        )

        assert "# Run_Name\ntest_name.out" in inp_data

        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            run_name="test_name",
        )

        assert "# Run_Name\ntest_name.out" in inp_data

        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
        )

        assert "# Run_Name\nnvt.out" in inp_data

        with pytest.raises(TypeError, match=r"must be a string"):
            generate_input(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=500,
                temperature=300.0 * u.K,
                run_name=1,
            )
        with pytest.raises(ValueError):
            generate_input(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=500,
                temperature=300.0 * u.K,
                run_name="test-fail",
            )

    def test_sim_type(self, onecomp_system):
        (system, moveset) = onecomp_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
        )

        assert "# Sim_Type\nnvt" in inp_data

        with pytest.raises(ValueError, match=r"Unsupported sim_type"):
            inp_data = mc.writers.inp_functions.get_sim_type("gccmc")

    def test_nbr_species(self, onecomp_system, twocomp_system):
        (system, moveset) = onecomp_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
        )
        assert "# Nbr_Species\n1" in inp_data
        (system, moveset) = twocomp_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
        )
        assert "# Nbr_Species\n2" in inp_data

    def test_vdw_style(self, twocomp_system, twobox_system):
        (system, moveset) = twocomp_system

        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
        )
        assert "# VDW_Style\nlj cut_tail 12.0" in inp_data

        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            vdw_style="none",
        )
        assert "# VDW_Style\nnone\n" in inp_data

        with pytest.raises(ValueError, match=r"Unsupported vdw_style"):
            inp_data = generate_input(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=500,
                temperature=300.0 * u.K,
                vdw_style="cutoff",
                vdw_cutoff=12.0 * u.angstrom,
            )

        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            cutoff_style="cut",
            vdw_cutoff=15.0 * u.angstrom,
        )
        assert "# VDW_Style\nlj cut 15.0" in inp_data

        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            cutoff_style="cut_shift",
            vdw_cutoff=15.0 * u.angstrom,
        )
        assert "# VDW_Style\nlj cut_shift 15.0" in inp_data

        with pytest.raises(ValueError, match=r"Only one box"):
            inp_data = generate_input(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=500,
                temperature=300.0 * u.K,
                vdw_cutoff_box2=10.0 * u.angstrom,
            )

        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            cutoff_style="cut_switch",
            vdw_cutoff=[12.0 * u.angstrom, 15.0 * u.angstrom],
        )
        assert "# VDW_Style\nlj cut_switch 12.0 15.0" in inp_data

        with pytest.raises(ValueError, match=r"requires an inner"):
            inp_data = generate_input(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=500,
                temperature=300.0 * u.K,
                cutoff_style="cut_switch",
                vdw_cutoff=12.0 * u.angstrom,
            )

        (system, moveset) = twobox_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
        )
        assert "# VDW_Style\nlj cut_tail 12.0\nlj cut_tail 12.0" in inp_data

        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            cutoff_style="cut_switch",
            vdw_cutoff_box1=[12.0 * u.angstrom, 15.0 * u.angstrom],
            vdw_cutoff_box2=[11.0 * u.angstrom, 13.0 * u.angstrom],
        )
        assert (
            "# VDW_Style\nlj cut_switch 12.0 15.0\nlj cut_switch 11.0 13.0"
            in inp_data
        )

        with pytest.raises(ValueError, match=r"Unsupported cutoff style"):
            inp_data = generate_input(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=500,
                temperature=300.0 * u.K,
                cutoff_style="cutoff",
                vdw_cutoff=12.0 * u.angstrom,
            )

    def test_charge_style(self, twocomp_system, twobox_system):
        (system, moveset) = twocomp_system

        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
        )
        assert "# Charge_Style\ncoul ewald 12.0 1e-05\n" in inp_data

        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            charge_style="cut",
        )
        assert "# Charge_Style\ncoul cut 12.0\n" in inp_data

        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            charge_style="dsf",
        )
        assert "# Charge_Style\ncoul dsf 12.0\n" in inp_data

        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            charge_style="dsf",
            dsf_damping=0.2,
        )
        assert "# Charge_Style\ncoul dsf 12.0 0.2\n" in inp_data

        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            charge_style="none",
        )
        assert "# Charge_Style\nnone\n" in inp_data

        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            charge_cutoff=15.0 * u.angstrom,
            ewald_accuracy=5e-6,
        )
        assert "# Charge_Style\ncoul ewald 15.0 5e-06\n" in inp_data

        with pytest.raises(ValueError, match=r"Only one box"):
            inp_data = generate_input(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=500,
                temperature=300.0 * u.K,
                charge_cutoff_box2=1.0 * u.angstrom,
            )

        (system, moveset) = twobox_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            charge_cutoff_box2=30.0 * u.angstrom,
            ewald_accuracy=5e-6,
        )
        assert (
            "# Charge_Style\ncoul ewald 12.0 5e-06\ncoul ewald 30.0 5e-06\n"
            in inp_data
        )

    def test_mixing_rule(self, onecomp_system):
        (system, moveset) = onecomp_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
        )
        assert "# Mixing_Rule\nlb\n" in inp_data

        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            mixing_rule="geometric",
        )
        assert "# Mixing_Rule\ngeometric\n" in inp_data

        mixing_dict = {"ls_138_s1 ls_140_s1": "1.0 1.0"}
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
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
                moveset=moveset,
                run_type="equilibration",
                run_length=500,
                temperature=300.0 * u.K,
                mixing_rule="custom",
            )

        with pytest.raises(ValueError, match=r"Unsupported mixing rule"):
            inp_data = generate_input(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=500,
                temperature=300.0 * u.K,
                mixing_rule="other",
            )

    def test_seeds(self, onecomp_system):
        (system, moveset) = onecomp_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
        )

        assert "# Seed_Info\n" in inp_data

        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            seeds=[1, 2],
        )

        assert "# Seed_Info\n1 2\n" in inp_data

        with pytest.raises(TypeError, match=r"argument should be a list"):
            inp_data = generate_input(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=500,
                temperature=300.0 * u.K,
                seeds=100,
            )

        with pytest.raises(ValueError, match=r"must be integers"):
            inp_data = generate_input(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=500,
                temperature=300.0 * u.K,
                seeds=[100, -1],
            )

    def test_rcut_min(self, onecomp_system):
        (system, moveset) = onecomp_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
        )

        assert "# Rcutoff_Low\n1.0\n" in inp_data

        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            rcut_min=10.0 * u.angstrom,
        )

        assert "# Rcutoff_Low\n10.0\n" in inp_data

        with pytest.raises(TypeError, match=r"unyt_array"):
            inp_data = generate_input(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=500,
                temperature=300.0 * u.K,
                rcut_min="hello",
            )

    def test_pair_energy(self, onecomp_system):
        (system, moveset) = onecomp_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            rcut_min=10.0 * u.angstrom,
        )

        assert "# Pair_Energy\ntrue\n" in inp_data

        with pytest.raises(TypeError, match=r"be of type boolean"):
            inp_data = generate_input(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=500,
                temperature=300.0 * u.K,
                pair_energy=1,
            )

    def test_max_molecules(self, twocomp_system, gcmc_system):
        (system, moveset) = twocomp_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
        )

        assert (
            "# Molecule_Files\nspecies1.mcf 10\nspecies2.mcf 100" in inp_data
        )

        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            max_molecules=[100, 1000],
        )

        assert (
            "# Molecule_Files\nspecies1.mcf 100\nspecies2.mcf 1000" in inp_data
        )

        (system, moveset) = gcmc_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            chemical_potentials=["none", 10.0 * (u.kJ / u.mol)],
        )

        assert (
            "# Molecule_Files\nspecies1.mcf 1\nspecies2.mcf 2010\n" in inp_data
        )

        (system, moveset) = twocomp_system
        with pytest.raises(TypeError, match=r"should be a list"):
            inp_data = generate_input(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=500,
                temperature=300.0 * u.K,
                max_molecules=100,
            )

        with pytest.raises(ValueError, match=r"Length of list specified"):
            inp_data = generate_input(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=500,
                temperature=300.0 * u.K,
                max_molecules=[100],
            )

    def test_boxes(self, onecomp_system, twobox_system, gcmc_system):
        (system, moveset) = onecomp_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
        )

        assert "# Box_Info\n1\ncubic\n50.0\n" in inp_data

        (system, moveset) = twobox_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
        )

        assert "# Box_Info\n2\ncubic\n50.0\n\ncubic\n50.0\n" in inp_data

        (system, moveset) = gcmc_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            chemical_potentials=["none", 10.0 * (u.kJ / u.mol)],
        )

        assert "# Box_Info\n1\ncubic\n29.84\n" in inp_data

    def test_temperature(self, onecomp_system, twobox_system):
        (system, moveset) = onecomp_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=200.0 * u.K,
        )

        assert "# Temperature_Info\n200.0\n" in inp_data

        (system, moveset) = twobox_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=200.0 * u.K,
        )

        assert "# Temperature_Info\n200.0\n200.0\n" in inp_data

        with pytest.raises(ValueError, match=r"less than zero"):
            inp_data = generate_input(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=500,
                temperature=-300.0 * u.K,
            )

        with pytest.raises(TypeError, match=r"unyt_array"):
            inp_data = generate_input(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=500,
                temperature="hi",
            )

    def test_pressure(self, twocomptwobox_system):
        (system, moveset) = twocomptwobox_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            pressure=2.0 * u.bar,
        )

        assert "# Pressure_Info\n2.0\n2.0\n" in inp_data

        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            pressure=2.0 * u.bar,
            pressure_box2=10.0 * u.bar,
        )

        assert "# Pressure_Info\n2.0\n10.0\n" in inp_data

        with pytest.raises(ValueError, match=r"Pressure must be specified"):
            inp_data = generate_input(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=500,
                temperature=300.0 * u.K,
            )

        with pytest.raises(TypeError, match=r"unyt_array"):
            inp_data = generate_input(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=500,
                temperature=300.0 * u.K,
                pressure="string",
            )

    def test_chempot(self, gcmc_system):
        (system, moveset) = gcmc_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            chemical_potentials=["none", 10.0 * (u.kJ / u.mol)],
        )

        assert "# Chemical_Potential_Info\nnone 10.0 \n" in inp_data

        with pytest.raises(
            ValueError, match=r"Chemical potential information"
        ):
            inp_data = generate_input(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=500,
                temperature=300.0 * u.K,
            )

        with pytest.raises(TypeError, match=r"unyt_array"):
            inp_data = generate_input(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=500,
                temperature=300.0 * u.K,
                chemical_potentials=["none", "string"],
            )

    def test_moveset_formatting(self, onecomp_system):
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

    def test_moveset_onecomp(self, onecomp_system):

        (system, moveset) = onecomp_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
        )

        assert "# Move_Probability_Info" in inp_data
        assert "# Done_Probability_Info" in inp_data
        assert "# Prob_Translation\n0.33\n2.0 \n" in inp_data
        assert "# Prob_Rotation\n0.33\n30.0 \n" in inp_data
        assert "# Prob_Angle" not in inp_data
        assert "# Prob_Dihedral" not in inp_data
        assert "# Prob_Regrowth\n0.34\n1.0 \n" in inp_data
        assert "# Prob_Volume" not in inp_data
        assert "# Prob_Insertion" not in inp_data
        assert "# Prob_Deletion" not in inp_data
        assert "# Prob_Swap" not in inp_data
        assert "# Prob_Ring" not in inp_data

        moveset.prob_angle = 0.1
        moveset.prob_translate = 0.3
        moveset.prob_rotate = 0.3
        moveset.prob_regrow = 0.3
        moveset.max_translate[0][0] = 10.0 * u.angstrom
        moveset.max_rotate[0][0] = 10.0 * u.degree

        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
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

    def test_moveset_twocomp(self, twocomp_system):

        (system, moveset) = twocomp_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
        )

        assert "# Move_Probability_Info" in inp_data
        assert "# Done_Probability_Info" in inp_data
        assert "# Prob_Translation\n0.33\n2.0 2.0 \n" in inp_data
        assert "# Prob_Rotation\n0.33\n30.0 30.0 \n" in inp_data
        assert "# Prob_Angle" not in inp_data
        assert "# Prob_Dihedral" not in inp_data
        assert "# Prob_Regrowth\n0.34\n0.5 0.5 \n" in inp_data
        assert "# Prob_Volume" not in inp_data
        assert "# Prob_Insertion" not in inp_data
        assert "# Prob_Deletion" not in inp_data
        assert "# Prob_Swap" not in inp_data
        assert "# Prob_Ring" not in inp_data

        moveset.prob_angle = 0.1
        moveset.prob_translate = 0.3
        moveset.prob_rotate = 0.3
        moveset.prob_regrow = 0.26
        moveset.max_translate[0][0] = 10.0 * u.angstrom
        moveset.max_rotate[0][0] = 10.0 * u.degree

        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
        )

        assert "# Move_Probability_Info" in inp_data
        assert "# Done_Probability_Info" in inp_data
        assert "# Prob_Translation\n0.3\n10.0 2.0 \n" in inp_data
        assert "# Prob_Rotation\n0.3\n10.0 30.0 \n" in inp_data
        assert "# Prob_Angle\n0.1\n" in inp_data
        assert "# Prob_Dihedral" not in inp_data
        assert "# Prob_Regrowth\n0.26\n0.5 0.5 \n" in inp_data
        assert "# Prob_Volume" not in inp_data
        assert "# Prob_Insertion" not in inp_data
        assert "# Prob_Deletion" not in inp_data
        assert "# Prob_Swap" not in inp_data
        assert "# Prob_Ring" not in inp_data

    def test_moveset_twobox(self, twobox_system):

        (system, moveset) = twobox_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
        )

        assert "# Move_Probability_Info" in inp_data
        assert "# Done_Probability_Info" in inp_data
        assert "# Prob_Translation\n0.3\n2.0 \n2.0 \n" in inp_data
        assert "# Prob_Rotation\n0.3\n30.0 \n30.0 \n" in inp_data
        assert "# Prob_Angle" not in inp_data
        assert "# Prob_Dihedral" not in inp_data
        assert "# Prob_Regrowth\n0.295\n1.0 \n" in inp_data
        assert "# Prob_Volume\n0.005\n500.0\n" in inp_data
        assert "# Prob_Insertion" not in inp_data
        assert "# Prob_Deletion" not in inp_data
        assert (
            "# Prob_Swap\n0.1\ncbmc \nprob_swap_species 1.0 \nprob_swap_from_box 0.5 0.5 \n"
            in inp_data
        )
        assert "# Prob_Ring" not in inp_data

    def test_moveset_twocomptwobox(self, twocomptwobox_system):

        (system, moveset) = twocomptwobox_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            pressure=1.0 * u.bar,
        )
        assert "# Move_Probability_Info" in inp_data
        assert "# Done_Probability_Info" in inp_data
        assert "# Prob_Translation\n0.3\n2.0 2.0 \n2.0 2.0 \n" in inp_data
        assert "# Prob_Rotation\n0.3\n30.0 30.0 \n30.0 30.0 \n" in inp_data
        assert "# Prob_Angle" not in inp_data
        assert "# Prob_Dihedral" not in inp_data
        assert "# Prob_Regrowth\n0.295\n0.5 0.5 \n" in inp_data
        assert "# Prob_Volume\n0.005\n500.0\n5000.0\n" in inp_data
        assert "# Prob_Insertion" not in inp_data
        assert "# Prob_Deletion" not in inp_data
        assert (
            "# Prob_Swap\n0.1\ncbmc cbmc \nprob_swap_species 0.5 0.5 \nprob_swap_from_box 0.5 0.5 \n"
            in inp_data
        )
        assert "# Prob_Ring" not in inp_data

    def test_moveset_gcmc(self, gcmc_system):

        (system, moveset) = gcmc_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            chemical_potentials=["none", 1.0 * (u.kJ / u.mol)],
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

        (system, moveset) = onecomp_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
        )

        assert "# Start_Type\nmake_config 10\n" in inp_data

        (system, moveset) = twocomp_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
        )

        assert "# Start_Type\nmake_config 10 100\n" in inp_data

        (system, moveset) = twobox_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
        )

        assert "# Start_Type\nmake_config 10\nmake_config 5\n" in inp_data

        (system, moveset) = twocomptwobox_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            pressure=1.0 * u.bar,
        )
        assert (
            "# Start_Type\nmake_config 10 100\nmake_config 1 5\n" in inp_data
        )

        (system, moveset) = gcmc_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            chemical_potentials=["none", 1.0 * (u.kJ / u.mol)],
        )
        assert "# Start_Type\nadd_to_config 1 0 box1.in.xyz 0 10\n" in inp_data

        # HACK to test read config
        system_copy = deepcopy(system)
        system_copy._mols_to_add = [[0, 0], [0, 0]]
        inp_data = generate_input(
            system=system_copy,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            chemical_potentials=["none", 1.0 * (u.kJ / u.mol)],
        )

        assert "# Start_Type\nread_config 1 0 box1.in.xyz\n" in inp_data

    def test_run_type(self, onecomp_system, twobox_system):
        (system, moveset) = onecomp_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
        )
        assert "# Run_Type\nequilibration 1000 \n" in inp_data

        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="production",
            run_length=500,
            temperature=300.0 * u.K,
        )
        assert "# Run_Type\nproduction 1000 \n" in inp_data
        with pytest.raises(ValueError, match=r"Invalid run type"):
            inp_data = generate_input(
                system=system,
                moveset=moveset,
                run_type="pro",
                run_length=500,
                temperature=300.0 * u.K,
            )

        (system, moveset) = twobox_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
        )

        assert "# Run_Type\nequilibration 1000 100\n" in inp_data

        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            thermal_stat_freq=100,
            vol_stat_freq=50,
        )

        assert "# Run_Type\nequilibration 100 50\n" in inp_data

        with pytest.raises(ValueError, match=r"must be an integer"):
            inp_data = generate_input(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=500,
                temperature=300.0 * u.K,
                thermal_stat_freq=10.2,
                vol_stat_freq=50,
            )

        with pytest.raises(ValueError, match=r"must be an integer"):
            inp_data = generate_input(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=500,
                temperature=300.0 * u.K,
                thermal_stat_freq=10,
                vol_stat_freq=1.2,
            )

    def test_length_info(self, onecomp_system):
        (system, moveset) = onecomp_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
        )

        assert (
            "# Simulation_Length_Info\nunits steps\nprop_freq 500\ncoord_freq 5000\nrun 500"
            in inp_data
        )

        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            steps_per_sweep=10,
            units="sweeps",
        )

        assert (
            "# Simulation_Length_Info\nunits sweeps\nprop_freq 500\ncoord_freq 5000\nrun 500\nsteps_per_sweep 10\n"
            in inp_data
        )
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            block_avg_freq=10,
        )
        assert (
            "# Simulation_Length_Info\nunits steps\nprop_freq 500\ncoord_freq 5000\nrun 500\nblock_averages 10\n"
            in inp_data
        )

        with pytest.raises(ValueError, match=r"Invalid units"):
            inp_data = generate_input(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=500,
                temperature=300.0 * u.K,
                units="stweeps",
            )
        with pytest.raises(ValueError, match=r"must be an integer"):
            inp_data = generate_input(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=500,
                temperature=300.0 * u.K,
                prop_freq=1.2,
            )
        with pytest.raises(ValueError, match=r"must be an integer"):
            inp_data = generate_input(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=500,
                temperature=300.0 * u.K,
                coord_freq=1.2,
            )
        with pytest.raises(ValueError, match=r"must be an integer"):
            inp_data = generate_input(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=5.2,
                temperature=300.0 * u.K,
            )
        with pytest.raises(ValueError, match=r"must be an integer"):
            inp_data = generate_input(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=500,
                temperature=300.0 * u.K,
                block_avg_freq=10.2,
            )
        with pytest.raises(ValueError, match=r"must be an integer"):
            inp_data = generate_input(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=500,
                temperature=300.0 * u.K,
                steps_per_sweep=10.2,
            )

    def test_property_info(self, onecomp_system, twobox_system):
        (system, moveset) = onecomp_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
        )

        assert (
            "# Property_Info 1\nenergy_total\nenergy_intra\nenergy_inter\nenthalpy\npressure\nvolume\nnmols\nmass_density\n"
            in inp_data
        )

        (system, moveset) = twobox_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
        )

        assert (
            "# Property_Info 1\nenergy_total\nenergy_intra\nenergy_inter\nenthalpy\npressure\nvolume\nnmols\nmass_density\n\n# Property_Info 2\nenergy_total\nenergy_intra\nenergy_inter\nenthalpy\npressure\nvolume\nnmols\nmass_density\n"
            in inp_data
        )

        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            properties=["energy_total", "enthalpy", "density"],
        )

        assert (
            "# Property_Info 1\nenergy_total\nenthalpy\ndensity\n\n# Property_Info 2\nenergy_total\nenthalpy\ndensity\n"
            in inp_data
        )

        with pytest.raises(ValueError, match=r"Invalid property"):
            inp_data = generate_input(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=500,
                temperature=300.0 * u.K,
                properties=["temperature"],
            )

    def test_fragment_files(self, onecomp_system):
        (system, moveset) = onecomp_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
        )

        assert "# Fragment_Files\n" in inp_data

    def test_verbose_log(self, onecomp_system):
        (system, moveset) = onecomp_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            verbose_log=True,
        )

        assert "# Verbose_Logfile\ntrue\n" in inp_data

        with pytest.raises(TypeError, match=r"Verbosity must be"):
            inp_data = generate_input(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=500,
                temperature=300.0 * u.K,
                verbose_log="true",
            )

    def test_cbmc_info(self, onecomp_system, twobox_system):
        (system, moveset) = onecomp_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
        )

        assert (
            "# CBMC_Info\nkappa_ins 10\nkappa_dih 10\nrcut_cbmc 6.0\n"
            in inp_data
        )

        (system, moveset) = twobox_system
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
        )

        assert (
            "# CBMC_Info\nkappa_ins 10\nkappa_dih 10\nrcut_cbmc 6.0 6.0\n"
            in inp_data
        )

        (system, moveset) = onecomp_system
        moveset.cbmc_rcut = [0.45 * u.nm]
        moveset.cbmc_n_insert = 2
        moveset.cbmc_n_dihed = 5
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
        )
        print(inp_data)

        assert (
            "# CBMC_Info\nkappa_ins 2\nkappa_dih 5\nrcut_cbmc 4.5\n"
            in inp_data
        )

    @pytest.mark.parametrize(
        "typ,value",
        [
            ("slitpore", 1.0 * u.angstrom),
            ("cylinder", 1.0 * u.angstrom),
            ("sphere", 1.0 * u.angstrom),
            ("interface", [1.0 * u.angstrom, 2.0 * u.angstrom]),
        ],
    )
    def test_write_restricted_gcmc(self, gcmc_system, typ, value):
        (system, moveset) = gcmc_system
        moveset.add_restricted_insertions(
            system.species_topologies, [[None, typ]], [[None, value]]
        )
        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            chemical_potentials=["none", 10.0 * (u.kJ / u.mol)],
        )

        if typ == "interface":
            assert (
                "\nrestricted_insertion {} {:0.1f} {:0.1f}\n".format(
                    typ, value[0].to_value(), value[1].to_value()
                )
                in inp_data
            )
        else:
            assert (
                "\nrestricted_insertion {} {:0.1f}\n".format(
                    typ, value.to_value()
                )
                in inp_data
            )

    @pytest.mark.parametrize(
        "typ,value",
        [
            ("slitpore", 30 * u.angstrom),
            ("cylinder", 30 * u.angstrom),
            ("sphere", 30 * u.angstrom),
            ("interface", [30 * u.angstrom, 50 * u.angstrom]),
        ],
    )
    def test_fail_restricted_gcmc(self, gcmc_system, typ, value):
        (system, moveset) = gcmc_system
        moveset.add_restricted_insertions(
            system.species_topologies, [[None, typ]], [[None, value]]
        )
        with pytest.raises(ValueError, match=r"Restricted insertion"):
            inp_data = generate_input(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=500,
                temperature=300.0 * u.K,
                chemical_potentials=["none", 10.0 * (u.kJ / u.mol)],
            )

    @pytest.mark.parametrize(
        "typ,value",
        [
            ("slitpore", 10.0 * u.angstrom),
            ("cylinder", 10.0 * u.angstrom),
            ("sphere", 10.0 * u.angstrom),
            ("interface", [10.0 * u.angstrom, 20.0 * u.angstrom]),
        ],
    )
    def test_write_restricted_gemc_npt(self, twocomptwobox_system, typ, value):
        (system, moveset) = twocomptwobox_system
        moveset.add_restricted_insertions(
            system.species_topologies,
            [[None, None], [None, typ]],
            [[None, None], [None, value]],
        )

        inp_data = generate_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=500,
            temperature=300.0 * u.K,
            pressure=1 * u.bar,
        )

        if typ == "interface":
            assert (
                "\nrestricted_insertion {} {:0.1f} {:0.1f}\n".format(
                    typ, value[0].to_value(), value[1].to_value()
                )
                in inp_data
            )
        else:
            assert (
                "\nrestricted_insertion {} {:0.1f}\n".format(
                    typ, value.to_value()
                )
                in inp_data
            )

    @pytest.mark.parametrize(
        "typ,value",
        [
            ("slitpore", 60 * u.angstrom),
            ("cylinder", 60 * u.angstrom),
            ("sphere", 60 * u.angstrom),
            ("interface", [10 * u.angstrom, 70 * u.angstrom]),
        ],
    )
    def test_fail_restricted_gemc_npt(self, twocomptwobox_system, typ, value):
        (system, moveset) = twocomptwobox_system
        moveset.add_restricted_insertions(
            system.species_topologies,
            [[None, None], [None, typ]],
            [[None, None], [None, value]],
        )
        with pytest.raises(ValueError, match=r"Restricted insertion"):
            inp_data = generate_input(
                system=system,
                moveset=moveset,
                run_type="equilibration",
                run_length=500,
                temperature=300.0 * u.K,
                pressure=1 * u.bar,
            )

    @pytest.mark.parametrize(
        "angle_style", [["fixed"], ["harmonic"], "fixed", "harmonic"]
    )
    def test_onecomp_angle_style(self, onecomp_system, angle_style):
        with temporary_directory() as tmp_dir:
            with temporary_cd(tmp_dir):
                (system, moveset) = onecomp_system
                write_mcfs(system, angle_style=angle_style)

    @pytest.mark.parametrize("angle_style", ["fixed", "harmonic"])
    def test_twocomp_angle_style(self, twocomp_system, angle_style):
        with temporary_directory() as tmp_dir:
            with temporary_cd(tmp_dir):
                (system, moveset) = twocomp_system
                write_mcfs(system, angle_style=[angle_style, angle_style])

    def test_angle_style_error(self, onecomp_system):
        (system, moveset) = onecomp_system
        with pytest.raises(ValueError, match="Invalid"):
            write_mcfs(system, angle_style=["charmm"])

    def test_rst_inp_not_exists(self):
        with pytest.raises(FileNotFoundError):
            _generate_restart_inp(
                restart_from="equil",
                run_name="equil.rst.001",
                run_type=None,
                run_length=None,
            )

    def test_rst_inp(self, onecomp_system):
        (system, moveset) = onecomp_system
        with temporary_directory() as tmp_dir:
            with temporary_cd(tmp_dir):
                write_input(
                    system=system,
                    moveset=moveset,
                    run_type="equilibration",
                    run_length=500,
                    temperature=300 * u.K,
                )
                write_restart_input(
                    restart_from="nvt",
                    run_name="nvt.rst.001",
                    run_type=None,
                    run_length=None,
                )
                assert Path("nvt.rst.001.inp").is_file()

    def test_rst_inp_invalid_run_length(self, onecomp_system):
        (system, moveset) = onecomp_system
        with temporary_directory() as tmp_dir:
            with temporary_cd(tmp_dir):
                write_input(
                    system=system,
                    moveset=moveset,
                    run_type="equilibration",
                    run_length=500,
                    temperature=300 * u.K,
                )
                with pytest.raises(ValueError):
                    write_restart_input(
                        restart_from="nvt",
                        run_name="nvt.rst.001",
                        run_type=None,
                        run_length=200,
                    )
                with pytest.warns(UserWarning):
                    write_restart_input(
                        restart_from="nvt",
                        run_name="nvt.rst.001",
                        run_type=None,
                        run_length=500,
                    )

    def test_rst_inp_switch_run_type(self, onecomp_system):
        (system, moveset) = onecomp_system
        with temporary_directory() as tmp_dir:
            with temporary_cd(tmp_dir):
                write_input(
                    system=system,
                    moveset=moveset,
                    run_type="equilibration",
                    run_length=500,
                    temperature=300 * u.K,
                )
                contents = _generate_restart_inp(
                    restart_from="nvt",
                    run_name="nvt.rst.001",
                    run_type="production",
                    run_length=1000,
                )
                assert "# Run_Type\nproduction" in contents
                assert "run 1000" in contents
                assert "# Start_Type\ncheckpoint nvt.out.chk" in contents
                assert "# Run_Name\nnvt.rst.001.out" in contents
                write_restart_input(
                    restart_from="nvt",
                    run_name="nvt.rst.001",
                    run_type="production",
                    run_length=1000,
                )
                contents = _generate_restart_inp(
                    restart_from="nvt.rst.001",
                    run_name="nvt.rst.002",
                    run_type=None,
                    run_length=None,
                )
                assert "# Run_Type\nproduction" in contents
                assert "run 1000" in contents
                assert (
                    "# Start_Type\ncheckpoint nvt.rst.001.out.chk" in contents
                )
                assert "# Run_Name\nnvt.rst.002.out" in contents

    def test_rst_twobox(self, twobox_system):
        (system, moveset) = twobox_system
        with temporary_directory() as tmp_dir:
            with temporary_cd(tmp_dir):
                write_input(
                    system=system,
                    moveset=moveset,
                    run_type="equilibration",
                    run_length=500,
                    temperature=300 * u.K,
                )
                contents = _generate_restart_inp(
                    restart_from="gemc",
                    run_name="gemc.rst.001",
                    run_type="production",
                    run_length=1000,
                )
                assert "# Run_Type\nproduction" in contents
                assert "run 1000" in contents
                assert "# Start_Type\ncheckpoint gemc.out.chk\n\n!" in contents
                assert "# Run_Name\ngemc.rst.001.out" in contents

    @pytest.mark.parametrize(
        "system_fixture, base_name",
        [
            ("onecomp_system", "nvt"),
            ("twobox_system", "gemc"),
        ],
    )
    def test_rst_multiple_rst(self, request, system_fixture, base_name):
        """
        Test the creation of a chain of input files, each of which restarts from
        the previous input file in the chain. This is useful within the context
        of an automatic equilibration detection loop, in which a simulation needs to
        be restarted if its not equilibrated.

        This test evaluates systems with one or two boxes. Two boxes might be 
        problematic because some start types require two lines in the # Start_Type
        section.
        
        This test ensures that the input files generated at each restart step:
        1. Contain a valid # Start_Type header.
        2. Follow the correct checkpoint line format (two entries, starting with "checkpoint").
        3. Include only comments or whitespace between the checkpoint line and the next # header.

        Parameters:
        - system_fixture: The fixture name of the system setup, allowing tests with 
                          both one-component and two-box systems.
        - base_name: The base name used for the generated files (e.g., "nvt" or "gemc").
        """

        # Access the system fixture based on parameter
        (system, moveset) = request.getfixturevalue(system_fixture)
        repeats = 3
        with temporary_directory() as tmp_dir:
            with temporary_cd(tmp_dir):
                for count in range(repeats):
                    if count == 0:
                        # Initial write of the input file
                        run_name = f"{base_name}.{count:03d}"
                        write_input(
                            run_name=run_name,
                            system=system,
                            moveset=moveset,
                            run_type="equilibration",
                            run_length=2,
                            temperature=300 * u.K,
                        )
                    else:
                        restart_from = f"{base_name}.{count - 1:03d}"
                        run_name = f"{base_name}.{count:03d}"
                        write_restart_input(
                            restart_from=restart_from,
                            run_name=run_name,
                            run_type="equilibration",
                            run_length=1000,
                        )

                        with open(f"{run_name}.inp", mode="r") as f:
                            inp_contents = f.readlines()

                        # Step 1: Check we have a # Start_Type header
                        start_idx = self.check_start_type_header(inp_contents)

                        # Step 2: Check the line after "# Start_Type" has two entries, with the first as "checkpoint"
                        self.check_checkpoint_line(inp_contents, start_idx)

                        # Step 3: Ensure only comments or whitespace are present until the next '#' header
                        self.check_only_comments_or_whitespace(inp_contents, start_idx)