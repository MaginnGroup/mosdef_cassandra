import pytest
from pathlib import Path

import mosdef_cassandra as mc
from mosdef_cassandra.tests.base_test import BaseTest
from mosdef_cassandra.runners.utils import get_restart_name
from mosdef_cassandra.utils.tempdir import temporary_directory, temporary_cd


class TestRunners(BaseTest):
    def test_mismatch_boxes(self, methane_oplsaa, box):
        with pytest.raises(ValueError, match=r"requires 1 simulation"):
            system = mc.System(
                [box, box], [methane_oplsaa], mols_to_add=[[10], [0]]
            )
            moveset = mc.MoveSet("nvt", [methane_oplsaa])
            mc.run(system, moveset, 300.0, "equilibration", 500)
        with pytest.raises(ValueError, match=r"requires 2 simulation"):
            system = mc.System([box], [methane_oplsaa], mols_to_add=[[10]])
            moveset = mc.MoveSet("gemc", [methane_oplsaa])
            mc.run(system, moveset, 300.0, "equilibration", 500)

    def test_corrupt_boxes(self, methane_oplsaa, box):
        with pytest.raises(TypeError, match=r"corrupted"):
            system = mc.System([box], [methane_oplsaa], mols_to_add=[[10]])
            moveset = mc.MoveSet("nvt", [methane_oplsaa])

            system.boxes[0] = 1

            mc.run(system, moveset, 300.0, "equilibration", 500)

    def test_corrupt_topologies(self, methane_oplsaa, box):
        with pytest.raises(TypeError, match=r"corrupted"):
            system = mc.System([box], [methane_oplsaa], mols_to_add=[[10]])
            moveset = mc.MoveSet("nvt", [methane_oplsaa])

            system._species_topologies = methane_oplsaa
            mc.run(system, moveset, 300.0, "equilibration", 500)

        with pytest.raises(TypeError, match=r"corrupted"):
            system = mc.System([box], [methane_oplsaa], mols_to_add=[[10]])
            moveset = mc.MoveSet("nvt", [methane_oplsaa])

            system.species_topologies[0] = 1
            mc.run(system, moveset, 300.0, "equilibration", 500)

    def test_corrupt_natoms(self, methane_oplsaa, box):
        with pytest.raises(ValueError, match=r"corrupted"):
            system = mc.System([box], [methane_oplsaa], mols_to_add=[[10]])
            moveset = mc.MoveSet("nvt", [methane_oplsaa])

            system.mols_in_boxes[0][0] = 10
            mc.run(system, moveset, 300.0, "equilibration", 500)

    def test_restart_run_name_simple(self):
        restart_from, run_name = get_restart_name("equil", "equil.rst")
        assert run_name == "equil.rst"

    def test_restart_run_name_equal(self):
        restart_from, run_name = get_restart_name("equil", "equil")
        assert run_name == "equil.rst.001"

    def test_restart_run_name_none(self):
        restart_from, run_name = get_restart_name("equil", None)
        assert run_name == "equil.rst.001"

    def test_restart_run_name_iterate(self):
        restart_from, run_name = get_restart_name("equil.rst.002", None)
        assert run_name == "equil.rst.003"

    def test_restart_run_name_max_iterations(self):
        with pytest.raises(ValueError, match=r"Maximum number"):
            get_restart_name("equil.rst.999", None)

    def test_restart_run_name_no_files(self):
        with pytest.raises(FileNotFoundError):
            get_restart_name(None, None)

    def test_restart_run_name_single_file(self):
        with temporary_directory() as tmp_dir:
            with temporary_cd(tmp_dir):
                Path("equil.inp").touch()
                restart_from, run_name = get_restart_name(None, None)
                assert run_name == "equil.rst.001"

        with temporary_directory() as tmp_dir:
            with temporary_cd(tmp_dir):
                Path("equil.rst.001.inp").touch()
                restart_from, run_name = get_restart_name(None, None)
                assert run_name == "equil.rst.002"

    def test_restart_run_name_multiple_files(self):
        with temporary_directory() as tmp_dir:
            with temporary_cd(tmp_dir):
                Path("equil.inp").touch()
                Path("equil.rst.001.inp").touch()
                Path("equil.rst.002.inp").touch()
                restart_from, run_name = get_restart_name(None, None)
                assert run_name == "equil.rst.003"

        with temporary_directory() as tmp_dir:
            with temporary_cd(tmp_dir):
                Path("equil.rst.001.inp").touch()
                restart_from, run_name = get_restart_name(None, None)
                assert run_name == "equil.rst.002"

    def test_restart_run_name_multiple_invalid_files(self):
        with temporary_directory() as tmp_dir:
            with temporary_cd(tmp_dir):
                Path("prod.inp").touch()
                Path("equil.rst.001.inp").touch()
                Path("equil.rst.002.inp").touch()
                with pytest.raises(ValueError, match=r"Multiple"):
                    get_restart_name(None, None)
