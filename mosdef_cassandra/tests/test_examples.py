import glob
import pytest

import mosdef_cassandra.examples as ex
import unyt as u
from mosdef_cassandra.utils.tempdir import *
from mosdef_cassandra.utils.exceptions import *
from mosdef_cassandra.tests.base_test import BaseTest


@pytest.mark.long
class TestExamples(BaseTest):
    def test_run_nvt(self):
        with temporary_directory() as tmp_dir:
            with temporary_cd(tmp_dir):
                ex.run_nvt()
                log_files = sorted(
                    glob.glob("./mosdef_cassandra*.log"), key=os.path.getmtime
                )
                log_file = log_files[-1]
                log_data = []
                save_data = False
                with open(log_file) as log:
                    for line in log:
                        if "CASSANDRA STANDARD" in line:
                            save_data = True
                        if save_data:
                            log_data.append(line)

                completed = False
                for line in log_data:
                    if "Cassandra simulation complete" in line:
                        completed = True
                assert completed

    @pytest.mark.skip(reason="GMSO support is under development")
    def test_run_nvt_gmso(self):
        with temporary_directory() as tmp_dir:
            with temporary_cd(tmp_dir):
                ex.run_nvt_gmso()
                log_files = sorted(
                    glob.glob("./mosdef_cassandra*.log"), key=os.path.getmtime
                )
                log_file = log_files[-1]
                log_data = []
                save_data = False
                with open(log_file) as log:
                    for line in log:
                        if "CASSANDRA STANDARD" in line:
                            save_data = True
                        if save_data:
                            log_data.append(line)

                completed = False
                for line in log_data:
                    if "Cassandra simulation complete" in line:
                        completed = True
                assert completed

    def test_run_nvt_custom_mixing(self):
        mixing_dict = {"opls_138_s1 opls_140_s1": "1.0 1.0"}
        custom_args = {
            "mixing_rule": "custom",
            "custom_mixing_dict": mixing_dict,
        }
        with temporary_directory() as tmp_dir:
            with temporary_cd(tmp_dir):
                ex.run_nvt(**custom_args)
                log_files = sorted(
                    glob.glob("./mosdef_cassandra*.log"), key=os.path.getmtime
                )
                log_file = log_files[-1]
                log_data = []
                save_data = False
                with open(log_file) as log:
                    for line in log:
                        if "CASSANDRA STANDARD" in line:
                            save_data = True
                        if save_data:
                            log_data.append(line)

                completed = False
                for line in log_data:
                    if "Cassandra simulation complete" in line:
                        completed = True
                assert completed

    @pytest.mark.parametrize("custom_args", [{"angle_style": ["fixed"]}, {}])
    def test_run_nvt_spce(self, custom_args):
        with temporary_directory() as tmp_dir:
            with temporary_cd(tmp_dir):
                ex.run_nvt_spce(**custom_args)
                log_files = sorted(
                    glob.glob("./mosdef_cassandra*.log"), key=os.path.getmtime
                )
                log_file = log_files[-1]
                log_data = []
                save_data = False
                with open(log_file) as log:
                    for line in log:
                        if "CASSANDRA STANDARD" in line:
                            save_data = True
                        if save_data:
                            log_data.append(line)

                completed = False
                for line in log_data:
                    if "Cassandra simulation complete" in line:
                        completed = True
                assert completed

    @pytest.mark.parametrize("fix_bonds", [True, False])
    def test_run_nvt_mbuild(self, fix_bonds):
        with temporary_directory() as tmp_dir:
            with temporary_cd(tmp_dir):
                if fix_bonds:
                    ex.run_nvt_mbuild(fix_bonds)
                    log_files = sorted(
                        glob.glob("./mosdef_cassandra*.log"),
                        key=os.path.getmtime,
                    )
                    log_file = log_files[-1]
                    log_data = []
                    save_data = False
                    with open(log_file) as log:
                        for line in log:
                            if "CASSANDRA STANDARD" in line:
                                save_data = True
                            if save_data:
                                log_data.append(line)

                    completed = False
                    for line in log_data:
                        if "Cassandra simulation complete" in line:
                            completed = True
                    assert completed
                else:
                    with pytest.raises(CassandraRuntimeError):
                        ex.run_nvt_mbuild(fix_bonds)

    def test_run_failure(self):
        custom_args = {"vdw_cutoff": 17.0 * u.angstrom}
        with temporary_directory() as tmp_dir:
            with temporary_cd(tmp_dir):
                with pytest.raises(
                    CassandraRuntimeError, match=r"Cassandra exited with"
                ):
                    ex.run_nvt(**custom_args)

    def test_run_npt(self):
        with temporary_directory() as tmp_dir:
            with temporary_cd(tmp_dir):
                ex.run_npt()
                log_files = sorted(
                    glob.glob("./mosdef_cassandra*.log"), key=os.path.getmtime
                )
                log_file = log_files[-1]
                log_data = []
                save_data = False
                with open(log_file) as log:
                    for line in log:
                        if "CASSANDRA STANDARD" in line:
                            save_data = True
                        if save_data:
                            log_data.append(line)

                completed = False
                for line in log_data:
                    if "Cassandra simulation complete" in line:
                        completed = True
                assert completed

    def test_run_gcmc(self):
        with temporary_directory() as tmp_dir:
            with temporary_cd(tmp_dir):
                ex.run_gcmc()
                log_files = sorted(
                    glob.glob("./mosdef_cassandra*.log"), key=os.path.getmtime
                )
                log_file = log_files[-1]
                log_data = []
                save_data = False
                with open(log_file) as log:
                    for line in log:
                        if "CASSANDRA STANDARD" in line:
                            save_data = True
                        if save_data:
                            log_data.append(line)

                completed = False
                for line in log_data:
                    if "Cassandra simulation complete" in line:
                        completed = True
                assert completed

    def test_run_gemc(self):
        with temporary_directory() as tmp_dir:
            with temporary_cd(tmp_dir):
                ex.run_gemc()
                log_files = sorted(
                    glob.glob("./mosdef_cassandra*.log"), key=os.path.getmtime
                )
                log_file = log_files[-1]
                log_data = []
                save_data = False
                with open(log_file) as log:
                    for line in log:
                        if "CASSANDRA STANDARD" in line:
                            save_data = True
                        if save_data:
                            log_data.append(line)

                completed = False
                for line in log_data:
                    if "Cassandra simulation complete" in line:
                        completed = True
                assert completed

    def test_run_nvt_mixture(self):
        with temporary_directory() as tmp_dir:
            with temporary_cd(tmp_dir):
                ex.run_nvt_mixture()
                log_files = sorted(
                    glob.glob("./mosdef_cassandra*.log"), key=os.path.getmtime
                )
                log_file = log_files[-1]
                log_data = []
                save_data = False
                with open(log_file) as log:
                    for line in log:
                        if "CASSANDRA STANDARD" in line:
                            save_data = True
                        if save_data:
                            log_data.append(line)

                completed = False
                for line in log_data:
                    if "Cassandra simulation complete" in line:
                        completed = True
                assert completed

    def test_run_gcmc_adsorption(self):
        with temporary_directory() as tmp_dir:
            with temporary_cd(tmp_dir):
                # For speed; plus the example uses UA methane
                ex.run_gcmc_adsorption(charge_style="none")
                log_files = sorted(
                    glob.glob("./mosdef_cassandra*.log"), key=os.path.getmtime
                )
                log_file = log_files[-1]
                log_data = []
                save_data = False
                with open(log_file) as log:
                    for line in log:
                        if "CASSANDRA STANDARD" in line:
                            save_data = True
                        if save_data:
                            log_data.append(line)

                completed = False
                for line in log_data:
                    if "Cassandra simulation complete" in line:
                        completed = True
                assert completed

    def test_run_gcmc_restricted(self):
        with temporary_directory() as tmp_dir:
            with temporary_cd(tmp_dir):
                ex.run_gcmc_restricted()
                log_files = sorted(
                    glob.glob("./mosdef_cassandra*.log"), key=os.path.getmtime
                )
                log_file = log_files[-1]
                log_data = []
                save_data = False
                with open(log_file) as log:
                    for line in log:
                        if "CASSANDRA STANDARD" in line:
                            save_data = True
                        if save_data:
                            log_data.append(line)

                completed = False
                for line in log_data:
                    if "Cassandra simulation complete" in line:
                        completed = True
                assert completed
