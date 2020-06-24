import glob
import os

import pytest

import mosdef_cassandra.examples as ex
import unyt as u
from mosdef_cassandra.utils.tempdir import *
from mosdef_cassandra.utils.exceptions import *
from mosdef_cassandra.tests.base_test import BaseTest


class TestExamples(BaseTest):
    @pytest.mark.parametrize("custom_args", [{"angle_style": ["fixed"]}, {}])
    def test_run_nvt(self, custom_args):
        with temporary_directory() as tmp_dir:
            with temporary_cd(tmp_dir):
                ex.run_nvt(custom_args=custom_args)
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
                        if save_data == True:
                            log_data.append(line)

                completed = False
                for line in log_data:
                    if "Cassandra simulation complete" in line:
                        completed = True
                assert completed

    def test_run_failure(self):
        kwargs = {"vdw_cutoff": 17.0 * u.angstrom}
        with temporary_directory() as tmp_dir:
            with temporary_cd(tmp_dir):
                with pytest.raises(
                    CassandraRuntimeError, match=r"Cassandra exited with"
                ):
                    ex.run_nvt(custom_args=kwargs)

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
                        if save_data == True:
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
                        if save_data == True:
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
                        if save_data == True:
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
                        if save_data == True:
                            log_data.append(line)

                completed = False
                for line in log_data:
                    if "Cassandra simulation complete" in line:
                        completed = True
                assert completed

    def test_run_gcmc_adsorption(self):
        with temporary_directory() as tmp_dir:
            with temporary_cd(tmp_dir):
                ex.run_gcmc_adsorption()
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
                        if save_data == True:
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
                        if save_data == True:
                            log_data.append(line)

                completed = False
                for line in log_data:
                    if "Cassandra simulation complete" in line:
                        completed = True
                assert completed
