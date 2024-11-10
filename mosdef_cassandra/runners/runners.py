import datetime
import subprocess
import os
import re


from mosdef_cassandra.runners.utils import check_system
from mosdef_cassandra.runners.utils import get_restart_name
from mosdef_cassandra.writers.writers import write_mcfs
from mosdef_cassandra.writers.writers import write_configs
from mosdef_cassandra.writers.writers import write_input
from mosdef_cassandra.writers.writers import write_pdb
from mosdef_cassandra.writers.writers import write_restart_input
from mosdef_cassandra.utils.detect import detect_cassandra_binaries
from mosdef_cassandra.utils.exceptions import CassandraRuntimeError


def run(system, moveset, run_type, run_length, temperature, run_dir=None, **kwargs):
    """Run the Monte Carlo simulation with Cassandra

    The following steps are performed: write the molecular connectivity
    files for each species to disk, write the starting structures
    (if any) to disk, generate and write the Cassandra input file to disk,
    call Cassandra to generate the required fragment libraries, and
    call Cassandra to run the MC simulation.

    Parameters
    ----------
    system : mosdef_cassandra.System
        the System to simulate
    moveset : mosdef_cassandra.MoveSet
        the MoveSet to simulate
    run_type : "equilibration" or "production"
        the type of run; in "equilibration" mode, Cassandra adaptively changes
        the maximum translation, rotation, and volume move sizes to achieve
        an acceptance ratio of 0.5
    run_length : int
        length of the MC simulation
    run_dir : str
        directory where the simulation will be run
    temperature : float
        temperature at which to perform the MC simulation
    **kwargs : keyword arguments
        any other valid keyword arguments, see
        ``mosdef_cassandra.print_valid_kwargs()`` for details
    """

    if run_dir is None:
        run_dir = os.getcwd()

    # Check that the user has the Cassandra binary on their PATH
    # Also need library_setup.py on the PATH and python2
    py, fraglib_setup, cassandra = detect_cassandra_binaries()

    # Sanity checks
    # TODO: Write more of these
    check_system(system, moveset)

    # Write MCF files
    if "angle_style" in kwargs:
        write_mcfs(system, angle_style=kwargs["angle_style"], run_dir=run_dir)
    else:
        write_mcfs(system, run_dir=run_dir)

    # Write starting configs (if needed)
    write_configs(system, run_dir=run_dir)

    # Write input file
    inp_file = write_input(
        system=system,
        moveset=moveset,
        run_type=run_type,
        run_length=run_length,
        temperature=temperature,
        run_dir=run_dir,
        **kwargs,
    )

    # Write pdb files (this step will be removed when frag generation
    # is incorporated into this workflow )
    for isp, top in enumerate(system.species_topologies):
        filename = os.path.join(run_dir, "species{}.pdb".format(isp + 1))
        write_pdb(top, filename)

    log_file = os.path.join(run_dir, "mosdef_cassandra_{}.log".format(
        datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f")
    ))

    # Run fragment generation
    print("Generating fragment libraries...")
    _run_fraglib_setup(
        py,
        fraglib_setup,
        cassandra,
        inp_file,
        log_file,
        len(system.species_topologies),
        run_dir,
    )

    # Run simulation
    print("Running Cassandra...")
    _run_cassandra(cassandra, inp_file, log_file, run_dir)


def restart(
    total_run_length=None, restart_from=None, run_name=None, run_type=None, run_dir=None,
):
    """Restart a Monte Carlo simulation from a checkpoint file with Cassandra

    The function requires the following in the working directory. These items
    would have all been generated for the original run:
        * Cassandra input (.inp) file named {restart_from}.inp
        * Cassandra checkpoint file (.chk) name {restart_from}.out.chk
        * MCF files for each species
        * Fragment libraries for each species

    The maximum translation, rotation, and volume move sizes are read from
    the checkpoint file. Similarly, the starting structure is taken from the
    checkpoint file. If the "restart_name" is not provided or if the
    "run_name" is the same as "restart_name", ".rst.N" will be appended to the
    "run_name".

    If you wish to extend a simulation you will need to specify the _total_
    number of simulation steps desired with the total_run_length option. For example,
    if your original run was 1e6 MC steps, but you wish to extend it by an
    additional 1e6 steps, use total_run_length=2000000.

    Parameters
    ----------
    total_run_length: int, optional, default=None
        total length of the MC simulation; if None, use original simulation length
    restart_from: str, optional, default=None
        name of run to restart from; if None, searches current
        directory for Cassandra inp files
    run_name: str, optional, default=None
        name of this run; if None, appends ".rst.NNN." to run_name,
        where "NNN" is the restart iteration "001", "002", ...,
    run_type : str, "equilibration" or "production", default=None
        the type of run; in "equilibration" mode, Cassandra adaptively changes
        the maximum translation, rotation, and volume move sizes to achieve
        an acceptance ratio of 0.5. If None, use the same choice as the
        previous run.
    run_dir : str, directory where the simulation was run
    """

    if run_dir is None:
        run_dir = os.getcwd()

    valid_run_types = ["equilibration", "equil", "production", "prod"]
    # Check that the user has the Cassandra binary on their PATH
    # Also need library_setup.py on the PATH and python2
    py, fraglib_setup, cassandra = detect_cassandra_binaries()

    # Parse the arguments
    if total_run_length is not None:
        if not isinstance(total_run_length, int):
            raise TypeError("`total_run_length` must be an integer")
    if run_name is not None:
        if not isinstance(run_name, str):
            raise TypeError("`run_name` must be a string")
    if run_type is not None:
        if (
            not isinstance(run_type, str)
            or run_type.lower() not in valid_run_types
        ):
            raise TypeError(f"`run_type` must be one of: {valid_run_types}")
        if run_type.lower() == "equil" or run_type.lower() == "equilibration":
            run_type = "equilibration"
        if run_type.lower() == "prod" or run_type.lower() == "production":
            run_type = "production"

    restart_from, run_name = get_restart_name(restart_from, run_name)
    checkpoint_name = os.path.join(run_dir, restart_from + ".out.chk")
    if not os.path.isfile(checkpoint_name):
        raise FileNotFoundError(
            f"Checkpoint file: {checkpoint_name} does not exist."
        )

    write_restart_input(restart_from, run_name, run_type, total_run_length, run_dir)

    log_file = os.path.join(run_dir, "mosdef_cassandra_{}.log".format(
        datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f")
    ))

    print("Running Cassandra...")
    _run_cassandra(cassandra, run_name + ".inp", log_file, run_dir)


def _run_fraglib_setup(
    py, fraglib_setup, cassandra, inp_file, log_file, nspecies, run_dir
):
    """Builds the fragment libraries required to run Cassandra.

    Requires python.
    """

    species_pdb_files = ""
    for isp in range(nspecies):
        species_pdb_files += "species{}.pdb ".format(isp + 1)

    fraglib_cmd = (
        "{py} {fraglib_setup} {cassandra} {inp_file} "
        "{species_pdb_files}".format(
            py=py,
            fraglib_setup=fraglib_setup,
            cassandra=cassandra,
            inp_file=inp_file,
            species_pdb_files=species_pdb_files,
        )
    )

    p = subprocess.Popen(
        fraglib_cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        cwd=run_dir,
    )
    out, err = p.communicate()

    with open(log_file, "a") as log:
        header = (
            "\n*************************************************\n"
            "******* CASSANDRA FRAGLIB STANDARD OUTPUT *******\n"
            "*************************************************\n\n"
        )
        log.write(header)
        log.write(_clean_cassandra_log(out))
        header = (
            "\n*************************************************\n"
            "******* CASSANDRA FRAGLIB STANDARD ERROR ********\n"
            "*************************************************\n\n"
        )
        log.write(header)
        log.write(_clean_cassandra_log(err))

    if p.returncode != 0 or "error" in err.lower() or "error" in out.lower():
        raise CassandraRuntimeError(
            "Cassandra fragment library generation failed, "
            "see {} for details".format(log_file)
        )


def _run_cassandra(cassandra, inp_file, log_file, run_dir):
    """Calls Cassandra"""
    cassandra_cmd = "{cassandra} {inp_file}".format(
        cassandra=cassandra, inp_file=inp_file
    )
    p = subprocess.Popen(
        cassandra_cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        cwd=run_dir,
    )
    out, err = p.communicate()

    with open(log_file, "a") as log:
        header = (
            "\n*************************************************\n"
            "*********** CASSANDRA STANDARD OUTPUT ***********\n"
            "*************************************************\n\n"
        )
        log.write(header)
        log.write(_clean_cassandra_log(out))
        header = (
            "\n*************************************************\n"
            "*********** CASSANDRA STANDARD ERROR ************\n"
            "*************************************************\n\n"
        )
        log.write(header)
        log.write(_clean_cassandra_log(err))

    if p.returncode != 0 or "error" in err.lower() or "error" in out.lower():
        raise CassandraRuntimeError(
            "Cassandra exited with an error, "
            "see {} for details.".format(log_file)
        )


def _clean_cassandra_log(string):
    """Strip the text formatting from Cassandra output"""
    string = re.sub(r"\^\[\[0m", "", string)
    string = re.sub(r"\^\[\[1m", "", string)
    return string
