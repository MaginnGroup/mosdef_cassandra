import datetime
import subprocess
import mbuild
import parmed

from mosdef_cassandra.writers.writers import write_mcfs
from mosdef_cassandra.writers.writers import write_configs
from mosdef_cassandra.writers.writers import write_input
from mosdef_cassandra.writers.writers import write_pdb
from mosdef_cassandra.utils.detect import detect_cassandra_binaries
from mosdef_cassandra.utils.exceptions import CassandraRuntimeError


def run(system, moveset, run_type, run_length, temperature, **kwargs):
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
    temperature : float
        temperature at which to perform the MC simulation
    **kwargs : keyword arguments
        any other valid keyword arguments, see
        ``mosdef_cassandra.print_valid_kwargs()`` for details
    """

    # Check that the user has the Cassandra binary on their PATH
    # Also need library_setup.py on the PATH and python2
    py, fraglib_setup, cassandra = detect_cassandra_binaries()

    # Sanity checks
    # TODO: Write more of these
    _check_system(system, moveset)

    # Write MCF files
    if "angle_style" in kwargs:
        write_mcfs(system, angle_style=kwargs["angle_style"])
    else:
        write_mcfs(system)

    # Write starting configs (if needed)
    write_configs(system)

    # Write input file
    inp_file = write_input(
        system=system,
        moveset=moveset,
        run_type=run_type,
        run_length=run_length,
        temperature=temperature,
        **kwargs
    )

    # Write pdb files (this step will be removed when frag generation
    # is incorporated into this workflow )
    for isp, top in enumerate(system.species_topologies):
        filename = "species{}.pdb".format(isp + 1)
        write_pdb(top, filename)

    log_file = "mosdef_cassandra_{}.log".format(
        datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f")
    )

    # Run fragment generation
    print("Generating fragment libraries...")
    _run_fraglib_setup(
        py,
        fraglib_setup,
        cassandra,
        inp_file,
        log_file,
        len(system.species_topologies),
    )

    # Run simulation
    print("Running Cassandra...")
    _run_cassandra(cassandra, inp_file, log_file)


def restart(system, moveset, run_type, run_length, temperature, **kwargs):
    """Restart a Monte Carlo simulation with Cassandra

    This function is used to restart a Cassandra simulation from a
    checkpoint file. No new MCF files are written to disk; the function
    assumes they already exist. No new fragment libraries are generated,
    again, these should already exist from the original run. The maximum
    translation, rotation, and volume move sizes are read from the checkpoint
    file so the values specified in mc.Moves are **not** used. Similarly,
    the starting structure is taken from the checkpoint file. The
    keyword argument "restart_name" should specify the old "run_name",
    i.e., the "run_name" that you are restarting from. If the "restart_name"
    is not provided or if the "run_name" is the same as "restart_name",
    "-rst" will be appended to the "run_name".

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
    temperature : float
        temperature at which to perform the MC simulation
    **kwargs : keyword arguments
        any other valid keyword arguments, see
        ``mosdef_cassandra.print_valid_kwargs()`` for details
    """

    # Check that the user has the Cassandra binary on their PATH
    # Also need library_setup.py on the PATH and python2
    py, fraglib_setup, cassandra = detect_cassandra_binaries()

    kwargs["restart"] = True

    # Sanity checks
    # TODO: Write more of these
    _check_system(system, moveset)

    log_file = "mosdef_cassandra_{}.log".format(
        datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f")
    )

    # Write input file
    inp_file = write_input(
        system=system,
        moveset=moveset,
        run_type=run_type,
        run_length=run_length,
        temperature=temperature,
        **kwargs
    )

    print("Running Cassandra...")
    _run_cassandra(cassandra, inp_file, log_file)


def _run_fraglib_setup(
    py, fraglib_setup, cassandra, inp_file, log_file, nspecies
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
    )
    out, err = p.communicate()

    with open(log_file, "a") as log:
        header = (
            "*************************************************\n"
            "******* CASSANDRA FRAGLIB STANDARD OUTPUT *******\n"
            "*************************************************\n"
        )
        log.write(header)
        log.write(out)
        header = (
            "*************************************************\n"
            "******* CASSANDRA FRAGLIB STANDARD ERROR ********\n"
            "*************************************************\n"
        )
        log.write(err)

    if p.returncode != 0 or "error" in err.lower() or "error" in out.lower():
        raise CassandraRuntimeError(
            "Cassandra fragment library generation failed, "
            "see {} for details".format(log_file)
        )


def _run_cassandra(cassandra, inp_file, log_file):
    """Calls Cassandra

    """
    cassandra_cmd = "{cassandra} {inp_file}".format(
        cassandra=cassandra, inp_file=inp_file
    )
    p = subprocess.Popen(
        cassandra_cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    out, err = p.communicate()

    with open(log_file, "a") as log:
        header = (
            "*************************************************\n"
            "*********** CASSANDRA STANDARD OUTPUT ***********\n"
            "*************************************************\n"
        )
        log.write(header)
        log.write(out)
        header = (
            "*************************************************\n"
            "*********** CASSANDRA STANDARD ERROR ************\n"
            "*************************************************\n"
        )
        log.write(err)

    if p.returncode != 0 or "error" in err.lower() or "error" in out.lower():
        raise CassandraRuntimeError(
            "Cassandra exited with an error, "
            "see {} for details.".format(log_file)
        )


def _check_system(system, moveset):
    """Run a series of sanity checks on the System and MoveSet

    """

    if moveset.ensemble == "gemc" or moveset.ensemble == "gemc_npt":
        if len(system.boxes) != 2:
            raise ValueError(
                "{} requested but {} simulation "
                "boxes provided as part of system. {} requires "
                "2 simulation boxes".format(
                    moveset.ensemble, len(system.boxes), moveset.ensemble
                )
            )
    else:
        if len(system.boxes) != 1:
            raise ValueError(
                "{} requested but {} simulation "
                "boxes provided as part of system. {} requires "
                "1 simulation box".format(
                    moveset.ensemble, len(system.boxes), moveset.ensemble
                )
            )

    for box in system.boxes:
        if not isinstance(box, mbuild.Box) and not isinstance(
            box, mbuild.Compound
        ):
            raise TypeError(
                "Not all System.boxes are mbuild.Box "
                "or mbuild.Compound objects. It appears "
                "your System object has been corrupted"
            )

    # TODO: Add check that species_topologies provided to the
    # System and MoveSet are the same

    if not isinstance(system.species_topologies, list):
        raise TypeError(
            "System.species_topologies should be a "
            "list. It appears your System object "
            "has been corrupted"
        )

    for species in system.species_topologies:
        if not isinstance(species, parmed.Structure):
            raise TypeError(
                "Each species should be a parmed.Structure. "
                "It appears your System object has been "
                "corrupted"
            )

    try:
        system.check_natoms()
    except:
        raise ValueError(
            "The number of atoms in one or more boxes "
            "does not match the number expected from "
            "System.species_topologies and system.mols_in_boxes. "
            "It appears your System object has been corrupted"
        )
