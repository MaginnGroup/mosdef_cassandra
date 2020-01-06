
import datetime
import subprocess
import mbuild
import parmed

from mosdef_cassandra.writers.writers import write_mcfs
from mosdef_cassandra.writers.writers import write_configs
from mosdef_cassandra.writers.writers import write_input
from mosdef_cassandra.writers.writers import write_pdb
from mosdef_cassandra.utils.detect import detect_cassandra_binaries

def run(system, moves, temperature, run_type, length, **kwargs):

    # Check that the user has the Cassandra binary on their PATH
    # Also need library_setup.py on the PATH and python2
    py2, fraglib_setup, cassandra = detect_cassandra_binaries()

    # Sanity checks
    # TODO: Write more of these
    _check_system(system, moves)

    # Write MCF files
    write_mcfs(system)

    # Write starting configs (if needed)
    write_configs(system)

    # Write input file
    inp_file = write_input(system, moves, temperature, run_type, length, **kwargs)

    # Write pdb files (this step will be removed when frag generation
    # is incorporated into this workflow )
    for isp, top in enumerate(system.species_topologies):
        filename = 'species{}.pdb'.format(isp+1)
        write_pdb(top,filename)

    log_file = 'mosdef_cassandra_{}.log'.format(
            datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))

    # Run fragment generation (again, will be removed...)
    print("Generating fragment libraries...")
    successful_fraglib = _run_fraglib_setup(fraglib_setup, cassandra,
                                            inp_file, log_file,
                                            len(system.species_topologies))

    # Run simulation
    if successful_fraglib:
        print("Running Cassandra...")
        _run_cassandra(cassandra, inp_file, log_file)
    else:
       raise ValueError("Cassandra failed due to unsuccessful "
                        "fragment generation")

def restart(system, moves, temperature, run_type, length, **kwargs):

    # Check that the user has the Cassandra binary on their PATH
    # Also need library_setup.py on the PATH and python2
    py2, fraglib_setup, cassandra = detect_cassandra_binaries()

    kwargs['restart'] = True

    # Sanity checks
    # TODO: Write more of these
    _check_system(system, moves)

    log_file = 'mosdef_cassandra_{}.log'.format(
            datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))

    # Write input file
    inp_file = write_input(system, moves, temperature, run_type, length, **kwargs)

    print("Running Cassandra...")
    _run_cassandra(cassandra, inp_file, log_file)

def _run_fraglib_setup(fraglib_setup, cassandra, inp_file, log_file,
                        nspecies):
    """Builds the fragment libraries required to run Cassandra.

    Requires python2.
    """

    species_pdb_files = ""
    for isp in range(nspecies):
        species_pdb_files += "species{}.pdb ".format(isp+1)

    fraglib_cmd = ( 'python2 {fraglib_setup} {cassandra} {inp_file} '
                    '{species_pdb_files}'.format(fraglib_setup=fraglib_setup,
                    cassandra=cassandra,inp_file=inp_file,
                    species_pdb_files=species_pdb_files)
                  )

    p = subprocess.Popen(fraglib_cmd,
            shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True)
    out,err = p.communicate()

    with open (log_file, 'a') as log:
        header = ("*************************************************\n"
                  "******* CASSANDRA FRAGLIB STANDARD OUTPUT *******\n"
                  "*************************************************\n"
                 )
        log.write(header)
        log.write(out)
        header = ("*************************************************\n"
                  "******* CASSANDRA FRAGLIB STANDARD ERROR ********\n"
                  "*************************************************\n"
                 )
        log.write(err)

    if p.returncode != 0:
        print('Cassandra fragment library generation failed, '
              'see {} for details'.format(log_file))
        return False
    return True

def _run_cassandra(cassandra, inp_file, log_file):
    """Calls Cassandra

    """
    cassandra_cmd = ('{cassandra} {inp_file}'.format(cassandra=cassandra,
                                                     inp_file=inp_file))
    p = subprocess.Popen(cassandra_cmd,
            shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True)
    out,err = p.communicate()

    with open (log_file, 'a') as log:
        header = ("*************************************************\n"
                  "*********** CASSANDRA STANDARD OUTPUT ***********\n"
                  "*************************************************\n"
                 )
        log.write(header)
        log.write(out)
        header = ("*************************************************\n"
                  "*********** CASSANDRA STANDARD ERROR ************\n"
                  "*************************************************\n"
                 )
        log.write(err)

    if p.returncode != 0 or "error" in err.lower():
        print('Cassandra error, see {}'.format(log_file))

def _check_system(system, moves):
    """Run a series of sanity checks on the System and Moves objects

    """

    if moves.ensemble == 'gemc' or moves.ensemble == 'gemc_npt':
        if len(system.boxes) != 2:
            raise ValueError('{} requested but {} simulation '
                    'boxes provided as part of system. {} requires '
                    '2 simulation boxes'.format(moves.ensemble,
                        len(system.boxes), moves.ensemble))
    else:
        if len(system.boxes) != 1:
            raise ValueError('{} requested but {} simulation '
                    'boxes provided as part of system. {} requires '
                    '1 simulation box'.format(moves.ensemble,
                        len(system.boxes), moves.ensemble))

    for box in system.boxes:
        if ( not isinstance(box,mbuild.Box) and
             not isinstance(box,mbuild.Compound) ):
            raise TypeError('Not all System.boxes are mbuild.Box '
                            'or mbuild.Compound objects. It appears '
                            'your System object has been corrupted')

    # TODO: Add check that species_topologies provided to System
    # and Moves objects are the same

    if not isinstance(system.species_topologies,list):
        raise TypeError('System.species_topologies should be a '
                        'list. It appears your System object '
                        'has been corrupted')

    for species in system.species_topologies:
        if not isinstance(species,parmed.Structure):
            raise TypeError('Each species should be a parmed.Structure. '
                        'It appears your System object has been '
                        'corrupted')

    try:
        system.check_natoms()
    except:
        raise ValueError('The number of atoms in one or more boxes '
                'does not match the number expected from '
                'System.species_topologies and system.species_in_boxes. '
                'It appears your System object has been corrupted')

