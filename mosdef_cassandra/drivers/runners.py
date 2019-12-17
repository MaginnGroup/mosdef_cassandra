
from mosdef_cassandra.drivers.writers import write_mcfs
from mosdef_cassandra.drivers.writers import write_configs
from mosdef_cassandra.drivers.writers import write_input

def run(system, moves, temperature, run_type, length, **kwargs):

    # Sanity checks
    # E.g., system has multiple boxes for GEMC

    ensemble = moves.ensemble
    if ensemble == 'gemc' or ensemble == 'gemc_npt':
        if len(system.boxes) != 2:
            raise ValueError('{} requested but {} simulation'
                    'boxes provided as part of system. {} requires'
                    '2 simulation boxes'.format(ensemble,
                        len(system.boxes), ensemble))

    # Write MCF files
    write_mcfs(system)

    # Write starting configs (if needed)
    write_configs(system)

    # Write input file
    write_input(system, moves, temperature, run_type, length, **kwargs)

    # Run fragment generation

    # Run simulation


