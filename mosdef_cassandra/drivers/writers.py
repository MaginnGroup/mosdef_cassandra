
import parmed
import mbuild
from mbuild.formats.cassandramcf import write_mcf

import mosdef_cassandra
from mosdef_cassandra.drivers.gen_inp import generate_input

def write_mcfs(system):

    if not isinstance(system, mosdef_cassandra.System):
        raise TypeError('"system" must be of type '
                'mosdef_cassandra.System')

    for species_count, species in enumerate(system.species_topologies):
        if not isinstance(species,parmed.Structure):
            raise TypeError('Your "system" object appears to have '
                    'been corrupted. Species {} is not a parmed'
                    '.Structure object'.format(species))
        if len(species.dihedrals) > 0 and len(species.rb_torsions) > 0:
            raise ValueError('Your species has both CHARMM style '
                    'dihedrals and Ryckaert-Bellemans style dihedrals. '
                    'Only a single dihedral style per species is '
                    'currently supported')
        elif len(species.dihedrals) > 0:
            dihedral_style = 'charmm'
        elif len(species.rb_torsions) > 0:
            dihedral_style = 'opls'
        else:
            dihedral_style = 'none'

        mcf_name = 'species{}.mcf'.format(species_count+1)
        write_mcf(species,mcf_name,angle_style='harmonic',
                dihedral_style=dihedral_style)

def write_configs(system):

    if not isinstance(system, mosdef_cassandra.System):
        raise TypeError('"system" must be of type '
                'mosdef_cassandra.System')

    for box_count, box in enumerate(system.boxes):
        if ( not isinstance(box,mbuild.Compound) and
             not isinstance(box,mbuild.Box) ):
            raise TypeError('Your "system" object appears to have '
                    'been corrupted. Box {} is not a mbuild'
                    '.Compound object'.format(species))

        # Only save if box has particles inside
        # This only occurs if box is an mbuild.Compound
        if isinstance(box,mbuild.Compound):
            xyz_name = 'box{}.in.xyz'.format(box_count+1)
            box.save(xyz_name)

def write_input(system, moves, temperature, run_type, length, **kwargs):

    inp_data = generate_input(system, moves, temperature, run_type,
                            length, **kwargs)

    if 'run_name' in kwargs:
        run_name = kwargs['run_name']
    else:
        run_name = moves.ensemble

    inp_name = run_name + '.inp'

    with open(inp_name,'w') as inp:
        inp.write(inp_data)


