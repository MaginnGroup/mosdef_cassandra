from __future__ import division

import os
import sys
import tempfile
import warnings
from distutils.spawn import find_executable
import subprocess
import contextlib
import shutil

import numpy as np

NM_TO_A = 10.0


def fill_box(compounds, n_compounds=None, box=None, density=None,
             temperature=300.0, compound_ratio=None, aspect_ratio=None,
             log_file=None, seed1=None, seed2=None):
    """Fill a box with one or more `mbuild.compound` using Cassandra.

    `fill_box` takes a single `mbuild.Compound` or a list of
    `mbuild.Compound`'s and returns an `mbuild.Compound` that has
    been filled to the user's specifications.

    When filling a system, two arguments of `n_compounds`, `box`,
    and `density` must be specified.

    If `n_compounds` and `box` are not None, `n_compounds` will be
    inserted into a box of dimensions specified by `box`.

    If `n_compounds` and `density` are not None, the required box size
    will be calculated. In this case, `n_compounds` must be an int and
    not a list of int.

    If `box` and `density` are not None, the required number of
    compounds will be calculated.

    If `box` is not specified, the default behavior is to generate a
    cubic box. Optionally, `aspect_ratio` can be passed to generate
    a non-cubic box with the desired aspect ratio.

    Parameters
    ----------
    compounds : pmd.Structure or list of multiple pmd.Structure objects
        Compound or list of compounds to place in the box.
    n_compounds : int or list of int
        Number of compounds to be placed in the box.
    box : mb.Box
        Box to be filled by compounds.
    density : float, units kg/m^3, default=None
        Target system density. If not None, either `n_compounds` or
        `box`, but not both, must be specified.
    seed1 : int, default=None
        Random seed to be passed to Cassandra.
        If none, a random number will be generated.
    seed2 : int, default=None
        Random seed to be passed to Cassandra.
        If none, a random number will be generated.
    compound_ratio : list, default=None
        Ratio of number of each compound to be put in box. Only used
        in the case more than one compound, with `n_compounds` not
        specified while `density` and `box` have been specified.
    aspect_ratio : list of float
        If a non-cubic box is desired, the ratio of box lengths in the x, y,
        and z directions.
    temp_file : str, default=None
        File name to write Cassandra's raw output to.
    update_port_locations : bool, default=False
        After packing, port locations can be updated, but since compounds
        can be rotated, port orientation may be incorrect.

    Returns
    -------
    filled : mb.Compound

    """
    # Check that the user has the Cassandra binary on their PATH
    # Also need library_setup.py on the PATH and python2
    py2, fraglib_setup, cassandra = _detect_cassandra_binaries()

    arg_count = 3 - [n_compounds, box, density].count(None)
    if arg_count != 2:
        msg = ("Exactly 2 of `n_compounds`, `box`, and `density` "
               "must be specified. {} were given.".format(arg_count))
        raise ValueError(msg)

    # Check box
    if box is not None:
        box = _validate_box(box)
    # Check compounds
    compounds, n_compounds = _validate_compounds(compounds, n_compounds)

    # Calculate n_compounds or box dimensions if required
    if density is not None:
        if n_compounds is not None:
            box = _calculate_box(compounds, n_compounds, density,
                    aspect_ratio)
        elif box is not None:
            n_compounds = _calculate_ncompounds(compounds, box, density,
                    compound_ratio)

    if log_file is None:
        log_file = 'cassandra-cbmc-log.txt'
        log_file = os.path.abspath(log_file)

    # Get random seeds
    seeds = _get_seeds(seed1, seed2)

    # Do everything with Cassandra in temporary directory
    with temporary_directory() as tmp_dir:
        with temporary_cd(tmp_dir):
            for i, compound in enumerate(compounds):
                mcf_file = 'species_{}.mcf'.format(i+1)
                pdb_file = 'species_{}.pdb'.format(i+1)
                # Infer dihedral style from parmed Structure
                # Assume harmonic angles
                if len(compound.rb_torsions) > 0:
                    dihedral_style = 'opls'
                elif len(compound.dihedrals) > 0:
                    dihedral_style = 'charmm'
                else:
                    dihedral_style = 'none'
                # Save MCF file for each mol
                write_mcf(compound, mcf_file, angle_style='harmonic',
                        dihedral_style=dihedral_style)
                # Save PDB file for each mol
                mb.load(compound).to_trajectory().save_pdb(pdb_file)

            # Generate a Cassandra input file
            inp_file = _write_cassandra_inp(n_compounds, box, temperature,
                                            seeds)

            # Print some info
            print("Placing n_compounds: {} in box dimensions "
                  " (nm, degrees): {}".format(n_compounds,box))

            # Call Cassandra
            print("Generating fragment libraries...")
            successful_fraglib = _run_fraglib_setup(fraglib_setup, cassandra,
                                 inp_file, len(compounds), log_file)

            if successful_fraglib:
                print("Running Cassandra to generate initial structure with CBMC...")
                _run_cassandra(cassandra, inp_file, log_file)
            else:
               raise MBuildError("Cassandra failed due to unsuccessful "
                                "fragment generation")

            # Extract new coords
            filled = _create_topology(compounds, n_compounds)
            # TODO: Look at what we can do about port locations
            filled.update_coordinates('pack.out.xyz',update_port_locations=False)
            filled.periodicity = np.asarray(box.lengths, dtype=np.float32)

    # Do everything with Cassandra in temporary directory

    return filled

def _validate_box(box):
    """Ensure that the box passed by the user can be formatted as an mbuild.Box

    Parameters
    ----------
    box : mbuild.Box or a tuple or list thereof
        Box or inputs to `mbuild.Box` to generate a `mbuild.Box`.

    Returns
    -------
    box : mbuild.Box
    """
    if isinstance(box, (list, tuple)):
        if len(box) == 3:
            box = Box(lengths=box)
        elif len(box) == 6:
            box = Box(mins=box[:3], maxs=box[3:])

    if not isinstance(box, Box):
        raise MBuildError('Unknown format for `box` parameter. Must pass a'
                          ' list/tuple of length 3 (box lengths) or length'
                          ' 6 (box mins and maxes) or an mbuild.Box object.')
    return box


def _new_xyz_file():
    """Generate PDB file using tempfile.NamedTemporaryFile.

    Return
    ------
    _ : file-object
        Temporary PDB file.
    """

    return tempfile.NamedTemporaryFile(suffix='.xyz', delete=False)


def _create_topology(compounds, n_compounds):
    """Return updated mBuild compound with new coordinates.

    Parameters
    ----------
    compounds : mb.Compound or list of mb.Compounds, required
        Compound(s) to add to the container.
    n_compounds : int or list of int, required
        Amount of compounds to container.

    Return
    ------
    container : mb.Compound
        Compound with added compounds from PACKMOL.
    """

    container = mb.Compound()
    for comp_structure, m_compound in zip(compounds, n_compounds):
        comp_compound = mb.load(comp_structure)
        for _ in range(m_compound):
            container.add(clone(comp_compound))
    return container


def _validate_compounds(compounds,n_compounds):

    if isinstance(compounds,pmd.Structure):
        compounds = [compounds]
    elif isinstance(compounds,(list,set)):
        for compound in compounds:
            if not isinstance(compound,pmd.Structure):
                msg = ('`Compound` argument must be a parmed.Structure '
                       'object or a list of multiple parmed.Structure '
                       'objects.')
                raise ValueError(msg)

    if n_compounds is not None:
        if isinstance(n_compounds,int):
            n_compounds = [n_compounds]
        elif isinstance(n_compounds,(list,set)):
            for n_compound in n_compounds:
                if not isinstance(n_compound,int):
                    msg = ('`n_compounds` argument must be of type '
                            '`int` or a list where each element is of '
                            'type `int`')
                    raise ValueError(msg)
        if len(compounds) != len(n_compounds):
            msg = ("`compounds` and `n_compounds` must be of equal length.")
            raise ValueError(msg)

    return compounds, n_compounds

def _calculate_box(compounds,n_compounds,density,aspect_ratio):
    """Calculate the required box dimensions given a specified density
    and number of compounds for each species"""

    total_mass = np.sum([n*np.sum([a.mass for a in c.atoms])
                        for c, n in zip(compounds, n_compounds)])
    # Conversion from (amu/(kg/m^3))**(1/3) to nm
    L = (total_mass/density)**(1/3)*1.1841763
    if aspect_ratio is None:
        box = _validate_box(Box(3*[L]))
    else:
        L *= np.prod(aspect_ratio) ** (-1/3)
        box = _validate_box(Box([val*L for val in aspect_ratio]))

    return box

def _calculate_ncompounds(compounds,box,density,compound_ratio):
    """Calculate the required number of compounds for each species
    given a specified density, box dimensions and (if more than one
    species), the desired ratio of the species"""

    if len(compounds) == 1:
        compound_mass = np.sum([a.mass for a in compounds[0].atoms])
        # Conversion from kg/m^3 / amu * nm^3 to dimensionless units
        n_compounds = [
                int(density/compound_mass*np.prod(box.lengths)*0.60224)]
    else:
        if compound_ratio is None:
            msg = ("Determing `n_compounds` from `density` and `box` "
                   "for systems with more than one compound type "
                   "requires `compound_ratio`")
            raise ValueError(msg)
        if len(compounds) != len(compound_ratio):
            msg = ("Length of `compound_ratio` must equal length of "
                   "`compounds`")
            raise ValueError(msg)
        prototype_mass = 0
        for c, r in zip(compounds, compound_ratio):
            prototype_mass += r * np.sum([a.mass for a in c.atoms])
        # Conversion from kg/m^3 / amu * nm^3 to dimensionless units
        n_prototypes = int(density/prototype_mass*np.prod(box.lengths)*0.60224)
        n_compounds = list()
        for c in compound_ratio:
            n_compounds.append(int(n_prototypes * c))

    return n_compounds


def _write_cassandra_inp(nbr_molecules,box,temperature,seeds):
    """Writes an input file for Cassandra"""

    # Cassandra works in units of Angstroms
    # But here we specify cutoffs in units of nm to
    # be consistent with mbuild
    lower_cutoff = 0.05
    upper_cutoff = 1.5
    cbmc_cutoff = 0.6
    upper_cutoff = _check_cutoff(upper_cutoff,box)
    cbmc_cutoff = _check_cutoff(cbmc_cutoff,box)
    nbr_species = len(nbr_molecules)
    molecule_files = ""
    start_type = "make_config "
    for imol in range(nbr_species):
        molecule_files += "species_{}.mcf {}\n".format(
                imol+1,nbr_molecules[imol])
        start_type += "{} ".format(nbr_molecules[imol])

    # Each basis vector is a COLUMN of the box_matrix
    formatted_box = "{} {} {}\n{} {} {}\n{} {} {}\n".format(
            box.lengths[0]*NM_TO_A,0.0,0.0,
            0.0,box.lengths[1]*NM_TO_A,0.0,
            0.0,0.0,box.lengths[2]*NM_TO_A)

    filename = 'pack.inp'
    with open(filename, 'w') as inpfile:
        inpfile.write("""! Input file for testing generating initial configuration

# Run_Name
{name}.out
!------------------------------------------------------------------------------

# Sim_Type
nvt_mc
!------------------------------------------------------------------------------

# Nbr_Species
{nbr_species}
!------------------------------------------------------------------------------

# VDW_Style
lj cut {upper_cutoff}
!------------------------------------------------------------------------------

# Charge_Style
coul ewald {upper_cutoff} 1e-3
!------------------------------------------------------------------------------

# Seed_Info
{seed1} {seed2}
!------------------------------------------------------------------------------

# Rcutoff_Low
{lower_cutoff}
!------------------------------------------------------------------------------

# Molecule_Files
{molecule_files}
!------------------------------------------------------------------------------

# Box_Info
1
cell_matrix
{formatted_box}
!------------------------------------------------------------------------------

# Temperature_Info
{temperature}
!------------------------------------------------------------------------------

# Move_Probability_Info
!------------------------------------------------------------------------------
!--Need to specify something here...won't affect CBMC used for make_config ----
!--Selected Prob_Angle bc input options not affected by nspecies---------------

# Prob_Angle
1.0

# Done_Probability_Info
!------------------------------------------------------------------------------

# Start_Type
{start_type}
!------------------------------------------------------------------------------

# Run_Type
production 1
!------------------------------------------------------------------------------

# Simulation_Length_Info
units        sweeps
prop_freq    1
coord_freq   1
run          0
!------------------------------------------------------------------------------

# Property_Info 1
energy_total
!------------------------------------------------------------------------------

# Fragment_Files
!------------------------------------------------------------------------------
!-----------------------------------library_setup.py will autofill this section

# CBMC_Info
kappa_ins 10
kappa_dih 10
rcut_cbmc {cbmc_cutoff}

END""".format(nbr_species=nbr_species,upper_cutoff=upper_cutoff*NM_TO_A,
        lower_cutoff=lower_cutoff*NM_TO_A,cbmc_cutoff=cbmc_cutoff*NM_TO_A,
        seed1=seeds[0],seed2=seeds[1],formatted_box=formatted_box,
        temperature=temperature,start_type=start_type,
        molecule_files=molecule_files))

    return filename

def _check_cutoff(cutoff,box):
    """ Makes sure that the cutoff is not more than one-half the
    box length in any dimension. If it is, reduces the cutoff
    accordingly.
    """
    min_box_length = min(box.lengths)
    if cutoff > 0.495*min_box_length:
        cutoff = 0.495*min_box_length

    return cutoff

def _run_fraglib_setup(fraglib_setup, cassandra, inp_file,
        nspecies, log_file):
    """Builds the fragment libraries required to run Cassandra.

    Requires python2.
    """

    species_pdb_files = ""
    for i in range(nspecies):
        species_pdb_files += "species_{}.pdb ".format(i+1)

    fraglib_cmd = ( 'python2 {fraglib_setup} {cassandra} {inp_file} '
                    '{species_pdb_files}'.format(fraglib_setup=fraglib_setup,
                    cassandra=cassandra,inp_file=inp_file,
                    species_pdb_files=species_pdb_files)
                  )

    p = subprocess.Popen(fraglib_cmd,
            shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True)
    out,err = p.communicate()

    with open (log_file, 'w') as log:
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
    """Calls Cassandra to generate an initial configuration

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


def run_name(name):
    # TODO: Check for spaces in name and remove
    inp_data = """
# Run_Name
{name}.out""".format(name=name)
    inp_data += """
!------------------------------------------------------------------------------
"""
    return inp_data

def sim_type(sim_type):
    sim_types = ['nvt_mc','npt_mc','gcmc','gemc']
    if sim_type.lower() not in sim_types:
        raise ValueError('Unsupported sim_type: {}. Supported options'
                'include {}'.format(sim_type,sim_types))

    inp_data = """
# Sim_Type
{sim_type}""".format(sim_type=sim_type)
    inp_data += """
!------------------------------------------------------------------------------
"""
    return inp_data

def nbr_species(nbr_species):
    inp_data = """
# Nbr_Species
{nbr_species}""".format(nbr_species=nbr_species)
    inp_data += """
!------------------------------------------------------------------------------
"""
    return inp_data

def vdw_style(vdw_styles,cut_styles,cutoffs):
    nbr_boxes = len(vdw_styles)
    assert len(vdw_styles) == len(cut_styles)
    assert len(vdw_styles) == len(cutoffs)
    valid_vdw_styles = ['lj','none']
    valid_cut_styles = { vstyle : [] for vstyle in valid_vdw_styles}
    valid_cut_styles['lj'].append('cut')
    valid_cut_styles['lj'].append('cut_tail')
    valid_cut_styles['lj'].append('cut_switch')
    valid_cut_styles['lj'].append('cut_shift')
    valid_cut_styles['none'].append(None)
    for vdw_style in vdw_styles:
        if vdw_style.lower() not in valid_vdw_styles:
            raise ValueError('Unsupported vdw_style: {}. Supported options '
                    'include {}'.format(vdw_style,vdw_styles))
    for cut_style,vdw_style in zip(cut_styles,vdw_styles):
        if cut_style.lower() not in valid_cut_styles[vdw_style]:
             raise ValueError('Unsupported cutoff style: {}. Supported '
                    'options for the selected vdw_style ({}) include '
                    '{}'.format(cut_style,vdw_style,
                        valid_cut_styles[vdw_style]))

    for cut_style,cutoff in zip(cut_styles,cutoffs):
        if cut_style.lower() == 'cut_switch':
            if len(cutoff) != 2:
                raise ValueError('Style "cut_switch" requires an inner '
                        'and outer cutoff. {} cutoffs were '
                        'specified'.format(len(cutoff)))

    inp_data = """
# VDW_Style
"""
    for vdw_style,cut_style,cutoff in zip(vdw_styles,cut_styles,cutoffs):
        if vdw_style.lower() == 'none':
            inp_data  += """
{vdw_style}""".format(vdw_style=vdw_style)
        else:
            if cut_style.lower() == 'cut_switch':
                inner_cutoff = cutoff[0]
                outer_cutoff = cutoff[1]
                inp_data  += """
{vdw_style} {cut_style} {inner_cutoff} {outer_cutoff}""".format(
                                          vdw_style=vdw_style,
                                          cut_style=cut_style,
                                          inner_cutoff=inner_cutoff,
                                          outer_cutoff=outer_cutoff)
            else:
                inp_data  += """
{vdw_style} {cut_style} {cutoff}""".format(vdw_style=vdw_style,
                                           cut_style=cut_style,
                                           cutoff=cutoff)
    inp_data += """
!------------------------------------------------------------------------------
"""

    return inp_data

def charge_style(charge_styles, cutoffs,
        ewald_accuracy=None, dsf_damping=None):
    nbr_boxes = len(charge_styles)
    assert len(charge_styles) == len(cutoffs)
    valid_charge_styles = ['none','cut','ewald','dsf']
    for charge_style in charge_styles:
        if charge_style.lower() not in valid_charge_styles:
            raise ValueError('Unsupported charge_style: {}. Supported options '
                    'include {}'.format(charge_style,charge_styles))
        if charge_style.lower() == 'ewald':
            if ewald_accuracy is None:
                raise ValueError('Ewald selected as the charge style but '
                        'no ewald accuracy provided')

    inp_data = """
# Charge_Style
"""
    for charge_style,cutoff in zip(charge_styles,cutoffs):
        if charge_style.lower() == 'none':
            inp_data  += """
{charge_style}""".format(charge_style=charge_style)

        elif charge_style.lower() == 'cut':
            inp_data += """
coul {charge_style} {cutoff}""".format(charge_style=charge_style,
                                       cutoff=cutoff)

        elif charge_style.lower() == 'ewald':
            inp_data += """
coul {charge_style} {cutoff} {accuracy}""".format(charge_style=charge_style,
                                                  cutoff=cutoff,
                                                  accuracy=ewald_accuracy)
        elif charge_style.lower() == 'dsf':
            if dsf_damping is not None:
                inp_data += """
coul {charge_style} {cutoff} {damping}""".format(charge_style=charge_style,
                                                 cutoff=cutoff,
                                                 damping=dsf_damping)
            else:
                inp_data += """
coul {charge_style} {cutoff}""".format(charge_style=charge_style,
                                       cutoff=cutoff)

    inp_data += """
!------------------------------------------------------------------------------
"""

    return inp_data


def mixing_rule(mixing_rule,custom_mixing_dict=None):
    valid_mixing_rules = ['lb','geometric','custom']
    if mixing_rule.lower() not in mixing_rules:
        raise ValueError('Unsupported mixing rule: {}. Supported options '
                    'include {}'.format(mixing_rule,valid_mixing_rules))
    if mixing_rule.lower() == 'custom' and custom_mixing_dict is None:
        raise ValueError('Custom mixing rule requested but no mixing '
                    'parmameters provided.')

    inp_data = """
# Mixing_Rule
{mixing_rule}""".format(mixing_rule=mixing_rule)

    if mixing_rule.lower() == 'custom':
        for pair, parms in mixing_dict.items():
            inp_data += """
{pair} {parms}
""".format(pair=pair,parms=parms)

    inp_data += """
!------------------------------------------------------------------------------
"""

    return inp_data

def seed_info(seed1=None,seed2=None):
    if seed1 is None:
        seed1 = np.random.randint(1,100000000)
    if seed2 is None:
        seed2 = np.random.randint(1,100000000)

    inp_data = """
# Seed_Info
{seed1} {seed2}""".format(seed1=seed1,seed2=seed2)

    inp_data += """
!------------------------------------------------------------------------------
"""

    return inp_data

def minimum_cutoff(cutoff):
     inp_data = """
# Rcutoff_Low
{cutoff}""".format(cutoff=cutoff)

    inp_data += """
!------------------------------------------------------------------------------
"""

    return inp_data

def pair_energy(save):
     inp_data = """
# Pair_Energy
{save}""".format(save=save)

    inp_data += """
!------------------------------------------------------------------------------
"""

    return inp_data

def molecule_files(max_molecule_dict):
     inp_data = """
# Molecule_Files"""

    for filename, max_mols in max_molecule_dict.items():
        inp_data = """
{filename} {max_mols}""".format(filename=filename,
                                max_mols=max_mols)

    inp_data += """
!------------------------------------------------------------------------------
"""

    return inp_data

def box_info(boxes):
"""Get the box info section of the input file

    Parameters
    ----------
    boxes : list
       list of box matrices with one box matrix
       per simulation box
"""
    nbr_boxes = len(boxes)
    for box in boxes:
        assert box.shape == (3,3)
    inp_data = """
# Box_Info
{nbr_boxes}""".format(nbr_boxes=nbr_boxes)

    box_types = []
    for box in boxes:
        if np.count_nonzero(box - np.diag(np.diagonal(box))) == 0:
            if np.all(np.diagonal(box) == box[0][0]):
                box_types.append('cubic')
            else:
                box_types.append('orthogonal')
        else:
            box_types.append('cell_matrix')

    for box,box_type in zip(boxes,box_types):
        inp_data += """
{box_type}
""".format(box_type=box_type)
        if box_type == 'cubic':
            inp_data += """
{dim}
""".format(dim=box[0][0])

        elif box_type == 'orthogonal':
            inp_data += """
{dim1} {dim2} {dim3}
""".format(dim1=box[0][0],dim2=box[1][1],dim3=box[2][2])

        else:
            inp_data += """
{ax} {bx} {cx}
{ay} {by} {cy}
{az} {bz} {cz}
""".format(ax=box[0][0],ay=box[0][1],az=box[0][2],
           bx=box[1][0],by=box[1][1],bz=box[1][2],
           cx=box[2][0],cy=box[2][1],cz=box[2][2])

    inp_data += """
!------------------------------------------------------------------------------
"""

    return inp_data

def temperature_info(temps):
"""Get the Temperature_Info section of the input file

   Parameters
   ----------
   temps : list
        list of temperatures with one for each box
"""
    nbr_boxes = len(temps)
    for temp in temps:
        if temp < 0.0:
            raise ValueError('Specified temperature ({}) is '
                    'less than zero'.format(temp))

    inp_data = """
# Temperature_Info"""

    for temp in temps:
        inp_data += """
{temperature}""".format(temperature=temp)

    inp_data += """
!------------------------------------------------------------------------------
"""

    return inp_data

def pressure_info(pressures):
"""Get the Pressure_Info section of the input file

   Parameters
   ----------
   pressures : list
        list of pressures with one for each box
"""

    inp_data = """
# Pressure_Info"""

    for press in pressures:
        inp_data += """
{pressure}""".format(pressure=press)

    inp_data += """
!------------------------------------------------------------------------------
"""

    return inp_data

def chemical_potential_info(chem_pots):
"""Get the Chemical_Potential_Info section of the input file

   Parameters
   ----------
   chem_pots : list
        list of chemical potentials with one for each species
        Non-insertable species should have None for the chemical potential
"""

    inp_data = """
# Chemical_Potential_Info
"""

    for chem_pot in chem_pots:
        inp_data += """{chem_pot} """.format(chem_pot=chem_pot)

    inp_data += """
!------------------------------------------------------------------------------
"""

    return inp_data


def move_probability_info(**kwargs):
"""Get the Move_Probability_Info section of the input file

   Parameters
   ----------
   kwargs : dict
        Dictionary of move probability information. Each valid keyword and
        associated information is described below

        'trans' :  [prob, box_i, box_j, ...]
                   where prob is the overall probability of
                   selecting a translation move and box_i/j are lists
                   containing the max displacement (angstroms) for
                   each species

        'rot'   :  [prob, box_i, box_j, ...]
                   where prob is the overall probability of
                   selecting a rotation move and box_i/j are lists
                   containing the max rotations (degrees) for
                   each species

        'angle' : prob
                   where prob is the overall probability of selecting
                   an angle move

        'dihed' : [prob, displacements]
                   where prob is the overall probability of selecting
                   a dihedral move and displacements is a list of the
                   maximum displacement (degrees) for dihedrals in
                   each species

        'regrow' : [prob, species_probs]
                   where prob is the overall probability of selecting
                   a regrowth move and species_probs is a list of the
                   probabilities of selecting a regrowth move for each
                   species

        'vol'    : [prob, displacements]
                   where prob is the overall probability of selecting
                   a volume move and displacements is a list of the
                   max volume change for each box

        'insert' : [prob, insertable]
                   where prob is the overall probability of selecting
                   a insertion/deletion move and insertable is a list
                   of booleans indicating whether each species is
                   insertable or not

        'swap'  : [prob, insertable, prob_species, prob_from_box]
                   where prob is the overall probability of selecting
                   a swap move and insertable is a list of booleans
                   indicating whether each species is insertable or not.
                   prob_species and prob_from_box are optional and
                   should be None if they are not to be specified.
                   prob_species is a list with one value for each
                   species that determines the probability of selecting
                   that species for the swap move. prob_from_box is a
                   list with one value for each box specifying the
                   probability of using the box as a donor box.

"""

    # First a sanity check on kwargs
    valid_args = [ 'trans', 'rot', 'angle', 'dihed', 'regrow', 'vol',
                   'insert', 'swap' ]

    for arg in kwargs:
        if arg not in valid_args:
            raise ValueError('Invaild probability info section {}. '
                    'Allowable options include {}'.format(
                     arg,valid_args))

    inp_data = """
# Move_Probability_Info
"""

    # Translation
    if 'trans' in kwargs:
        trans = kwargs['trans']
        if not isinstance(trans,list):
            raise ValueError('Translate probability information not '
                    'formatted properly')
        if not isinstance(trans[0],float):
            raise ValueError('Probability of translation move must be '
                    'a floating point value')
        for sp_displacements in trans[1:]:
            if not isinstance(sp_displacements,list):
                raise ValueError('Translate probability information not '
                    'formatted properly')
        inp_data += """
# Prob_Translation
{prob_trans}
""".format(prob_trans=trans[0])

        for sp_displacements in trans[1:]:
            for max_displace in sp_displacements:
                inp_data += """{} """.format(max_displace)
            inp_data += "\n"
        inp_data += """
!------------------------------------------------------------------------------
"""

    # Rotation
    if 'rot' in kwargs:
        rotate = kwargs['rot']
        if not isinstance(rotate,list):
            raise ValueError('Rotation probability information not '
                    'formatted properly')
        if not isinstance(rotate[0],float):
            raise ValueError('Probability of rotation move must be '
                    'a floating point value')
        for sp_displacements in rotate[1:]:
            if not isinstance(sp_displacements,list):
                raise ValueError('Rotation probability information not '
                    'formatted properly')
        inp_data += """
# Prob_Rotation
{prob_rotate}
""".format(prob_rotate=rotate[0])

        for sp_displacements in rotate[1:]:
            for max_displace in sp_displacements:
                inp_data += """{} """.format(max_displace)
            inp_data += "\n"

        inp_data += """
!------------------------------------------------------------------------------
"""

    # Angle
    if 'angle' in kwargs:
        angle = kwargs['angle']
        if not isinstance(angle,float):
            raise ValueError('Angle probability information not '
                    'formatted properly')
 
        inp_data += """
# Prob_Angle
{prob_angle}
""".format(prob_angle=angle)
        inp_data += """
!------------------------------------------------------------------------------
"""

    # Dihedral
    if 'dihed' in kwargs:
        dihed = kwargs['dihed']
        if not isinstance(dihed,list):
            raise ValueError('Dihedral probability information not '
                    'formatted properly')
        if not isinstance(dihed[0],float):
            raise ValueError('Probability of dihedral move must be '
                    'a floating point value')
        for sp_displacements in dihed[1:]:
            if not isinstance(sp_displacements,list):
                raise ValueError('Dihedral probability information not '
                    'formatted properly')
        inp_data += """
# Prob_Dihedral
{prob_dihed}
""".format(prob_dihed=dihed[0])

        for sp_displacements in dihed[1:]:
            for max_displace in sp_displacements:
                inp_data += """{} """.format(max_displace)
            inp_data += "\n"

        inp_data += """
!------------------------------------------------------------------------------
"""

    # Regrowth
    if 'regrow' in kwargs:
        regrow = kwargs['regrow']
        if not isinstance(regrow,list):
            raise ValueError('Regrowth probability information not '
                    'formatted properly')
        if len(regrow) != 2:
            raise ValueError('Regrowth probability information not '
                    'formatted properly')
        if not isinstance(regrow[0],float):
            raise ValueError('Probability of regrowth move must be '
                    'a floating point value')
        if not isinstance(regrow[1],list):
            raise ValueError('Regrowth probability information not '
                    'formatted properly')
        for sp_probs in regrow[1]:
            if not isinstance(sp_probs,float):
                raise ValueError('Probability of selecting each '
                    'species for a regrowth move must be a '
                    'floating point value')

        inp_data += """
# Prob_Regrowth
{prob_regrow}
""".format(prob_regrow=regrow[0])

        for sp_prob in regrow[1]:
            inp_data += """{} """.format(max_displace)
        
        inp_data += """
!------------------------------------------------------------------------------
"""

    # Volume
    if 'vol' in kwargs:
        vol = kwargs['vol']
        if not isinstance(vol,list):
            raise ValueError('Volume probability information not '
                    'formatted properly')
        if len(vol) != 2:
            raise ValueError('Volume probability information not '
                    'formatted properly')
        if not isinstance(vol[0],float):
            raise ValueError('Probability of volume move must be '
                    'a floating point value')
        if not isinstance(vol[1],list):
            raise ValueError('Volume probability information not '
                    'formatted properly')
        for max_displace in vol[1]:
            if not isinstance(max_displace,float):
                raise ValueError('Max displacement for volume move '
                    'must be a floating point value')

        inp_data += """
# Prob_Volume
{prob_vol}""".format(prob_vol=vol[0])
        
        for max_displace in vol[1]:
            inp_data += """
{max_displace}""".format(max_displace=max_displace)

        inp_data += """
!------------------------------------------------------------------------------
"""

    # Insert/delete
    if 'insert' in kwargs:
        insert = kwargs['insert']
        if not isinstance(insert,list):
            raise ValueError('Insertion probability information not '
                    'formatted properly')
        if len(insert) != 2:
            raise ValueError('Insertion probability information not '
                    'formatted properly')
        if not isinstance(insert[0],float):
            raise ValueError('Probability of insertion move must be '
                    'a floating point value')
        if not isinstance(insert[1],list):
            raise ValueError('Insertion probability information not '
                    'formatted properly')
        for insertable in insert[1]:
            if not isinstance(insertable,bool):
                raise ValueError('Whether or not a species is insertable'
                    'must be a boolean value')

        inp_data += """
# Prob_Insertion
{prob_insert}
""".format(prob_insert=insert[0])
        
        for insertable in insert[1]:
            if insertable:
                inp_data += """cbmc """
            else:
                inp_data += """none """

        inp_data += """
!------------------------------------------------------------------------------
"""
        inp_data += """
# Prob_Deletion
{prob_insert}
""".format(prob_insert=insert[0])
 
        inp_data += """
!------------------------------------------------------------------------------
"""
 
    # Swap
    if 'swap' in kwargs:
        swap = kwargs['swap']
        if not isinstance(swap,list):
            raise ValueError('Swap probability information not '
                    'formatted properly')
        if len(swap) != 4:
            raise ValueError('Swap probability information not '
                    'formatted properly')
        if not isinstance(swap[0],float):
            raise ValueError('Probability of swap move must be '
                    'a floating point value')
        if not isinstance(swap[1],list):
            raise ValueError('Insertion probability information not '
                    'formatted properly')
        for insertable in swap[1]:
            if not isinstance(insertable,bool):
                raise ValueError('Whether or not a species is insertable'
                    'must be a boolean value')
        if swap[2] is not None:
            if not isinstance(swap[2],list):
                raise ValueError('Swap probability information not '
                    'formatted properly')
            for prob in swap[2]:
                if not isinstance(prob,float):
                    raise ValueError('Probability of selecting species '
                            'for a swap move must be a floating point '
                            'value')
        if swap[3] is not None:
            if not isinstance(swap[3],list):
                raise ValueError('Swap probability information not '
                    'formatted properly')
            for prob in swap[3]:
                if not isinstance(prob,float):
                    raise ValueError('Probability of selecting box '
                            'as donor for a swap move must be a '
                            'floating point value')

        inp_data += """
# Prob_Swap
{prob_swap}
""".format(prob_swap=swap[0])
        
        for insertable in swap[1]:
            if insertable:
                inp_data += """cbmc """
            else:
                inp_data += """none """

        if swap[2] is not None:
            inp_data += """
prob_swap_species """
            for prob in swap[2]:
                inp_data += '{} '.format(prob)

        if swap[3] is not None:
            inp_data += """
prob_swap_from_box """
            for prob in swap[3]:
                inp_data += '{} '.format(prob)
        inp_data += """
!------------------------------------------------------------------------------
"""
 
    inp_data += """
# Done_Probability_Info
!------------------------------------------------------------------------------
"""

    return inp_data


def start_type(start_types):
"""Get the Start_Type section of the input file

   Parameters
   ----------
   start_type : list
        list of start_type with one for each box
"""

    inp_data = """
# Start_Type"""

    for start_type in start_types:
        inp_data += """
{start_type}""".format(start_type=start_type)

    inp_data += """
!------------------------------------------------------------------------------
"""

    return inp_data

def run_type(run_type,thermal_freq,vol_freq=None):
"""Get the Run_Type section of the input file

   Parameters
   ----------
   run_type : string
        'equilibration' or 'production'
   thermal_freq : int
        frequency of updating thermal move displacement 
        widths/output statistics
   vol_freq : int
        frequency of updating volume move displacement 
        widths/output statistics
"""

    valid_run_types = ['equilibration','production']
    if run_type not in valid_run_types:
        raise ValueError('Invalid run type specified {} '
                'Allowable options include {}'.format(run_type,
                    valid_run_types))
    if not isinstance(thermal_freq,int):
        raise ValueError('thermal_freq must be an integer')
    if vol_freq is not None and not isinstance(vol_freq,int):
        raise ValueError('thermal_freq must be an integer')

    inp_data = """
# Run_Type
{run_type} {thermal_freq} """.format(run_type=run_type)

    if vol_freq is not None:
        inp_data += '{vol_freq}'.format(vol_freq=vol_freq)

    inp_data += """
!------------------------------------------------------------------------------
"""

    return inp_data

def simulation_length_info(units,prop_freq,coord_freq,length
        steps_per_sweep=None,block_average_freq=None):
"""Get the Simulation_Length_Info section of the input file

   Parameters
   ----------
   units : string
        'minutes', 'steps', or 'sweeps'
   prop_freq : int
        frequency of writing property info
   coord_freq : int
        frequency of writing coordinates to file
   length : int
        number of (units) to run the simulation
   steps_per_sweep : int, optional
        number of steps in a MC sweep
   block_avg_freq : int, optional
        write properties as block averages, averaged over
        block_avg_freq (units)
"""

    valid_units = ['minutes','steps','sweeps']
    if units not in valid_units:
        raise ValueError('Invalid units specified {} Allowable options '
                'include {}'.format(units,valid_units))
    if not isinstance(prop_freq,int):
        raise ValueError('prop_freq must be an integer')
    if not isinstance(coord_freq,int):
        raise ValueError('coord_freq must be an integer')
    if not isinstance(length,int):
        raise ValueError('length must be an integer')
    if steps_per_sweep is not None:
        if not isinstance(steps_per_sweep,int):
            raise ValueError('steps_per_sweep must be an integer')
    if block_avg_freq is not None:
        if not isinstance(block_average_freq,int):
            raise ValueError('block_avg_freq must be an integer')


    inp_data = """
# Simulation_Length_Info
units {units}
prop_freq {prop_freq}
coord_freq {coord_freq}
run {length}""".format(units=units,prop_freq=prop_freq,
                       coord_freq=coord_freq,length=length)

    if steps_per_sweep is not None:
        inp_data += """
steps_per_sweep {steps_per_sweep}
""".format(steps_per_sweep=steps_per_sweep)
    if block_avg_freq is not None:
        inp_data += """
block_avg_freq {block_avg_freq}
""".format(block_avg_freq=block_avg_freq)

    inp_data += """
!------------------------------------------------------------------------------
"""

    return inp_data

def get_fragment_files():

    inp_data = """
# Fragment_Files
!------------------------------------------------------------------------------
"""

    return inp_data

def get_verbose_log(verbose):

    if not isinstance(verbose,bool):
        raise ValueError('Verbosity must be a boolean')

    inp_data = """
# Verbose_Logfile
{verbose}
!------------------------------------------------------------------------------
""".format(verbose)

    return inp_data

def get_cbmc_info(n_insert, n_dihed, cutoffs):
    """Get the CBMC_Info section of the input file

    Parameters
    ----------
    n_insert : int
        number of insertion sites to attempt for CBMC
    n_dihed : int
        number of dihedral angles to attempt for CBMC
    cutoffs : list
        list containing CBMC cutoff values for each box
    """

    if not isinstance(n_insert,int):
        raise ValueError('Number of CBMC insertion attempts must be'
                    'an integer')
    if not isinstance(n_dihed,int):
        raise ValueError('Number of CBMC dihedral angle attempts must be'
                    'an integer')
    if not isinstance(cutoffs,list):
        raise ValueError('Cutoff information improperly specified')
    for cutoff in cutoffs:
        if not isinstance(cutoff,float):
            raise ValueError('CBMC cutoff must be a float')


    inp_data = """
# CBMC_Info
kappa_ins {n_insert}
kappa_dih {n_dihed}
rcut_cbmc""".format(n_insert=n_insert,n_dihed=n_dihed)

    for cutoff in cutoffs:
        inp_data += ' {cutoff}'.format(cutoff=cutoff)

    inp_data += """
!------------------------------------------------------------------------------
"""

    return inp_data












