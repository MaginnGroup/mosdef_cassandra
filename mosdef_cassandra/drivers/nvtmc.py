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


