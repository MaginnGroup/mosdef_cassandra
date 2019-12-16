from __future__ import division

import numpy as np

def run_name(name):
    """Get the Run_Name section of the input file"""

    if not isinstance(name,str):
        raise TypeError('name: {} must be a string'.format(name))

    name = name.replace(' ','-')
    inp_data = """
# Run_Name
{name}.out
!------------------------------------------------------------------------------
""".format(name=name)

    return inp_data

def sim_type(sim_type):
    """Get the Sim_Type section of the input file"""

    sim_types = ['nvt_mc','npt_mc','gcmc','gemc', 'gemc_npt']
    if sim_type not in sim_types:
        raise ValueError('Unsupported sim_type: {}. Supported options'
                'include {}'.format(sim_type,sim_types))

    inp_data = """
# Sim_Type
{sim_type}
!------------------------------------------------------------------------------
""".format(sim_type=sim_type)

    return inp_data

def nbr_species(nbr_species):
    """Get the Nbr_Species section of the input file"""

    if not isinstance(nbr_species,int):
        raise TypeError('nbr_species must be an int')

    inp_data = """
# Nbr_Species
{nbr_species}
!------------------------------------------------------------------------------
""".format(nbr_species=nbr_species)

    return inp_data

def vdw_style(vdw_styles,cut_styles,cutoffs):
    """Get the VDW_Style section of the input file

    Parameters
    ----------
    vdw_styles : list
        list of vdw_style for each box, one entry per box
    cut_styles : list
        list of cutoff_style for each box, one entry per box. For a
        box with vdw_style == 'none', the cutoff style is None
    cutoffs : list
        list with cutoffs for each box, one entry per box For a
        box with vdw_style == 'none', the cutoff is None
    """

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
        if vdw_style not in valid_vdw_styles:
            raise ValueError('Unsupported vdw_style: {}. Supported options '
                    'include {}'.format(vdw_style,vdw_styles))
    for cut_style,vdw_style in zip(cut_styles,vdw_styles):
        if cut_style not in valid_cut_styles[vdw_style]:
             raise ValueError('Unsupported cutoff style: {}. Supported '
                    'options for the selected vdw_style ({}) include '
                    '{}'.format(cut_style,vdw_style,
                        valid_cut_styles[vdw_style]))

    for cut_style,cutoff in zip(cut_styles,cutoffs):
        if cut_style == 'cut_switch':
            if len(cutoff) != 2:
                raise ValueError('Style "cut_switch" requires an inner '
                        'and outer cutoff. {} cutoffs were '
                        'specified'.format(len(cutoff)))

    inp_data = """
# VDW_Style
"""
    for vdw_style,cut_style,cutoff in zip(vdw_styles,cut_styles,cutoffs):
        if vdw_style == 'none':
            inp_data  += """
{vdw_style}""".format(vdw_style=vdw_style)
        else:
            if cut_style == 'cut_switch':
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
    """Get the Charge_Style section of the input file

    Parameters
    ----------
    charge_styles : list
        list of charge styles, one for each box
    cutoffs :
        list of coulombic cutoffs, one for each box. For a box with
        charge style 'none', the cutoff should be None
    ewald_accuracy : float, optional
        accuracy of ewald sum. Required if charge_style == ewald
    dsf_damping : float, optional
        value for dsf damping.
    """
    assert len(charge_styles) == len(cutoffs)
    valid_charge_styles = ['none','cut','ewald','dsf']
    for charge_style in charge_styles:
        if charge_style not in valid_charge_styles:
            raise ValueError('Unsupported charge_style: {}. Supported options '
                    'include {}'.format(charge_style,charge_styles))
        if charge_style == 'ewald':
            if ewald_accuracy is None:
                raise ValueError('Ewald selected as the charge style but '
                        'no ewald accuracy provided')

    inp_data = """
# Charge_Style
"""
    for charge_style,cutoff in zip(charge_styles,cutoffs):
        if charge_style == 'none':
            inp_data  += """
{charge_style}""".format(charge_style=charge_style)

        elif charge_style == 'cut':
            inp_data += """
coul {charge_style} {cutoff}""".format(charge_style=charge_style,
                                       cutoff=cutoff)

        elif charge_style == 'ewald':
            inp_data += """
coul {charge_style} {cutoff} {accuracy}""".format(charge_style=charge_style,
                                                  cutoff=cutoff,
                                                  accuracy=ewald_accuracy)
        elif charge_style == 'dsf':
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
    if mixing_rule not in mixing_rules:
        raise ValueError('Unsupported mixing rule: {}. Supported options '
                    'include {}'.format(mixing_rule,valid_mixing_rules))
    if mixing_rule == 'custom' and custom_mixing_dict is None:
        raise ValueError('Custom mixing rule requested but no mixing '
                    'parmameters provided.')

    inp_data = """
# Mixing_Rule
{mixing_rule}""".format(mixing_rule=mixing_rule)

    if mixing_rule == 'custom':
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
            raise TypeError('Translate probability information not '
                    'formatted properly')
        if not isinstance(trans[0],float):
            raise TypeError('Probability of translation move must be '
                    'a floating point value')
        for sp_displacements in trans[1:]:
            if not isinstance(sp_displacements,list):
                raise TypeError('Translate probability information not '
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
            raise TypeError('Rotation probability information not '
                    'formatted properly')
        if not isinstance(rotate[0],float):
            raise TypeError('Probability of rotation move must be '
                    'a floating point value')
        for sp_displacements in rotate[1:]:
            if not isinstance(sp_displacements,list):
                raise TypeError('Rotation probability information not '
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
            raise TypeError('Angle probability information not '
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
            raise TypeError('Dihedral probability information not '
                    'formatted properly')
        if not isinstance(dihed[0],float):
            raise TypeError('Probability of dihedral move must be '
                    'a floating point value')
        for sp_displacements in dihed[1:]:
            if not isinstance(sp_displacements,list):
                raise TypeError('Dihedral probability information not '
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
            raise TypeError('Regrowth probability information not '
                    'formatted properly')
        if len(regrow) != 2:
            raise TypeError('Regrowth probability information not '
                    'formatted properly')
        if not isinstance(regrow[0],float):
            raise TypeError('Probability of regrowth move must be '
                    'a floating point value')
        if not isinstance(regrow[1],list):
            raise TypeError('Regrowth probability information not '
                    'formatted properly')
        for sp_probs in regrow[1]:
            if not isinstance(sp_probs,float):
                raise TypeError('Probability of selecting each '
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
            raise TypeError('Volume probability information not '
                    'formatted properly')
        if len(vol) != 2:
            raise TypeError('Volume probability information not '
                    'formatted properly')
        if not isinstance(vol[0],float):
            raise TypeError('Probability of volume move must be '
                    'a floating point value')
        if not isinstance(vol[1],list):
            raise TypeError('Volume probability information not '
                    'formatted properly')
        for max_displace in vol[1]:
            if not isinstance(max_displace,float):
                raise TypeError('Max displacement for volume move '
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
            raise TypeError('Insertion probability information not '
                    'formatted properly')
        if len(insert) != 2:
            raise TypeError('Insertion probability information not '
                    'formatted properly')
        if not isinstance(insert[0],float):
            raise TypeError('Probability of insertion move must be '
                    'a floating point value')
        if not isinstance(insert[1],list):
            raise TypeError('Insertion probability information not '
                    'formatted properly')
        for insertable in insert[1]:
            if not isinstance(insertable,bool):
                raise TypeError('Whether or not a species is insertable'
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
            raise TypeError('Swap probability information not '
                    'formatted properly')
        if len(swap) != 4:
            raise TypeError('Swap probability information not '
                    'formatted properly')
        if not isinstance(swap[0],float):
            raise TypeError('Probability of swap move must be '
                    'a floating point value')
        if not isinstance(swap[1],list):
            raise TypeError('Insertion probability information not '
                    'formatted properly')
        for insertable in swap[1]:
            if not isinstance(insertable,bool):
                raise TypeError('Whether or not a species is insertable'
                    'must be a boolean value')
        if swap[2] is not None:
            if not isinstance(swap[2],list):
                raise TypeError('Swap probability information not '
                    'formatted properly')
            for prob in swap[2]:
                if not isinstance(prob,float):
                    raise TypeError('Probability of selecting species '
                            'for a swap move must be a floating point '
                            'value')
        if swap[3] is not None:
            if not isinstance(swap[3],list):
                raise TypeError('Swap probability information not '
                    'formatted properly')
            for prob in swap[3]:
                if not isinstance(prob,float):
                    raise TypeError('Probability of selecting box '
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

def simulation_length_info(units,prop_freq,coord_freq,length,
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


