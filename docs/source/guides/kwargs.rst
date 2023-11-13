Keyword Arguments
=================

Nearly all options of MoSDeF Cassandra can be controlled through
the use of keyword arguments to the ``run``/``restart`` functions.
These arguments can be specified individually or provided to the
``run``/``restart`` functions via a dictionary. The
dictionary-based approach is preferred if there are a large
number of keyword arguments to keep the number of explicit
arguments to the ``run``/``restart`` functions manageable.

Usage
+++++

Below is an example of providing the ``vdw_cutoff`` option
to ``run`` as an extra keyword argument.

.. code:: python

  mc.run(
    system=system,
    moveset=moveset,
    run_type="equilibration",
    run_length=1000,
    temperature=300.0 * u.K,
    vdw_cutoff=9.0 * u.angstroms
  )

or as a dictionary, where the ``**`` operator is used to
expand the dictionary.

.. code:: python

  custom_args = {
    'vdw_cutoff': 9.0 * u.angstroms,
    'charge_cutoff': 9.0 * u.angstroms,
  }

  mc.run(
    system=system,
    moveset=moveset,
    run_type="equilibration",
    run_length=1000,
    temperature=300.0 * u.K,
    **custom_args
  )

Valid arguments
+++++++++++++++

A list of the valid keyword arguments is provided below
with a brief explanation. If more detail is required, please
consult the Cassandra user manual. Most arguments below have
nearly a one-to-one mapping with options of the Cassandra
input file.

``run_name``
~~~~~~~~~~~~
| **Type:** ``str``
| **Description:** run name to prepend to output files
| **Default:** ensemble of simulation (i.e., ``"nvt"``, ``"npt"``, etc.)


``restart``
~~~~~~~~~~~
| **Type:** ``bool``
| **Description:** if ``True``, restart from a Cassandra ``.chk`` file
| **Default:** ``False``


``restart_name``
~~~~~~~~~~~~~~~~
| **Type:** ``str``
| **Description:** name of the checkpoint file without the extension, i.e.,
  the ``run_name`` for the simulation from which you wish to restart
| **Default:**
| **Notes:** only relevant if ``restart=True``


``verbose_log``
~~~~~~~~~~~~~~~
| **Type:** ``bool``
| **Description:** if ``True``, print the Cassandra log file with additional verbosity
| **Default:** ``False``


``vdw_style``
~~~~~~~~~~~~~
| **Type:** ``str``
| **Description:** type of van der Waals interactions. Valid options include ``lj`` or ``none``
| **Default:** ``lj``


``cutoff_style``
~~~~~~~~~~~~~~~~
| **Type:** ``str``
| **Description:** method of handling the cutoff for the van der Waals interactions. Valid
   options include ``cut_tail``, ``cut_switch``, ``cut_shift``
| **Default:** ``cut_tail``

``vdw_cutoff``
~~~~~~~~~~~~~~
| **Type:** ``unyt_quantity, dimensions=length``, except for ``cutoff_style="cut_switch"``, which
  requires a list of, ``[inner_cutoff, outer_cutoff]``.
| **Description:** cutoff distance for van der Waals interactions
| **Default:** ``12.0 * u.angstroms``
| **Notes:** In a system with multiple boxes, per
  box values can be specified with ``vdw_cutoff_box1`` and ``vdw_cutoff_box2``
  keywords. If provided, these will override the ``vdw_cutoff``.

``charge_style``
~~~~~~~~~~~~~~~~
| **Type:** ``str``
| **Description:** method of computing electrostatic energy, options include ``none``,
    ``ewald``, or ``dsf``
| **Default:** ``ewald``


``charge_cutoff``
~~~~~~~~~~~~~~~~~
| **Type:** ``unyt_quantity, dimensions=length``
| **Description:** cutoff distance for short-range portion of charged interactions
| **Default:** ``12.0 * u.angstroms``
| **Notes:** In a system with multiple boxes, per
  box values can be specified with ``charge_cutoff_box1`` and ``charge_cutoff_box2``
  keywords. If provided, these will override the ``charge_cutoff``. In GEMC simulations
  where the vapor box is much larger than the liquid box, it may be necessary to increase
  the charge cutoff of the vapor box to maintain the desired ``ewald_accuracy`` without
  exceeding the maximum number of k-space vectors.

``ewald_accuracy``
~~~~~~~~~~~~~~~~~~
| **Type:** ``float``
| **Description:** relative accuracy of ewald summation
| **Default:** ``1.0e-5``
| **Notes:** Only relevant if ``charge_style="ewald"``

``dsf_damping``
~~~~~~~~~~~~~~~
| **Type:** ``float``
| **Description:** damping parameter for ``dsf`` charge style
| **Default:** ``None``
| **Notes:** Only relevant if ``charge_style="dsf"``


``mixing_rule``
~~~~~~~~~~~~~~~
| **Type:** ``str``
| **Description:** the type of mixing rule to apply to van der Waals interactions. Options include
  ``lb`` (Lorentz-Berthelot), ``geometric`` or ``custom``
| **Default:** ``lb``

``custom_mixing_dict``
~~~~~~~~~~~~~~~~~~~~~~
| **Type:** ``dict``
| **Description:** dictionary specifying the custom mixing rules. One key-value pair
 is specified per pair of atomtypes. The key is a string of the species combination,
 and the value is a list of the relevant parameters. For example, the two atom types
 are ``opls_140`` and ``opls_141`` and the mixed epsilon and sigma are ``10.0 * u.Unit('kJ/mol')``
 and ``3.0 * u.angstrom``, then the ``dict`` would be:

 .. code-block:: python
   
    { 'opls_140 opls_141': [10.0 * u.Unit('kJ/mol'), 3.0 * u.angstrom] }

| **Default:** ``None``

``seeds``
~~~~~~~~~
| **Type:** ``list`` of two ``ints``
| **Description:** the starting seeds for the random number generator.
| **Default:** selected at random

``rcut_min``
~~~~~~~~~~~~
| **Type:** ``unyt_quantity, dimensions = length``
| **Description:** minimum distance to calculate interaction energy. If particles are
  closer than this distance the energy is taken as infinity and the move is automatically
  rejected. If the value is too large moves that might possibly be accepted will be
  unecessarily rejected.
| **Default:** ``1.0 * u.angstrom``

``pair_energy``
~~~~~~~~~~~~~~~
| **Type:** ``bool``
| **Description:** store pair interactions energies (requires more memory but may be faster)
| **Default:** ``True``

``max_molecules``
~~~~~~~~~~~~~~~~~
| **Type:** ``list`` of ``ints``, ``len=n_species``
| **Description:** maximum number of molecules of each species. Cassandra will
  exit if the number of molecules of a species exceeds this number at any point
  during a simulation. 
| **Default:** Number of molecules in the ``System`` for ``nvt``, ``npt``,
  ``gemc``, ``gemc_npt``, and non-insertable species in ``gcmc``. Number of
  molecules in the ``System`` plus 500 for insertable molecules in ``gcmc``.
| **Notes:** The default may need to be overridden in GCMC if the
  initial configuration has many fewer molecules than at equilibrium.


``pressure``
~~~~~~~~~~~~
| **Type:** ``unyt_quantity``, valid units of pressure 
| **Description:** desired pressure (NPT or GEMC-NPT) ensembles
| **Default:** ``None``
| **Notes:** in GEMC-NPT, different pressures for ``box1`` and ``box2`` can be specified
  with the ``pressure_box1`` and ``pressure_box2``. If specified, these values will override
  the value in ``pressure``.


``chemical_potentials``
~~~~~~~~~~~~~~~~~~~~~~~
| **Type:** ``list`` of ``unyt_array/unyt_quantity`` with units of ``energy/mol``, or ``"none"``
  for species that are not insertable
| **Description:** specify the desired chemical potential for each species (``gcmc``)
| **Default:** ``None``

``thermal_stat_freq``
~~~~~~~~~~~~~~~~~~~~~
| **Type:** ``int``
| **Description:** frequency, in number of thermal moves, of printing statistics and (if ``run_type="equilibration"``),
    updating the maximum translation and rotation sizes
| **Default:** ``1000``
| **Notes:** in ``equilibration`` mode, the maximum translation and rotation move sizes
  are continuously adjusted to target 50% of moves accepted.


``vol_stat_freq``
~~~~~~~~~~~~~~~~~
| **Type:** ``int``
| **Description:** frequency, in number of volume moves, of printing statistics and (if ``run_type="equilibration"``),
    updating the maximum volume move size
| **Default:** ``100``
| **Notes:** in ``equilibration`` mode, the maximum volume move size
  is continuously adjusted to target 50% of moves accepted.


``units``
~~~~~~~~~
| **Type:** ``str``
| **Description:** units for measuring simulation length, valid options include ``minutes``, ``steps``, or ``sweeps``
| **Default:** ``steps``

``steps_per_sweep``
~~~~~~~~~~~~~~~~~~~
| **Type:** ``int``
| **Description:** the number of MC steps in one MC sweep
| **Default:** ``None``
| **Notes:** required if ``units="steps"``. A standard choice is one sweep is one attempted move per molecule in the system.

``prop_freq``
~~~~~~~~~~~~~
| **Type:** ``int``
| **Description:** frequency of writing thermo properties to the ``.prp`` file
| **Default:** 500
| **Notes:** units determined by the ``units`` argument

``coord_freq``
~~~~~~~~~~~~~~
| **Type:** ``int``
| **Description:** frequency of writing system coordinates to the ``.xyz`` file
| **Default:** 5000
| **Notes:** units determined by the ``units`` argument

``block_avg_freq``
~~~~~~~~~~~~~~~~~~
| **Type:** ``int``
| **Description:** block average size 
| **Default:** ``None``
| **Notes:** units determined by the ``units`` argument

``properties``
~~~~~~~~~~~~~~
| **Type:** ``list`` of ``str``
| **Description:** list of properties to write to the ``.prp`` file. Valid options include: ``energy_total``, ``energy_intra``, ``energy_bond``, ``energy_angle``, ``energy_diheral``, ``energy_improper``, ``energy_intravdw``, ``energy_intraq``, ``energy_inter``, ``energy_intervdw``, ``energy_lrc``, ``energy_interq``, ``energy_recip``, ``energy_self``, ``enthalpy``, ``pressure``, ``pressure_xx``, ``pressure_yy``, ``pressure_zz``, ``volume``, ``nmols``, ``density``, ``mass_density``.
| **Default:** ``["energy_total", "energy_intra", "energy_inter", "enthalpy", "pressure", "volume", "nmols", "mass_density"]``

``widom_insertions``
~~~~~~~~~~~~~~~~~~~~
| **Type:** ``list`` of ``dicts``
| **Description:** One ``dict`` per box.  The dictionary keys are the species numbers of the Widom test particle species, and each dictionary entry is a list of two ``ints``: ``[n_ins, widom_freq, n_subgroups]``, where ``n_ins`` is the number of Widom insertions to be performed after every ``widom_freq`` MC steps (or MC ``sweeps`` if ``units="sweeps"``) and ``n_subgroups`` is the number of Widom insertion subgroups per Widom insertion frame.
| **Default:** ``None``
| **Notes:** units of ``widom_freq`` cannot be time units, so they default to ``steps`` if ``units="minutes"``.


``cell_list``
~~~~~~~~~~~~~~~~~~~~
| **Type:** ``bool`` or ``str``
| **Description:** ``True`` if cell list overlap detection is to be used for Widom insertions.
| **Default:** ``False``

``adaptive_rmin``
~~~~~~~~~~~~~~~~~~~~
| **Type:** ``bool``, ``int``, or ``float``
| **Description:** Maximum desired intermolecular nonbonded atom pair energy, normalized by kBT.  Setting as ``True`` sets no value in the input file, leaving it as the Cassandra's default value (708.0).
| **Default:** ``False``

