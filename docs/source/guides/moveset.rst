MoveSet
=======

The ``MoveSet`` contains all the information related to the
Monte Carlo moves that will be attempted during the simulation.
This includes the types of moves, the probability of selecting
each move type, and the other related choices, such as the
maximum translation distance, maximum volume move size,
configurational biasing options, etc. The desired ensemble and
the species in the system are used to assign default values to
all of the attributes of the ``MoveSet``. Nonetheless, most
attributes can be edited once the object has been created.

The ``MoveSet`` is one of two objects that must be created prior to
running a MC simulation in MoSDeF Cassandra. Creating the ``MoveSet``
requires specification of the following:

* The desired ensemble
* A list of the unique chemical species in the system (species_list)

Instantiating the ``MoveSet`` normally appears as follows:

.. code-block:: python

  import mosdef_cassandra as mc
  moveset = mc.MoveSet("nvt", species_list)

Attributes
++++++++++

The ``MoveSet`` contains attributes that can be grouped into
the following four categories.

**Overall attributes**, specified as a single number for the entire system:

* ``ensemble`` - ensemble of the MC simulation (``nvt``, ``npt``, ``gcmc``, ``gemc``, ``gemc_npt``)
* ``prob_translate`` - probability of attempting a translation move
* ``prob_rotate`` - probability of attempting a rotation move
* ``prob_angle`` - probability of attempting an angle change move
* ``prob_dihedral`` - probability of attempting a dihedral change move
* ``prob_regrow`` - probability of attempting a regrowth move
* ``prob_volume`` - probability of attempting a volume change move
* ``prob_insert`` - probability of attempting a molecule insertion move
* ``prob_swap`` - probability of attempting a molecule swap move
* ``max_volume`` - maximum volume move size (except for ``gemc_npt``, where it is optionally per-box)
* ``cbmc_n_insert`` - number of locations to attempt a CBMC insertion
* ``cbmc_n_dihed`` - number of dihedral angles to attempt when regrowing a molecule with CBMC
* ``cbmc_rcut`` - cutoff to use when calculating energies during CBMC trials (optionally specified per-box)

**Attributes specified per-species:**

* ``insertable`` - boolean, is species insertable
* ``prob_regrow_species`` - probability of attempting a regrowth move with each species
* ``prob_swap_species`` - probability of attempting a swap move with each species
* ``max_dihedral`` - maximum dihedral angle change for dihedral change move

**Attributes specified per-box:**

* ``prob_swap_from_box`` - probability of selecting each box as donor for a swap move

**Attributes specified per-box-per-species:**

* ``max_translate`` - maximum translation distance
* ``max_rotate`` - maximum rotation angle


Printing the contents of the MoveSet
++++++++++++++++++++++++++++++++++++

Imagine we have created a ``MoveSet`` as follows:

.. code-block:: python

  moveset = mc.MoveSet('nvt', species_list)

We can then print the current contents with:

.. code-block:: python

  moveset.print()

Example output for a single species (OPLS-AA methane):

.. code-block::


  Ensemble:  nvt
  
  Probability of selecting each move type:
  
      Translate: 0.33
      Rotate:    0.33
      Regrow:    0.34
      Volume:    0.0
      Insert:    0.0
      Delete:    0.0
      Swap:      0.0
      Angle:     0.0
      Dihedral:  0.0
  
  CBMC selections:
  
      Number of trial positions: 10
      Number of trial dihedral angles: 10
      CBMC cutoff(s):
          Box 1: 6.0
  
  
  Per species quantities:
  
                               species1
                               ========
      Max translate (Ang):     2.00          (Box 1)
      Max rotate (deg):        30.00         (Box 1)
      Insertable:              False
      Max dihedral:            0.00
      Prob swap:               0.00
      Prob regrow:             1.00
  
  
  Max volume (Ang^3):
      Box 1: 0.0


Default values for attempting each move type
++++++++++++++++++++++++++++++++++++++++++++

``prob_translate``, ``prob_rotate``, ``prob_angle``, ``prob_dihedral``,
``prob_regrow``, ``prob_volume``, ``prob_insert``, and ``prob_swap`` are the
probabilities of selecting each of those respective move types. The default
move probabilities are as follows for each ensemble. Move probabilities that are
not explicitly defined have a default probability of 0.0 for that ensemble.


NVT:
~~~~

* ``prob_translate = 0.33``
* ``prob_rotate = 0.33``
* ``prob_regrow = 0.34``

NPT:
~~~~

* ``prob_translate = 0.33``
* ``prob_rotate = 0.33``
* ``prob_regrow = 0.335``
* ``prob_volume = 0.005``

GCMC:
~~~~~

* ``prob_translate = 0.25``
* ``prob_rotate = 0.25``
* ``prob_regrow = 0.30``
* ``prob_insert = 0.1``

.. note::
    In GCMC the deletion probability is set equal to the insertion
    probability, making the sum of the move probabilities 1.0

GEMC:
~~~~~

* ``prob_translate = 0.30``
* ``prob_rotate = 0.30``
* ``prob_regrow = 0.295``
* ``prob_swap = 0.1``
* ``prob_volume = 0.005``

GEMC-NPT:
~~~~~~~~~

* ``prob_translate = 0.30``
* ``prob_rotate = 0.30``
* ``prob_regrow = 0.295``
* ``prob_swap = 0.1``
* ``prob_volume = 0.005``


Default values for other quantities
+++++++++++++++++++++++++++++++++++

* ``max_translate``: 2.0 Angstroms
* ``max_rotate`` : 30.0 degrees
* ``max_volume`` : 500 Angstroms\ :sup:`3` for Box 1, 5000 Angstroms\ :sup:`3` for Box 2
* ``max_dihedral`` : 0.0 degrees
* ``cbmc_n_insert`` : 10
* ``cbmc_n_dihed`` : 10
* ``cbmc_rcut`` : 6.0 Angstroms


``max_translate`` and ``max_rotate`` are specified per-box-per-species.
For example, if the system contained two species and the ensemble
was GEMC (a two-box ensemble), then the default max translate would be
``[[2.0,2.0],[2.0,2.0]]``. To set the max translation distance of species 1 in
box 2 to 30.0 Angstroms, set ``max_translate = [[2.0,2.0],[30.0,2.0]]``.

.. note::
    Exceptions to the above values are implemented based upon the topologies
    provided in ``species_list``. The maximum rotation of single particle
    species is set to ``0.0`` degrees. Species that are multi-particle but
    contain zero bonds are considered fixed and not insertable; the maximum translation
    and rotation are set to ``0.0`` Angstroms and ``0.0`` degrees, respectively.
