The Moves object
================

The ``mc.Moves`` object contains attributes that can be grouped into
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
* ``max_volume`` - maximum volume move size (except for ``gemc_npt``, where it is per-box)

**Attributes specified per-species:**

* ``insertable`` - boolean, is species insertable
* ``prob_regrow_species`` - probability of attempting a regrowth move with each species
* ``prob_swap_species`` - probability of attempting a swap move with each species

**Attributes specified per-box:**

* ``prob_swap_from_box`` - probability of selecting each box as donor for a swap move

**Attributes specified per-box-per-species:**

* ``max_translate`` - maximum translation distance
* ``max_rotate`` - maximum rotation angle

The ``ensemble`` and ``species_list`` are used to determine the default
values assigned to the other attributes of the ``mc.Moves`` object.

Printing the contents of the Moves object
+++++++++++++++++++++++++++++++++++++++++

Imagine we have created a moves object as follows:

.. code-block:: python

  moves = mc.Moves('nvt', species_list)

We can then print the current contents with:

.. code-block:: python

  moves.print()

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
  
  
  Per species quantities:
  
                           species1     
                           ========     
  Max translate (Ang):     2.00          (box 1)
  Max rotate (deg):        30.00         (box 1)
  Insertable:              True          
  Max dihedral:            0.00          
  Prob swap:               0.00          
  Prob regrow:             1.00          
  
  
  Max volume (Ang^3):
  Box 1: 0.0


Default move probabilities:
+++++++++++++++++++++++++++

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
~~~~~

* ``prob_translate = 0.30``
* ``prob_rotate = 0.30``
* ``prob_regrow = 0.295``
* ``prob_swap = 0.1``
* ``prob_volume = 0.005``


Default move sizes:
+++++++++++++++++++

``max_translate`` and ``max_rotate`` are the per-box-per-species maximum
translation distances (in Angstroms) and maximum rotation angles (in degrees).
The default maximum translation and rotation are 2.0 Angstrom and 30.0 degrees,
respectively. For example, if the system contained two species and the ensemble
was GEMC (a two-box ensemble), then the default max translate would be
``[[2.0,2.0],[2.0,2.0]]``. To set the max translation distance of species 1 in
box 2 to 30.0 Angstroms, set ``max_translate = [[2.0,2.0],[30.0,2.0]]``.

.. note::
    Exceptions to the above values are implemented based upon the topologies
    provided in ``species_list``. The maximum rotation of single particle
    species is set to ``0.0`` degrees. Species that are multi-particle but
    contain zero bonds are considered fixed; the maximum translation
    and rotation are set to ``0.0`` Angstroms and ``0.0`` degrees, respectively.
