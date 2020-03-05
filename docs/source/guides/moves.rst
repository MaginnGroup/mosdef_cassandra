The Moves object
================

The ``mc.Moves`` object contains the following:

* ``ensemble``
* ``prob_translate``
* ``prob_rotate``
* ``prob_angle``
* ``prob_dihedral``
* ``prob_regrow``
* ``prob_volume``
* ``prob_insert``
* ``prob_swap``
* ``max_translate``
* ``max_rotate``
* ``max_volume``
* ``sp_insertable``
* ``sp_prob_regrow``
* ``sp_prob_swap``
* ``prob_swap_from_box``

``ensemble`` is a string describing the desired ensemble for the simulation. The
ensemble and ``species_list`` determine the default values assigned to all the
other attributes of the ``mc.Moves`` object.

Printing the contents of the Moves object
+++++++++++++++++++++++++++++++++++++++++

Imagine we have created a moves object as follows:

.. code-block:: python

  moves = mc.Moves('nvt', species_list)

We can then print the current contents with:

.. code-block:: python

  moves.print()

Default move probabilities:
+++++++++++++++++++++++++++

``prob_translate``, ``prob_rotate``, ``prob_angle``, ``prob_dihedral``,
``prob_regrow``, ``prob_volume``, ``prob_insert``, and ``prob_swap`` are the
probabilities of selecting each of those respective move types. The default
move probabilities are as follows for each ensemble. Move probabilities that are
not explicitly defined have a default probability of 0.0 for that ensemble.


NVT:
~~~~

* ``prob_translate = 0.35``
* ``prob_rotate = 0.35``
* ``prob_regrow = 0.30``

NPT:
~~~~

* ``prob_translate = 0.35``
* ``prob_rotate = 0.34``
* ``prob_regrow = 0.30``
* ``prob_volume = 0.02``

GCMC:
~~~~~

* ``prob_translate = 0.25``
* ``prob_rotate = 0.25``
* ``prob_regrow = 0.30``
* ``prob_insert = 0.1``

.. note::
    In GCMC the deletion probability is set equal to the insertion
    probability, making the sum of the move probabilities to 1.0

GEMC:
~~~~~

* ``prob_translate = 0.29``
* ``prob_rotate = 0.29``
* ``prob_regrow = 0.30``
* ``prob_swap = 0.1``
* ``prob_volume = 0.02``

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
