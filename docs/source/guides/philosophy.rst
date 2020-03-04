
Concepts
========

The goal of mosdef_cassandra is to integrate the `Cassandra Monte Carlo (MC)
code <https://cassandra.nd.edu>`_ with the `Molecular Simulation Design
Framework (MoSDeF) <https://mosdef.org>`_. The Cassandra
`documentation <https://cassandra.nd.edu/index.php/documentation>`_
contains a fairly comprehensive background on MC simulations. Canonical
textbooks by `Frenkel and Smit
<https://www.sciencedirect.com/book/9780122673511/understanding-molecular-simulation>`_,
`Allen and Tidesley
<https://www.oxfordscholarship.com/view/10.1093/oso/9780198803195.001.0001/oso-9780198803195>`_,
and `Tuckerman <https://onlinelibrary.wiley.com/doi/10.1002/anie.201105752>`_ also provide
useful reference materials. Here however, we focus on the mosdef_cassandra package,
proceeding with the assumption that the reader has a basic understanding of MC methods as
applied in molecular simulations.

mosdef_cassandra divides the process of setting up a MC simulation into the
creation of two objects, the ``mosdef_cassandra.System`` object
and the ``mosdef_cassandra.Moves`` object. Each is described in detail
below. First, we provide a bit of context for the mosdef_cassandra package.

Philosophy
----------

Setting up a Monte Carlo simulation should not be difficult. The simulation
setup procedure should not be prone to error. The process should be easily
repeatible and extensible. And none of the prior goals should sacrifice the
complete flexibility required by the expert simulator. These simple ideas
underpin the development of mosdef_cassandra.

What are the components of a Monte Carlo simulation?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

What is the simplest organization of those components?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Where do we go from here?
~~~~~~~~~~~~~~~~~~~~~~~~~

The System object
-----------------

The ``mc.System`` object contains the following:

* ``box_list``
* ``species_list``
* ``mols_in_boxes``
* ``mols_to_add``

These items comprise a complete description of the system to be simulated in
Cassandra. The box information, forcefield information, number of species,
and coordinates of any initial structure are contained within this object.

``box_list`` is a ``list`` of the simulation boxes. For boxes that are initially
empty, the element of the list is a ``mbuild.Box``. For boxes that are already
occupied, the list element is a ``mbuild.Compound``. The length of ``box_list`` will
normally be one (NVT, NPT, GCMC) or two (GEMC).

``species_list`` is a ``list`` of the unique chemical species in the system. For
example, if simulating a mixture of methane and ethane, there are two species.
Therefore, the length of ``species_list`` would be two. Each element of ``species_list`` is
a parameterized ``parmed.Structure``. This structure contains all the forcefield details for
a given species.

``mols_in_boxes`` is a ``list`` containing the number of molecules of each
species currently existing in each box. ``len(mols_in_boxes)`` must equal
``n_boxes = len(box_list)``. In other words, ``n_boxes`` is the number of
simulation boxes. For each box ``ibox``, there is a list
with ``len(mols_in_boxes[ibox]) = n_species``, where ``n_species`` is the number
of unique species in the system. If setting up a system with a single
simulation box containing 10 methane molecules and 5 ethane molecules,
``mols_in_boxes = [[10,5]]``.

``mols_to_add`` is a ``list`` containing the number of molecules of each
species to add each box before beginning the simulation. The format of
``mols_to_add`` is analogous to the format of ``mols_in_boxes``.
``len(mols_to_add) = n_boxes``, and for each box, there is a list
with ``len(mols_to_add[ibox]) = n_species``. If setting up a system with a
single simulation box to which we wish to add 10 methane molecules and 0 ethane
molecules, ``mols_to_add = [[10,0]]``.


The Moves object
----------------

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


