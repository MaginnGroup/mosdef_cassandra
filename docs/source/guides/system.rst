
The System object
=================

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
occupied, the list element is a ``mbuild.Compound``. The length of ``box_list``
will normally be one (NVT, NPT, GCMC) or two (GEMC).

``species_list`` is a ``list`` of the unique chemical species in the system. For
example, if simulating a mixture of methane and ethane, there are two species.
Therefore, the length of ``species_list`` would be two. Each element of
``species_list`` is a parameterized ``parmed.Structure``. This structure
contains all the forcefield details for a given species.

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

