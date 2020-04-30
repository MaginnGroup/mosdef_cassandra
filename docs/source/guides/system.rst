
The System object
=================

The ``System`` contains all the details of the system (i.e.,
*what* is being simulated). This includes the simulation box(es), any initial
structure(s) in the simulation box(es), and the force field parameters
describing the interactions between particles in the system.

The ``System`` is one of two objects that must be created prior to running
a MC simulation. When the system object is created, mosdef_cassandra performs
several consistency checks and ensures all the required details are
provided.

Creating the ``mc.System`` object requires specifying the following:

* A list of the simulation boxes in the system (``box_list``)
* A list of the unique chemical species in the system (``species_list``)

and perhaps,

* The number of molecules in each box (``mols_in_boxes``)
* The number of molecules for Cassandra to add (``mols_to_add``)

These items comprise a complete description of the system to be simulated in
Cassandra. The box information, forcefield information, number of species,
and coordinates of any initial structure are contained within this object.

box_list
~~~~~~~~

The ``box_list`` is a Python ``list`` of the simulation boxes in the system.
It should contain a single item in the case of Monte Carlo simulations
performed in the NVT, NpT, or GCMC ensembles, and two items for simulations in
the GEMC or GEMC-NpT ensembles.

The simulation boxes can be empty or contain some initial structure. If the
simulation box is empty, then the list element should be a ``mbuild.Box``.
An ``mbuild.Box`` can be created as follows:

.. code-block:: python

  box = mbuild.Box(lengths=[3.0, 3.0, 3.0], angles=[90., 90., 90.])


where the lengths are specified in units of nanometers and the box angles
are specified in degrees. If the angles are not specified they are assumed to
be 90 degrees.

If the simulation box is not empty, but rather contains some initial structure,
then the list element for that box should be a ``mbuild.Compound`` object.

For a single-box simulation with an initially empty simulation box:

.. code-block:: Python

  box = mbuild.Box([3.,3.,3.])
  box_list = [box]

For a two-box simulation where one box is initially empty and the other
contains a structure from a PDB file and loaded with ``mbuild``.

.. code-block:: Python

  zeolite = mbuild.load("zeolite.pdb")
  vapor_box = mbuild.Box([3.,3.,3.])

  box_list = [zeolite, vapor_box]

  .. warning::

    If an initial structure (i.e., an ``mbuild.Compound``) is provided, the
    order of the atoms must follow a *very* specific order. Each complete
    molecule must appear one at a time. The order of the atoms in each molecule
    must *exactly* match the order of the atoms in the relevant species provided
    in the ``species_list``. If there are multiple different species, then all
    the molecules of species1 must be provided before any molecules of species2,
    and so on. 

.. note::

  The box lengths are taken from the ``periodicity`` attribute of the
  ``mbuild.Compound`` object. The box angles are taken from
  ``mbuild.Compound.boundingbox.angles``. This is temporary solution due to the
  fact the ``mbuild.Compound.boundingbox.lengths`` attribute is calculated
  on-the-fly from the extent of the particles in the ``Compound`` rather than
  storing periodic box information, while the ``Compound.perodocitiy`` attribute
  contains no information regarding the box angles.

species_list
~~~~~~~~~~~~

The ``species_list`` is a Python ``list`` of the unique chemical species in the
system. For example, if simulating a mixture of methane and ethane, there are
two species; if simulating pure methane, there is one species, regardless of
the number of methane molecules in the system. Therefore, the
length of ``species_list`` would be two in the first example and one in the
second example. Each element of ``species_list`` is a parameterized
``parmed.Structure``. This structure contains all the forcefield details for a
given species.

For example, if we are simulating a mixture of methane and ethane with the
OPLS-AA force field, we could use the following sequence of steps to generate
the species list. Note that with mbuild and foyer we are able to generate a
molecule with force field parameters starting from a SMILES string with only
a few lines of Python code.

.. code-block:: python

  import mbuild
  import foyer

  methane = mbuild.load("C", smiles=True)
  ethane = mbuild.load("CC", smiles=True)

  ff = foyer.forcefields.load_OPLSAA()

  methane_ff = ff.apply(methane)
  ethane_ff = ff.apply(ethane)

  species_list = [methane_ff, ethane_ff]

.. note::

  The order of the items in species list determines the labeling of
  the species. The first is considered species1, the second species2, and
  so forth.

mols_in_boxes
~~~~~~~~~~~~~

The ``mols_in_boxes`` is a ``list`` containing the number of molecules of each
species currently in each box specified in ``box_list``. If the simulation
box(es) are empty, then ``mols_in_boxes`` does not need to be specified. If
specified, it is provided as a nested list with ``shape=(n_boxes, n_species)``.

For example, for a system with one box and one species where the
initial structure contains 100 molecules:

.. code-block:: Python

  mols_in_boxes = [[100]]

For a single box with two species, 25 molecules of the first species and 75
molecules of the second species:

.. code-block:: Python

  mols_in_boxes = [[25, 75]]

For two boxes with a single species, where the first box has 100 molecules and
the second box is empty:

.. code-block:: Python

  mols_in_boxes = [[100], [0]]

For two boxes with two species, where the first box has 300 molecules of
the first species and 50 molecules of the second species, and the second box
has 30 molecules of the first species and 100 molecules of the second:

.. code-block:: Python

  mols_in_boxes = [[300, 50], [30, 100]]

When the ``System`` object is created, it verifies that the number of atoms
provided in each box match the number of atoms specified by ``mols_in_boxes``.
The number of atoms per molecule are determined from the species provided
in the ``species_list``.

``mols_in_boxes = [[100]]``. For a
single box with two species, 25 molecules of the first species and 75 molecules
of the second species ``mols_in_boxes = [[25, 75]]``

``len(mols_in_boxes)`` must equal
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
