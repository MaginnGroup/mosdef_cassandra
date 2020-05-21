
System
======

The ``System`` contains all the details of the system (i.e.,
*what* is being simulated). This includes the simulation box(es), any initial
structure(s) in the simulation box(es), and the force field parameters
describing the interactions between particles in the system.

The ``System`` is one of two objects that must be created prior to running
a MC simulation. Creating the ``System`` requires specification of the following:

* A list of the simulation boxes (``box_list``)
* A list of the unique chemical species (``species_list``)

and perhaps,

* The number of molecules already present in each box
  in the ``box_list`` (``mols_in_boxes``)
* The number of molecules for Cassandra to add to each box
  prior to beginning the MC simulation (``mols_to_add``)

Instantiating the ``System`` normally appears as follows:

.. code-block:: python

  import mosdef_cassandra as mc
  system = mc.System(
      box_list,
      species_list,
      mols_in_boxes=mols_in_boxes,
      mols_to_add=mols_to_add
  )

These items comprise a complete description of the system to be simulated in
Cassandra. The box information, forcefield information, number of species,
and coordinates of any initial structure are contained within this object.

box_list
~~~~~~~~

The ``box_list`` is a Python ``list`` of the simulation boxes in the system.
It should contain a single item in the case of simulations
performed in the NVT, NPT, or GCMC ensembles, and two items for simulations in
the GEMC or GEMC-NPT ensembles.

Each simulation box can be empty or contain an initial structure. If a
simulation box is empty, then the list element should be an ``mbuild.Box``.
An ``mbuild.Box`` can be created as follows:

.. code-block:: python

  box = mbuild.Box(lengths=[3.0, 3.0, 3.0], angles=[90., 90., 90.])


where the lengths are specified in units of nanometers and the box angles
are specified in degrees. If the angles are not specified they are taken as
90 degrees.

If a simulation box contains an initial structure, then the list
element should be an ``mbuild.Compound`` object. mBuild supports reading
many common simulation file formats via ``mbuild.load``. See the
`mBuild documentation <https://mbuild.mosdef.org/en/stable/>`_
for more details.

.. warning::

  If an initial structure (i.e., an ``mbuild.Compound``) is provided, the
  order of atoms is *very* important. Each complete molecule
  must appear one after another. Within each molecule, the order of 
  atoms in *must match* the order of atoms in the
  relevant species provided in the ``species_list``. If there are
  multiple different species, then all molecules of species1 must
  be provided before any molecules of species2, and so on. We hope
  to relax these restrictions in future releases.

For a single-box simulation with an initially empty simulation box:

.. code-block:: Python

  box = mbuild.Box([3.,3.,3.])
  box_list = [box]

For a two-box simulation where one box is initially empty and initial
structure for the other simulation box is loaded from a PDB file
with ``mbuild``:

.. code-block:: Python

  zeolite_box = mbuild.load("zeolite.pdb")
  vapor_box = mbuild.Box([3.,3.,3.])

  box_list = [zeolite_box, vapor_box]

In this case, the initial structure and box dimensions for the box
containing the zeolite are taken from the PDB file. Note that
the box dimensions can be manually edited by changing the
``mbuild.Compound.periodicity`` attribute.

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
system. For example, a simulation of pure methane contains one unique chemical
species (methane), regardless of the number of methane molecules in the
simulation. A simulation containing a mixture of methane and ethane has two
unique chemical species. Therefore, in the first example, the
``species_list`` contains a single item and in the second example the
``species_list`` contains two items. Each item in the ``species_list`` is
a ``parmed.Structure``. All the forcefield required force field
parameters for each species must be in their respective
``parmed.Structure``.

.. note::
  
  The ``parmed.Structure`` will be replaced with a
  ``gmso.Topology`` as the GMSO package matures.

For example, to simulate a mixture of methane and ethane with the
OPLS-AA force field, we could use the following sequence of steps to generate
the species list. Note that mbuild and foyer allow us to generate a
molecule with force field parameters from a SMILES string and a few lines
of Python code.

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

  The order of items in species list determines the labeling of
  the species. The first is considered species1, the second species2, and
  so forth.

mols_in_boxes
~~~~~~~~~~~~~

The ``mols_in_boxes`` is a ``list`` containing the number of molecules of each
species currently in each box specified in ``box_list``. If all simulation
box(es) are empty, ``mols_in_boxes`` does not need to be specified. When
specified, it is a nested list with ``shape=(n_boxes, n_species)``.
This is perhaps easier to explain with a few examples.

Consider a system with one simulation box and one species.
If the initial structure provided in ``box_list`` contains
100 molecules of that species, then:

.. code-block:: Python

  mols_in_boxes = [[100]]

For a system with one simulation box and two species, where there
are 25 molecules of the first species and 75 molecules of the
second species:

.. code-block:: Python

  mols_in_boxes = [[25, 75]]

For a system with two simulation boxes and one species, where the first box
contains 100 molecules and the second box is empty:

.. code-block:: Python

  mols_in_boxes = [[100], [0]]

For a system with two boxes and two species; the first box has 300 molecules of
the first species and 50 molecules of the second species, the second box
has 30 molecules of the first species and 100 molecules of the second:

.. code-block:: Python

  mols_in_boxes = [[300, 50], [30, 100]]

When the ``System`` object is created, it verifies that the number of atoms
provided in each box match the number of atoms specified by ``mols_in_boxes``.
The number of atoms per molecule are determined from the species provided
in the ``species_list``.

mols_to_add
~~~~~~~~~~~~~
Cassandra can insert molecules in a simulation box prior to
starting an MC simulation. Therefore, you can provide an empty simulation
box and request Cassandra to add some number of molecules before beginning the
simulation. This capability is controlled through the ``mols_to_add`` option.
The format of ``mols_to_add`` is analogous to ``mols_in_boxes``. If
specified, it is provided as a nested list with ``shape=(n_boxes, n_species)``.


For example, consider a system with a single simulation box and two species.
If we wish to add 10 molecules of the first species and 0 molecules
of the second species, we could use:

.. code-block:: Python

  mols_to_add = [[10,0]]

.. warning::
  If ``mols_to_add`` is too large for the given box/species, the MC simulation
  may never begin. Cassandra will be stuck attempting (and failing) to insert
  the requested number of molecules.
