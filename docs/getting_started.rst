Getting Started
===============

The following assumes you have already installed mosdef_cassandra
as detailed in our installation instructions. More details on the
``System`` and ``Moves`` objects can be found under Concepts.

Let's start by setting up an NVT Monte Carlo simulation of OPLS-AA
methane. We will use the ``mbuild`` and ``foyer`` packages to create
the methane molecule, and ``mosdef_cassandra`` to run the MC simulation.
We begin with the required imports:

.. code-block:: python

    import mbuild
    import foyer
    import mosdef_cassandra as mc

Next, we create an all-atom methane molecule from a `SMILES
<https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html>`_ string:

.. code-block:: python

    methane = mbuild.load("C", smiles=True)

``methane`` is a single all-atom methane molecule. It is an
``mbuild.Compound`` object. The ``mbuild.Compound`` object contains particles
for each element (C, H, H, H) in the molecule and bond information that
describes the particle connectivity. However, there are no forcefield parameters
associated with ``methane``.

To add forcefield parameters to our ``methane``, we load the OPLS-AA forcefield
from foyer. The OPLS-AA forcefield is distributed with foyer. Note that not
all atomtypes are currently defined.

.. code-block:: python

   oplsaa = foyer.forcefields.load_OPLSAA()

We then apply the forcefield using foyer:

.. code-block:: python

    methane_ff = oplsaa.apply(methane)

``methane_ff`` is a ``parmed.Structure`` object that contains forcefield
parameters.

Now that we have a molecule with forcefield parameters, the next step is
to define our simulation box. Since Cassandra can add molecules to a
simulation box at the start of the simulation we can begin with an
empty simulation box. We will define an ``mbuild.Box`` with the box
lengths specified in nanometers:

.. code-block:: python

    box = mbuild.Box([3.0, 3.0, 3.0])

.. warning::
    Even though the default units of Cassandra are Angstroms, the
    ``mbuild.Box`` object should be specified in nanometers. This is
    to maintain consistency with the units of mbuild and foyer.

Next, we create the four primary components of the ``System`` object:
the ``box_list``, ``species_list``, ``mols_in_boxes``, and ``mols_to_add``.
The ``box_list`` is simply a list of the simulation boxes in the system.
In this case, since we are performing an NVT simulation there is only our
single ``box``. The ``species_list`` is a list of the unique chemical
species in our system. Here we only have methane molecules. ``mols_in_boxes``
describes the number of molecules in each simulation box. Here we have
one box, with a single component. Finally, ``mols_to_add`` describes the
number of molecules that we wish to add to each box prior to beginning
the simulation in Cassandra. We will add 50 methane molecules for this example.

.. code-block:: python

    box_list = [box]
    species_list = [methane_ff]
    mols_in_boxes = [[0]]
    mols_to_add = [[50]]


We now combine the four components created above into a single
``System`` object.

.. code-block:: python

    system = mc.System(box_list, species_list, mols_to_add=mols_to_add)

.. note::
    ``mols_in_boxes`` and ``mols_to_add`` are optional arguments when creating
    the ``System`` object. If they are not provided it is assumed that the
    values are zero for all species in all boxes.

.. note::
    Each item in the ``species_list`` must be a ``parmed.Structure`` object with
    the associated forcefield parameters. For example, ``species_list =
    [methane]`` would not work because the ``mbuild.Compound`` object does not
    contain any forcefield parameters.


Now we create a ``Moves`` object. This object contains all selections related to
the ``# Move_Probabilities`` section of the Cassandra input file. In addition
to the probability of performing different types of MC moves, the ``Moves``
object also contains the maximum move sizes (e.g., maximum translation distance),
whether each species is insertable, and so on. To create the moves object, we
must specify the ensemble we wish to perform our MC simulation in and the 
``species_list``.

.. code-block:: python

    ensemble = 'nvt'
    moves = mc.Moves(ensemble, species_list)

Many attributes of the moves object can be edited after it is created. This allows
complete control over all the move-related selections in Cassandra.

The only remaining step is to run the simulation. The ``mc.run`` function requires
five arguments: the ``System`` object, the ``Moves`` object, the temperature,
a selection of ``equilibration`` or ``production``, and the number of steps to run.

.. code-block:: python

    mc.run(system, moves, 300.0, 'equilibration', 10000)

A large number of additional keyword arguments can be provided inline or as part
of a keyword dictionary. See XX for details.



