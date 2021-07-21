Quickstart Guide
================

The following assumes you have MoSDeF Cassandra installed. If not, please
refer to our :doc:`installation guide <install>`. More details of the core
MoSDeF Cassandra functionality can be found under the Guides section of
the documentation.

Let's start by setting up an NVT Monte Carlo simulation of OPLS-AA
methane. We will use the ``mbuild`` and ``foyer`` packages to create
the methane molecule, and ``mosdef_cassandra`` to run the Monte Carlo
simulation. We begin with the required imports:

.. code-block:: python

    import mbuild
    import foyer
    import unyt as u
    import mosdef_cassandra as mc

Next, we create an all-atom methane molecule from a `SMILES
<https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html>`_ string:

.. code-block:: python

    methane = mbuild.load("C", smiles=True)

``methane`` is a single all-atom methane molecule. It is an
``mbuild.Compound``. ``methane`` contains particles for each
element (C, H, H, H) in the molecule, coordinates associated
with each particle, and the bonds that describe the particle
connectivity. However, there are no forcefield parameters
associated with ``methane``.

To add forcefield parameters to ``methane``, we first load the OPLS-AA
forcefield from foyer. The OPLS-AA forcefield is distributed with foyer.
Be aware that not all atomtypes are currently defined.

.. code-block:: python

   oplsaa = foyer.forcefields.load_OPLSAA()

We then apply the forcefield using foyer:

.. code-block:: python

    methane_ff = oplsaa.apply(methane)

``methane_ff`` is a ``parmed.Structure`` that contains all the
forcefield parameters for our methane molecule.

Now that we have a molecule with forcefield parameters, the next step is
to define our simulation box. Since Cassandra can add molecules to a
simulation box before the start of a simulation, we can begin with an
empty simulation box. We will define an ``mbuild.Box`` with the box
lengths specified in nanometers:

.. code-block:: python

    box = mbuild.Box([3.0, 3.0, 3.0])

.. warning::
    Even though most quantities in MoSDeF Cassandra must be
    :doc:`specified with the unyt package <../guides/unyts>`,
    the ``mbuild.Box`` object is specified
    in nanometers without using ``unyt``. This is because
    ``mbuild`` does not currently support ``unyt``.


Next, we create the ``System`` object. It has two required arguments and
two optional arguments, depending on your system. The ``box_list`` and
``species_list`` are always specified. The ``box_list`` is simply a list
of the simulation boxes in the system. In this case, since we are performing
an NVT simulation there is only our single ``box``. The ``species_list`` is a
list of the unique chemical species in our system. Here we only have methane.

The two system-dependent arguments are ``mols_in_boxes`` and ``mols_to_add``.
Here we have an empty initial box, so we don't need to specify
``mols_in_boxes``. Finally, ``mols_to_add`` specifies the
number of molecules that we wish to add to each box prior to beginning
the simulation in Cassandra. We will add 50 methane molecules for this example.

.. code-block:: python

    box_list = [box]
    species_list = [methane_ff]
    mols_to_add = [[50]]

.. note::
    ``mols_in_boxes`` and ``mols_to_add`` are lists with one entry
    for each box. Each entry is itself a list, with one entry for
    each species in the ``species_list``.

We now combine the four components created above into a ``System``:

.. code-block:: python

    system = mc.System(box_list, species_list, mols_to_add=mols_to_add)

.. note::
    ``mols_in_boxes`` and ``mols_to_add`` are optional arguments when creating
    the ``System`` object. If not provided, the values are taken as zero for
    all species in all boxes.

.. note::
    Each item in the ``species_list`` must be a ``parmed.Structure`` object with
    the associated forcefield parameters. For example, ``species_list =
    [methane]`` would not work because unlike ``methane_ff``, ``methane`` is 
    a ``mbuild.Compound`` and does not contain forcefield parameters.

Now we create our ``MoveSet``. The ``MoveSet`` contains all selections
related to the MC moves that will be performed during the simulation.
In addition to the probability of performing different types of MC moves,
the ``MoveSet`` contains the maximum move sizes (e.g., maximum translation distance),
whether each species is insertable, and more. To create the ``MoveSet``, we
specify the ensemble in which we wish to perform the MC simulation and provide
the ``species_list``.

.. code-block:: python

    ensemble = 'nvt'
    moveset = mc.MoveSet(ensemble, species_list)

Some attributes of the ``MoveSet`` can be edited after it is created. This
allows complete control over all the move-related selections in Cassandra. To
view the current selections, use ``moveset.print()``.

The final step is to run the simulation. The ``run`` function requires
five arguments: the ``System``, ``MoveSet`` object, a selection of
``"equilibration"`` or ``"production"`` (``run_type``), the simulation length
(``run_length``), and the desired temperature. Note that since the temperature
is a physical quantity it must be specified with
:doc:`units attached <../guides/unyts>`.

.. code-block:: python

    mc.run(
        system=system,
        moveset=moveset,
        run_type="equilibration",
        run_length=10000,
        temperature=300.0 * u.K
    )

A large number of additional keyword arguments can be provided inline or as part
of a keyword dictionary. See ``mc.print_valid_kwargs()`` for a complete list of
the available keyword arguments.
