
Examples
========

Below we provide a few simple examples of short Monte Carlo simulations with
MoSDeF Cassandra.

.. note::
  Many of these examples require the openbabel package to create molecules from a
  SMILES string. Though openbabel is not a required dependency of
  MoSDeF Casssandra, we strongly recommend that users install it
  (``conda install -c conda-forge openbabel``) if they
  would like that functionality.

NVT simulation of methane
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    import mbuild
    import foyer
    import mosdef_cassandra as mc
    import unyt as u

    # Use mBuild to create a methane molecule
    methane = mbuild.load("C", smiles=True)

    # Create an empty mbuild.Box
    box = mbuild.Box(lengths=[3.0, 3.0, 3.0])

    # Load force field
    oplsaa = foyer.forcefields.load_OPLSAA()

    # Use foyer to apply force field to methane
    methane_ff = oplsaa.apply(methane)

    # Create box and species list
    box_list = [box]
    species_list = [methane_ff]

    # Use Cassandra to insert some initial number of methane molecules
    mols_to_add = [[50]]

    # Define the System
    system = mc.System(box_list, species_list, mols_to_add=mols_to_add)
    # Define the MoveSet
    moveset = mc.MoveSet("nvt", species_list)

    # Run a simulation at 300 K for 10000 MC moves
    mc.run(
        system=system,
        moveset=moveset,
        run_type="equilibration",
        run_length=10000,
        temperature=300.0 * u.K,
    )


NPT simulation of methane
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    import mbuild
    import foyer
    import mosdef_cassandra as mc
    import unyt as u

    # Use mbuild to create molecules
    methane = mbuild.load("C", smiles=True)

    # Create an empty mbuild.Box
    box = mbuild.Box(lengths=[3.0, 3.0, 3.0])

    # Load force field
    oplsaa = foyer.forcefields.load_OPLSAA()

    # Use foyer to apply force field
    methane_ff = oplsaa.apply(methane)

    # Create box and species list
    box_list = [box]
    species_list = [methane_ff]

    # Use Cassandra to insert some initial number of species
    mols_to_add = [[5]]

    # Define the System
    system = mc.System(box_list, species_list, mols_to_add=mols_to_add)
    # Define the MoveSet
    moveset = mc.MoveSet("npt", species_list)

    # Here we must specify the pressure since we are performing a
    # NpT simulation. It can be provided in the custom_args dictionary
    # or as a keyword argument to the "run" function.
    custom_args = {
        "pressure": 1.0 * u.bar,
    }

    # Run a simulation with at 300 K with 10000 MC moves
    mc.run(
        system=system,
        moveset=moveset,
        run_type="equilibration",
        run_length=10000,
        temperature=300.0 * u.K,
        **custom_args,
    )

NVT simulation of methane and propane mixture
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    import mbuild
    import foyer
    import mosdef_cassandra as mc
    import unyt as u

    # Use mbuild to create methane and propane molecules
    methane = mbuild.load("C", smiles=True)
    propane = mbuild.load("CCC", smiles=True)

    # Create an empty mbuild.Box
    box = mbuild.Box(lengths=[3.0, 3.0, 3.0])

    # Load force field
    oplsaa = foyer.forcefields.load_OPLSAA()

    # Use foyer to apply the force field
    typed_methane = oplsaa.apply(methane)
    typed_propane = oplsaa.apply(propane)

    # Create box and species list
    box_list = [box]
    species_list = [typed_methane, typed_propane]

    # Use Cassandra to insert some initial number of species
    mols_to_add = [[100, 50]]

    system = mc.System(box_list, species_list, mols_to_add=mols_to_add)
    moveset = mc.MoveSet("nvt", species_list)

    mc.run(
        system=system,
        moveset=moveset,
        run_type="equilibration",
        run_length=10000,
        temperature=200.0 * u.K,
    )

GEMC simulation of methane (united atom)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    import mbuild
    import foyer
    import mosdef_cassandra as mc
    import unyt as u

    # Use mbuild to create a coarse-grained CH4 bead
    methane = mbuild.Compound(name="_CH4")

    # Create two empty mbuild.Box
    # (vapor = larger, liquid = smaller)
    liquid_box = mbuild.Box(lengths=[3.0, 3.0, 3.0])
    vapor_box = mbuild.Box(lengths=[4.0, 4.0, 4.0])

    # Load force field
    trappe = foyer.forcefields.load_TRAPPE_UA()

    # Use foyer to apply force field
    typed_methane = trappe.apply(methane)

    # Create box and species list
    box_list = [liquid_box, vapor_box]
    species_list = [typed_methane]

    mols_to_add = [[350], [100]]

    system = mc.System(box_list, species_list, mols_to_add=mols_to_add)
    moveset = mc.MoveSet("gemc", species_list)

    moveset.prob_volume = 0.010
    moveset.prob_swap = 0.11

    thermo_props = [
        "energy_total",
        "energy_intervdw",
        "pressure",
        "volume",
        "nmols",
        "mass_density",
    ]

    custom_args = {
        "run_name": "equil",
        "charge_style": "none",
        "rcut_min": 2.0 * u.angstrom,
        "vdw_cutoff": 14.0 * u.angstrom,
        "units": "sweeps",
        "steps_per_sweep": 450,
        "coord_freq": 50,
        "prop_freq": 10,
        "properties": thermo_props,
    }

    mc.run(
        system=system,
        moveset=moveset,
        run_type="equilibration",
        run_length=250,
        temperature=151.0 * u.K,
        **custom_args,
    )

    # Update run_name and restart_name
    custom_args["run_name"] = "prod"
    custom_args["restart_name"] = "equil"

    mc.restart(
        system=system,
        moveset=moveset,
        run_type="production",
        run_length=750,
        temperature=151.0 * u.K,
        **custom_args,
    )

GCMC simulation of methane
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    import mbuild
    import foyer
    import mosdef_cassandra as mc
    import unyt as u

    # Use mbuild to create a methane
    methane = mbuild.load("C", smiles=True)

    # Create an empty mbuild.Box
    box = mbuild.Box(lengths=[10.0, 10.0, 10.0])

    # Load force field
    oplsaa = foyer.forcefields.load_OPLSAA()

    # Use foyer to apply the force field
    methane_ff = oplsaa.apply(methane)

    # Create box and species list
    box_list = [box]
    species_list = [methane_ff]

    mols_to_add = [[100]]

    system = mc.System(box_list, species_list, mols_to_add=mols_to_add)
    moveset = mc.MoveSet("gcmc", species_list)

    custom_args = {
        "chemical_potentials": [-35.0 * (u.kJ / u.mol)],
        "prop_freq": 100,
    }

    mc.run(
        system=system,
        moveset=moveset,
        run_type="equilibration",
        run_length=1000,
        temperature=300.0 * u.K,
        **custom_args,
    )

GCMC simulation of methane adsorption in a solid framework
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    import mbuild
    import foyer
    import mosdef_cassandra as mc
    import unyt as u

    from mosdef_cassandra.examples.structures import carbon_lattice


    # Load a structure created with mbuild
    lattice = carbon_lattice()
    # Use mbuild to create a methane
    methane = mbuild.load("C", smiles=True)

    # Load force field
    trappe = foyer.forcefields.load_TRAPPE_UA()
    oplsaa = foyer.forcefields.load_OPLSAA()

    # Use foyer to apply the force fields
    typed_lattice = trappe.apply(lattice)
    methane_ff = oplsaa.apply(methane)

    # Create box and species list
    box_list = [lattice]
    species_list = [typed_lattice, methane_ff]

    # Since we have an occupied box we need to specify
    # the number of each species present in the initial config
    mols_in_boxes = [[1, 0]]

    system = mc.System(box_list, species_list, mols_in_boxes=mols_in_boxes)
    moveset = mc.MoveSet("gcmc", species_list)

    custom_args = {
        "chemical_potentials": ["none", -30.0 * (u.kJ / u.mol)],
        "rcut_min": 0.5 * u.angstrom,
        "vdw_cutoff": 14.0 * u.angstrom,
        "charge_cutoff": 14.0 * u.angstrom,
        "coord_freq": 100,
        "prop_freq": 10,
    }

    mc.run(
        system=system,
        moveset=moveset,
        run_type="equilibration",
        run_length=10000,
        temperature=300.0 * u.K,
        **custom_args,
    )


NVT simulation of SPC/E water
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    import mbuild
    import foyer
    import mosdef_cassandra as mc
    import unyt as u
    from mosdef_cassandra.utils.get_files import get_example_ff_path, get_example_mol2_path

    # Load water with SPC/E geometry from mol2 file
    molecule = mbuild.load(get_example_mol2_path("spce"))

    # Create an empty mbuild.Box
    box = mbuild.Box(lengths=[3.0, 3.0, 3.0])

    # Load force field
    spce = foyer.Forcefield(get_example_ff_path("spce"))

    # Use foyer to apply force field
    molecule_ff = spce.apply(molecule)

    # Create box and species list
    box_list = [box]
    species_list = [molecule_ff]

    # Use Cassandra to insert some initial number of species
    mols_to_add = [[50]]

    # Define the System
    system = mc.System(box_list, species_list, mols_to_add=mols_to_add)
    # Define the MoveSet
    moveset = mc.MoveSet("nvt", species_list)

    # Note here we need to use the angle_style="fixed" keyword argument
    # SPC/E geometry is rigid; default angle style is "harmonic"
    custom_args = {"angle_style": ["fixed"]}

    # Run a simulation with at 300 K with 10000 MC moveset
    mc.run(
        system=system,
        moveset=moveset,
        run_type="equilibration",
        run_length=10000,
        temperature=300.0 * u.K,
        **custom_args,
    )
