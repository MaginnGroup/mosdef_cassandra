
MoSDeF Cassandra
================

|Citing|
|CodeCov|
|Azure|

.. |Citing| image:: https://img.shields.io/badge/cite-mosdef__cassandra-blue
   :target: https://githib.com/rsdefever/mosdef_cassandra
.. |Codecov| image:: https://codecov.io/gh/rsdefever/mosdef_cassandra/branch/master/graph/badge.svg?token=xsgvdWprHp
.. |Azure| image:: https://img.shields.io/azure-devops/build/rdefever/aa33040a-0348-429f-a1ed-c9a8c69222c1/1

Overview
~~~~~~~~

**MoSDeF Cassandra** is a Python wrapper for the **Cassandra** Monte Carlo code.
The wrapper interfaces with the **MoSDeF** tools and provides a user-friendly
interface to Cassandra without sacrificing any capabilities of Cassandra.

.. warning::
  **MoSDeF Cassandra** is still in early development (0.x releases). The API may
  change unexpectedly.

Resources
~~~~~~~~~

* Reference Documentation: Examples, tutorials, guides, and API documentation
* :doc:`Installation guide <install>`: Instructions for installing MoSDeF Cassandra
* `GitHub repository <https://github.com/rsdefever/mosdef_cassandra>`_: View the source code, contribute, and raise issues
* `Cassandra <https://cassandra.nd.edu>`_: Home of the Cassandra Monte Carlo package
* `MoSDeF tools <https://mosdef.org>`_: A generic collection of tools for constructing systems and applying forcefield parameters for particle-based simulations

Citation
~~~~~~~~

Please cite **MoSDeF Cassandra**, **Cassandra**, and the **MoSDeF** suite of
tools if you use this tool in your research. See :doc:`here <citing>` for details.

Installation
~~~~~~~~~~~~

Installation instructions are :doc:`here <install>`. A conda installation will
be added in the near future.

Examples
~~~~~~~~

**MoSDeF Cassandra** provides a Python interface to **Cassandra**. The workflow
consists of first constructing a system and move set. These two objects are
passed to the runner that performs the Monte Carlo simulation with
**Cassandra**. We use classes from the **MoSDeF** tools to structure some of the
simulation inputs. The example below demonstrates an NVT Monte Carlo simulation
of OPLS methane. No input files are required. Everything required to run the
Monte Carlo calculation is contained in the script below.

.. code-block:: python

  import mbuild
  import foyer
  import mosdef_cassandra as mc

  # Create a methane molecule from a SMILES string
  methane = mbuild.load("C", smiles=True)

  # Load the forcefield via foyer
  ff = foyer.forcefields.load_OPLSAA()

  # Apply the forcefield parameters to methane with foyer
  methane_ff = ff.apply(methane)

  # Define an empty simulation box
  box = mbuild.Box([3.0, 3.0, 3.0])

  # Define the boxes, species in the system, molecules in the box
  ensemble = 'nvt'
  box_list = [box]
  species_list = [methane_ff]
  molecules_to_add = [[100]]

  # Create the System object
  system = mc.System(box_list, species_list, mols_to_add=molecules_to_add)

  # Create the Moves object
  moves = mc.Moves(ensemble, species_list)

  # Run a Monte Carlo simulation!
  mc.run(
      system=system,
      moves=moves,
      run_type="equilibration",
      run_length=1000,
      temperature=300.0
  )



.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install
   getting_started
   concepts
   api



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
