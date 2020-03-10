
MoSDeF Cassandra
================ 
|Citing|
|License|
|CodeCov|
|Azure|

.. |Citing| image:: https://img.shields.io/badge/cite-mosdef__cassandra-blue
   :target: https://mosdef-cassandra.readthedocs.io/en/latest/reference/citing.html
.. |Codecov| image:: https://codecov.io/gh/MaginnGroup/mosdef_cassandra/branch/master/graph/badge.svg
.. |Azure| image:: https://dev.azure.com/MaginnGroup/mosdef_cassandra/_apis/build/status/MaginnGroup.mosdef_cassandra?branchName=master
.. |License| image:: https://img.shields.io/github/license/maginngroup/mosdef_cassandra

Overview
~~~~~~~~

**MoSDeF Cassandra** is a Python wrapper for the
`Cassandra <https://cassandra.nd.edu>`_ Monte Carlo code.
The wrapper interfaces with the `MoSDeF <https://mosdef.org>`_
tools and provides a user-friendly interface to Cassandra
without sacrificing any capabilities of Cassandra.

Warning
~~~~~~~

**MoSDeF Cassandra** is still in early development (0.x releases). The API may
change unexpectedly.

Resources
~~~~~~~~~

* `Reference documentation <https://mosdef-cassandra.readthedocs.io>`_:
  Examples, tutorials, guides, and API documentation
* `Installation guide <https://mosdef-cassandra.readthedocs.io/en/latest/getting_started/install.html>`_:
  Instructions for installing MoSDeF Cassandra
* `GitHub repository <https://github.com/MaginnGroup/mosdef_cassandra>`_:
  View the source code, contribute, and raise issues
* `Cassandra <https://cassandra.nd.edu>`_:
  Home of the Cassandra Monte Carlo package
* `MoSDeF tools <https://mosdef.org>`_:
  A generic collection of tools for constructing systems and applying
  forcefield parameters for particle-based simulations

Citation
~~~~~~~~

Please cite **MoSDeF Cassandra**, **Cassandra**, and the **MoSDeF** suite of
tools if you use this tool in your research. Details 
`here <https://mosdef-cassandra.readthedocs.io/en/latest/reference/citing.html>`_.

Installation
~~~~~~~~~~~~

Installation instructions are
`here <https://mosdef-cassandra.readthedocs.io/en/latest/getting_started/install.html>`_.
A conda installation will be added in the near future.

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

Credits
~~~~~~~

Development of MoSDeF Cassandra was supported by the National Science Foundation
under grant NSF Grant Number 1835874. Any opinions, findings, and conclusions or
recommendations expressed in this material are those of the author(s) and do
not necessarily reflect the views of the National Science Foundation.


