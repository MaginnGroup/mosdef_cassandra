
Examples
========

Below we provide a few simple examples of short Monte Carlo simulations with
MoSDeF Cassandra.

NVT simulation of methane
~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../../mosdef_cassandra/examples/nvt.py
  :language: python

NPT simulation of methane
~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../../mosdef_cassandra/examples/npt.py
  :language: python
  
NVT simulation of methane and propane mixture
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../../mosdef_cassandra/examples/nvt_mixture.py
  :language: python

GEMC simulation of methane
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../../mosdef_cassandra/examples/gemc.py
  :language: python

GCMC simulation of methane
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../../mosdef_cassandra/examples/gcmc.py
  :language: python

GCMC simulation of methane adsorption in a solid framework
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../../mosdef_cassandra/examples/gcmc_adsorption.py
  :language: python

.. note::
  The above examples require the openbabel package to create molecules from a
  SMILES string. Though openbabel is not a required dependency of
  MoSDeF Casssandra, we strongly recommend that users install it
  (``conda install -c conda-forge openbabel``) if they
  would like that functionality.

