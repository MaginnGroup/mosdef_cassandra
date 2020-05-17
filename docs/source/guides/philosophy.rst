
Philosophy
==========

Performing a Monte Carlo simulation should be straightforward
and intuitive. The simulation setup procedure should not be prone to error.
The process should be easily repeatible and extensible by others. And none of
the prior goals can sacrifice the complete flexibility required by the expert
simulator.

.. _keyconcepts:

Key concepts
============

MoSDeF Cassandra integrates the `Cassandra Monte Carlo (MC)
code <https://cassandra.nd.edu>`_ with the `Molecular Simulation Design
Framework (MoSDeF) <https://mosdef.org>`_. This integration enables
users to build systems, apply force field parameters, and setup
and run MC simulations from within a single Python script. In
principle, MoSDeF Cassandra should make it easier to create `TRUE
simulations <https://www.tandfonline.com/doi/full/10.1080/00268976.2020.1742938>`_
-- simulations that are transparent, reproducible,
usable by others, and extensible. In addition to improving
simulation workflow reproducibility, MoSDeF Cassandra also provides
a user-friendly interface to Cassandra to improve the Cassandra
user experience for beginning and expert simulators alike.

The Cassandra
`documentation <https://cassandra.nd.edu/index.php/documentation>`_
contains a fairly comprehensive detail on the theory of MC simulations,
with a particular focus on the algorithms used in Cassandra. Canonical
textbooks by `Frenkel and Smit
<https://www.sciencedirect.com/book/9780122673511/understanding-molecular-simulation>`_,
`Allen and Tidesley
<https://www.oxfordscholarship.com/view/10.1093/oso/9780198803195.001.0001/oso-9780198803195>`_,
and `Tuckerman <https://onlinelibrary.wiley.com/doi/10.1002/anie.201105752>`_ also provide
useful reference materials for background on MC simulations. Here
however, we focus on the MoSDeF Cassandra package, proceeding with
the assumption that the reader has a basic understanding of MC methods as
applied in molecular simulations.

Organization
~~~~~~~~~~~~

.. The development of the MoSDeF Cassandra interface is motivated by a pair of
   simple questions: *What are the components of a Monte Carlo simulation?* and
   *What is the simplest logical organization of these components?*

The MoSDeF Cassandra interface is motivated by a simple
question: **What is the simplest logical organization of the components of a
Monte Carlo simulation?**

We organize the components of a Monte Carlo simulation into two primary groups: the
*system* and the *move set*. The *system* is what you are simulating; the simulation
boxes, any initial structures in the simulation boxes, and the forcefield parameters
used to describe the physics of the simulation. The *move set* are what you do during
the simulation; the types of MC moves that are attempted, the probabilities of
each, and any other parameters associated with the attempted moves.

Implementation
~~~~~~~~~~~~~~

MoSDeF Cassandra implements this organization by mapping the process of
setting up and running an MC simulation into three discrete steps:

1. Create the ``System``. This defines the simulation
   boxes, any initial configuration(s), and the forcefield parameters for all
   molecules in the simulation. Further detail of the :doc:`System <system>`.
2. Create the ``MoveSet``. This defines the all the selections related to
   different possible Monte Carlo moves. Further detail of the
   :doc:`MoveSet <moveset>`.
3. Pass the ``System`` and ``MoveSet`` to the ``run()`` function, along
   with a few remaining selections required for MC. Further detail of the
   :doc:`run function <runners>`.



