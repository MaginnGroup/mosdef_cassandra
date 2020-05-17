Run Monte Carlo
===============

To run a Monte Carlo simulation, use:

.. code-block:: python

  mc.run(
    system=system,
    moveset=moveset,
    run_type="equilibration",
    run_length=1000,
    temperature=300.0
  )

The ``run`` function has five required arguments: a ``System``,
``MoveSet``, a choice of ``run_type``, the ``run_length``,
and the ``temperature``. Other optional arguments can be specified
individually or with a dictionary. For example, if we were performing
and NPT simulation and needed to specify the pressure, we could do the
following:

.. code-block:: python

  mc.run(
    system=system,
    moveset=moveset,
    run_type="equilibration",
    run_length=1000,
    temperature=300.0,
    pressure=1.0
  )

or, if we wished to use a dictionary:

.. code-block:: python

  custom_args = {
    'pressure' : 1.0
  }

  mc.run(
    system=system,
    moveset=moveset,
    run_type="equilibration",
    run_length=1000,
    temperature=300.0,
    **custom_args
  )

The dictionary-based approach is easier to read when
specifying a larger number of custom options. For example:

.. code-block:: python

  custom_args = {
    'pressure' : 1.0,
    'cutoff_style' : 'cut_shift',
    'vdw_cutoff' : 14.0,
    'units' : 'sweeps',
    'prop_freq' : 10,
    'coord_freq' : 100
  }

  mc.run(
    system=system,
    moveset=moveset,
    run_type="equilibration",
    run_length=1000,
    temperature=300.0,
    **custom_args
  )

Restarting a simulation
=======================

MoSDeF Cassandra also supports restarting from a checkpoint file.
This is particularly helpful when switching from an equilibration
to production simulation.

The procedure follows:

.. code-block:: python

  mc.run(
    system=system,
    moveset=moveset,
    run_type="equilibration",
    run_length=1000,
    temperature=300.0,
    run_name="equil"
  )

  mc.restart(
    system=system,
    moveset=moveset,
    run_type="production",
    run_length=1000,
    temperature=300.0,
    restart_name="equil",
    run_name="prod"
  )

Notice the usage of ``run_name`` in both commands and ``restart_name`` in the
call to ``restart``. The output from the equilibration is named ``"equil"``.
Therefore, when we use ``restart``, we specify that it should restart
from the output files named ``"equil"``, while the new ``run_name`` is
``"prod"``.

.. note::
  If the ``run_type`` is ``"equilibration"``, Cassandra adjusts the
  maximum translation, rotation, and volume move sizes to achieve a
  50% acceptance ratio. If the ``run_type`` is ``"production"``, the
  maximum move sizes are fixed to the specified values.

.. warning::
  When using ``restart``, the maximum translation, rotation and volume
  move sizes are read from the checkpoint file and the values in the
  ``MoveSet`` are ignored.
