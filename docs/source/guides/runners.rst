The run function
================

To run a Monte Carlo simulation, we use:

.. code-block:: python

  mc.run(
    system=system,
    moves=moves,
    run_type="equilibration",
    run_length=1000,
    temperature=300.0
  )

The function has five required arguments: an ``mc.System`` object,
an ``mc.Moves`` object, a choice of ``run_type``, the ``run_length``,
and the ``temperature``. Other optional arguments can be specified
individually or with a dictionary. For example, if we were performing
and NPT simulation and wished to specify the pressure, we could do the
following:

.. code-block:: python

  mc.run(
    system=system,
    moves=moves,
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
    moves=moves,
    run_type="equilibration",
    run_length=1000,
    temperature=300.0,
    **custom_args
  )

The dictionary-based approach is cleaner when specify a larger number of custom
options. For example:

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
    moves=moves,
    run_type="equilibration",
    run_length=1000,
    temperature=300.0,
    **custom_args
  )

The restart function
====================

MoSDeF Cassandra also supports restarting from a checkpoint file.
This is particularly helpful when switching from an equilibration
to production simulation. 

The procedure follows:

.. code-block:: python

  mc.run(
    system=system,
    moves=moves,
    run_type="equilibration",
    run_length=1000,
    temperature=300.0,
    run_name="equil"
  )

  mc.restart(
    system=system,
    moves=moves,
    run_type="production",
    run_length=1000,
    temperature=300.0,
    restart_name="equil",
    run_name="prod"
  )

Notice the usage of ``run_name`` in both commands and ``restart_name`` in the
call to ``mc.restart``. The output from the equilibration is named ``"equil"``.
Therefore, when we use ``mc.restart``, we specify that it should restart
from the output files named ``"equil"``.

.. note:: 
  In Cassandra, during an "equilibration",
  the move sizes are adjusted to achieve a 50% acceptance ratio. In
  a "production" run the move sizes are fixed.

.. warning::
  If using ``mc.restart()``, the move sizes are read from the
  checkpoint file and therefore the move sizes in the ``Moves``
  object are NOT used.




