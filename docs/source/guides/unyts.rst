Unyts
=====

`Unyt <https://unyt.readthedocs.io/en/stable/>`_
is a Python library for working with physical units.
In MoSDeF Cassandra, all quantities that have physical units
associated with them must be specified as a ``unyt_quantity``.
This approach yields several benefits. Users can specify
quantities in any (dimensionally valid) units they desire,
and thus do not need to dig through the reference manual
to determine the correct units for each quantity.
Possible errors in unit conversions are mitigated, and
we remove any possible ambiguity with regards to the
units of physical quantities in MoSDeF Cassandra scripts.

Basic usage
+++++++++++

Adding units to quantities is as easy as:

.. code-block:: python

  import unyt as u
  temperature = 300 * u.K

Compound units can be specified as:

.. code-block:: python
  
  import unyt as u
  energy = 100 * u.Unit('kJ/mol')

If a quantity or array is specified as a ``unyt_quantity``
or ``unyt_array``, then performing a unit conversion
is as simple as:

.. code-block:: python

  import unyt as u
  energy = 100 * u.Unit('kJ/mol')
  energy.in_units('kcal/mol')
 
The value (without units) can be extracted as:

.. code-block:: python

  energy.in_units('kcal/mol').value

Unyts in MoSDeF Cassandra
+++++++++++++++++++++++++

The base data structure of ``unyt`` is the ``unyt_array`` or ``unyt_quantity``
(a subclass of numpy ndarray) which carries both a value and a unit.
One of the main functionalities of **unyt** is the ability to convert units.
In **MoSDeF Cassandra**, a user can pass in a ``unyt_quantity`` of any valid
unit type which will get then get converted into the standard unit
specified by **Cassandra**. Unyt arrays are expected for values with
units, such as cutoffs, angles, volumes, pressures, and temperatures.
Unyt arrays are **not** expected for dimensionless values such as probabilities.  A
list of arguments and their required type can be viewed by running
``mosdef_cassandra.print_valid_kwargs``.


Important Cavaets
+++++++++++++++++

``mBuild`` does not use the ``unyt`` package. The distance units in
``mBuild`` are **nanometers**.



