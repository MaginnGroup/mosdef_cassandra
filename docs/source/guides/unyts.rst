Unyts
=====

Unyt is a Python library for working with physical units.
The package documentation can be found
`here <https://unyt.readthedocs.io/en/stable/>`_. In MoSDeF
Cassandra, all quantities that have physical units
associated with them must be specified as a ``unyt_quantity``.
The benefit of this approach is that users can specify
quantities in any units they desire. Users are not forced
to reference the manual to determine the default units of
Cassandra. Using the ``unyt`` package also removes
any ambiguity with regards to the physical units of
any quantity in a MoSDeF Cassandra script.

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



Cavaets
+++++++





