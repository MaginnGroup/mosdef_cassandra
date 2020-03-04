Installation
============

At this stage, the package must be installed from source.  A conda
installation will be added in the near future.

First, make sure you have the required dependencies. These include:

* Cassandra >= 1.2
* mbuild
* foyer

mbuild and foyer can be installed via conda:

.. code-block:: bash

    conda install -c omnia -c mosdef -c conda-forge mbuild
    conda install -c omnia -c mosdef -c conda-forge foyer

Cassandra must be installed from source. Once you have downloaded the tarball
(available `here <https://cassandra.nd.edu/index.php/download>`_):

.. code-block:: bash

    tar -xzvf Cassandra_V1.2.tar.gz
    cd Cassandra_V1.2/Src
    make -f Makefile.gfortran
    cd ../
    mkdir bin/
    mv Src/cassandra_gfortran.exe ./bin/.
    cp Script/Frag_Library_Setup/library_setup.py ./bin/.

.. warning::
  gfortran >= 9.0 causes problems with compilation. Please use gfortran < 9.0
  to compile Cassandra.

And finally add ``Cassandra_V1.2/bin`` to your ``PATH``:

.. code-block:: bash

    export PATH=path_to_install/Cassandra_V1.2/bin:${PATH}

Unless you add the preceding line to your ``.bashrc`` you will need to
run it every time you open a new terminal window.

.. warning::
    Though a conda installation of Cassandra is available,
    it is not currently compatibile with MoSDeF Cassandra.

Run the following commands to install mosdef_cassandra:

.. code-block:: bash

    git clone git@github.com:rsdefever/mosdef_cassandra.git
    cd mosdef_cassandra/
    pip install .

You can test your installation by opening up a Python interpreter and typing:

.. code-block:: python

    import mosdef_cassandra as mc
    (python2, fraglib_setup, cassandra) = mc.utils.detect_cassandra_binaries()

If the module is imported without error and is able to find the required
binaries (``python2``, ``cassandra_gfortran.exe``, and ``library_setup.py``,
you have successfully installed the package. Example output from the second
line is:

.. code-block:: text

    Using the following executables for Cassandra:
    Python: /usr/bin/python2
    library_setup:
    /afs/crc.nd.edu/user/r/rdefever/software/Cassandra_V1.2/bin/library_setup.py
    Cassandra:
    /afs/crc.nd.edu/user/r/rdefever/software/Cassandra_V1.2/bin/cassandra_gfortran_openMP.exe

.. warning::
  openbabel is not a required dependency for mbuild or MoSDeF Cassandra.
  However, most of the examples require openbabel to generate
  an ``mbuild.Compound`` from a SMILES string
  (e.g., ``molecule = mbuild.load("CC", smiles=True)``). Therefore, we highly
  recommend installing openbabel via conda:
  ``conda install -c conda-forge openbabel``.
