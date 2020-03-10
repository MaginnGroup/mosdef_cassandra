Installation
============

At this stage, the package must be installed from source.  A conda
installation will be added in the near future.

First, clone MoSDeF Cassandra from GitHub to a location of your choosing:

.. code-block:: bash

    git clone git@github.com:maginngroup/mosdef_cassandra.git

Next, install the required dependencies. These include:

* Cassandra >= 1.2
* mbuild
* foyer
* numpy

mbuild, foyer, and numpy can be installed via conda:

.. code-block:: bash

    conda install -c conda-forge -c mosdef -c omnia --file mosdef_cassandra/requirements.txt
    conda install -c conda-forge openbabel

The second line installs openbabel. It is not a strict requirement,
but is **highly recommended** and necessary for many of the examples.

Cassandra must be installed from source. Once you have downloaded the tarball
(available `here <https://cassandra.nd.edu/index.php/download>`_):

.. code-block:: bash

    tar -xzvf Cassandra_V1.2.tar.gz
    cd Cassandra_V1.2/Src
    make -f Makefile.gfortran
    cd ../
    mkdir bin/
    mv Src/cassandra_gfortran.exe ./bin/.
    cp Scripts/Frag_Library_Setup/library_setup.py ./bin/.

.. warning::
  gfortran >= 9.0 causes problems with compilation. Please use gfortran < 9.0
  to compile Cassandra.

.. note::
    You may also wish to use the openMP version. In that case use the
    ``Makefile.gfortran.openMP`` and move the relevant executable to
    ``bin/``. Depending on system size, Cassandra the openMP version
    may offer speedups for up to ~8 cores. The number of OMP threads
    can be controlled by setting the ``OMP_NUM_THREADS`` environment
    variable, e.g., ``export OMP_NUM_THREADS=8``.


Add ``Cassandra_V1.2/bin`` to your ``PATH``:

.. code-block:: bash

    export PATH=path_to_install/Cassandra_V1.2/bin:${PATH}

Unless you add the preceding line to your ``.bashrc`` you will need to
run it every time you open a new terminal window.

.. warning::
    Though a conda installation of Cassandra is available,
    it is not currently compatibile with MoSDeF Cassandra.

Finally, run the following commands to complete the installation of
mosdef_cassandra:

.. code-block:: bash

    cd PATH_TO_REPO/mosdef_cassandra/
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
    /user/username/software/Cassandra_V1.2/bin/library_setup.py
    Cassandra:
    /user/username/software/Cassandra_V1.2/bin/cassandra_gfortran_openMP.exe

