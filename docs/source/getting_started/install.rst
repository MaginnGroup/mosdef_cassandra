Installation
============

We recommend the conda installation for most users.

Installing from conda
~~~~~~~~~~~~~~~~~~~~~

Assuming you already have
`conda <https://docs.conda.io/en/latest/miniconda.html>`_ installed,
you can create a new conda environment and install MoSDeF
Cassandra with a single command:

.. code-block:: bash

    conda create --name mc mosdef_cassandra foyer -c conda-forge -c mosdef -c omnia

The command creates a new conda environment (``mc``) and installs
``mosdef_cassandra`` and ``foyer``. The ``-c`` options specify the conda channels
that are searched, with priority decreasing from left to right. To use the
environment, run ``conda activate mc``.

You can test your installation by opening up a Python interpreter and typing:

.. code-block:: python

    import mosdef_cassandra as mc
    (py, fraglib_setup, cassandra) = mc.utils.detect_cassandra_binaries()

If the module is imported without error and is able to find the required
binaries (``python``, ``cassandra.exe``, and ``library_setup.py``,
you have successfully installed the package. Example output from the second
line is:

.. code-block:: text

    Using the following executables for Cassandra:
    Python: /Users/username/anaconda3/envs/mc-prod/bin/python
    library_setup: /Users/username/anaconda3/envs/mc-prod/bin/library_setup.py
    Cassandra: /Users/ryandefever/anaconda3/envs/mc-prod/bin/cassandra.exe

Installing from source
~~~~~~~~~~~~~~~~~~~~~~

MoSDeF Cassandra may alternatively be installed from source. First, clone
MoSDeF Cassandra from GitHub to a location of your choosing:

.. code-block:: bash

    git clone git@github.com:maginngroup/mosdef_cassandra.git

Next, install the required dependencies. You can use the dependencies listed
in ``requirements.txt`` or ``requirements-dev.txt``. However, if you are
installing from source we recommend the latter:

.. code-block:: bash

    conda install -c conda-forge -c mosdef -c omnia --file mosdef_cassandra/requirements-dev.txt

Finally, run the following commands to complete the installation of
MoSDeF Cassandra:

.. code-block:: bash

    cd mosdef_cassandra/
    pip install .


Installing Cassandra from source
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::
    Installing Cassandra from source is unnecessary unless you wish to modify
    the source code of Cassandra or use a hardware specific (e.g., intel) compiler.

Once you have downloaded the tarball (available
`here <https://github.com/MaginnGroup/Cassandra/releases>`_):

.. code-block:: bash

    tar -xzvf Cassandra-1.2.2.gz
    cd Cassandra-1.2.2/Src
    make -f Makefile.gfortran
    cd ../
    mkdir bin/
    mv Src/cassandra_gfortran.exe ./bin/.
    cp Scripts/Frag_Library_Setup/library_setup.py ./bin/.

.. note::
    You may also wish to use the openMP version. In that case use the
    ``Makefile.gfortran.openMP`` and move the relevant executable to
    ``bin/``. Depending on system size, Cassandra the openMP version
    may offer speedups for up to ~8 cores. The number of OMP threads
    can be controlled by setting the ``OMP_NUM_THREADS`` environment
    variable, e.g., ``export OMP_NUM_THREADS=8``.


Add ``Cassandra-1.2.2/bin`` to your ``PATH``:

.. code-block:: bash

    export PATH=path_to_install/Cassandra-1.2.2/bin:${PATH}

Unless you add the preceding line to your ``.bashrc`` you will need to
run it every time you open a new terminal window.

