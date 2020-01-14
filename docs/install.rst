Installation
============

The package must be installed from source.  A conda installation will be
added in the future.

First, make sure you have the required dependencies. These include:

* Cassandra V1.2
* mbuild
* foyer

mbuild and foyer can be installed via conda::

    conda install -c conda-forge -c omnia -c mosdef mbuild
    conda install -c conda-forge -c omnia -c mosdef foyer

Cassandra must be installed from source. Once you have downloaded the tarball
(available `here <https://cassandra.nd.edu/index.php/download>`_)::

    tar -xzvf Cassandra_V1.2.tar.gz
    cd Cassandra_V1.2/Src
    make -f Makefile.gfortran
    cd ../
    mkdir bin/
    mv Src/cassandra_gfortran.exe ./bin/.
    cp Script/Frag_Library_Setup/library_setup.py ./bin/.

Run the following commands to install mosdef_cassandra::

    git clone git@github.com:rsdefever/mosdef_cassandra.git
    cd mosdef_cassandra/
    pip install .

You can test your installation by opening up a Python interpreter and typing::

    import mosdef_cassandra as mc
    mc.detect_cassandra_binaries()

If the module is imported without error and is able to find the required
binaries (``python2``, ``cassandra_gfortran.exe``, and ``library_setup.py``,
you have successfully installed the package.


