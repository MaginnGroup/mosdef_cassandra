# mosdef_cassandra
MoSDeF compatible wrapper for Cassandra Monte Carlo code

WARNING: This package is under development and should not be considered ready for production use.

### Installation

To install this package, run the following commands (recommended to do from within a conda environment)

    git clone git@github.com:rsdefever/mosdef_cassandra.git
    cd mosdef_cassandra/
    pip install .

### Dependencies

`mosdef_cassandra` requires `mbuild` and `foyer`. These can be installed via conda:

    conda install mbuild
    conda install foyer

`mosdef_cassandra` also requires the Cassandra Monte Carlo code.
It can be found [here](https://cassandra.nd.edu/). Once you have installed
Cassandra (see installation instructions at prior link), you will need to
make sure that the Cassandra executable (e.g., `cassandra_gfortran.exe`)
and `library_setup.py` are somewhere in your `PATH`. As an example of
compiling Cassandra and then adding the appropriate location to your `PATH`:

    tar -xzvf Cassandra_V1.2.tar.gz
    cd Cassandra_V1.2/Src
    make -f Makefile.gfortran
    cd ../
    mkdir bin/
    mv Src/cassandra_gfortran.exe ./bin/.
    cp Script/Frag_Library_Setup/library_setup.py ./bin/.

    # Execute the following everytime you wish to use Cassandra or add to your .bashrc
    export PATH=/home/USER/software/Cassandra_V1.2/bin:${PATH}

You will also need a `python2` installation and on your `PATH`.

### Future plans
* Addition of fragment library generation to `mosdef_cassandra` so that a call
to a separate python script is not required.


### List of Cassandra capabilities that are not (yet) supported
* Mie potential (support coming with migration to topology.Topology object)
* Fixed angles (support coming with migration to topology.Topology object)

