# mosdef_cassandra
MoSDeF compatible wrapper for Cassandra Monte Carlo code

WARNING: This package is under development and should not be considered ready for production use.

### Example usage

First we import the required packages. `mbuild` for system construction, `foyer` for atomtyping
and forcefield application, and `mosdef_cassandra` for Monte Carlo simulation.

    import mbuild
    import foyer
    import mosdef_cassandra as mc

Next, create an all-atom methane molcule from a smiles string (`C`). Use `foyer` to
load and apply the OPLS all atom forcefield.

    methane = mbuild.load('C',smiles=True)
    methane_typed = foyer.forcefields.load_OPLSAA().apply(methane)

Create a few variables to describe the system you wish to simulate. Here we are
going to simulate in the `nvt` ensemble. We create an empty simulation box with
`mbuild`.

    ensemble = 'nvt'
    box = mbuild.Box([30.,30.,30.])

Next, we create a list of all the boxes (there is only one), and a list of all
the unique species (again, there is only the methane). Finally, we specify that
we want to add 100 methane molecules to box 1. The format of `species_to_add`
is `[[list_for_box_1], [list_for_box_2]]` where `list_for_box_1` and `list_for_box_2`
are `[number_species_1_to_add, number_species_2_to_add, ...]`.

    box_list = [box]
    species_list = [methane_typed]
    species_to_add = [[100]]

We now combine all of the above information to define our system:

    # Define system
    system = mc.System(box_list,species_list,species_to_add=species_to_add)

Next we get a default set of move probabilities:

    # Define move set
    moves = mc.Moves(ensemble,species_list)

The default move probabilities can be edited by editing the `moves` object. Finally
we run the simulation. The required arguments are `system`, `moves`, the temperature
(`300`), whether the simulation is an equilibration or production (`equil` or `prod`),
and the simulation length (by default, number of MC moves). Other arguments can be
added with keywords or by specifying a dictionary of keywords and arguments.

    # Run simulation
    mc.run(system, moves, 300., 'equil', 10000)

Assuming all runs smoothly, your current directory will now contain a number of
new files. The standard output and standard error from the simulation can be found
in `mosdef_cassandra.log`.

### Installation

To install this package, run the following commands (recommended to do from within a conda environment)

    git clone git@github.com:rsdefever/mosdef_cassandra.git
    cd mosdef_cassandra/
    pip install .

### Dependencies

`mosdef_cassandra` requires `mbuild` and `foyer`. These can be installed via conda:

    conda install -c conda-forge -c omnia -c mosdef mbuild
    conda install -c conda-forge -c omnia -c mosdef foyer

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

