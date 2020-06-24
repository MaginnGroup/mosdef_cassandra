import parmed
import mbuild
from mbuild.formats.cassandramcf import write_mcf

import mosdef_cassandra
from mosdef_cassandra.writers.inp_functions import generate_input


def write_mcfs(system, angle_style="harmonic"):

    if angle_style == "harmonic":
        angle_style = [angle_style] * len(system.species_topologies)

    if not isinstance(system, mosdef_cassandra.System):
        raise TypeError('"system" must be of type ' "mosdef_cassandra.System")

    for species_count, species in enumerate(system.species_topologies):
        if not isinstance(species, parmed.Structure):
            raise TypeError(
                'Your "system" object appears to have '
                "been corrupted. Species {} is not a parmed"
                ".Structure object".format(species)
            )
        if len(species.dihedrals) > 0 and len(species.rb_torsions) > 0:
            raise ValueError(
                "Your species has both CHARMM style "
                "dihedrals and Ryckaert-Bellemans style dihedrals. "
                "Only a single dihedral style per species is "
                "currently supported"
            )
        elif len(species.dihedrals) > 0:
            dihedral_style = "charmm"
        elif len(species.rb_torsions) > 0:
            dihedral_style = "opls"
        else:
            dihedral_style = "none"

        mcf_name = "species{}.mcf".format(species_count + 1)
        write_mcf(
            species,
            mcf_name,
            angle_style=angle_style[species_count],
            dihedral_style=dihedral_style,
        )


def write_configs(system):

    if not isinstance(system, mosdef_cassandra.System):
        raise TypeError('"system" must be of type ' "mosdef_cassandra.System")

    for box_count, box in enumerate(system.boxes):
        if not isinstance(box, mbuild.Compound) and not isinstance(
            box, mbuild.Box
        ):
            raise TypeError(
                'Your "system" object appears to have '
                "been corrupted. Box {} is not a mbuild"
                ".Compound object".format(species)
            )

        # Only save if box has particles inside
        # This only occurs if box is an mbuild.Compound
        if isinstance(box, mbuild.Compound):
            xyz_name = "box{}.in.xyz".format(box_count + 1)
            box.save(xyz_name, overwrite=True)


def write_input(system, moveset, run_type, run_length, temperature, **kwargs):

    if "run_name" not in kwargs:
        kwargs["run_name"] = moveset.ensemble

    if "restart" in kwargs and kwargs["restart"]:
        if "restart_name" not in kwargs:
            kwargs["restart_name"] = kwargs["run_name"]
        if kwargs["restart_name"] == kwargs["run_name"]:
            kwargs["run_name"] = kwargs["run_name"] + "-rst"

    inp_data = generate_input(
        system=system,
        moveset=moveset,
        run_type=run_type,
        run_length=run_length,
        temperature=temperature,
        **kwargs
    )

    inp_name = kwargs["run_name"] + ".inp"

    with open(inp_name, "w") as inp:
        inp.write(inp_data)

    return inp_name


def print_inputfile(
    system, moveset, run_type, run_length, temperature, **kwargs
):
    """Print an example Cassandra input file to screen

    This function allows one to look at the Cassandra input file that
    will be generated without running the MC simulation. The arguments
    are identical mosdef_cassandra.run

    Parameters
    ----------
    system : mosdef_cassandra.System
        the System to simulate
    moves : mosdef_cassandra.Moves
        the Move set to simulate
    run_type : "equilibration" or "production"
        the type of run; in "equilibration" mode, Cassandra adaptively changes
        the maximum translation, rotation, and volume move sizes to achieve
        an acceptance ratio of 0.5
    run_length : int
        length of the MC simulation
    temperature : float
        temperature at which to perform the MC simulation
    **kwargs : keyword arguments
        any other valid keyword arguments, see
        ``mosdef_cassandra.print_valid_kwargs()`` for details
    """

    if "run_name" not in kwargs:
        kwargs["run_name"] = moveset.ensemble

    if "restart" in kwargs and kwargs["restart"]:
        if "restart_name" not in kwargs:
            kwargs["restart_name"] = kwargs["run_name"]
        if kwargs["restart_name"] == kwargs["run_name"]:
            kwargs["run_name"] = kwargs["run_name"] + "-rst"

    inp_data = generate_input(
        system=system,
        moveset=moveset,
        run_type=run_type,
        run_length=run_length,
        temperature=temperature,
        **kwargs
    )

    print(inp_data)


def write_pdb(molecule, filename):

    # Generate CONECT records
    conect = {atom.idx: [] for atom in molecule.atoms}
    for atom in molecule.atoms:
        atidx = atom.idx
        for bond in molecule.bonds:
            if bond.atom1.idx == atidx:
                conect[atidx].append(bond.atom2.idx)
            elif bond.atom2.idx == atidx:
                conect[atidx].append(bond.atom1.idx)

    with open(filename, "w") as pdb:
        pdb.write("REMARK 1   Created by mosdef_cassandra\n")
        pdb.write(
            "CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f}"
            " {:9s}{:3d}\n".format(
                molecule.box[0],
                molecule.box[1],
                molecule.box[2],
                molecule.box[3],
                molecule.box[4],
                molecule.box[5],
                molecule.space_group,
                1,
            )
        )
        for atom in molecule.atoms:
            pdb.write(
                "ATOM  {:5d} {:4s} RES A{:4d}    "
                "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}"
                "          {:>2s}  \n".format(
                    atom.idx + 1,
                    atom.name,
                    0,
                    atom.xx,
                    atom.xy,
                    atom.xz,
                    1.0,
                    0.0,
                    atom.element_name,
                )
            )
        for atidx, atomlist in conect.items():
            pdb.write("CONECT{:5d}".format(atidx + 1))
            for at2idx in atomlist:
                pdb.write("{:5d}".format(at2idx + 1))
            pdb.write("\n")
