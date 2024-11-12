import parmed
import mbuild
import gmso
from gmso.formats.mcf import write_mcf as gmso_write_mcf
from mbuild.formats.cassandramcf import write_mcf
from pathlib import Path
from warnings import warn

from mosdef_cassandra import System, MoveSet
from mosdef_cassandra.writers.inp_functions import generate_input


def write_mcfs(system, angle_style="harmonic"):
    """Write a MCF file for a given mosdef_cassandra.System
    Parameters
    ----------
    system : mosdef_cassandra.System
        System to simulate in Cassandra
    angle_style : str, default="harmonic"
        Angle style for the system, valid arguments: "harmonic", "fixed"
    """
    if type(angle_style) == str:
        angle_style = [angle_style] * len(system.species_topologies)

    for astyle in angle_style:
        if astyle not in ["harmonic", "fixed"]:
            raise ValueError(
                'Invalid "angle_style" {} given.'.format(angle_style)
            )

    if not isinstance(system, System):
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

        if all(
            isinstance(top, parmed.Structure) for top in system.original_tops
        ):
            write_mcf(
                species,
                mcf_name,
                angle_style=angle_style[species_count],
                dihedral_style=dihedral_style,
            )
        elif all(
            isinstance(top, gmso.Topology) for top in system.original_tops
        ):
            gmso_write_mcf(system.original_tops[species_count], mcf_name)


def write_configs(system):

    if not isinstance(system, System):
        raise TypeError('"system" must be of type ' "mosdef_cassandra.System")

    for box_count, box in enumerate(system.boxes):
        if not isinstance(box, mbuild.Compound) and not isinstance(
            box, mbuild.Box
        ):
            raise TypeError(
                'Your "system" object appears to have '
                "been corrupted. Box {} is not a mbuild"
                ".Compound or mbuild.Box object".format(box)
            )

        # Only save if box has particles inside
        # This only occurs if box is an mbuild.Compound
        if isinstance(box, mbuild.Compound):
            xyz_name = "box{}.in.xyz".format(box_count + 1)
            box.save(xyz_name, overwrite=True)


def write_input(system, moveset, run_type, run_length, temperature, **kwargs):

    if "run_name" not in kwargs:
        kwargs["run_name"] = moveset.ensemble

    inp_data = generate_input(
        system=system,
        moveset=moveset,
        run_type=run_type,
        run_length=run_length,
        temperature=temperature,
        **kwargs,
    )

    inp_name = kwargs["run_name"] + ".inp"

    with open(inp_name, "w") as inp:
        inp.write(inp_data)

    return inp_name


def write_restart_input(restart_from, run_name, run_type, run_length):
    """Write an input file for a restart run"""
    input_contents = _generate_restart_inp(
        restart_from, run_name, run_type, run_length
    )
    with open(run_name + ".inp", "w") as f:
        f.write(input_contents)


def _generate_restart_inp(restart_from, run_name, run_type, run_length):
    """Create the input file for a restart"""
    # Extract contents of old input file
    old_inpfile_name = restart_from + ".inp"
    if not Path(old_inpfile_name).is_file():
        raise FileNotFoundError(
            f"Input file {old_inpfile_name} does not exist."
        )

    inp_contents = []
    with open(old_inpfile_name) as f:
        for line in f:
            inp_contents.append(line.strip())

    # Edit sections run_name, run_type, run_length
    for idx, line in enumerate(inp_contents):
        if "# Run_Name" in line:
            inp_contents[idx + 1] = run_name + ".out"
        if "# Start_Type" in line:
            inp_contents[idx + 1] = "checkpoint " + restart_from + ".out.chk"
            i = idx + 2
            while i < len(inp_contents) and not inp_contents[
                i
            ].strip().startswith("#"):
                if not inp_contents[i].strip().startswith("!"):
                    inp_contents[i] = (
                        ""  # Replace non-comment lines with an empty string
                    )
                i += 1
        if run_type is not None:
            if "# Run_Type" in line:
                old_contents = inp_contents[idx + 1].split()
                inp_contents[idx + 1] = (
                    run_type + " " + " ".join(old_contents[1:])
                )
        if run_length is not None:
            if "# Simulation_Length_Info" in line:
                # Verify new run length is >= original
                if run_length < int(inp_contents[idx + 4].split()[1]):
                    raise ValueError(
                        "Total run length on restart cannot be less than "
                        "the original run length. Please see the mc.restart "
                        "documentation for more details."
                    )
                if run_length == int(inp_contents[idx + 4].split()[1]):
                    warn(
                        "Total run length on restart is equal to the "
                        "original run length. This will not extend your "
                        " simulation. Please see the mc.restart "
                        "documentation for more details."
                    )
                inp_contents[idx + 4] = "run " + str(run_length)

    new_inp_contents = ""
    for line in inp_contents:
        new_inp_contents += line + "\n"

    return new_inp_contents


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
    moveset : mosdef_cassandra.MoveSet
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

    inp_data = generate_input(
        system=system,
        moveset=moveset,
        run_type=run_type,
        run_length=run_length,
        temperature=temperature,
        **kwargs,
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
