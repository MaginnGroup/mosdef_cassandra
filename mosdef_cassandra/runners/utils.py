import mbuild
import parmed
import glob
import re


def check_system(system, moveset):
    """Run a series of sanity checks on the System and MoveSet"""

    if moveset.ensemble == "gemc" or moveset.ensemble == "gemc_npt":
        if len(system.boxes) != 2:
            raise ValueError(
                "{} requested but {} simulation "
                "boxes provided as part of system. {} requires "
                "2 simulation boxes".format(
                    moveset.ensemble, len(system.boxes), moveset.ensemble
                )
            )
    else:
        if len(system.boxes) != 1:
            raise ValueError(
                "{} requested but {} simulation "
                "boxes provided as part of system. {} requires "
                "1 simulation box".format(
                    moveset.ensemble, len(system.boxes), moveset.ensemble
                )
            )

    for box in system.boxes:
        if not isinstance(box, mbuild.Box) and not isinstance(
            box, mbuild.Compound
        ):
            raise TypeError(
                "Not all System.boxes are mbuild.Box "
                "or mbuild.Compound objects. It appears "
                "your System object has been corrupted"
            )

    # TODO: Add check that species_topologies provided to the
    # System and MoveSet are the same

    if not isinstance(system.species_topologies, list):
        raise TypeError(
            "System.species_topologies should be a "
            "list. It appears your System object "
            "has been corrupted"
        )

    for species in system.species_topologies:
        if not isinstance(species, parmed.Structure):
            raise TypeError(
                "Each species should be a parmed.Structure. "
                "It appears your System object has been "
                "corrupted"
            )

    try:
        system.check_natoms()
    except:
        raise ValueError(
            "The number of atoms in one or more boxes "
            "does not match the number expected from "
            "System.species_topologies and system.mols_in_boxes. "
            "It appears your System object has been corrupted"
        )


def get_restart_name(restart_from, run_name):
    """Get the run name for a restart"""
    if restart_from is None:
        # Search for inp files
        inp_files = glob.glob("*.inp")
        if len(inp_files) == 0:
            raise FileNotFoundError("No previous input files found!")
        else:
            match = re.search(r"(.*)(\.rst\.)(\d{3})\.inp$", inp_files[0])
            if match is None:
                match = re.search(r"(.*)\.inp$", inp_files[0])
            base_name = match.group(1)
            if all([base_name in inp_file for inp_file in inp_files]):
                match = re.search(
                    r"(.*)(\.rst\.)(\d{3})\.inp$", sorted(inp_files)[-1]
                )
                if match is None:
                    match = re.search(r"(.*)\.inp$", sorted(inp_files)[-1])
                    restart_from = match.group(1)
                else:
                    restart_from = (
                        match.group(1) + match.group(2) + match.group(3)
                    )
            else:
                raise ValueError(
                    f"Multiple inp files: {inp_files} in working directory. "
                    "Please specify the file you wish to restart from with "
                    "the `restart_from` argument."
                )

    if run_name is None or run_name == restart_from:
        if ".rst." in restart_from:
            match = re.search(r"(.*)(.rst.)(\d{3})$", restart_from)
            base_name = match.group(1)
            restart_iter = int(match.group(3))
            if restart_iter >= 999:
                raise ValueError(
                    "Maximum number of restart iterations exceeded."
                )
            else:
                run_name = base_name + f".rst.{restart_iter + 1:03d}"
        else:
            run_name = restart_from + ".rst.001"

    return restart_from, run_name
