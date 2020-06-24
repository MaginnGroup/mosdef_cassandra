import mbuild
import foyer
import numpy as np
import mosdef_cassandra as mc
import unyt as u
from scipy.spatial.transform import Rotation
from mosdef_cassandra.utils.get_ff import get_ff_path


def run_nvt(custom_args={}):

    # Use mbuild to create molecules
    molecule = mbuild.load("O", smiles=True)
    bond_length = 0.1 * u.nm
    angle = 109.5 * u.degree

    molecule = mbuild.Compound(name="SOL")
    pos = np.asarray([0.5, 0.0, 0.0])
    o = mbuild.Compound(name="O", pos=pos)
    molecule.add(o)

    pos = np.asarray([0.5 + bond_length.to_value(), 0.0, 0.0])
    h1 = mbuild.Compound(name="H", pos=pos)
    molecule.add(h1)
    molecule.add_bond([o, h1])

    pos = o.pos - bond_length.to_value() * rotate2d(o.pos, angle.to_value())
    h2 = mbuild.Compound(name="H", pos=pos)
    molecule.add(h2)
    molecule.add_bond([o, h2])

    # Create an empty mbuild.Box
    box = mbuild.Box(lengths=[3.0, 3.0, 3.0])

    # Load forcefields
    spce = foyer.Forcefield(get_ff_path("spce"))

    # Use foyer to apply forcefields
    molecule_ff = spce.apply(molecule)

    # Create box and species list
    box_list = [box]
    species_list = [molecule_ff]

    # Use Cassandra to insert some initial number of species
    mols_to_add = [[50]]

    # Define the system object
    system = mc.System(box_list, species_list, mols_to_add=mols_to_add)
    # Get the move probabilities
    moveset = mc.MoveSet("nvt", species_list)

    # Run a simulation with at 300 K with 10000 MC moveset
    mc.run(
        system=system,
        moveset=moveset,
        run_type="equilibration",
        run_length=10000,
        temperature=300.0 * u.K,
        **custom_args,
    )


def rotate2d(vec, angle):
    """Rotate a vector `vec` anticlockwise by angle `angle` degrees"""

    angle = np.radians(180.0 - angle)

    xx = np.cos(angle) * vec[0] - np.sin(angle) * vec[1]
    yy = np.sin(angle) * vec[0] + np.cos(angle) * vec[1]

    vec = np.asarray([xx, yy, 0.0])
    vec /= np.linalg.norm(vec)

    return vec


if __name__ == "__main__":
    run_nvt()
