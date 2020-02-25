import mbuild


def carbon_lattice():
    carbon = mbuild.Compound(name="_CH4")
    angles = [90.0, 90.0, 90.0]
    carbon_locations = [
        [0.0, 0.0, 0.0],
        [1.0 / 2.0, 1.0 / 2.0, 0.0],
        [0.0, 1.0 / 2.0, 1.0 / 2.0],
        [1.0 / 2.0, 0.0, 0.0],
        [0.0, 0.0, 1.0 / 2.0],
    ]
    basis = {"C": carbon_locations}
    carbon_dict = {"C": carbon}
    lattice = mbuild.Lattice(
        lattice_spacing=[0.746, 0.746, 0.746],
        angles=angles,
        lattice_points=basis,
    )

    system = lattice.populate(compound_dict=carbon_dict, x=4, y=4, z=4)

    return system
