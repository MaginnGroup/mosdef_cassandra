from copy import deepcopy
from constrainmol import ConstrainedMolecule
from gmso.external.convert_parmed import to_parmed
import gmso
import numpy as np
import mbuild
import parmed


class System(object):
    def __init__(
        self,
        boxes,
        species_topologies,
        mols_in_boxes=None,
        mols_to_add=None,
        fix_bonds=True,
    ):
        """A class to contain the system to simulate in Cassandra

        A System comprises the initial simulation box(es) (empty or
        occupied), the topologies of each species to be simulated,
        and the number of each species to be added to the simulation
        box(es) prior to the start of the simulation. These three
        items are represented by ``boxes``, ``species_topologies``,
        and ``mols_to_add``. If providing a box with existing
        species, you are required to specify ``mols_in_boxes``,
        the number of each species that already exists.

        Each argument is specified as a list, with either one
        element for each box or one element for each species.
        Arguments must be provided as a list even in the case
        of a single species or single box.

        Parameters
        ----------
        boxes : list
            one element per box. Each element should be a
            mbuild.Compound or mbuild.Box
        species_topologies : list
            list of parmed.Structures or gmso.Topology, with one species per element
        mols_in_boxes: list, optional
            one element per box. Each element is a list of length
            n_species, specifying the number of each species that
            are currently in each box
        mols_to_add : list, optional
            one element per box. Each element is a list of length
            n_species, specifying the number of each species that
            should be added to each box
        fix_bonds : boolean, optional, default=True
            update the bond lengths in any initial structure
            (i.e., boxes) to match the values specified in
            the species_topologies

        Returns
        -------
        mosdef_cassandra.System
        """
        self._boxes = None
        self._species_topologies = None
        self._mols_in_boxes = None
        self._mols_to_add = None
        self.original_tops = None

        # @setter decorators used to protect boxes, species
        # topologies, and mols_in_boxes from modification.
        # Error handling also occurs there
        self.boxes = boxes
        self.species_topologies = species_topologies
        self.mols_in_boxes = mols_in_boxes
        self.mols_to_add = mols_to_add

        # Error checking on number of particles in box(es)
        # vs. number from self.mols_in_boxes and
        # self.species_topologies
        self.check_natoms()

        # Fix the coordinates if the user provides a starting structure
        if fix_bonds:
            self.fix_bonds()

    # TODO: one possibility is to return list(self._boxes)
    # rather than self._boxes --> this prevents list items from
    # being edited ¯\_(ツ)_/¯
    @property
    def boxes(self):
        return self._boxes

    @boxes.setter
    def boxes(self, boxes):
        if self._boxes is None:
            self._boxes = []
            if not isinstance(boxes, list):
                raise TypeError(
                    '"boxes" should be a list. See '
                    "help(mosdef_Cassandra.System) for details."
                )
            for box in boxes:
                if isinstance(box, mbuild.Compound):
                    self._boxes.append(mbuild.clone(box))
                elif isinstance(box, mbuild.Box):
                    self._boxes.append(deepcopy(box))
                else:
                    raise TypeError(
                        "Each box should be an "
                        "mbuild.Compound or mbuild.Box object"
                    )
        else:
            raise AttributeError(
                "Box(es) cannot be modified after "
                "System object is created. Create a new System "
                "object if you wish to change the box(es)"
            )

    @property
    def original_tops(self):
        return self._original_topology

    @original_tops.setter
    def original_tops(self, original_topology):
        self._original_topology = original_topology

    @property
    def species_topologies(self):
        return self._species_topologies

    @species_topologies.setter
    def species_topologies(self, species_topologies):
        self._constrained_species = []
        if self._species_topologies is None:
            if not isinstance(species_topologies, list):
                raise TypeError(
                    '"species_topologies" should be a list. '
                    "See help(mosdef_Cassandra.System) for details."
                )

            if not (
                all(
                    isinstance(top, gmso.Topology)
                    for top in species_topologies
                )
                or all(
                    isinstance(top, parmed.Structure)
                    for top in species_topologies
                )
            ):

                raise TypeError(
                    "Each species should be a "
                    "parmed.Structure or gmso.Topology"
                    "and must be of the same type"
                )

            self._species_topologies = []
            self.original_tops = species_topologies

            for top in species_topologies:
                subtops = []
                if isinstance(top, gmso.Topology):
                    for molecule in top.unique_site_labels(name_only=True):
                        subtops.append(
                            top.create_subtop("molecule", (molecule, 1))
                        )

                    if len(subtops) > 1:
                        raise ValueError(
                            "GMSO Topology must contain only one molecule type. For example, "
                            "if you have a box of water, you must have a single water molecule "
                            "type in the topology. If you have a box of water and methane, you "
                            "must have two separate single molecule topologies, one for water "
                            " and one for methane."
                        )
                    subtops[0].box = top.box
                    top = to_parmed(subtops[0])
                    # top = to_parmed(top)

                # If no bonds in topology don't try to apply constraints
                # Store "None" in _constrained_species instead
                if len(top.bonds) > 0:
                    constrain = ConstrainedMolecule(top)
                    constrain.solve()
                    top.coordinates = constrain.xyz
                    self._constrained_species.append(constrain)
                else:
                    self._constrained_species.append(None)

                self._species_topologies.append(parmed.structure.copy(top))

        else:
            raise AttributeError(
                "species_topologies cannot be "
                "modified after System object is created. "
                "Create a new System object if you wish to "
                "edit the species_topolgies"
            )

    @property
    def mols_in_boxes(self):
        return list(self._mols_in_boxes)

    @mols_in_boxes.setter
    def mols_in_boxes(self, mols_in_boxes):
        if self._mols_in_boxes is None:
            # Make sure mols_in_boxes and mols_to_add are
            # properly formatted even if use specified none
            n_species = len(self.species_topologies)
            n_boxes = len(self.boxes)
            if mols_in_boxes is None:
                mols_in_boxes = [[0] * n_species] * n_boxes
            # Error checking first
            if not isinstance(mols_in_boxes, list):
                raise TypeError(
                    '"mols_in_boxes" should be a list. '
                    "See help(mosdef_Cassandra.System) for details."
                )
            if len(mols_in_boxes) != n_boxes:
                raise ValueError(
                    "The number of boxes inferred from the "
                    'length of "mols_in_boxes" must match the '
                    "number of boxes inferred from the length of "
                    '"boxes"'
                )
            for species_in_box in mols_in_boxes:
                if not isinstance(species_in_box, list):
                    raise TypeError(
                        '"mols_in_boxes" should be a list '
                        "with one list for each box."
                        "See help(mosdef_Cassandra.System) for details."
                    )
                if len(species_in_box) != n_species:
                    raise ValueError(
                        "The number of each species existing "
                        "in each box must be specified (even if = 0)"
                    )
                for current_n in species_in_box:
                    if not isinstance(current_n, int):
                        raise TypeError(
                            "The number of each species existing "
                            "in each box must be specified "
                            "as an integer"
                        )
            # Save data
            self._mols_in_boxes = deepcopy(mols_in_boxes)
        else:
            raise AttributeError(
                "mols_in_boxes cannot be "
                "modified after System object is created. "
                "Create a new System object if you wish to "
                "edit mols_in_boxes"
            )

    @property
    def mols_to_add(self):
        return self._mols_to_add

    @mols_to_add.setter
    def mols_to_add(self, mols_to_add):
        n_species = len(self.species_topologies)
        n_boxes = len(self.boxes)
        if mols_to_add is None:
            mols_to_add = [[0] * n_species] * n_boxes
        if not isinstance(mols_to_add, list):
            raise TypeError(
                '"mols_to_add" should be a list. '
                "See help(mosdef_Cassandra.System) for details."
            )
        if len(mols_to_add) != n_boxes:
            raise ValueError(
                "The number of boxes inferred from the "
                'length of "mols_to_add" must match the '
                "number of boxes inferred from the length of "
                '"boxes"'
            )
        for add_to_box in mols_to_add:
            if not isinstance(add_to_box, list):
                raise TypeError(
                    '"mols_to_add" should be a list '
                    "with one list for each box. "
                    "See help(mosdef_Cassandra.System) for details."
                )
            if len(add_to_box) != n_species:
                raise ValueError(
                    "The number of each species to "
                    "be added to each box must be specified "
                    "(even if it is = 0)"
                )
            for add_n in add_to_box:
                if not isinstance(add_n, int):
                    raise TypeError(
                        "The number of each species to "
                        "be added to each box must be specified "
                        "as an integer"
                    )

        self._mols_to_add = deepcopy(mols_to_add)

    def check_natoms(self):
        """Confirm that the number of existing atoms in each box
        agrees with the number of atoms specified from the combination
        of the number of atoms in each species and the number of each
        species in the box.
        """
        n_species = len(self.species_topologies)
        n_boxes = len(self.boxes)
        atoms_per_species = [len(top.atoms) for top in self.species_topologies]
        atoms_in_box = [
            np.sum(np.multiply(atoms_per_species, self.mols_in_boxes[ibox]))
            for ibox in range(n_boxes)
        ]

        # If the box is empty it should be an mbuild.Box object. If occupied
        # it should be an mbuild.Compound object.
        for ibox, box in enumerate(self.boxes):
            if isinstance(box, mbuild.Compound):
                if box.n_particles != atoms_in_box[ibox]:
                    err_msg = (
                        "The number of atoms in box {} ({}) "
                        "does not match the number of atoms "
                        "calculated to be in the box ({}) from "
                        "the number of atoms per species ({}) "
                        "and the number of species specified "
                        "in each box (mols_in_boxes "
                        "= {})".format(
                            ibox + 1,
                            box.n_particles,
                            atoms_in_box[ibox],
                            atoms_per_species,
                            self.mols_in_boxes[ibox],
                        )
                    )

                    if box.n_particles == 1:
                        addtl_msg = (
                            "NOTE: mbuild.Compound objects "
                            "cannot contain zero particles. "
                            "If you wish to specify an empty "
                            "box please use an mbuild.Box "
                            "object instead."
                        )
                        raise ValueError(addtl_msg + "\n" + err_msg)
                    raise ValueError(err_msg)
            elif isinstance(box, mbuild.Box):
                if sum(self.mols_in_boxes[ibox]) > 0:
                    raise ValueError(
                        "Box {} is an mbuild.Box object "
                        "but species_in_box ({}) indicates that "
                        "molecules should already be present. If you "
                        "wish to provide a starting structure for "
                        "Box {} then it must be a mbuild.Compound "
                        "object".format(
                            ibox + 1, self.mols_in_boxes[ibox], ibox + 1
                        )
                    )

    def fix_bonds(self):
        """Apply the bond length constraints to each molecule in the system"""
        for ibox, box in enumerate(self.boxes):
            if isinstance(box, mbuild.Box):
                continue
            unconstrained_coordinates = box.xyz
            constrained_coordinates = np.zeros(box.xyz.shape)
            n_species = len(self.species_topologies)
            idx_offset = 0
            for isp in range(n_species):
                constrain = self._constrained_species[isp]
                n_mols = self.mols_in_boxes[ibox][isp]
                n_atoms = self._species_topologies[isp].coordinates.shape[0]
                # If no constrained molecule grab all the coordinates
                # in one big block
                if constrain is None:
                    start_idx = idx_offset
                    end_idx = idx_offset + n_mols * n_atoms
                    constrained_coordinates[start_idx:end_idx] = (
                        unconstrained_coordinates[start_idx:end_idx]
                    )
                # Else we apply the constraints one molecule
                # at a time
                else:
                    for imol in range(n_mols):
                        start_idx = idx_offset + imol * n_atoms
                        end_idx = idx_offset + (imol + 1) * n_atoms
                        constrain.update_xyz(
                            unconstrained_coordinates[start_idx:end_idx]
                            * 10.0  # nm to Angstrom
                        )
                        constrain.solve()
                        constrained_coordinates[start_idx:end_idx] = (
                            constrain.xyz / 10.0
                        )  # Angstrom to nm
                # Now we're done with isp; update idx_offset
                idx_offset += n_mols * n_atoms

            box.xyz = constrained_coordinates
