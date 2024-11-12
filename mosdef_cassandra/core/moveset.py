from copy import deepcopy
from unyt import dimensions
from mosdef_cassandra.utils.units import validate_unit, validate_unit_list
from gmso.external.convert_parmed import to_parmed
import gmso
import parmed
import warnings
import unyt as u


class MoveSet(object):
    def __init__(self, ensemble, species_topologies):
        """A class to contain all the move probabilities and related
        values required to perform a simulation in ``Cassandra``.

        A MoveSet contains the move probabilities
        and other related quantities (e.g., max translation/rotation)
        that are required to run Cassandra. When the MoveSet
        is created the specified ``ensemble`` and ``species_topologies``
        are used to generate initial guesses for all required values.
        Depending upon the specifics of your system, these guesses may
        be very reasonable or downright terrible. Use the same
        ``species_topologies`` for your call to ``mosdef_cassandra.System()``
        and ``mosdef_cassandra.MoveSet()``.

        Parameters
        ----------
        ensemble : str
            string describing the desired ensembled. Supported
            values include ``'nvt'``, ``'npt'``, ``'gcmc'``,
            ``'gemc'``, ``'gemc_npt'``
        species_topologies : list
            list of ``parmed.Structures``, with one species per element

        Returns
        -------
        ``mosdef_cassandra.MoveSet``

        """

        if not isinstance(species_topologies, list):
            raise TypeError(
                "species_topologies should be a " "list of species"
            )

        if not (
            all(isinstance(top, gmso.Topology) for top in species_topologies)
            or all(
                isinstance(top, parmed.Structure) for top in species_topologies
            )
        ):

            raise TypeError(
                "Each species should be a "
                "parmed.Structure or gmso.Topology"
                "and must be of the same type"
            )

        # Extract self._n_species

        self._n_species = len(species_topologies)

        # Set the ensemble
        self.ensemble = ensemble

        # Infer the number of boxes
        if (
            self.ensemble == "nvt"
            or self.ensemble == "npt"
            or self.ensemble == "gcmc"
        ):
            self._n_boxes = 1
        else:
            self._n_boxes = 2

        # Set '_restricted_typed' and '_restricted_value'
        self._restricted_type = None
        self._restricted_value = None

        # Define default probabilities
        # Most are ensemble-dependent
        self.prob_angle = 0.0
        self.prob_dihedral = 0.0

        if self.ensemble == "nvt":
            self.prob_translate = 0.33
            self.prob_rotate = 0.33
            self.prob_regrow = 0.34
            self.prob_volume = 0.0
            self.prob_insert = 0.0
            self.prob_swap = 0.0
        elif self.ensemble == "npt":
            self.prob_translate = 0.33
            self.prob_rotate = 0.33
            self.prob_regrow = 0.335
            self.prob_volume = 0.005
            self.prob_insert = 0.0
            self.prob_swap = 0.0
        # GCMC sums to 0.9 b/c symmetric prob_delete
        elif self.ensemble == "gcmc":
            self.prob_translate = 0.25
            self.prob_rotate = 0.25
            self.prob_regrow = 0.30
            self.prob_volume = 0.0
            self.prob_insert = 0.1
            self.prob_swap = 0.0
        elif self.ensemble == "gemc":
            self.prob_translate = 0.30
            self.prob_rotate = 0.30
            self.prob_regrow = 0.295
            self.prob_volume = 0.005
            self.prob_insert = 0.0
            self.prob_swap = 0.1
        elif self.ensemble == "gemc_npt":
            self.prob_translate = 0.30
            self.prob_rotate = 0.30
            self.prob_regrow = 0.295
            self.prob_volume = 0.005
            self.prob_insert = 0.0
            self.prob_swap = 0.1
        else:
            raise ValueError("Uh oh, how did we end up here?")

        # Max translation and rotations specified per-species-per-box
        self.max_translate = [
            [2.00 * u.angstrom] * self._n_species
        ] * self._n_boxes
        self.max_rotate = [[30.0 * u.degree] * self._n_species] * self._n_boxes

        # Prob swap and max vol are per-box
        self.prob_swap_from_box = [1.0 / self._n_boxes] * self._n_boxes

        # Default max deltas for volume moves
        if self.ensemble == "npt" or self.ensemble == "gemc":
            self.max_volume = [500.0 * (u.angstrom**3)]
        elif self.ensemble == "gemc_npt":
            self.max_volume = [
                500.0 * (u.angstrom**3),
                5000.0 * (u.angstrom**3),
            ]
        else:
            self.max_volume = [0.0 * (u.angstrom**3)]

        # Set the default CBMC options
        self.cbmc_n_insert = 10
        self.cbmc_n_dihed = 10
        self.cbmc_rcut = 6.0 * u.angstrom

        # Remaining options are per-species
        self.max_dihedral = [0.0 * u.degree] * self._n_species
        self.prob_regrow_species = [1.0] * self._n_species
        if self.ensemble in ["gcmc", "gemc", "gemc_npt"]:
            self.insertable = [True] * self._n_species
        else:
            self.insertable = [False] * self._n_species
        if self.ensemble in ["gemc", "gemc_npt"]:
            self.prob_swap_species = [1.0] * self._n_species
        else:
            self.prob_swap_species = [0.0] * self._n_species

        # Here we handle species-wise exceptions
        for ispec, species in enumerate(species_topologies):

            if isinstance(species, gmso.Topology):
                species = to_parmed(species)

            if len(species.atoms) == 1:
                for ibox in range(self._n_boxes):
                    self.max_rotate[ibox][ispec] = 0.0 * u.degree
                self.prob_regrow_species[ispec] = 0.0
            elif len(species.bonds) == 0:
                print(
                    "Treating {} as a non-insertable rigid species "
                    "since it has no bonds".format(species)
                )
                for ibox in range(self._n_boxes):
                    self.max_translate[ibox][ispec] = 0.0 * u.angstrom
                    self.max_rotate[ibox][ispec] = 0.0 * u.degree
                self.prob_regrow_species[ispec] = 0.0
                self.insertable[ispec] = False
                self.prob_swap_species[ispec] = 0.0

        # Correct species_prob_regrow
        if sum(self.prob_regrow_species) > 0:
            sp_regrowth_prob = 1.0 / sum(self.prob_regrow_species)
            for i, prob in enumerate(self.prob_regrow_species):
                if prob > 0.0:
                    self.prob_regrow_species[i] = sp_regrowth_prob

        if sum(self.prob_swap_species) > 0:
            # Correct species_prob_swap
            prob_swap_species = 1.0 / sum(self.prob_swap_species)
            for idx, insert in enumerate(self.insertable):
                if insert:
                    self.prob_swap_species[idx] = prob_swap_species

        # If all species have no prob regrowth, set prob_regrow to
        # zero and redistribute prob to translate/rotate
        if sum(self.prob_regrow_species) == 0.0:
            self.prob_translate += self.prob_regrow / 2.0
            self.prob_rotate += self.prob_regrow / 2.0
            self.prob_regrow = 0.0

        # If all species are not rotatable change prob rotation
        # move to zero. Redistribute prob to translate
        if self.ensemble == "gemc" or self.ensemble == "gemc_npt":
            if (
                sum(self.max_rotate[0]).to_value()
                + sum(self.max_rotate[1]).to_value()
                == 0.0
            ):
                self.prob_translate += self.prob_rotate
                self.prob_rotate = 0.0
        else:
            if sum(self.max_rotate[0]).to_value() == 0.0:
                self.prob_translate += self.prob_rotate
                self.prob_rotate = 0.0

    def add_restricted_insertions(
        self, species_topologies, restricted_type, restricted_value
    ):
        """Add restricted insertions for specific species and boxes

        Parameters
        ----------
        species_topologies : list
            list of ``parmed.Structures`` containing one list per box of species
        restricted_type : list
            list of restricted insertion types containing one list per box of species
        restricted_value : list
            list of restricted insertion values (unyt arrays) containing one list per box of species
        """
        if self._restricted_type and self._restricted_value:
            warnings.warn(
                "Restricted insertion has been previously"
                " added and will be replaced."
            )
        if self.ensemble not in ["gcmc", "gemc", "gemc_npt"]:
            raise ValueError(
                "Restricted insertions are only valid for"
                " 'gcmc', 'gemc', and 'gemc_npt' ensembles."
            )
        if len(restricted_type) != len(restricted_value):
            raise ValueError(
                "Length of 'restricted_type' and "
                " 'restricted_value' must match."
            )
        for box in restricted_type:
            if isinstance(box, (str, int, float)):
                raise TypeError(
                    "Restricted type must be passed as a list"
                    " of lists corresponding to each box."
                )
            if len(box) != len(species_topologies):
                raise ValueError(
                    "Length of 'species' and "
                    " length of box list in 'restricted_type'"
                    " must match.  `species` has a length of {}"
                    " and the box list in 'restricted_type' has a "
                    " length of {}".format(len(species_topologies), len(box))
                )
        for box in restricted_value:
            if isinstance(box, (str, int, float)):
                raise TypeError(
                    "Restricted value must be passed as a list"
                    " of lists corresponding to each box."
                )
            if len(box) != len(species_topologies):
                raise ValueError(
                    "Length of 'species' and "
                    " length of species list in 'restricted_value'"
                    " must match.  `species` has a length of {}"
                    " and the box list in 'restricted_value' has a "
                    " length of {}".format(len(species_topologies), len(box))
                )
        if self.ensemble == "gcmc" and len(restricted_type) != 1:
            raise ValueError(
                "GCMC ensemble contains 1 box but"
                " `restricted_type` of length {}"
                " was passed.".format(len(restricted_type))
            )
        if self.ensemble in ["gemc", "gemc_npt"] and len(restricted_type) != 2:
            raise ValueError(
                "GEMC ensembles contain 2 boxes but"
                " `restricted_type` of length {}"
                " was passed.".format(len(restricted_type))
            )

        for types, values in zip(restricted_type, restricted_value):
            for typ, val in zip(types, values):
                if not typ and not val:
                    pass
                elif typ and not val:
                    raise ValueError(
                        "`restricted_type` {} was passed"
                        " but `restricted_value` is None.".format(typ, val)
                    )
                elif val and not typ:
                    raise ValueError(
                        "`restricted_value` {} was passed"
                        " but `restricted_type` is None.".format(val, typ)
                    )
                else:
                    _check_restriction_type(typ, val)
                    # Check units of restricted value
                    if typ == "interface":
                        [validate_unit(i, dimensions.length) for i in val]
                    else:
                        validate_unit(val, dimensions.length)

        self._restricted_type = restricted_type
        self._restricted_value = restricted_value

    @property
    def ensemble(self):
        return self._ensemble

    @ensemble.setter
    def ensemble(self, ensemble):
        if hasattr(self, "_ensemble"):
            raise AttributeError(
                "Ensemble cannot be changed. Please create a new MoveSet instead."
            )
        valid_ensembles = ["nvt", "npt", "gcmc", "gemc", "gemc_npt"]
        if ensemble not in valid_ensembles:
            raise ValueError(
                'Invalid ensemble "{}" Supported '
                "ensembles include {}".format(ensemble, valid_ensembles)
            )
        self._ensemble = ensemble

    @property
    def prob_translate(self):
        return self._prob_translate

    @prob_translate.setter
    def prob_translate(self, prob_translate):
        prob_translate = self._validate_probability(
            prob_translate,
            "prob_translate",
        )
        self._prob_translate = prob_translate

    @property
    def prob_rotate(self):
        return self._prob_rotate

    @prob_rotate.setter
    def prob_rotate(self, prob_rotate):
        prob_rotate = self._validate_probability(
            prob_rotate,
            "prob_rotate",
        )
        self._prob_rotate = prob_rotate

    @property
    def prob_angle(self):
        return self._prob_angle

    @prob_angle.setter
    def prob_angle(self, prob_angle):
        prob_angle = self._validate_probability(
            prob_angle,
            "prob_angle",
        )
        self._prob_angle = prob_angle

    @property
    def prob_dihedral(self):
        return self._prob_dihedral

    @prob_dihedral.setter
    def prob_dihedral(self, prob_dihedral):
        prob_dihedral = self._validate_probability(
            prob_dihedral,
            "prob_dihedral",
        )
        self._prob_dihedral = prob_dihedral

    @property
    def prob_regrow(self):
        return self._prob_regrow

    @prob_regrow.setter
    def prob_regrow(self, prob_regrow):
        prob_regrow = self._validate_probability(
            prob_regrow,
            "prob_regrow",
        )
        self._prob_regrow = prob_regrow

    @property
    def prob_volume(self):
        return self._prob_volume

    @prob_volume.setter
    def prob_volume(self, prob_volume):
        prob_volume = self._validate_probability(
            prob_volume,
            "prob_volume",
        )
        if prob_volume > 0.0:
            if self.ensemble == "nvt" or self.ensemble == "gcmc":
                raise ValueError(
                    "Ensemble is {}. prob_volume cannot be "
                    "non-zero in the {} ensemble".format(
                        self._ensemble, self.ensemble
                    )
                )
        elif prob_volume == 0.0:
            if (
                self.ensemble == "npt"
                or self.ensemble == "gemc"
                or self.ensemble == "gemc_npt"
            ):
                raise ValueError(
                    "Ensemble is {}. prob_volume must be "
                    "> 0.0 in this ensemble".format(self.ensemble)
                )
        # Pass all checks. Update prob_volume.
        self._prob_volume = prob_volume

    @property
    def prob_insert(self):
        return self._prob_insert

    @prob_insert.setter
    def prob_insert(self, prob_insert):
        prob_insert = self._validate_probability(
            prob_insert,
            "prob_insert",
        )
        if self.ensemble != "gcmc" and prob_insert != 0.0:
            raise ValueError(
                "Ensemble is {}. Insertion probability "
                "must be = 0.0".format(self.ensemble)
            )
        if self.ensemble == "gcmc" and prob_insert == 0.0:
            raise ValueError(
                "Ensemble is {}. Insertion probability "
                "must be > 0.0".format(self.ensemble)
            )
        self._prob_insert = prob_insert

    @property
    def prob_swap(self):
        return self._prob_swap

    @prob_swap.setter
    def prob_swap(self, prob_swap):
        prob_swap = self._validate_probability(
            prob_swap,
            "prob_swap",
        )
        if self.ensemble != "gemc" and self.ensemble != "gemc_npt":
            if prob_swap != 0.0:
                raise ValueError(
                    "Ensemble is {}. Swapping probability "
                    "must be = 0.0".format(self.ensemble)
                )
        if self.ensemble == "gemc" or self.ensemble == "gemc_npt":
            if prob_swap == 0.0:
                raise ValueError(
                    "Ensemble is {}. Swapping probability "
                    "must be > 0.0".format(self.ensemble)
                )

        self._prob_swap = prob_swap

    @property
    def max_translate(self):
        return self._max_translate

    @max_translate.setter
    def max_translate(self, max_translate):
        max_translate = validate_unit_list(
            max_translate,
            (self._n_boxes, self._n_species),
            dimensions.length,
            "max_translate",
        )
        for max_val in max_translate.flatten():
            if max_val.to_value() < 0.0:
                raise ValueError(
                    "Max translation values cannot be less than zero"
                )
        self._max_translate = max_translate

    @property
    def max_rotate(self):
        return self._max_rotate

    @max_rotate.setter
    def max_rotate(self, max_rotate):
        max_rotate = validate_unit_list(
            max_rotate,
            (self._n_boxes, self._n_species),
            dimensions.angle,
            "max_rotate",
        )
        for max_val in max_rotate.flatten():
            if (
                max_val.to_value("degree") < 0.0
                or max_val.to_value("degree") > 360.0
            ):
                raise ValueError(
                    "Max rotation values must be between 0.0 and 360.0 degrees."
                )
        self._max_rotate = max_rotate

    @property
    def max_dihedral(self):
        return self._max_dihedral

    @max_dihedral.setter
    def max_dihedral(self, max_dihedral):
        max_dihedral = validate_unit_list(
            max_dihedral,
            (self._n_species,),
            dimensions.angle,
            "max_dihedral",
        )
        for max_val in max_dihedral:
            if (
                max_val.to_value("degree") < 0.0
                or max_val.to_value("degree") > 360.0
            ):
                raise ValueError(
                    "Max dihedral rotation values must be between 0.0 and 360.0 degrees."
                )
        self._max_dihedral = max_dihedral

    @property
    def prob_swap_from_box(self):
        return self._prob_swap_from_box

    @prob_swap_from_box.setter
    def prob_swap_from_box(self, prob_swap_from_box):

        if (
            not isinstance(prob_swap_from_box, list)
            or len(prob_swap_from_box) != self._n_boxes
        ):
            raise TypeError(
                "prob_swap_from_box must be a list with length "
                "(number of boxes)"
            )
        validated_prob_swap_from_box = []
        for prob_swap in prob_swap_from_box:
            prob_swap = self._validate_probability(
                prob_swap,
                "prob_swap_from_box",
            )
            validated_prob_swap_from_box.append(prob_swap)
        self._prob_swap_from_box = validated_prob_swap_from_box

    @property
    def max_volume(self):
        return self._max_volume

    @max_volume.setter
    def max_volume(self, max_volume):
        if type(max_volume) not in (list, u.unyt_array):
            if self.ensemble == "gemc_npt":
                max_volume = [max_volume] * self._n_boxes
            else:
                max_volume = [max_volume]

        if self.ensemble == "gemc_npt":
            shape = (self._n_boxes,)
        else:
            shape = (1,)

        max_volume = validate_unit_list(
            max_volume,
            shape,
            dimensions.length**3,
            "max_volume",
        )
        for max_vol in max_volume.flatten():
            if max_vol < 0.0:
                raise ValueError("max_volume cannot be less than zero.")
        self._max_volume = max_volume

    @property
    def insertable(self):
        return self._insertable

    @insertable.setter
    def insertable(self, insertable):

        if (
            not isinstance(insertable, list)
            or len(insertable) != self._n_species
        ):
            raise TypeError(
                "insertable must be a list with length " "(number of species)"
            )
        for insert in insertable:
            if not isinstance(insert, bool):
                raise TypeError(
                    "The insertability of each species "
                    "must be provided as a boolean type."
                )

        self._insertable = insertable

    @property
    def prob_swap_species(self):
        return self._prob_swap_species

    @prob_swap_species.setter
    def prob_swap_species(self, prob_swap_species):

        if (
            not isinstance(prob_swap_species, list)
            or len(prob_swap_species) != self._n_species
        ):
            raise TypeError(
                "prob_swap_species must be a list with length "
                "(number of species)"
            )
        validated_prob_swap_species = []
        for prob_swap in prob_swap_species:
            prob_swap = self._validate_probability(
                prob_swap,
                "prob_swap_species",
            )
            validated_prob_swap_species.append(prob_swap)
        self._prob_swap_species = validated_prob_swap_species

    @property
    def prob_regrow_species(self):
        return self._prob_regrow_species

    @prob_regrow_species.setter
    def prob_regrow_species(self, prob_regrow_species):
        if (
            not isinstance(prob_regrow_species, list)
            or len(prob_regrow_species) != self._n_species
        ):
            raise TypeError(
                "prob_regrow_species must be a list with length "
                "(number of species)"
            )
        validated_prob_regrow_species = []
        for prob_regrow in prob_regrow_species:
            prob_regrow = self._validate_probability(
                prob_regrow, "prob_regrow"
            )
            validated_prob_regrow_species.append(prob_regrow)
        self._prob_regrow_species = validated_prob_regrow_species

    @property
    def cbmc_n_insert(self):
        return self._cbmc_n_insert

    @cbmc_n_insert.setter
    def cbmc_n_insert(self, cbmc_n_insert):
        if type(cbmc_n_insert) != int:
            raise TypeError("cbmc_n_insert must be of type int")
        if cbmc_n_insert <= 0:
            raise ValueError("cbmc_n_insert must be greater than zero")
        self._cbmc_n_insert = cbmc_n_insert

    @property
    def cbmc_n_dihed(self):
        return self._cbmc_n_dihed

    @cbmc_n_dihed.setter
    def cbmc_n_dihed(self, cbmc_n_dihed):
        if type(cbmc_n_dihed) != int:
            raise TypeError("cbmc_n_dihed must be of type int")
        if cbmc_n_dihed <= 0:
            raise ValueError("cbmc_n_dihed must be greater than zero")
        self._cbmc_n_dihed = cbmc_n_dihed

    @property
    def cbmc_rcut(self):
        return self._cbmc_rcut

    @cbmc_rcut.setter
    def cbmc_rcut(self, cbmc_rcut):
        if type(cbmc_rcut) not in (list, u.unyt_array):
            cbmc_rcut = [cbmc_rcut] * self._n_boxes
        cbmc_rcut = validate_unit_list(
            cbmc_rcut,
            (self._n_boxes,),
            dimensions.length,
            "cbmc_rcut",
        )

        for rcut in cbmc_rcut.flatten():
            if rcut.to_value() < 0.0:
                raise ValueError("cbmc_rcut cannot be less than zero.")

        self._cbmc_rcut = cbmc_rcut

    def print(self):
        """Print the current contents of the MoveSet"""

        contents = """
Ensemble:  {ensemble}

Probability of selecting each move type:

    Translate: {prob_translate}
    Rotate:    {prob_rotate}
    Regrow:    {prob_regrow}
    Volume:    {prob_volume}
    Insert:    {prob_insert}
    Delete:    {prob_delete}
    Swap:      {prob_swap}
    Angle:     {prob_angle}
    Dihedral:  {prob_dihedral}
""".format(
            ensemble=self.ensemble,
            prob_translate=self.prob_translate,
            prob_rotate=self.prob_rotate,
            prob_regrow=self.prob_regrow,
            prob_volume=self.prob_volume,
            prob_insert=self.prob_insert,
            prob_delete=self.prob_insert,
            prob_swap=self.prob_swap,
            prob_angle=self.prob_angle,
            prob_dihedral=self.prob_dihedral,
        )

        contents += """
CBMC selections:

    Number of trial positions: {n_insert}
    Number of trial dihedral angles: {n_dihed}
    CBMC cutoff(s): 
""".format(
            n_insert=self.cbmc_n_insert,
            n_dihed=self.cbmc_n_dihed,
        )

        for idx, value in enumerate(self.cbmc_rcut):
            contents += "        Box {}: {}\n".format(idx + 1, value)

        contents += "\n\nPer species quantities:\n\n"
        contents += "                             "
        for idx in range(self._n_species):
            contents += "species{idx}     ".format(idx=idx + 1)
        contents += "\n"
        contents += "                             "
        for idx in range(self._n_species):
            contents += "========     ".format(idx=idx + 1)
        contents += "\n"
        contents += "    Max translate (Ang):     "
        for box, max_translate_box in enumerate(self.max_translate):
            if box > 0:
                contents += "                             "
            for idx, max_translate in enumerate(max_translate_box):
                contents += "{max_trans:4.2f}          ".format(
                    max_trans=max_translate
                )
            contents += "(Box {box})".format(box=box + 1)
            contents += "\n"
        contents += "    Max rotate (deg):        "
        for box, max_rotate_box in enumerate(self.max_rotate):
            if box > 0:
                contents += "                             "
            for idx, max_rotate in enumerate(max_rotate_box):
                contents += "{max_rot:4.2f}         ".format(
                    max_rot=max_rotate
                )
            contents += "(Box {box})".format(box=box + 1)
            contents += "\n"
        contents += "    Insertable:              "
        for idx, insert in enumerate(self.insertable):
            contents += "{insert}          ".format(insert=insert)
        contents += "\n"
        contents += "    Max dihedral:            "
        for idx, max_dih in enumerate(self.max_dihedral):
            contents += "{max_dih:4.2f}          ".format(max_dih=max_dih)
        contents += "\n"
        contents += "    Prob swap:               "
        for idx, prob_swap in enumerate(self.prob_swap_species):
            contents += "{prob_swap:4.2f}          ".format(
                prob_swap=prob_swap
            )
        contents += "\n"
        contents += "    Prob regrow:             "
        for idx, prob_regrow in enumerate(self.prob_regrow_species):
            contents += "{regrow:4.2f}          ".format(regrow=prob_regrow)
        contents += "\n"

        contents += "\n\nMax volume (Ang^3):\n"
        for box, max_vol in enumerate(self.max_volume):
            contents += "    Box {box}: {max_vol}\n".format(
                box=box + 1, max_vol=max_vol
            )

        if self._restricted_type != None:
            contents += "\nRestricted Insertions (Ang):\n"
            for box in range(self._n_boxes):
                for species, (typ, value) in enumerate(
                    zip(
                        self._restricted_type[box], self._restricted_value[box]
                    )
                ):
                    if typ == "sphere":
                        contents += "Box {box}, Species {species}: sphere, R = {r_value}\n".format(
                            box=box + 1, species=species + 1, r_value=value
                        )
                    elif typ == "cylinder":
                        contents += "Box {box}, Species {species}: cylinder, R = {r_value}\n".format(
                            box=box + 1, species=species + 1, r_value=value
                        )
                    elif typ == "slitpore":
                        contents += "Box {box}, Species {species}: slitpore, z_max = {z_max}\n".format(
                            box=box + 1, species=species + 1, z_max=value
                        )
                    elif typ == "interface":
                        contents += "Box {box}, Species {species}: interface, z_min = {z_min}, z_max = {z_max}\n".format(
                            box=box + 1,
                            species=species + 1,
                            z_min=value[0],
                            z_max=value[1],
                        )
                    else:
                        contents += (
                            "Box {box}, Species {species}: None\n".format(
                                box=box + 1, species=species + 1
                            )
                        )

        print(contents)

    def _validate_probability(self, probability, name):
        if type(probability) not in (float, int):
            raise TypeError(f"{name} must be of type float")
        else:
            probability = float(probability)
        if probability < 0.0 or probability > 1.0:
            raise ValueError(f"{name} must be between 0.0 and 1.0.")
        return probability


def _check_restriction_type(restriction_type, restriction_value):
    valid_restrict_types = ["sphere", "cylinder", "slitpore", "interface"]
    # Check restriction insertion type
    if restriction_type not in valid_restrict_types:
        raise ValueError(
            'Invalid restriction type "{}".  Supported '
            "restriction types include {}".format(
                restriction_type, valid_restrict_types
            )
        )
    # Check if correct number of arguments passed
    if restriction_type == "interface":
        if len(restriction_value) != 2:
            raise ValueError(
                "Invalid number of arguments passed."
                "{} arguments for restriction type {}"
                "were passed.  2 are required".format(
                    len(restriction_value), restriction_type
                )
            )
    else:
        if not isinstance(restriction_value, u.unyt_array):
            raise TypeError(
                "Invalid type for `restriction_value` passed. A"
                " single argument of type `unyt_array"
                " should be passed".format(restriction_type)
            )
