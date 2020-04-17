from copy import deepcopy

import parmed
import warnings


class Moves(object):
    def __init__(self, ensemble, species_topologies):
        """A class to contain all the move probabilities and related
        values required to perform a simulation in ``Cassandra``.

        A Moves object contains a variety of move probabilities
        and other related quantities (e.g., max translation/rotation)
        that are required to run Cassandra. When the moves object
        is created the specified ``ensemble`` and ``species_topologies``
        are used to generate guesses for all required values.
        Depending upon the specifics of your system, these guesses may
        be very reasonable or downright terrible. Use the same
        ``species_topologies`` for your call to ``mosdef_cassandra.System()``
        and ``mosdef_cassandra.Moves()``. Consult the Cassandra user
        manual for more details on the meaning of different move
        probabilities.

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
        ``mosdef_cassandra.Moves``

        """

        if not isinstance(species_topologies, list):
            raise TypeError(
                "species_topologies should be a " "list of species"
            )
        for species in species_topologies:
            if not isinstance(species, parmed.Structure):
                raise TypeError("each species should be a " "parmed.Structure")
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
            self.prob_translate = 0.35
            self.prob_rotate = 0.35
            self.prob_regrow = 0.30
            self.prob_volume = 0.0
            self.prob_insert = 0.0
            self.prob_swap = 0.0
        elif self.ensemble == "npt":
            self.prob_translate = 0.34
            self.prob_rotate = 0.34
            self.prob_regrow = 0.30
            self.prob_volume = 0.02
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
            self.prob_translate = 0.29
            self.prob_rotate = 0.29
            self.prob_regrow = 0.30
            self.prob_volume = 0.02
            self.prob_insert = 0.0
            self.prob_swap = 0.1
        elif self.ensemble == "gemc_npt":
            self.prob_translate = 0.29
            self.prob_rotate = 0.29
            self.prob_regrow = 0.30
            self.prob_volume = 0.02
            self.prob_insert = 0.0
            self.prob_swap = 0.1
        else:
            raise ValueError("Uh oh, how did we end up here?")

        # Max translation and rotations specified per-species-per-box
        self.max_translate = [[2.00] * self._n_species] * self._n_boxes
        self.max_rotate = [[30.0] * self._n_species] * self._n_boxes

        # Prob swap and max vol are per-box
        self.prob_swap_from_box = [1.0 / self._n_boxes] * self._n_boxes

        # Default max deltas for volume moves
        if self.ensemble == "npt" or self.ensemble == "gemc":
            self.max_volume = [500.0]
        elif self.ensemble == "gemc_npt":
            self.max_volume = [500.0, 5000.0]
        else:
            self.max_volume = [0.0]

        # Remaining options are per-species
        self.max_dihedral = [0.0] * self._n_species
        self.sp_insertable = [True] * self._n_species
        self.sp_prob_swap = [1.0] * self._n_species
        self.sp_prob_regrow = [1.0] * self._n_species

        # Here we handle species-wise exceptions
        for ispec, species in enumerate(species_topologies):
            if len(species.atoms) == 1:
                for ibox in range(self._n_boxes):
                    self.max_rotate[ibox][ispec] = 0.0
                self.sp_prob_regrow[ispec] = 0.0
            elif len(species.bonds) == 0:
                print(
                    "Treating {} as a non-insertable rigid species "
                    "since it has no bonds".format(species)
                )
                for ibox in range(self._n_boxes):
                    self.max_translate[ibox][ispec] = 0.0
                    self.max_rotate[ibox][ispec] = 0.0
                self.sp_prob_regrow[ispec] = 0.0
                self.sp_insertable[ispec] = False
                self.sp_prob_swap[ispec] = 0.0

        # Correct species_prob_regrow
        if sum(self.sp_prob_regrow) > 0:
            sp_regrowth_prob = 1.0 / sum(self.sp_prob_regrow)
            for i, prob in enumerate(self.sp_prob_regrow):
                if prob > 0.0:
                    self.sp_prob_regrow[i] = sp_regrowth_prob

        if sum(self.sp_prob_swap) > 0:
            # Correct species_prob_swap
            sp_prob_swap = 1.0 / sum(self.sp_prob_swap)
            for i, insertable in enumerate(self.sp_insertable):
                if insertable:
                    self.sp_prob_swap[i] = sp_prob_swap

        # If all species have no prob regrowth, set prob_regrow to
        # zero and redistribute prob to translate/rotate
        if sum(self.sp_prob_regrow) == 0.0:
            self.prob_translate += self.prob_regrow / 2.0
            self.prob_rotate += self.prob_regrow / 2.0
            self.prob_regrow = 0.0

        # If all species are not rotatable change prob rotation
        # move to zero. Redistribute prob to translate
        if self.ensemble == "gemc" or self.ensemble == "gemc_npt":
            if sum(self.max_rotate[0]) + sum(self.max_rotate[1]) == 0.0:
                self.prob_translate += self.prob_rotate
                self.prob_rotate = 0.0
        else:
            if sum(self.max_rotate[0]) == 0.0:
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
            list of restricted insertion values containing one list per box of species
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

        self._restricted_type = restricted_type
        self._restricted_value = restricted_value

    @property
    def ensemble(self):
        return self._ensemble

    @ensemble.setter
    def ensemble(self, ensemble):
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
        if type(prob_translate) not in (float, int):
            raise TypeError("prob_translate must be of type float")
        else:
            prob_translate = float(prob_translate)
        if prob_translate < 0.0 or prob_translate > 1.0:
            raise ValueError("Probability must be between 0.0 and 1.0.")
        self._prob_translate = prob_translate

    @property
    def prob_rotate(self):
        return self._prob_rotate

    @prob_rotate.setter
    def prob_rotate(self, prob_rotate):
        if type(prob_rotate) not in (float, int):
            raise TypeError("prob_rotate must be of type float")
        else:
            prob_rotate = float(prob_rotate)
        if prob_rotate < 0.0 or prob_rotate > 1.0:
            raise ValueError("Probability must be between 0.0 and 1.0.")
        self._prob_rotate = prob_rotate

    @property
    def prob_angle(self):
        return self._prob_angle

    @prob_angle.setter
    def prob_angle(self, prob_angle):
        if type(prob_angle) not in (float, int):
            raise TypeError("prob_angle must be of type float")
        else:
            prob_angle = float(prob_angle)
        if prob_angle < 0.0 or prob_angle > 1.0:
            raise ValueError("Probability must be between 0.0 and 1.0.")
        self._prob_angle = prob_angle

    @property
    def prob_dihedral(self):
        return self._prob_dihedral

    @prob_dihedral.setter
    def prob_dihedral(self, prob_dihedral):
        if type(prob_dihedral) not in (float, int):
            raise TypeError("prob_dihedral must be of type float")
        else:
            prob_dihedral = float(prob_dihedral)
        if prob_dihedral < 0.0 or prob_dihedral > 1.0:
            raise ValueError("Probability must be between 0.0 and 1.0.")
        self._prob_dihedral = prob_dihedral

    @property
    def prob_regrow(self):
        return self._prob_regrow

    @prob_regrow.setter
    def prob_regrow(self, prob_regrow):
        if type(prob_regrow) not in (float, int):
            raise TypeError("prob_regrow must be of type float")
        else:
            prob_regrow = float(prob_regrow)
        if prob_regrow < 0.0 or prob_regrow > 1.0:
            raise ValueError("Probability must be between 0.0 and 1.0.")
        self._prob_regrow = prob_regrow

    @property
    def prob_volume(self):
        return self._prob_volume

    @prob_volume.setter
    def prob_volume(self, prob_volume):
        if type(prob_volume) not in (float, int):
            raise TypeError("prob_volume must be of type float")
        else:
            prob_volume = float(prob_volume)
        if prob_volume < 0.0 or prob_volume > 1.0:
            raise ValueError("Probability must be between 0.0 and 1.0.")
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
        if type(prob_insert) not in (float, int):
            raise TypeError("prob_insert must be of type float")
        else:
            prob_insert = float(prob_insert)
        if prob_insert < 0.0 or prob_insert > 1.0:
            raise ValueError("Probability must be between 0.0 and 1.0.")
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
        if type(prob_swap) not in (float, int):
            raise TypeError("prob_swap must be of type float")
        else:
            prob_swap = float(prob_swap)
        if prob_swap < 0.0 or prob_swap > 1.0:
            raise ValueError("Probability must be between 0.0 and 1.0.")
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

        if (
            not isinstance(max_translate, list)
            or len(max_translate) != self._n_boxes
        ):
            raise ValueError(
                "max_translate must be a list with shape "
                "(number of boxes, number of species)"
            )
        for max_translate_box in max_translate:
            if (
                not isinstance(max_translate_box, list)
                or len(max_translate_box) != self._n_species
            ):
                raise ValueError(
                    "max_translate must be a list with "
                    "shape (number of boxes, number of species)"
                )
            for max_val in max_translate_box:
                if type(max_val) not in (float, int):
                    raise TypeError(
                        "Max translation values must be " "of type float"
                    )
                else:
                    max_val = float(max_val)
                if max_val < 0.0:
                    raise ValueError(
                        "Max translation values cannot " "be less than zero"
                    )

        self._max_translate = max_translate

    @property
    def max_rotate(self):
        return self._max_rotate

    @max_rotate.setter
    def max_rotate(self, max_rotate):

        if (
            not isinstance(max_rotate, list)
            or len(max_rotate) != self._n_boxes
        ):
            raise ValueError(
                "max_rotate must be a list with shape "
                "(number of boxes, number of species)"
            )
        for max_rotate_box in max_rotate:
            if (
                not isinstance(max_rotate_box, list)
                or len(max_rotate_box) != self._n_species
            ):
                raise ValueError(
                    "max_rotate must be a list with "
                    "shape (number of boxes, number of species)"
                )
            for max_val in max_rotate_box:
                if type(max_val) not in (float, int):
                    raise TypeError(
                        "Max rotation values must be " "of type float"
                    )
                else:
                    max_val = float(max_val)
                if max_val < 0.0:
                    raise ValueError(
                        "Max rotation values cannot " "be less than zero"
                    )

        self._max_rotate = max_rotate

    @property
    def max_dihedral(self):
        return self._max_dihedral

    @max_dihedral.setter
    def max_dihedral(self, max_dihedral):

        if (
            not isinstance(max_dihedral, list)
            or len(max_dihedral) != self._n_species
        ):
            raise ValueError(
                "max_dihedral must be a list with length "
                "(number of species)"
            )
        for max_val in max_dihedral:
            if type(max_val) not in (float, int):
                raise TypeError("Max dihedral values must be " "of type float")
            else:
                max_val = float(max_val)
            if max_val < 0.0:
                raise ValueError(
                    "Max dihedral values cannot " "be less than zero"
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
            raise ValueError(
                "prob_swap_from_box must be a list with length "
                "(number of boxes)"
            )
        for prob_swap in prob_swap_from_box:
            if type(prob_swap) not in (float, int):
                raise TypeError(
                    "Probability of swapping from a box "
                    "must be of type float"
                )
            else:
                prob_swap = float(prob_swap)
            if prob_swap < 0.0:
                raise ValueError(
                    "Probability of swapping from a box "
                    "cannot be less than zero"
                )

        self._prob_swap_from_box = prob_swap_from_box

    @property
    def max_volume(self):
        return self._max_volume

    @max_volume.setter
    def max_volume(self, max_volume):
        if not isinstance(max_volume, list):
            raise ValueError(
                "max_volume must be a list with " "length (number of boxes)"
            )
        if self.ensemble != "gemc":
            if len(max_volume) != self._n_boxes:
                raise ValueError(
                    "max_volume must be a list with "
                    "length (number of boxes)"
                )
        else:
            if len(max_volume) != 1:
                raise ValueError(
                    "max_volume must be a list of " "length (1) for gemc"
                )
        for max_vol in max_volume:
            if type(max_vol) not in (float, int):
                raise TypeError(
                    "Maximum volume change for a box " "must be of type float"
                )
            else:
                max_vol = float(max_vol)
            if max_vol < 0.0:
                raise ValueError(
                    "Maximum volume change for a box "
                    "cannot be less than zero"
                )

        self._max_volume = max_volume

    @property
    def sp_insertable(self):
        return self._sp_insertable

    @sp_insertable.setter
    def sp_insertable(self, sp_insertable):

        if (
            not isinstance(sp_insertable, list)
            or len(sp_insertable) != self._n_species
        ):
            raise ValueError(
                "sp_insertable must be a list with length "
                "(number of species)"
            )
        for insertable in sp_insertable:
            if not isinstance(insertable, bool):
                raise TypeError(
                    "The insertability of each species "
                    "must be provided as a boolean type."
                )

        self._sp_insertable = sp_insertable

    @property
    def sp_prob_swap(self):
        return self._sp_prob_swap

    @sp_prob_swap.setter
    def sp_prob_swap(self, sp_prob_swap):

        if (
            not isinstance(sp_prob_swap, list)
            or len(sp_prob_swap) != self._n_species
        ):
            raise ValueError(
                "sp_prob_swap must be a list with length "
                "(number of species)"
            )
        for prob_swap in sp_prob_swap:
            if type(prob_swap) not in (float, int):
                raise TypeError(
                    "Probability of swapping a species "
                    "must be of type float"
                )
            else:
                prob_swap = float(prob_swap)
            if prob_swap < 0.0:
                raise ValueError(
                    "Probability of swapping a species "
                    "cannot be less than zero"
                )

        self._sp_prob_swap = sp_prob_swap

    @property
    def sp_prob_regrow(self):
        return self._sp_prob_regrow

    @sp_prob_regrow.setter
    def sp_prob_regrow(self, sp_prob_regrow):

        if (
            not isinstance(sp_prob_regrow, list)
            or len(sp_prob_regrow) != self._n_species
        ):
            raise ValueError(
                "sp_prob_regrow must be a list with length "
                "(number of species)"
            )
        for prob_regrow in sp_prob_regrow:
            if type(prob_regrow) not in (float, int):
                raise TypeError(
                    "Probability of regrowing a species "
                    "must be of type float"
                )
            else:
                prob_regrow = float(prob_regrow)
            if prob_regrow < 0.0:
                raise ValueError(
                    "Probability of regrowing a species "
                    "cannot be less than zero"
                )

        self._sp_prob_regrow = sp_prob_regrow

    def print(self):
        """Print the current contents of Moves"""

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

        contents += "\n\nPer species quantities:\n\n"
        contents += "                         "
        for idx in range(self._n_species):
            contents += "species{idx}     ".format(idx=idx + 1)
        contents += "\n"
        contents += "                         "
        for idx in range(self._n_species):
            contents += "========     ".format(idx=idx + 1)
        contents += "\n"
        contents += "Max translate (Ang):     "
        for (box, max_translate_box) in enumerate(self.max_translate):
            for (idx, max_translate) in enumerate(max_translate_box):
                contents += "{max_trans:4.2f}          ".format(
                    max_trans=max_translate
                )
            contents += "(box {box})".format(box=box + 1)
            contents += "\n"
        contents += "Max rotate (deg):        "
        for (box, max_rotate_box) in enumerate(self.max_rotate):
            for (idx, max_rotate) in enumerate(max_rotate_box):
                contents += "{max_rot:4.2f}         ".format(
                    max_rot=max_rotate
                )
            contents += "(box {box})".format(box=box + 1)
            contents += "\n"
        contents += "Insertable:              "
        for (idx, insert) in enumerate(self.sp_insertable):
            contents += "{insert}          ".format(insert=insert)
        contents += "\n"
        contents += "Max dihedral:            "
        for (idx, max_dih) in enumerate(self.max_dihedral):
            contents += "{max_dih:4.2f}          ".format(max_dih=max_dih)
        contents += "\n"
        contents += "Prob swap:               "
        for (idx, prob_swap) in enumerate(self.sp_prob_swap):
            contents += "{prob_swap:4.2f}          ".format(
                prob_swap=prob_swap
            )
        contents += "\n"
        contents += "Prob regrow:             "
        for (idx, prob_regrow) in enumerate(self.sp_prob_regrow):
            contents += "{regrow:4.2f}          ".format(regrow=prob_regrow)
        contents += "\n"

        contents += "\n\nMax volume (Ang^3):\n"
        for (box, max_vol) in enumerate(self.max_volume):
            contents += "Box {box}: {max_vol}\n".format(
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
                        contents += "Box {box}, Species {species}: None\n".format(
                            box=box + 1, species=species + 1
                        )

        print(contents)


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
        if not isinstance(restriction_value, (float, int)):
            raise TypeError(
                "Restriction type is {}. A"
                ' single argument of type "int"'
                'or "float" should be passed'.format(restriction_type)
            )
