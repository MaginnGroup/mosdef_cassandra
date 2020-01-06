
import parmed

class Moves(object):

    def __init__(self,ensemble,species_topologies):
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

        valid_ensembles = ['nvt','npt','gcmc','gemc','gemc_npt']
        if ensemble not in valid_ensembles:
            raise ValueError('Invalid ensemble "{}" Supported '
                'ensembles include {}'.format(
                    ensemble,valid_ensembles))
        self.ensemble = ensemble

        if not isinstance(species_topologies,list):
            raise TypeError('species_topologies should be a '
                            'list of species')
        for species in species_topologies:
            if not isinstance(species,parmed.Structure):
                raise TypeError('each species should be a '
                                'parmed.Structure')

        # Start by setting all move probabilities to zero
        self.prob_translate = 0.0
        self.prob_rotate = 0.0
        self.prob_angle = 0.0
        self.prob_dihedral = 0.0
        self.prob_regrow = 0.0
        self.prob_volume = 0.0
        self.prob_insert = 0.0
        self.prob_swap = 0.0

        # Define default probabilities for each ensemble
        if self.ensemble == 'nvt':
            self.prob_translate = 0.35
            self.prob_rotate = 0.35
            self.prob_regrow = 0.30
        elif self.ensemble == 'npt':
            self.prob_translate = 0.34
            self.prob_rotate = 0.34
            self.prob_regrow = 0.30
            self.prob_volume = 0.02
        # GCMC sums to 0.9 b/c symmetric prob_delete
        elif self.ensemble == 'gcmc':
            self.prob_translate = 0.25
            self.prob_rotate = 0.25
            self.prob_regrow = 0.30
            self.prob_insert = 0.1
        elif self.ensemble == 'gemc':
            self.prob_translate = 0.29
            self.prob_rotate = 0.29
            self.prob_regrow = 0.30
            self.prob_swap = 0.1
            self.prob_volume = 0.02
        elif self.ensemble == 'gemc_npt':
            self.prob_translate = 0.29
            self.prob_rotate = 0.29
            self.prob_regrow = 0.30
            self.prob_swap = 0.1
            self.prob_volume = 0.02
        else:
            raise ValueError('Uh oh, how did we end up here?')

        # Max translation and rotations specified per-box
        self.max_translate = []
        self.max_rotate = []
        self.max_translate.append([])
        self.max_rotate.append([])
        self.prob_swap_from_box = [1.0]
        self.max_volume = []
        # Other options not
        self.sp_insertable = []
        self.sp_prob_swap = []
        self.sp_prob_regrow = []

        # Here we handle species-wise selections
        # including per-box selections for box1
        for species in species_topologies:
            if len(species.atoms) == 1:
                self.max_translate[0].append(2.0)
                self.max_rotate[0].append(0.0)
                self.sp_prob_regrow.append(0.0)
                self.sp_insertable.append(True)
                self.sp_prob_swap.append(1.0)
            elif len(species.bonds) == 0:
                print('Treating {} as a non-insertable rigid species '
                        'since it has no bonds'.format(species))
                self.max_translate[0].append(0.0)
                self.max_rotate[0].append(0.0)
                self.sp_prob_regrow.append(0.0)
                self.sp_insertable.append(False)
                self.sp_prob_swap.append(0.0)
            else:
                self.max_translate[0].append(2.0)
                self.max_rotate[0].append(30.0)
                self.sp_prob_regrow.append(1.0)
                self.sp_insertable.append(True)
                self.sp_prob_swap.append(1.0)

        # Add for box two
        if self.ensemble  == 'gemc' or self.ensemble == 'gemc_npt':
            self.max_translate.append([])
            self.max_rotate.append([])
            self.prob_swap_from_box = [0.5,0.5]
            for species in species_topologies:
                if len(species.atoms) == 1:
                    self.max_translate[1].append(2.0)
                    self.max_rotate[1].append(0.0)
                elif len(species.bonds) == 0:
                    self.max_translate[1].append(0.0)
                    self.max_rotate[1].append(0.0)
                else:
                    self.max_translate[1].append(2.0)
                    self.max_rotate[1].append(30.0)

        # Default max deltas for volume moves
        if ( self.ensemble == 'npt' or
             self.ensemble == 'gemc'
             or self.ensemble == 'gemc_npt' ):
            self.max_volume.append(500.)
            if self.ensemble == 'gemc_npt':
                self.max_volume.append(5000.)
        else:
            self.max_volume.append(0.0)

        # Correct species_prob_regrow
        if sum(self.sp_prob_regrow) > 0:
            sp_regrowth_prob = 1.0/sum(self.sp_prob_regrow)
            for i,prob in enumerate(self.sp_prob_regrow):
                if prob > 0.0:
                    self.sp_prob_regrow[i] = sp_regrowth_prob

        if sum(self.sp_prob_swap) > 0:
            # Correct species_prob_swap
            sp_prob_swap = 1.0/sum(self.sp_prob_swap)
            for i,insertable in enumerate(self.sp_insertable):
                if insertable:
                    self.sp_prob_swap[i] = sp_prob_swap

        # If all species have no prob regrowth, set prob_regrow to
        # zero and redistribute prob to translate/rotate
        if sum(self.sp_prob_regrow) == 0.0:
            self.prob_translate += self.prob_regrow/2.
            self.prob_rotate += self.prob_regrow/2.
            self.prob_regrow = 0.0

        # If all species are not rotatable change prob rotation
        # move to zero. Redistribute prob to translate
        if self.ensemble == 'gemc' or self.ensemble == 'gemc_npt':
            if sum(self.max_rotate[0]) + sum(self.max_rotate[1]) == 0.0:
                self.prob_translate += self.prob_rotate
                self.prob_rotate = 0.0
        else:
            if sum(self.max_rotate[0]) == 0.0:
                self.prob_translate += self.prob_rotate
                self.prob_rotate = 0.0


