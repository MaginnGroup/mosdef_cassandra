
from copy import deepcopy
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

        if not isinstance(species_topologies,list):
            raise TypeError('species_topologies should be a '
                            'list of species')
        for species in species_topologies:
            if not isinstance(species,parmed.Structure):
                raise TypeError('each species should be a '
                                'parmed.Structure')
        # Extract self._n_species
        self._n_species = len(species_topologies)

        # Set the ensemble
        self.ensemble = ensemble

        # Infer the number of boxes
        if (self.ensemble == 'nvt' or
            self.ensemble == 'npt' or
            self.ensemble == 'gcmc'):
            self._n_boxes = 1
        else:
            self._n_boxes = 2

        # Define default probabilities
        # Most are ensemble-dependent
        self.prob_angle = 0.0
        self.prob_dihedral = 0.0

        if self.ensemble == 'nvt':
            self.prob_translate = 0.35
            self.prob_rotate = 0.35
            self.prob_regrow = 0.30
            self.prob_volume = 0.0
            self.prob_insert = 0.0
            self.prob_swap = 0.0
        elif self.ensemble == 'npt':
            self.prob_translate = 0.34
            self.prob_rotate = 0.34
            self.prob_regrow = 0.30
            self.prob_volume = 0.02
            self.prob_insert = 0.0
            self.prob_swap = 0.0
        # GCMC sums to 0.9 b/c symmetric prob_delete
        elif self.ensemble == 'gcmc':
            self.prob_translate = 0.25
            self.prob_rotate = 0.25
            self.prob_regrow = 0.30
            self.prob_volume = 0.0
            self.prob_insert = 0.1
            self.prob_swap = 0.0
        elif self.ensemble == 'gemc':
            self.prob_translate = 0.29
            self.prob_rotate = 0.29
            self.prob_regrow = 0.30
            self.prob_volume = 0.02
            self.prob_insert = 0.0
            self.prob_swap = 0.1
        elif self.ensemble == 'gemc_npt':
            self.prob_translate = 0.29
            self.prob_rotate = 0.29
            self.prob_regrow = 0.30
            self.prob_volume = 0.02
            self.prob_insert = 0.0
            self.prob_swap = 0.1
        else:
            raise ValueError('Uh oh, how did we end up here?')

        # Max translation and rotations specified per-species-per-box
        self.max_translate = [ [2.00] * self._n_species ] * self._n_boxes
        self.max_rotate    = [ [30.0] * self._n_species ] * self._n_boxes

        # Prob swap and max vol are per-box
        self.prob_swap_from_box = [1.0/self._n_boxes] * self._n_boxes

        # Default max deltas for volume moves
        if ( self.ensemble == 'npt' or
             self.ensemble == 'gemc' ):
            self.max_volume = [500.]
        elif self.ensemble == 'gemc_npt':
            self.max_volume = [500., 5000.]
        else:
            self.max_volume = [0.0]

        # Remaining options are per-species
        self.sp_insertable = [True] * self._n_species
        self.sp_prob_swap = [1.0/self._n_species] * self._n_species
        self.sp_prob_regrow = [1.0/self._n_species] * self._n_species

        # Here we handle species-wise exceptions
        for ispec,species in enumerate(species_topologies):
            if len(species.atoms) == 1:
                for ibox in range(self._n_boxes):
                    self.max_rotate[ibox][ispec] = 0.0
                self.sp_prob_regrow[ispec] = 0.0
            elif len(species.bonds) == 0:
                print('Treating {} as a non-insertable rigid species '
                        'since it has no bonds'.format(species))
                for ibox in range(self._n_boxes):
                    self.max_translate[ibox][ispec] = 0.0
                    self.max_rotate[ibox][ispec] = 0.0
                self.sp_prob_regrow[ispec] = 0.0
                self.sp_insertable[ispec] = False
                self.sp_prob_swap[ispec] = 0.0

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

    @property
    def ensemble(self):
        return self._ensemble

    @ensemble.setter
    def ensemble(self,ensemble):
        valid_ensembles = ['nvt','npt','gcmc','gemc','gemc_npt']
        if ensemble not in valid_ensembles:
            raise ValueError('Invalid ensemble "{}" Supported '
                'ensembles include {}'.format(
                    ensemble,valid_ensembles))
        self._ensemble = ensemble

    @property
    def prob_translate(self):
        return self._prob_translate

    @prob_translate.setter
    def prob_translate(self,prob_translate):
        try:
            prob_translate = float(prob_translate)
        except (ValueError,TypeError):
            raise TypeError('prob_translate must be of type float')
        if prob_translate < 0.0 or prob_translate > 1.0:
            raise ValueError('Probability must be between 0.0 and 1.0.')
        self._prob_translate = prob_translate

    @property
    def prob_rotate(self):
        return self._prob_rotate

    @prob_rotate.setter
    def prob_rotate(self,prob_rotate):
        try:
            prob_rotate = float(prob_rotate)
        except (ValueError,TypeError):
            raise TypeError('prob_rotate must be of type float')
        if prob_rotate < 0.0 or prob_rotate > 1.0:
            raise ValueError('Probability must be between 0.0 and 1.0.')
        self._prob_rotate = prob_rotate

    @property
    def prob_angle(self):
        return self._prob_angle

    @prob_angle.setter
    def prob_angle(self,prob_angle):
        try:
            prob_angle = float(prob_angle)
        except (ValueError,TypeError):
            raise TypeError('prob_angle must be of type float')
        if prob_angle < 0.0 or prob_angle > 1.0:
            raise ValueError('Probability must be between 0.0 and 1.0.')
        self._prob_angle = prob_angle

    @property
    def prob_dihedral(self):
        return self._prob_dihedral

    @prob_dihedral.setter
    def prob_dihedral(self,prob_dihedral):
        try:
            prob_dihedral = float(prob_dihedral)
        except (ValueError,TypeError):
            raise TypeError('prob_dihedral must be of type float')
        if prob_dihedral < 0.0 or prob_dihedral > 1.0:
            raise ValueError('Probability must be between 0.0 and 1.0.')
        self._prob_dihedral = prob_dihedral

    @property
    def prob_angle(self):
        return self._prob_angle

    @prob_angle.setter
    def prob_angle(self,prob_angle):
        try:
            prob_angle = float(prob_angle)
        except (ValueError,TypeError):
            raise TypeError('prob_angle must be of type float')
        if prob_angle < 0.0 or prob_angle > 1.0:
            raise ValueError('Probability must be between 0.0 and 1.0.')
        self._prob_angle = prob_angle

    @property
    def prob_regrow(self):
        return self._prob_regrow

    @prob_regrow.setter
    def prob_regrow(self,prob_regrow):
        try:
            prob_regrow = float(prob_regrow)
        except (ValueError,TypeError):
            raise TypeError('prob_regrow must be of type float')
        if prob_regrow < 0.0 or prob_regrow > 1.0:
            raise ValueError('Probability must be between 0.0 and 1.0.')
        self._prob_regrow = prob_regrow

    @property
    def prob_volume(self):
        return self._prob_volume

    @prob_volume.setter
    def prob_volume(self,prob_volume):
        try:
            prob_volume = float(prob_volume)
        except (ValueError,TypeError):
            raise TypeError('prob_volume must be of type float')
        if prob_volume < 0.0 or prob_volume > 1.0:
            raise ValueError('Probability must be between 0.0 and 1.0.')
        if prob_volume > 0.0:
            if self._ensemble == 'nvt' or self._ensemble == 'gcmc':
                raise ValueError('Ensemble is {}. prob_volume cannot be '
                        'non-zero in the {} ensemble'.format(
                         self._ensemble,self.ensemble))
        elif prob_volume == 0.0:
            if ( self._ensemble == 'npt' or
                 self._ensemble == 'gemc' or
                 self._ensemble == 'gemc_npt' ):
                raise ValueError('Ensemble is {}. prob_volume must be '
                        '> 0.0 in this ensemble'.format(self._ensemble))
        # Pass all checks. Update prob_volume.
        self._prob_volume = prob_volume

    @property
    def prob_insert(self):
        return self._prob_insert

    @prob_insert.setter
    def prob_insert(self,prob_insert):
        try:
            prob_insert = float(prob_insert)
        except (ValueError,TypeError):
            raise TypeError('prob_insert must be of type float')
        if prob_insert < 0.0 or prob_insert > 1.0:
            raise ValueError('Probability must be between 0.0 and 1.0.')
        self._prob_insert = prob_insert

    @property
    def prob_swap(self):
        return self._prob_swap

    @prob_swap.setter
    def prob_swap(self,prob_swap):
        try:
            prob_swap = float(prob_swap)
        except (ValueError,TypeError):
            raise TypeError('prob_swap must be of type float')
        if prob_swap < 0.0 or prob_swap > 1.0:
            raise ValueError('Probability must be between 0.0 and 1.0.')
        self._prob_swap = prob_swap

    @property
    def max_translate(self):
        return self._max_translate

    @max_translate.setter
    def max_translate(self,max_translate):

        if hasattr(self,'_max_translate'):
            saved = deepcopy(self._max_translate)
        else:
            saved = None
        self._max_translate = max_translate

        if ( not isinstance(self._max_translate,list) or
             len(self._max_translate) != self._n_boxes  ):
            self._max_translate = saved
            raise ValueError('max_translate must be a list with shape '
                '(number of boxes, number of species)')
        for max_translate_box in self._max_translate:
            if ( not isinstance(max_translate_box,list) or
                 len(max_translate_box) != self._n_species  ):
                self._max_translate = saved
                raise ValueError('max_translate must be a list with '
                    'shape (number of boxes, number of species)')
            for max_val in max_translate_box:
                try:
                    max_val = float(max_val)
                except (ValueError,TypeError):
                    self._max_translate = saved
                    raise TypeError('Max translation values must be '
                        'of type float')
                if max_val < 0.0:
                    self._max_translate = saved
                    raise ValueError('Max translation values cannot '
                        'be less than zero')


    def _check_box_species(self):
        """Check validity of ibox and ispecies for current moves object"""

        # Infer the number of boxes
        if (self.ensemble == 'nvt' or
            self.ensemble == 'npt' or
            self.ensemble == 'gcmc'):
            self._n_boxes = 1
        else:
            self._n_boxes = 2

        if not isinstance(ibox,int):
            raise TypeError('The box index ("ibox") must be of type int')
        if ibox < 0:
            raise ValueError('The box index ("ibox") cannot be negative')
        if ibox > 1:
            raise ValueError('A maximum of two boxes is currently '
                    'supported. The box index ("ibox") must be 0 or 1')
        if ibox > 0 and self._n_boxes == 1:
            raise ValueError('The {} ensemble only supports a single '
                    'simulation box. The box index ("ibox") must be 0. '
                    'If you wish to simulate in a multi-box ensemble, '
                    'you must create a new moves object with the desired '
                    'ensemble.'.format(self._ensemble))

        if not isinstance(ispecies,int):
            raise TypeError('The species index ("ispecies") must be '
                'of type int')
        if ispecies < 0:
            raise ValueError('The species index ("ispecies") cannot '
                'be negative')
        if ispecies >= self._n_species:
            raise ValueError('The Moves object was originally created '
                'for {} species. If you wish to simulate '
                'with a different number of species, please '
                'create a new Moves object.')





