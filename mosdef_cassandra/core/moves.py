
import parmed

class MoveProbabilities:

    def __init__(self,ensemble,species_topologies):

        valid_ensembles = ['nvt','npt','gcmc','gemc','gemc_npt']
        if ensemble not in valid_ensembles:
            raise ValueError('Invalid ensemble "{}" Supported '
                'ensembles include {}'.format(
                    ensemble,valid_ensembles))
        self.ensemble = ensemble

        if not isinstance(species_topologies,list):
            raise TypeError('species should be a list of species')
        for species in species_topologies:
            if not isinstance(species,parmed.Structure):
                raise TypeError('species not of type parmed.Structure')

        # Start by setting all move probabilities to zero
        self.prob_translate = 0.0
        self.prob_rotate = 0.0
        self.prob_angle = 0.0
        self.prob_dihedral = 0.0
        self.prob_regrow = 0.0
        self.prob_vol = 0.0
        self.prob_insert = 0.0
        self.prob_swap = 0.0

        # Define default probabilities for each ensemble
        if self.ensemble == 'nvt':
            self.prob_translate = 0.35
            self.prob_rotate = 0.35
            self.prob_regrow = 0.30
        elif self.ensemble == 'npt':
            self.prob_translate = 0.345
            self.prob_rotate = 0.345
            self.prob_regrow = 0.30
            self.prob_vol = 0.01
        elif self.ensemble == 'gcmc':
            self.prob_translate = 0.30
            self.prob_rotate = 0.30
            self.prob_regrow = 0.30
            self.prob_insert = 0.1
        elif self.ensemble == 'gemc':
            self.prob_translate = 0.30
            self.prob_rotate = 0.30
            self.prob_regrow = 0.30
            self.prob_swap = 0.1
        elif self.ensemble == 'gemc_npt':
            self.prob_translate = 0.295
            self.prob_rotate = 0.295
            self.prob_regrow = 0.30
            self.prob_swap = 0.1
            self.prob_vol = 0.01
        else:
            raise ValueError('Uh oh, how did we end up here?')

        # Max translation and rotations specified per-box
        self.max_translate = []
        self.max_rotate = []
        self.max_translate.append([])
        self.max_rotate.append([])
        self.prob_swap_from_box = [1.0]
        # Other options not
        self.insertable = []
        self.sp_prob_regrow = []

        # Here we handle species-wise selections
        # including per-box selections for box1
        for species in species_topologies:
            if len(species.atoms) == 1:
                self.max_translate[0].append(2.0)
                self.max_rotate[0].append(0.0)
                self.sp_prob_regrow.append(0.0)
                self.sp_insertable.append(True)
            elif len(species.bonds) == 0:
                print('Treating {} as a non-insertable rigid species '
                        'since it has no bonds'.format(species))
                self.max_translate[0].append(0.0)
                self.max_rotate[0].append(0.0)
                self.sp_prob_regrow.append(0.0)
                self.sp_insertable.append(False)
            else:
                self.max_translate[0].append(2.0)
                self.max_rotate[0].append(30.0)
                self.sp_prob_regrow.append(1.0)
                self.sp_insertable.append(True)

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
        if self.ensemble == 'npt' or self.ensemble == 'gemc_npt':
            self.max_volume = []
            self.max_volume.append(300.)
            if self.ensemble == 'gemc_npt':
                self.max_volume.append(3000.)

        # Correct species_prob_regrow
        sp_regrowth_prob = 1.0/sum(self.sp_prob_regrow)
        for i,prob in enumerate(self.sp_prob_regrow):
            if prob > 0.0:
                self.sp_prob_regrow[i] = sp_regrowth_prob
        
        # Correct species_prob_swap
        sp_swap_prob = 1.0/sum(self.sp_swap_prob)
        for i,insertable in enumerate(self.sp_insertable):
            if insertable:
                self.sp_swap_prob[i] = sp_swap_prob

        # If all species have no prob regrowth, set prob_regrow to
        # zero and redistribute prob to translate/rotate
        if sum(self.sp_prob_regrow) == 0.0:
            self.prob_translate += self.prob_regrow/2.
            self.prob_rotate += self.prob_regrow/2.
            self.prob_regrow = 0.0

        # If all species are not rotatable change prob rotation
        # move to zero. Redistribute prob to translate
        if sum(self.max_rotate[0]) == 0.0:
            self.prob_translate += self.prob_rotate
            prob_rotate = 0.0
        

