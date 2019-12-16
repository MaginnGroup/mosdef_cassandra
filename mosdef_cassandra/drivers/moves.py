
import parmed

class MoveProbabilities:

    def __init__(self,ensemble,species):

        valid_ensembles = ['nvt','npt','gcmc','gemc','gemc_npt']
        if ensemble not in valid_ensembles:
            raise ValueError('Invalid ensemble "{}" Supported '
                'ensembles include {}'.format(
                    ensemble,valid_ensembles))
        self.ensemble = ensemble

        if not isinstance(species,list):
            raise TypeError('species should be a list of species')
        for spec in species:
            if not isinstance(species,parmed.Structure):
                raise TypeError('species not of type parmed.Structure')

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
            self.prob_translate = 0.25
            self.prob_rotate = 0.25
            self.prob_regrow = 0.30
            self.prob_insert = 0.1
        elif self.ensemble == 'gemc':
            self.prob_translate = 0.35
            self.prob_rotate = 0.35
            self.prob_regrow = 0.30
        elif self.ensemble == 'gemc_npt':
            self.prob_translate = 0.345
            self.prob_rotate = 0.345
            self.prob_regrow = 0.30
            self.prob_vol = 0.01
        else:
            raise ValueError('Uh oh, how did we end up here?')

        # Max translation and rotations specified per-box
        self.sp_max_translate.box1 = []
        self.sp_max_rotate.box1 = []
        # Other options not
        self.sp_regrow_prob = []
        for spec in species:
            if spec is rigid:
                self.max_translate.box1.append(0.0)
                self.max_rotate.box1.append(0.0)
                self.sp_regrow_prob.append(0.0)
            elif spec.natoms == 1:
                self.max_translate.box1.append(2.0)
                self.max_rotate.box1.append(0.0)
                self.sp_regrow_prob.append(0.0)
            else:
                self.max_translate.box1.append(2.0)
                self.max_rotate.box1.append(30.0)
                self.sp_regrow_prob.append(1.0)

        # Add for box two
        if self.ensemble  == 'gemc' or self.ensemble == 'gemc_npt':
            self.sp_max_translate.box1 = []
            self.sp_max_rotate.box1 = []
            for spec in species:
                if spec is rigid:
                    self.max_translate.box2.append(0.0)
                    self.max_rotate.box2.append(0.0)
                elif spec.natoms == 1:
                    self.max_translate.box2.append(2.0)
                    self.max_rotate.box2.append(0.0)
                else:
                    self.max_translate.box2.append(2.0)
                    self.max_rotate.box2.append(30.0)


