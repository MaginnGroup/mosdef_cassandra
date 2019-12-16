
import mbuild
import parmed

class System:
    def __init__(self,boxes,species_topologies,species_in_boxes,species_to_add):
        """A class to contain the system to simulate in Cassandra

        boxes : list
            one element per box. Each element should be a
            mbuild.Compound with box information (and maybe atoms)
        species_topologies : list
            list of pmd.Structures, with one species per element
        species_in_boxes: list
            one element per box. Each element is a list of length
            n_species, specifying the number of each species that
            are currently in each box
        species_to_add : list
            one element per box. Each element is a list of length
            n_species, specifying the number of each species that
            should be added to each box
        """
        if not isinstance(boxes,list):
            raise TypeError('"boxes" should be a list. See'
                    'help(mosdef_Cassandra.System) for details.')
        if not isinstance(species_topologies,list):
            raise TypeError('"species_topologies" should be a list.'
                    'See help(mosdef_Cassandra.System) for details.')
        if not isinstance(species_to_add,list):
            raise TypeError('"species_to_add" should be a list.'
                    'See help(mosdef_Cassandra.System) for details.')
        if len(boxes) != len(species_to_add):
            raise ValueError('The number of boxes inferred from the'
                    'length of "species_to_add" must match the'
                    'number of boxes inferred from the length of'
                    '"boxes"')
        if len(boxes) != len(species_in_boxes):
            raise ValueError('The number of boxes inferred from the'
                    'length of "species_in_boxes" must match the'
                    'number of boxes inferred from the length of'
                    '"boxes"')
        for box in boxes:
            if not isinstance(box,mbuild.Compound):
                raise TypeError('Each box should be an'
                        'mbuild.Compound object')
        for topology in species_topologies:
            if not isinstance(topology,parmed.Structure):
                raise TypeError('Each species should be a'
                        'parmed.Structure')
        n_species = len(species_topologies)
        for add_to_box in species_to_add:
            if len(add_to_box) != n_species:
                raise ValueError('The number of each species to'
                        'be added to each box must be specified'
                        '(even if it is = 0)')

        self.boxes = copy.deepcopy(boxes)
        self.species_topologies = copy.deepcopy(species_topologies)
        self.species_in_boxes = copy.deepcopy(species_in_boxes)
        self.species_to_add = copy.deepcopy(species_to_add)


