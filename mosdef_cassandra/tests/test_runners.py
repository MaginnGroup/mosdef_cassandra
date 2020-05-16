import pytest

import mosdef_cassandra as mc
from mosdef_cassandra.tests.base_test import BaseTest


class TestRunners(BaseTest):
    def test_mismatch_boxes(self, methane_oplsaa, box):
        with pytest.raises(ValueError, match=r"requires 1 simulation"):
            system = mc.System(
                [box, box], [methane_oplsaa], mols_to_add=[[10], [0]]
            )
            moves = mc.MoveSet("nvt", [methane_oplsaa])
            mc.run(system, moves, 300.0, "equilibration", 500)
        with pytest.raises(ValueError, match=r"requires 2 simulation"):
            system = mc.System([box], [methane_oplsaa], mols_to_add=[[10]])
            moves = mc.MoveSet("gemc", [methane_oplsaa])
            mc.run(system, moves, 300.0, "equilibration", 500)

    def test_corrupt_boxes(self, methane_oplsaa, box):
        with pytest.raises(TypeError, match=r"corrupted"):
            system = mc.System([box], [methane_oplsaa], mols_to_add=[[10]])
            moves = mc.MoveSet("nvt", [methane_oplsaa])

            system.boxes[0] = 1

            mc.run(system, moves, 300.0, "equilibration", 500)

    def test_corrupt_topologies(self, methane_oplsaa, box):
        with pytest.raises(TypeError, match=r"corrupted"):
            system = mc.System([box], [methane_oplsaa], mols_to_add=[[10]])
            moves = mc.MoveSet("nvt", [methane_oplsaa])

            system._species_topologies = methane_oplsaa
            mc.run(system, moves, 300.0, "equilibration", 500)

        with pytest.raises(TypeError, match=r"corrupted"):
            system = mc.System([box], [methane_oplsaa], mols_to_add=[[10]])
            moves = mc.MoveSet("nvt", [methane_oplsaa])

            system.species_topologies[0] = 1
            mc.run(system, moves, 300.0, "equilibration", 500)

    def test_corrupt_natoms(self, methane_oplsaa, box):
        with pytest.raises(ValueError, match=r"corrupted"):
            system = mc.System([box], [methane_oplsaa], mols_to_add=[[10]])
            moves = mc.MoveSet("nvt", [methane_oplsaa])

            system.mols_in_boxes[0][0] = 10
            mc.run(system, moves, 300.0, "equilibration", 500)
