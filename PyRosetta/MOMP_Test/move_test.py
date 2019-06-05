from rosetta import *
from pyrosetta import *
from pyrosetta.toolbox import *
from pyrosetta.teaching import *
import random

init()

momp = pose_from_pdb("MOMP_subunit_VDs.pdb")

movemap = MoveMap()

movemap.set_bb_true_range(5,10)

movemap.show(10)

pmm = PyMOLMover()
pmm.apply(momp)

#Small Movers are similar to pertubring the backbone for monte carlo simulations
#They move phi(i) and psi(i)
# >>movemap = MoveMap()
# >>movemap.set_bb(True)
# >>n_moves = 5
# >>kT = 1.0
# >>smallmover = SmallMover(movemap, n_moves, kT)
# >>smallmover.angle_max(10)
# >>smallmover.angle_max(‘E’, 5)
# >>smallmover.angle_max(‘H’, 10)
# >>smallmover.angle_max(‘L’, 20)
# >>smallmover.apply(pose)

#Shear movers are similar 
#They move phi(i) and psi(i-1) to prevent large changes in the whole protein