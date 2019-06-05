from pyrosetta import *
from pyrosetta.rosetta import core
from pyrosetta.toolbox import *
from pyrosetta.teaching import PyMOLMover
from pyrosetta.rosetta.protocols.simple_moves import *
from pyrosetta.rosetta.protocols.moves import *
from pyrosetta.rosetta.protocols.minimization_packing import *


# Initialize
init() 
start = pose_from_pdb("original_subunit.pdb")
test = Pose()
test.assign(start)

# Set PyMOL names and apply
start.pdb_info().name("start")
test.pdb_info().name("test")

pmm = PyMOLMover()
pmm.keep_history(True)
pmm.apply(start)
# pmm.apply(test)

#Score Function
E5_SF = create_score_function("ref2015_cst")


# Metropolis criterion
# Mover set up
kT = 1.0
n_moves = 1
# movemap = MoveMap()
# movemap.set_bb(True)
# movemap.set_bb(False)
# movemap.set_bb(50, True)
# movemap.set_bb(51, True)


# small_mover = SmallMover(movemap, kT, n_moves)


small_mover.angle_max("H", 25)
small_mover.angle_max("E", 25)
small_mover.angle_max("L", 25)

# Applying movers
# small_mover.apply(test)
# pmm.apply(test)


#Min Mover

test2 = Pose()


test2.assign(start)
test2.pdb_info().name("test2")
observer = AddPyMOLObserver(test2, True)



mm4060 = MoveMap()
mm4060.set_bb_true_range(40,60)

shear_mover = ShearMover(mm4060, kT, n_moves)

shear_mover.apply(test2)

min_mover = MinMover()
min_mover.score_function(E5_SF)

min_mover.apply(test2)






