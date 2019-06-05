from pyrosetta import *
from pyrosetta.rosetta import core
from pyrosetta.toolbox import *
from pyrosetta.teaching import PyMOLMover
from pyrosetta.rosetta.protocols.simple_moves import *
from pyrosetta.rosetta.protocols.moves import *
from pyrosetta.rosetta.protocols.minimization_packing import *


# Initialize
init() 

pose = pose_from_pdb("2ARC_chainA.pdb")
pose.pdb_info().name("2ARC")
start = Pose()
start.assign(pose)

ft = FoldTree()
ft.add_edge(161, 1, -1)

pose.fold_tree(ft)


start.pdb_info().name("start")
pmm = PyMOLMover()
pmm.keep_history(True)

# pmm = PyMOLMover()
# pmm.apply(pose)
observer = AddPyMOLObserver(pose, True)

scorefxn = create_score_function("ref2015_cst")


kT = 1
mc = MonteCarlo(pose, scorefxn, kT)
mc.boltzmann(pose)

movemap = MoveMap()
movemap.set_bb(False)
movemap.set_bb_true_range(1,16)

n_moves = 10

small_mover = SmallMover(movemap, kT, n_moves)
# small_trial_mc = TrialMover(small_mover, mc)
shear_mover = ShearMover(movemap, kT, n_moves)

# reapeat_small_trial_mc = RepeatMover(small_trial_mc, 5)

minmover = MinMover()
minmover.score_function(scorefxn)

seq_mover = SequenceMover()
seq_mover.add_mover(small_mover)
seq_mover.add_mover(minmover)
seq_mover.add_mover(shear_mover)
seq_mover.add_mover(minmover)

single_mc_move = TrialMover(seq_mover, mc) 


pmm.apply(pose)
pmm.apply(start)
pmm.send_movemap(pose, movemap)

relax = pyrosetta.rosetta.protocols.relax.FastRelax()

relax.set_scorefxn(scorefxn)
relax.set_movemap(movemap)

relax.apply(pose)

# ang_max = 25
# small_mover.angle_max("H", ang_max)
# small_mover.angle_max("E", ang_max)
# small_mover.angle_max("L", ang_max)

# shear_mover.angle_max("H", ang_max)
# shear_mover.angle_max("E", ang_max)
# shear_mover.angle_max("L", ang_max)

# ###### Main Loop ######
# for j in range(5):

# 	for i in range(200):

# 		single_mc_move.apply(pose)
# 		if i%3 == 0:
# 			pmm.apply(pose)

# 	ang_max = ang_max/2

# 	small_mover.angle_max("H", ang_max)
# 	small_mover.angle_max("E", ang_max)
# 	small_mover.angle_max("L", ang_max)

# 	shear_mover.angle_max("H", ang_max)
# 	shear_mover.angle_max("E", ang_max)
# 	shear_mover.angle_max("L", ang_max)
	





