from rosetta import *
from pyrosetta import *

from rosetta.protocols.loops.loop_closure.kinematic_closure import *

import rosetta.protocols.loops.loop_mover.refine

from sys import exit
from random import randrange
import math

pyrosetta.init()

MAX_KIC_BUILD_ATTEMPTS = 10000

p = pose_from_pdb( "../inputs/result.pdb" )

starting_p = Pose()
# starting_p.assign( protocols )

# scorefxn_low  = protocols.loops.get_cen_scorefxn() #  create_score_function( 'cen_std' )
# scorefxn_high = protocols.loops.get_fa_scorefxn()  #  create_score_function_ws_patch( 'standard', 'score12' )

pymol = PyMOLMover() # If Pymol server is running, centroid stage will display

loops_begin = [1, 18, 41, 49, 104, 115, 172, 186, 200, 216, 262, 274, 331, 340, 355, 369]
loops_end = [8, 33, 41, 94, 105, 164, 174, 190, 208, 252, 264, 323, 332, 347, 359, 375]
loops_cut = [math.ceil(sum(dat)/2) for dat in zip(loops_begin,loops_end)] #Average of the loop start and end

strands_begin = [9, 34, 42, 95, 106, 165, 175, 191, 209, 253, 265, 324, 333, 348, 360, 376]
strands_end = [17, 40, 48, 103, 114, 165, 171, 185, 199, 215, 261, 273, 330, 339, 354, 368, 382]

#Not sure what a Loops() object is good for
# my_loops = protocols.loops.Loops()
# for i in range(len(loops_begin)):
# 	my_loop = protocols.loops.Loop(loop_begin[i], loop_end[i], loop_cut[i])
# 	my_loops.add_loop( my_loop )

ft = FoldTree() #Define the loops so there's a cut point
mm = MoveMap() #only allows the loops to move

mm.set_chi( True )

for i in range(len(loops_begin)):
	ft.add_edge(loops_begin[i], loops_cut[i]-1, -1) #1 to 4
	ft.add_edge(loops_begin[i], loops_end[i], 1) #1 to 8
	ft.add_edge(loops_end[i], loops_cut[i], -1) #8 to 5

	mm.set_bb_true_range(loops_begin[i], loops_end[i] )
	
for i in range(len(strands_begin)):
	ft.add_edge(strands_begin[i], strands_end[i], -1)

print(ft)
print(ft.check_fold_tree())

# if(ft.check_fold_tree()):
# 	p.fold_tree(ft)










