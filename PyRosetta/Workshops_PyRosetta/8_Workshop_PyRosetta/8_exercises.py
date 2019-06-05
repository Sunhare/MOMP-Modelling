from rosetta import *
from pyrosetta import *
from pyrosetta.teaching import *

from rosetta.protocols.loops.loop_closure.ccd import *
from rosetta.protocols.loops.loop_mover.refine import *
from rosetta.protocols.loops.loop_mover.perturb import *
from rosetta.protocols.loops.loop_closure.kinematic_closure import *

init()

pose = pose_from_file("./test_in.pdb")

# Fold Tree
ft = FoldTree()
ft.add_edge(1, 13, -1)
ft.add_edge(13, 19, -1)
ft.add_edge(13, 26, 1)
ft.add_edge(26, 20, -1)
ft.add_edge(26, 116, -1)

print( ft )
ft.check_fold_tree()

pose.fold_tree(ft)

for res in (10, 13, 16, 23, 26, 30):
    pose.set_phi(res, 180)
    pose.dump_pdb("loop" + str(res) + ".pdb")

pose.pdb_info().name("start")

observer =AddPyMOLObserver(pose,True) 

# pmm = PyMOLMover()
# pmm.apply(pose)
# no longer supported: pmm.send_foldtree(pose)
# no longer supported: pmm.view_foldtree_diagram(pose, ft)

ft.clear()
ft.simple_tree(116)
ft.new_jump(76, 85, 80)

# Cyclic Coordination Descent (CCD) Loop Closure
movemap = MoveMap()
movemap.set_bb(True)
movemap.set_chi(True)

#Loop(start, end, cut), idk why I couldn't find this simple definition anywhere
loop1 = pyrosetta.rosetta.protocols.loops.Loop(15, 24, 19)
add_single_cutpoint_variant(pose, loop1)
ccd = pyrosetta.rosetta.protocols.loops.loop_closure.ccd.CCDLoopClosureMover(loop1, movemap)
ccd.apply(pose)
kic_mover = KinematicMover()
kic_mover.set_pivots(16, 20, 24)
kic_mover.apply(pose)

# sw = SwitchResidueTypeSetMover("centroid")
# sw.apply(pose)
# kic_perturb = LoopMover_Perturb_KIC(loop1)
# kic_perturb.apply(pose)
# sw = SwitchResidueTypeSetMover("fa_standard")
# sw.apply(pose)
# kic_refine = LoopMover_Refine_KIC(loop1)
# kic_refine.apply(pose)






