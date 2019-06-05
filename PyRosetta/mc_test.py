from rosetta import *
from pyrosetta import *
from pyrosetta.toolbox import *
from pyrosetta.teaching import *
import random

init()

#Set up input poly alanine
pose=Pose()
polyA = "A"*30
pose = pose_from_sequence(polyA, "fa_standard")
pose.pdb_info().name("polyA")

starting_pose= Pose()
starting_pose.assign(pose)

def make_perfect_helix(in_pose):
	n = in_pose.total_residue()
	for i in range(1,n-1):
		in_pose.set_phi(i, -64)
		in_pose.set_psi(i,-41)
	return in_pose

# native_pose = Pose()
# native_pose.assign(make_perfect_helix(pose))


#set up score function
scorefxn = ScoreFunction()
scorefxn.set_weight(hbond_sr_bb,1.0)
scorefxn.set_weight(vdw, 1.0)


#set up mover
def perturb_bb(pose):
	resnum = random.randint(1, pose.total_residue())
	pose.set_phi(resnum, pose.phi(resnum)-25+random.random()*50)
	pose.set_psi(resnum, pose.psi(resnum)-25+random.random()*50)
	return pose

#set up visualization
pmm = PyMOLMover()
pmm.keep_history(True)

kT = 1
mc = MonteCarlo(pose, scorefxn, kT)
def mc_fold(mc, pose, visualization=False):
	#set up MonteCarlo object

	global pmm #used for visualization
	mc.reset(pose)
	for i in range(1,2000): 
		perturb_bb(pose)
		mc.boltzmann(pose)
		if (i%1000 == 0):
			 mc.recover_low(pose)
			 if visualization == True:
				 pmm.apply(pose)
	#output lowest-energy structure
	# mc.recover_low(pose)

	return pose

#PyMOL Align different frames 
# align  "polyA", "perfectH", mobile_state=99, target_state=2
# jd = PyJobDistributor("polyA_output", 3, scorefxn)

# perfectH = pose_from_pdb("perfectH.pdb")
jd = PyJobDistributor("helix_fold", 3, scorefxn)
# jd.native_pose = perfectH

index = 0
# base_name = pose.pdb_info().name()
while (jd.job_complete == False): 
	pose.assign(starting_pose) 
	
	# new_name = base_name+"_"+str(index)
	# pose.pdb_info().name(new_name)
	# index+1
	mc_fold(mc, pose)	
	# jd.additional_decoy_info("loop_RMSD: " + str(loop_rmsd) + " res49A_res20B" + str(res_dist))
	pmm.apply(pose)
	# print(pose.pdb_info().name())
	jd.output_decoy(pose)

# for temp in range(4,0,-1): #Go from kT = 4 to kT = 1
# 	print("Folding for temp="+str(temp))
# 	mc_fold(pose, scorefxn, temp, True)


# dump_pdb(p, “mc_final.pdb”)




