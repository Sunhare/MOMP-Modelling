from pyrosetta import *
import pyrosetta.rosetta as rosetta
from pyrosetta.rosetta import *

from pyrosetta import init, Pose, pose_from_file, Vector1, create_score_function, PyJobDistributor, MoveMap, SwitchResidueTypeSetMover

# from pyrosetta.rosetta import core, protocols
from rosetta.protocols.rigid import (
    RigidBodyRandomizeMover, RigidBodyPerturbMover, RigidBodySpinMover,
    partner_upstream, partner_downstream
)
from rosetta.protocols import *
from pyrosetta.toolbox import *

pyrosetta.init()

import math
import random

#Returns true or false based on probability given by temperature
def MetropolisCriterion(energy, temp):
	#P is the probability to accept a move
	kT = temp; # Depends on temperature
	if energy >= 0: #If unfavorable
		P = math.e ** (-1*energy/kT)
	else: #If favorable
		P = 1
	return True if random.random() < P else False 

#Simulation Constants
iters = 10000


#Make initial pose
# pose = pose_from_sequence('A'*10) #Initial pose
pose = pose_from_sequence('CMTTWNPTISGSGIC')
n = pose.total_residue()

set_constraints = protocols.constraint_movers.ConstraintSetMover()
set_constraints.constraint_file("hevnp_vd4_inputs/hevnp_vd4_constraints.cst")
set_constraints.apply(pose)

ft = FoldTree()                                                             
ft.add_edge(1, 10, -1)                                                       
# ft.add_edge(5, 10, -1)                                                      
# pose.fold_tree(ft) 

#Pymol visualization

pymol = PyMOLMover()
pymol.keep_history(True)
pymol.apply(pose)


#Score function
mc_scorefxn = create_score_function("ref2015_cst")
mc_scorefxn.set_weight(core.scoring.ScoreType.atom_pair_constraint, 1000.0)


#Initial Simulation Poses
curr_pose = Pose()
curr_pose.assign(pose)
set_constraints.apply(curr_pose)

low_pose = Pose()
prev_pose = Pose()

#Main Loop
currEnergy = mc_scorefxn(curr_pose)
lowEnergy = currEnergy

temp_factor = 10

sd = 25 # standard deviation for angle


for i in range(1,iters+1):
	print("I: {0} Current E:\t{1}\tLowest E:\t{2}".format(i, currEnergy, lowEnergy))
	prev_pose.assign(curr_pose) #Save previous pose
	
	#Make a Monte Carlo move
	res = random.randint(1,n) #Choose a residue


	phi_psi_coinflip = random.randint(1,2) #Choose phi or psi to change
	if phi_psi_coinflip == 1:
		curr_angle = curr_pose.phi(res)
		new_angle = random.gauss(curr_angle, sd)
		curr_pose.set_phi(res, new_angle)
	else:
		curr_angle = curr_pose.psi(res)
		new_angle = random.gauss(curr_angle, sd)
		curr_pose.set_psi(res, new_angle)

	newEnergy = mc_scorefxn(curr_pose)
	energyDiff = newEnergy - currEnergy

	#Evaluate
	#If the move is accepted

	if(MetropolisCriterion(energyDiff, temp_factor) == True): 
		if(newEnergy < lowEnergy): #if lowest energy we've seen
			lowEnergy = newEnergy
			low_pose.assign(curr_pose)

		currEnergy = newEnergy
		pymol.apply(curr_pose)
		pymol.send_energy(pose, "atom_pair_constraint")
	#If rejected
	else:
		curr_pose.assign(prev_pose) #Revert

	
	if iters <= 10:
		mc_scorefxn.show(pose)

	# Simulated Annealing?
	if math.floor(i%(iters/5)) == 0:
		temp_factor -= 2
		curr_pose.assign(low_pose) #Don't waste too much time
		# sd = 5*temp_factor # standard deviation for angle
	if temp_factor < 1:
		temp_factor = 1
	

low_pose.pdb_info().name("low_pose")
pymol.apply(low_pose)



