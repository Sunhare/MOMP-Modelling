from pyrosetta import *
from pyrosetta.rosetta import core
from pyrosetta.toolbox import pose_from_rcsb
from pyrosetta.rosetta.core.scoring import CA_rmsd
from pymol import cmd

import random
import math
import os
import numpy as np

init()

def MetropolisCriterion(energy, temp):
	#P is the probability to accept a move
	kT = temp; # Depends on temperature
	if energy >= 0: #If unfavorable
		P = math.e ** (-1*energy/kT)
	else: #If favorable
		P = 1
	return True if random.random() < P else False 

def makeRama(pose):
	n = pose.total_residue()

	
	filename = str(pose.pdb_info().name())+"_rama.txt"

	if os.path.exists(filename):
	    append_write = 'a' # append if already exists
	else:
	    append_write = 'w' # make a new file if not

	ramachandran = open(filename, append_write)

	# Phi \t Psi \t Residue

	for i in range(1, n+1):
		curr_phi = pose.phi(i)%(180 * np.sign(pose.phi(i)))
		curr_psi = pose.psi(i)%(180 * np.sign(pose.psi(i)))

		rama_point = str(curr_phi)+ "\t" +str(curr_psi) + "\n"
		ramachandran.write(rama_point)

	return

def new_select(self): 
	import pymol 
	pymol.cmd.select("all")


def new_align(self, mobile, target):
	from pymol import cmd
	cmd.align(mobile, target)


def EvsRMSD(pose, native, energy):
	filename = str(pose.pdb_info().name())+"_E_RMSD.txt"

	if os.path.exists(filename):
	    append_write = 'a' # append if already exists
	else:
	    append_write = 'w' # make a new file if not

	E_RMSD = open(filename, append_write)

	RMSD = CA_rmsd(pose, native)

	to_write = str(RMSD) + "\t" + str(energy) + "\n"
	E_RMSD.write(to_write)
	return

def myMC(sequence, iters, visualize, temp_factor, native=None, calculations=False):
	if native == None:
		pose = pose_from_sequence(sequence)
	else:
		pose = pose_from_sequence(native.sequence())
	n = pose.total_residue()
	sd = 25

	if(visualize == True):
		PyMover = PyMOLMover()
		PyMover.keep_history(True)
		PyMover.apply(pose)

	MC_SF = create_score_function("ref2015_cst")
	# MC_SF = ScoreFunction()
	# MC_SF.set_weight(core.scoring.ScoreType.fa_atr, 1.0)
	# MC_SF.set_weight(core.scoring.ScoreType.fa_rep, 1.0)
	# MC_SF.set_weight(core.scoring.ScoreType.hbond_sr_bb, 1.0)
	# MC_SF.set_weight(core.scoring.ScoreType.hbond_lr_bb, 1.0)
	# MC_SF.set_weight(core.scoring.ScoreType.hbond_bb_sc, 1.0)
	# MC_SF.set_weight(core.scoring.ScoreType.hbond_sc, 1.0)


	# pose = Pose()
	# pose.assign(pose)
	low_pose = Pose()
	prev_pose = Pose()

	currEnergy = MC_SF(pose)
	lowEnergy = currEnergy

	for i in range(1,iters+1):
		# print("I: {0} Current E:\t{1}\tLowest E:\t{2}".format(i, currEnergy, lowEnergy))
		prev_pose.assign(pose) #Save previous pose
	
		#Make a Monte Carlo move
		res = random.randint(1,n) #Choose a residue
		phi_psi_coinflip = random.randint(1,2) #Choose phi or psi to change
		if phi_psi_coinflip == 1:
			curr_angle = pose.phi(res)
			new_angle = random.gauss(curr_angle, sd)
			pose.set_phi(res, new_angle)
		else:
			curr_angle = pose.psi(res)
			new_angle = random.gauss(curr_angle, sd)
			pose.set_psi(res, new_angle)

		newEnergy = MC_SF(pose)
		energyDiff = newEnergy - currEnergy

		if(MetropolisCriterion(energyDiff, temp_factor) == True): 
			if(newEnergy < lowEnergy): #if lowest energy we've seen
				lowEnergy = newEnergy
				low_pose.assign(pose)

			currEnergy = newEnergy
			if(visualize == True):
				PyMover.apply(pose)
		#If rejected
		else:
			pose.assign(prev_pose) #Revert

	#Post code processing
	print("Iters: {0} Current E:\t{1}\tLowest E:\t{2}".format(iters, currEnergy, lowEnergy))
	if(calculations):
		EvsRMSD(low_pose, native, lowEnergy)
		makeRama(low_pose)

# seq = 'AAAAAAAAAA'
# myMC(seq, 5000, True, 1)


