from pyrosetta import *
init()

polyA = "A"*20

polyA_extended = pose_from_sequence(polyA, "fa_standard")

pymol_extended = PyMOLMover()
pymol_extended.apply(polyA_extended)

pymol_helix = PyMOLMover()
polyA_helix = polyA_extended
for i in range(1,polyA_helix.total_residue()+1):
	polyA_helix.set_phi(i, -64)
	polyA_helix.set_psi(i, -47)

pymol_helix.apply(polyA_helix)
