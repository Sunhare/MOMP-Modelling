from pyrosetta import *
from pyrosetta.toolbox import cleanATOM
from pyrosetta.toolbox import pose_from_rcsb

#Fetch protein using command line
# wget http://www.rcsb.org/pdb/files/1YY8.pdb.gz

init()

cleanATOM("1YY8.pdb")

pose = pose_from_pdb("1YY8.clean.pdb")
# pose = pose_from_rcsb("1YY8")

print(pose)
print(pose.sequence())
print("Protein has", pose.total_residue(), "residues.")
print(pose.residue(500).name())

print(pose.pdb_info().chain(500))
print(pose.pdb_info().number(500))

print(pose.pdb_info().pdb2pose('A', 100))    
print(pose.pdb_info().pose2pdb(25))

# get_secstruct(pose)
#Why.
set_ss_from_phipsi(pose) 
pose.secstruct()
