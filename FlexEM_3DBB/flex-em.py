#  ============================================================================
#
#  a. Define the type of oprimization: conjugate gradients minimization (CG)
#     or molecular dynamics (MD)
#  b. Edit input parameters: density map parameters and probe structure
#
#  =======================  Maya Topf, 4 Dec 2007 =============================
#  =======================  Latest update: 6/8/09 ============================= 

from modeller import *
import shutil
import sys, os, os.path
import string
import math
from CG import opt_cg
from MD import opt_md

env = environ()
#Things to try
"""
Fit SvF5LDV in 1ake mrc - if fail, problem in rigid or pdb
Fit mdl1_fit in 30_monomer - if fail, problem in mrc

"""

############### INPUT PARAMETERS ##################
optimization = 'MD'                     # type of optimization: CG / MD
input_pdb_file = '3dBB_mod10_trimer_fit.pdb'         #5LDV
# input_pdb_file = 'MOMP_SvF_5ldv_wrongfit.pdb'         #AKE
# input_pdb_file = 'mdl1_fit.pdb'          # input model for optimization

code = 'CTMC'    #5LDV                       # 4 letter code of the structure
# code = '5ake' #AKE
# em_map_file = 'wei_momp.mrc'  #5LDV               # name of EM density map (mrc)
em_map_file = 'wei_momp.mrc' #AKE
# em_map_file = '1akeA_10A.mrc' #AKE
format='MRC'                            # map format: MRC or XPLOR
# apix=2.7                                  # voxel size: A/pixel #5LDV
apix = 2.7 #AKE
# box_size=128                             # size of the density map (cubic) #5LDV
box_size = 128 #AKE
resolution=8.0                         # resolution (A)
# x= -6.494; y=13.381; z=-5.410 #AKE		                # origin of the map (A)
# x=-54*(apix*-1); y=-57*(apix*-1); z=-46*(apix*-1) #5LDV
x= 0; y=0; z=0
path = './'	                        # path to work directory
init_dir = 1                            # number of the initial directory 
num_of_runs = 1                         # number of flex-em CG runs 
num_of_iter = 5                         # number of iterations in each MD run 
          				# (default=4) 
# rigid_filename = 'rigid_1ake.txt'            # rigid bodies file name #AKE
# rigid_filename = 'rigid.txt'            # rigid bodies file name #5LDV
rigid_filename = '3dBB_rigid.txt'            # rigid bodies file name #AKE
############### RUN OPTIMIZATION ##################
path = os.path.abspath(path)

#  CG
# ----
if optimization == 'CG':
    for i in range(init_dir,init_dir+num_of_runs):
        scratch = path + '/' + str(i) + '_cg/'
        os.system("mkdir -p " + scratch)
        os.system("cp " + path + '/'+ input_pdb_file + " " + scratch)
        os.chdir(scratch)
        opt_cg(path, code, str(i), 55*i, em_map_file, input_pdb_file,
                format, apix, box_size, resolution, x, y, z, rigid_filename)

#  MD
# ----
elif optimization == 'MD':
        scratch = path + '/' + str(init_dir) + '_md/'
        os.system("mkdir -p " + scratch)
        os.system("cp " + path + '/'+ input_pdb_file + " " + scratch)
        os.chdir(scratch)
        opt_md(path, code, str(init_dir), 10*init_dir, em_map_file, input_pdb_file,
                format, apix, box_size, resolution, x, y, z, rigid_filename, num_of_iter)
