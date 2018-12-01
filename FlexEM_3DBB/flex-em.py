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

############### INPUT PARAMETERS ##################
optimization = 'MD'                                # type of optimization: CG / MD
input_pdb_file = '3dBB_trimer_fit_308.pdb'         # input model for optimization
code = '3DBB'                                      # 4 letter code of the structure          
em_map_file = '308.mrc'                            # name of EM density map (mrc)
format='MRC'                                       # map format: MRC or XPLOR

apix = 2.2                                         # voxel size: A/pixel 
box_size = 128                                     # size of the density map (cubic)
resolution=12.0                                    # resolution (A)
x= 0; y=0; z=0                                     # coordinates of map/model (make sure they're fit)

path = './'	                                       # path to work directory
init_dir = 2                                       # number of the initial directory 
num_of_runs = 1                                    # number of flex-em CG runs 
num_of_iter = 5                                    # number of iterations in each MD run (default=4) 
          				
#MOST OF THE PROBLEMS ARE WITH THIS FILE
rigid_filename = '3dBB_rigid_308.txt'            # rigid bodies file name #AKE








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
