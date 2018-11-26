#  ============================================================================
#
#  Conjugate gradients minimization with cross-correlation, non-bonded
#  interaction and stereochemical restraints.
#
#  =======================  Maya Topf, 4 Dec 2007 =============================
#  =======================  Latest update: 6/8/09 =============================

from modeller import *
from modeller.automodel import refine
from modeller.scripts import complete_pdb
from modeller.optimizers import conjugate_gradients,actions
from modeller import schedule
from rigid import load_rigid
from random import *
import shutil
import sys, os, os.path
import string

class opt_cg:

    def __init__(self, path, code, run_num, rand, em_map_file, input_pdb_file, format, apix, box_size, res, x, y, z, rigid_filename):

        env = environ(rand_seed=-rand)
        env.libs.topology.read(file='$(LIB)/top_heav.lib')
        env.libs.parameters.read(file='$(LIB)/par.lib')
        env.io.hetatm = True

        log.verbose()

# reading the density
        den = density(env, file = path + '/' + em_map_file,
                     em_density_format=format, em_map_size=box_size,
                     density_type='GAUSS', voxel_size=apix,
                     resolution=res,px=x,py=y,pz=z)

        env.edat.density = den
        env.edat.dynamic_sphere = True

# read pdb file
        aln = alignment(env)
        mdl2 = model(env, file=input_pdb_file)
        aln.append_model(mdl2, align_codes=input_pdb_file,
                        atom_files=input_pdb_file)
        mdl = model(env, file=input_pdb_file)
        aln.append_model(mdl, align_codes=code, atom_files=code)
        aln.align(gap_penalties_1d=(-600, -400))
        mdl.clear_topology()
        mdl.generate_topology(aln[input_pdb_file])
        mdl.transfer_xyz(aln)
        mdl.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')

# Remove chain names
        for c in mdl.chains:
            c.name = ' '
        mdl.write(file=code+'_ini.pdb')

# Build restraints
        sel_all = selection(mdl)
        mdl.restraints.make(sel_all, restraint_type='stereo',
                         spline_on_site=False)

        mdl.restraints.make(sel_all, aln=aln, restraint_type='phi-psi_binormal',
                   spline_on_site=True,
                   residue_span_range=(0, 99999))

# define rigid and flexible regions
        sel_rigid=[]
        rand_rigid=[]
        load_rigid(path,mdl,sel_rigid,rand_rigid,rigid_filename)
        for n in sel_rigid:
            mdl.restraints.rigid_bodies.append(rigid_body(n))

        list = []
        for n in sel_all:
            for m in sel_rigid:
                if n in m:
                    a=1
                    break
                else:
                    a=0
            if a==0: list.append(n)
        sel_flex=selection(list)

        print sel_all
        print sel_rigid
        print sel_flex

# randomization of the rigid bodies
        if len(sel_flex) > 0:
           sel_flex.randomize_xyz(deviation=1.0)

        seed(rand)
        max_trans = 10
        max_ang = 30
        for num in sel_rigid:
            if rand_rigid[sel_rigid.index(num)]==1:

                rand_rot_x = uniform(-1,1)
                rand_rot_y = uniform(-1,1)
                rand_rot_z = uniform(-1,1)
                rand_ang = uniform(0,max_ang)
                rand_trans_x = uniform(-max_trans,max_trans)
                rand_trans_y = uniform(-max_trans,max_trans)
                rand_trans_z = uniform(-max_trans,max_trans)

                num.rotate_mass_center([rand_rot_x,rand_rot_y,rand_rot_z],rand_ang)

                num.translate([rand_trans_x,rand_trans_y,rand_trans_z])

        mdl.write(file=code+'_rand.pdb')

# CG minimization
        print "conjugate_gradients"
        print "directory num %s" %run_num
        CG = conjugate_gradients()

        sched = schedule.schedule(2,
            [ schedule.step(CG, 2, physical.values(default=0.01)),
              schedule.step(CG, 5, physical.values(default=0.1)),
              schedule.step(CG, 10, physical.values(default=0.2)),
              schedule.step(CG, 50, physical.values(default=0.5, soft_sphere=0.1, em_density=100)),
              schedule.step(CG, 9999, physical.values(default=0.5, soft_sphere=1.0, em_density=10000)),
              schedule.step(CG, 9999, physical.values(default=1.0, soft_sphere=1.0, em_density=10000))])

        trc_step=5
        trc_file = open('CG'+run_num+'.trc', "a")

        i=0
        for step in sched:
            i+=1
            step.optimize(sel_all, output='REPORT', max_iterations=200,
                          actions=actions.trace(trc_step,trc_file))
            pdb_file='cg'+run_num+'_'+str(i)+'.pdb'
            sel_all.write(file=pdb_file)

        trc_file.close()

# print final energy, CC and final model
        print "final energy all "
        scal = physical.values(default=1.0, em_density=10000)
        eval = sel_all.energy(schedule_scale=scal)

        print "final cc "
        scal = physical.values(default=0.0, em_density=1.0)
        eval = sel_all.energy(schedule_scale=scal)

        sel_all.write(file='final'+run_num+'_cg.pdb')

        os.system("rm -f *.MRC")
