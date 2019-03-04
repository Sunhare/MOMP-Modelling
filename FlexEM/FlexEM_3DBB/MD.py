#  ============================================================================
#
#  Simulated annealing molecular dynamics optimization with cross-correlation,
#  non-bonded interaction, and stereochemical restraints.
#
#  =======================  Maya Topf, 4 Dec 2007 =============================
#  =======================  Latest update: 6/8/09 =============================

from modeller import *
from modeller.automodel import refine
from modeller.scripts import complete_pdb
from modeller.optimizers import conjugate_gradients,molecular_dynamics,actions
from modeller import schedule
from rigid import load_rigid
from random import *
import shutil
import sys, os, os.path
import string

class opt_md:

    def __init__(self, path, code, run_num, rand, em_map_file, input_pdb_file, format, apix, box_size, res, x, y, z, rigid_filename, num_of_iter):

        env = environ(rand_seed=-rand)
        env.libs.topology.read(file='$(LIB)/top_heav.lib')
        env.libs.parameters.read(file='$(LIB)/par.lib')
        env.io.hetatm = True

        log.verbose()

        den = density(env, file=path + '/' + em_map_file,
                     em_density_format=format, em_map_size=box_size,
                     density_type='GAUSS', voxel_size=apix,
                     resolution=res,px=x,py=y,pz=z)

        env.edat.density = den
        env.edat.dynamic_sphere = True

# read pdb
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

# renumber chains
        for c in mdl.chains:
            c.name = ' '
        mdl.write(file=code+'_ini.pdb')

# Build restraints
        sel_all = selection(mdl)
        mdl.restraints.make(sel_all, restraint_type='stereo', spline_on_site=False)

        mdl.restraints.make(sel_all, aln=aln, restraint_type='phi-psi_binormal',
                   spline_on_site=True,
                   residue_span_range=(0, 99999))

# define rigid bodies
        sel_rigid=[]
        rand_rigid=[]
        load_rigid(path,mdl,sel_rigid,rand_rigid,rigid_filename)
        for n in sel_rigid:
            mdl.restraints.rigid_bodies.append(rigid_body(n))

        mdl.restraints.write(file=code+'.rsr')

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

# print selections
        print sel_all
        print sel_rigid
        print sel_flex

# start simulated annealing molecular dynamics
        print "MD annealing"
        scal = physical.values(default=1.0, em_density=10000)

        cap=0.39
        timestep=5.0
        icount=0

        MD =  molecular_dynamics(cap_atom_shift=cap, md_time_step=timestep,
                                    md_return='FINAL', output='REPORT',
                                    schedule_scale=scal)

        trc_file = open('MD'+run_num+'.trc', "a")

        a = 1
        b = num_of_iter + 1
        while  a < b:

# heating the system
            equil_its=100
            equil_equil=20
            equil_temps=(150.0, 250.0, 500.0, 1000.0)
            trc_step=5
            init_vel = True

            for (its, equil, temps) in [(equil_its, equil_equil, equil_temps)]:
                for temp in temps:
                    MD.optimize(sel_all, max_iterations=its,
                       temperature=temp, init_velocities=init_vel, equilibrate=equil,
                       actions=[actions.trace(trc_step,trc_file)])

                    scal = physical.values(default=0.0, em_density=1.0)
                    (molpdf, terms) = sel_all.energy(schedule_scale=scal)
                    print "iteration number= %s  step= %d %d  temp= %d  EM score= %.3f" %(a,icount,its,int(temp),-molpdf)

                    scal = physical.values(default=1.0, em_density=10000)
                    icount=icount+equil_its
                init_vel=False

# cooling the system
            equil_its=200
            equil_temps=(800.0, 500.0, 250.0, 150.0, 50.0, 0.0)

            MD =  molecular_dynamics(cap_atom_shift=cap, md_time_step=timestep,
                                      md_return='FINAL', output='REPORT',
                                      schedule_scale=scal)

            for (its, equil, temps) in [(equil_its, equil_equil, equil_temps)]:
                for temp in temps:

                    MD.optimize(sel_all, max_iterations=its,
                       temperature=temp, init_velocities=init_vel, equilibrate=equil,
                       actions=[actions.trace(trc_step,trc_file)])

                    scal = physical.values(default=0.0, em_density=1.0)
                    (molpdf, terms) = sel_all.energy(schedule_scale=scal)
                    print "step= %d %d   temp= %d    EM score= %.3f" %(icount,its,int(temp),-molpdf)
                    scal = physical.values(default=1.0, em_density=10000)
                    icount=icount+equil_its

            filename='md'+run_num+'_'+str(a)+'.pdb'
            sel_all.write(file=filename)
            a+=1

        trc_file.close()

        print "MD step %d: energy all (scaling 1:10000)" % icount
        eval = sel_all.energy(schedule_scale=scal)

# final minimization with and without CC restraints
        CG = conjugate_gradients()
        print " final conjugate_gradients"

        CG.optimize(sel_all, output='REPORT', max_iterations=200,
                    schedule_scale=scal)
        eval = sel_all.energy(schedule_scale=scal)

        scal = physical.values(default=1.0, em_density=0.0)
        CG.optimize(sel_all, output='REPORT', max_iterations=200,
                    schedule_scale=scal)
        eval = sel_all.energy(schedule_scale=scal)

        sel_all.write(file='final'+run_num+'_mdcg.pdb')

        os.system("rm -f *.MRC")
