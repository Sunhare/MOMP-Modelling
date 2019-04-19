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

def dump_args(func):
    # "This decorator dumps out the arguments passed to a function before calling it"
    argnames = func.func_code.co_varnames[:func.func_code.co_argcount]
    fname = func.func_name
    def echo_func(*args,**kwargs):
        print(fname, "(", ', '.join(
            '%s=%r' % entry
            for entry in zip(argnames,args[:len(argnames)])+[("args",list(args[len(argnames):]))]+[("kwargs",kwargs)]) +")")
    return echo_func


test_cst = dump_args(pyrosetta.rosetta.core.scoring.constraints.Constraint.read_def())