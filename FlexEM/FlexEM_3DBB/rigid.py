#  ============================================================================
#
#                           reading rigid bodies   
#
#  =======================  Maya Topf, 4 Dec 2007 =============================
#  =======================  Latest update: 6/8/09 =============================

from modeller import *

def load_rigid(path,mdl,sel_rigid,rand_rigid,rigid_filename):
    in_file = open(path + '/' + rigid_filename,"r")
    flag = 0
    while True:
        in_line = in_file.readline()
        if len(in_line) == 0:
            break
        if in_line[0] == '#':
            continue
        elif in_line[0] == ' ':
            break
        elif in_line == '\n':
            break
        else:
            resnum = in_line.split()
            i=0
            list=selection()
            if resnum[-1] == 'nr':
                randflag = 0
                length = len(resnum) - 1
            else:
                randflag = 1
                length = len(resnum)
            while i < length:
                list.add(mdl.residue_range(resnum[i], resnum[i+1]))
                i=i+2
            sel_rigid.append(list)
            rand_rigid.append(randflag)
    in_file.close()
