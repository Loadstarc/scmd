#!/usr/bin/python3
'''
Recalculate & check the trajectories by gaussian
'''
import numpy as np
import sys,os
import lib.softapi.gaucal 
from lib.constant import autoan
from lib.constant import au2kcal

def readtraj(filename):
    trajin = open(filename,'r').readlines()
    energy = []
    coord = []
    clines = []
    symb = []
    read_symb = True
    for c,line in enumerate(trajin):
        if ('the total potential energy' in line):
            energy.append(float(line.split()[5]))
        if ('The coordinate' in line):
            clines.append(c+4)
    for i in clines:
        co = []
        for l in trajin[i:i+4]:
            co.append(np.array(l.split()[1:],dtype=float))
            if (read_symb == True):
                symb.append(l.split()[0])
                read_symb == False
        co = np.array(co,dtype=float)
        coord.append(co)
    return energy,coord,symb

def initgjf(co,symb,name):
    title = ['%chk=ene.chk','%mem=32gb','%nproc=16','#p ub3lyp/6-31g(d) \n','Calculating energy\n','0 3']
    jobname = "%s_1.gjf"%(name)
    out = open(jobname,'w')
    for line in title:
        out.write(line+"\n")
    for a in range(4):
        out.write('{0:>2s} {1:>16.9f}  {2:>16.9f}  {3:>16.9f}\n'\
                .format(symb[a],co[a][0],co[a][1],co[a][2])),
    out.write('\n'),
    out.close()


def recal_traj(filename):
    gaussian = lib.softapi.gaucal.gaussian()
    samname = filename.split('.')[0]
    Enet,coords,symb = readtraj(filename)
    os.system('mkdir %s_recal'%(samname))
    os.chdir('%s_recal'%(samname))
    for num,co in enumerate(coords):
        nsam = num + 1
        if (nsam == 1):
            initgjf (co,symb,samname)
        else:
            jobname = gaussian.reinpotfile(samname, co/autoan, symb, nsam, 4, True)
            os.system('rm %s_%i.gjf'%(samname,num))
            os.system('rm %s_%i.log'%(samname,num))
        logname = gaussian.runsoft('%s_%i.gjf'%(samname,nsam),'r')
        Eref = float(os.popen("grep 'SCF Done' %s"%(logname)).readlines()[0].split()[4])
        print ('Step %i> NN energy: %.9f hartree; DFT energy: %.9f; Diff: %.9f kcal/mol'\
                %(nsam,Enet[num],Eref,abs(Enet[num]-Eref)*au2kcal))
