#!/usr/bin/python3
import numpy as np
import sys,os
import lib.softapi.tfcal 
from lib.constant import autoan


def irc_analy(filename):
    tfcal = lib.softapi.tfcal.nnpes(4,['C','H','H','O'])
    irclog = os.popen('readirc traj %s'%(filename)).readlines()
    irclabel = os.popen("grep 'Point Number:' %s"%(filename)).readlines()
    ircene = os.popen("grep 'SCF Done' %s"%(filename)).readlines()
    logfile = open('irc_ana.log','w')

    for c,stru in enumerate(irclog):
        #get info
        label = irclabel[c].split()[5]
        logfile.write(irclabel[c])
        #get NN message
        coord = np.array(stru.strip().split(),dtype=float).reshape(4,3) / autoan
        tfcal.getco(coord)
        F_nn,E_nn = tfcal.cal(4,check=True)
        #get gaussian energy
        refene = float(ircene[c].split()[4])
        dele = abs((E_nn-refene)*627.5)
        logfile.write('Gaussian energy is %.9f hartree, NNPES energy is %.9f hartree\nDifference is %.9f kcal/mol\n' %(refene,E_nn,dele))

    print ('IRC Analysis complete.')
    logfile.close()
    os.system('drawirc')
