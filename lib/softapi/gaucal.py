#!/usr/bin/env python
'''
set up the API for use quantum chemistry self
'''

import os
import re
import shutil
import datetime
import time
import numpy as np

import lib.iofile as iofile
from lib.constant import evtoau, autoan


#####################################################################################
#                            class gaussian                                         #
#####################################################################################

class gaussian():
        '''
        use the gaussian soft
        '''
        def __init__(self):
                self = self

#-------------------------------------------------------------------------------------

        def getname(self, filename):
                outname = filename + '.gjf'

                return outname

#-------------------------------------------------------------------------------------

        def reinpotfile(self, filename, newpos, symb, intcycle, numele,ifgau):
                '''
                use the new position to create a new gjf file
                '''
                oldname = filename + "_" + str(intcycle - 1) + '.gjf'
                oldlog = filename + "_" + str(intcycle - 1) + '.log'
                newname = filename + "_" + str(intcycle) + '.gjf'

                oldfile = open(oldname,'r').readlines()
                newfile = open(newname, 'w')

                newpos = np.array(newpos, dtype = float) * autoan
                linenu = 4 # ensure the positine line number.
                for i, line in enumerate(oldfile):
                        if line[0] in ['%', '#']:
                                linenu += 1
                        if ('#p' in line and intcycle == 2 and ifgau == True):
                                oldfile[i] = oldfile[i].strip() + ' guess=read\n'

                for i in range(linenu): 
                        newfile.write(oldfile[i])

                for i in range(linenu, linenu + numele): 
                        newfile.write(' {0:<5}{1[0]:15.8f}{1[1]:15.8f}{1[2]:15.8f}\n'.format(symb[i - linenu], newpos[i - linenu]))

                for i in range(linenu + numele, len(oldfile)):
                        newfile.write(oldfile[i])

                newfile.close()
                os.system('rm %s'%oldname)
                os.system('rm %s'%oldlog)
                return newname

#-------------------------------------------------------------------------------------

        def getpos(self, filename, numele):
                '''
                from the first gjf file to get the first position
                '''
                oldfile = open(filename).readlines()
                filelen = len(oldfile)

                linenu = 4 # ensure the positine line number.
                for i in range(filelen):
                        if oldfile[i][0] in ['%', '#']:
                                linenu += 1

                numpos  = oldfile[linenu: linenu + numele]

                pos = []
                sym = []
                for i in range(numele):
                        posarray = numpos[i].split()
                        symbol = posarray[0].title()

                        pos.append(posarray[1:])
                        sym.append(symbol)

                pos = np.array(pos, dtype = float) / autoan
                sym = tuple(sym)
                pos = tuple(map(tuple, pos))

                return pos, sym
        
#-------------------------------------------------------------------------------------

        def runsoft(self, runfile, keysoft):
                '''
                for run the gaussian progam, and return the log file
                '''
                outname = runfile.split('.')[0] + '.log'

                if keysoft == 'r':
                        os.system('g16 %s' %runfile)

                return outname

#------------------------------------------------------------------------------------

        def parseout(self, resultname, numele, outname, starttime, if_nbo=None):
                resultfile = open(resultname).readlines()

                signf = 0
                signe = []
                for i, line in enumerate(resultfile):
                        if 'Forces ' in line:
                                signf = i

                        if 'SCF Done' in line:
                                genergy = float(line.split()[4])
                                signe.append(line)

                        if 'Total Energy' in line:
                                signe.append(line)

                        if (if_nbo == True):
                                if ('Summary of Natural Population Analysis' in line):
                                        nboline = i+6

                if 'Forces ' in resultfile[signf]:
                        numf = resultfile[signf + 3: signf + 3 + numele]
                else:
                        iofile.werror(outname, 'forces', starttime)
                        print('no force in outfile')
                        exit(1)

                sfar = []
                for i in range(numele):
                        eachforce = numf[i].split()
                        sfar.append(eachforce[2:])

                sfar = np.array(sfar, dtype = float)
                far = tuple(map(tuple, sfar))
                
                ii = signe[-1]
                if 'SCF Done' in ii:
                        energyrow = ii.split()[4]

                elif 'Total Energy' in ii:
                        energyrow = ii.split()[-1]

                else:
                        iofile.werror(outname, 'energy', starttime)
                        print('no energy in outfile')
                        exit(1)

                if (if_nbo == True):
                        nbores = []
                        for line in resultfile[nboline:nboline+numele]:
                                nbores.append(line.split())
                        return far,genergy,nbores

                return far,genergy

