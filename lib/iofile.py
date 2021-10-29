#!/usr/bin/env python
'''
for create output file
'''
import numpy as np
import datetime
from lib.constant import *

def initdata(Samname):
    print ("===============Initializing Data================")
    gjfin = "initial/" + Samname + ".gjf"
    login = "initial/" + Samname + ".log"
    momin = "initial/" + Samname + ".dat"
    #Read momentum
    with open(momin) as momf:
        moml = momf.readlines()
    momr = []
    symbol = []
    for line in moml[2:]:
        symbol.append(line.split()[0])
        atmas = float(mas[line.split()[0]])
        momr.append(np.array(line.split()[1:],dtype=float)*atmas)
    momr = np.array(momr, dtype = float)
    print ("> Momentums Get,echo here")
    print_info(momr,symbol)
    #get coordinates
    gjfl = open(gjfin,'r').readlines()
    coord = []
    for line in gjfl[8: -1]:
        coord.append(line.split()[1:])
    coord = np.array(coord, dtype = float) / autoan
    atnum = len(gjfl[8: -1])
    print (">Coordinate Get, echo here(unit:Bohr)")
    print_info(coord,symbol)
    #read force
    logl = open(login,'r').readlines()
    for i, line in enumerate(logl):
        if 'Forces ' in line: fline = (i + 3)
        if 'SCF Done' in line: eline = i
    force = []
    for line in logl[fline: fline + atnum]:
        force.append(line.split()[2:])
    force = np.array(force, dtype = float)
    eneini = float(logl[eline].split()[4])
    print ("> Force Get, echo here(unit:hartree/ang)")
    print_info(force,symbol)
    pos = [coord]
    mom = [momr]
    far = [force]
    print ("==============Initialization Over===============")

    return pos,far,mom,symbol,atnum

def restart(Samname):
    traj = open(Samname+".zzip",'r').readlines()
    posa,fara,moma,symbol = [],[],[],[]

    bl = []
    for l in traj[-14:-10]:
        symbol.append(l.split()[0])
        bl.append(l.split()[1:])
    bl = np.array(bl,dtype=float)
    moma.append(bl)

    bl = []
    for l in traj[-25:-21]:
        bl.append(l.split()[1:])
    bl = np.array(bl,dtype=float)
    fara.append(bl)

    bl = []
    for l in traj[-36:-32]:
        bl.append(l.split()[1:])
    bl = np.array(bl,dtype=float) / autoan
    posa.append(bl)
    step = int(traj[-7].split()[0])
    return posa,fara,moma,symbol,step


def print_info(msg,sbl):
    for c,at in enumerate(msg):
        print ('%3s %16.9f %16.9f %16.9f'%(sbl[c],at[0],at[1],at[2]))

def ptmat(a):  # print any matrix with any dimention
    n = len(a[0])
    rows = n // 5
    rmnd = n % 5
    for i in range(rows):
        for ii in range(5):
            print("    ", end = "")
            print("{0:9d} ".format(i * 5 + ii + 1), end = "")
        print('\n')

        for ij, array in enumerate(a):
            print("{0:3d} ".format(ij + 1), end = "")
            for aa in array[i * 5: (i + 1) * 5]:
                print(" {:15.10f} ".format(aa), end = "")
            print()
        print()

    if rmnd:
        for j in range(rmnd):
            print("    ", end = "")
            jl = rows * 5 + j + 1
            print("{0:9d} ".format(jl), end = "")
        print('\n')

        for jj, array in enumerate(a):
            print("{0:3d} ".format(jj + 1), end = "")
            for aa in array[rows * 5:]:
                print(" {:15.10f} ".format(aa), end = "")
            print()
        print()

def writecon(outname, sym, con, group, numele):
        '''
        writecon means write content.
        because the out file has to many similarity content about force, momentum and position
        so use this function to write down the similarity content in every cycle
        '''
        con = np.array(con, dtype = float)
        con = con * autoan if group == 'coordinate' else con

        outfile = open(outname, 'a+')
        outfile.write('\n' + 'The {0}'.format(group).center(75) + '\n')
        outfile.write('-'*75 + '\n')
        outfile.write(' {0:<12}{1:>16}{2:>16}{3:>16}\n'.format('Symbol', 'x', 'y', 'z'))
        outfile.write('-'*75 + '\n')

        for i in range(numele):
                outfile.write('  {0:<16}{1[0]:>16.9f}{1[1]:>16.9f}{1[2]:>16.9f}\n'.format(sym[i], con[i]))

        outfile.write('-'*75 + '\n\n')
        outfile.close()

        return

#------------------------------------------------------------------------------------

def werror(outname, keyword, starttime):
        '''
        for write down the error infomation
        '''

        endtime = datetime.datetime.now()
        usetime = endtime - starttime

        outfile = open(outname, 'a+')
        outfile.write('*'*99 + '\n')
        outfile.write('ERROR %s' %keyword + '\n')
        outfile.write('Grad'*24 + '\n')
        outfile.write('*'*99 + '\n')
        outfile.write('no {0} in outfile\n'.format(keyword))
        outfile.write('start program at {0}\n'.format(starttime))
        outfile.write('use time {0}\n'.format(usetime))
        outfile.close()

        return

#------------------------------------------------------------------------------------

def endtime(outname, keyword, intcycle, starttime):
        endtime = datetime.datetime.now()
        usetime = endtime - starttime

        outfile = open(outname, 'a+')
        outfile.write('\n' + '*'*99 + '\n')
        outfile.write('cycle %d normal exit' %intcycle + '\n')
        outfile.write('Start program at  ' + str(starttime) + '\n')
        outfile.write('Use time  ' + str(usetime) + '\n')
        outfile.write('%s normal exit\n' %keyword)
        outfile.close()
        
        return

#------------------------------------------------------------------------------------

def outstate(outname, state):
        outfile = open(outname, 'a+')
        outfile.write('the calculate state is %s' %state)
        outfile.close()
        return

#------------------------------------------------------------------------------------

def outnbo(outname, nbores):
        outfile = open(outname, 'a+')
        outfile.write('Summary of Natural Population Analysis:\n')
        outfile.write('%6s %4s %8s %12s %12s %12s %12s\n'%('Atom','No','Charge','Core','Valence' ,'Rydberg','Total'))
        for  nboat in nbores:
                outfile.write('%6s %4s %8s %12s %12s %12s %12s\n'%(nboat[0],nboat[1],nboat[2],nboat[3],nboat[4],nboat[5],nboat[6]))
        outfile.close()
        return
#------------------------------------------------------------------------------------

def outenergy(outname, energy, keyword):
        outfile = open(outname, 'a+')
        outfile.write('the total %s energy is %f a.u. \n' %(keyword, energy))
        outfile.close()

        return

#------------------------------------------------------------------------------------

def wstate(outname, state, deltav, flip):
        outfile = open(outname, 'a+')
        outfile.write('the state in %d \n' %state)

        if deltav:
                gap = np.array(deltav[-1], dtype = float)
                for i, line in enumerate(gap):
                        outfile.write('the {0} gap is {:.4f}\t'.format(i + 1, line))

        else:
                outfile.write('adiabatic molecular dynamic')

        outfile.write('\n')
        outfile.close()

        return

#------------------------------------------------------------------------------------

def whop(outname, intcycle, state):
        outfile = open(outname, 'a+')
        outfile.write('\nthe state at %d hop to state %d \n' %(intcycle, state))
        outfile.close()

        return

#------------------------------------------------------------------------------------

def writeab(outname, asquare, bsquare):
        outfile = open(outname, 'a+')
        outfile.write('\neffective nonadiabatic coupling asquare is {:.6f}\n'.format(asquare))
        outfile.write('effective collision energy bsquare is {:.6f}\n'.format(bsquare))
        outfile.close()

        return

#------------------------------------------------------------------------------------

def wpes(outname, allpes, genergy, spin):
        outfile = open(outname, 'a+')

        for i, energy in enumerate(allpes):
                outfile.write('the %s%d state pes energy is %f a.u.\n to ground delta energy is %.4f (eV)\n' %(spin, i + 1, energy + genergy, energy / evtoau))

        outfile.close()

        return

#------------------------------------------------------------------------------------

def ncycle(outname, intcycle):

        outfile = open(outname, 'a+')
        outfile.write('\n' + '='*84 + '\n')
        outfile.write('Cycle {0:d}'.format(intcycle).center(84) + '\n')
        outfile.write('='*84 + '\n'*2)
        outfile.close()

        return

#------------------------------------------------------------------------------------

def cyclend(outname, intcycle):
        outfile = open(outname, 'a+')
        outfile.write('\n %d cycle end \n' %intcycle)
        outfile.close()

        return

#------------------------------------------------------------------------------------

def spindr(outname, spin, state):
        outfile = open(outname, 'a+')
        outfile.write('\nthe spin state is {0}{1}\n'.format(spin, state))
        outfile.close()

        return

#------------------------------------------------------------------------------------

def trunline(outname):
        outfile = open(outname, 'rb+')
        allfile = outfile.readlines()
        endfile = len(allfile[-1])
        outfile.seek(-(endfile + 1), 1)
        outfile.truncate()
        outfile.close()

        return

