#!/usr/bin/python3 -u
'''
NN-PES MD program beta
Switch back process was added.
Writen by: Shichen Lin
e-mail: linsc97@163.com
'''
import sys,os
import datetime
import numpy as np
import lib.hop as hop
import lib.iofile as iofile
import lib.constant as constant
import lib.softapi.gaucal as gaucal

from copy import copy
from lib.iofile import *
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

def get_divnum(path):
        name = os.listdir(path)
        for f in name:
            if ('.gjf' in f): return int(f.split('_')[1].split('.')[0])


def check_bdlen(pos0):
        bdinfo = [[1,2],[1,3],[1,4]]
        for bond in bdinfo:
                if (np.linalg.norm(pos0[bond[0]-1]-pos0[bond[1]-1]) * constant.autoan > 2.5):
                        print ('Msg:Atom %i - Atom %i bond cleveraged, MD terminated.' % (bond[0],bond[1]))
                        return True
        return False

def calc_torsion_angle(coordinates,at):
        """
        Calculate single torsion angle
        :param coordinates: Multidimensional array of coordinates
        :return: Torsion angle
        """
        """
        a, b, c, d are the points that make up a torsion angle.
        """
        a = coordinates[at[0]-1]
        b = coordinates[at[1]-1]
        c = coordinates[at[2]-1]
        d = coordinates[at[3]-1]
        vector_1, vector_2, vector_3 = b - a, c - b, d - c
        norm_vector_1 = np.cross(vector_1, vector_2)   # a-b-c plain
        norm_vector_2 = np.cross(vector_2, vector_3)   # b-c-d plain
        norm_vector_1x2 = np.cross(norm_vector_1, norm_vector_2)
        x = np.dot(norm_vector_1, norm_vector_2)
        y = np.dot(norm_vector_1x2, vector_2) / np.linalg.norm(vector_2)
        radian = np.arctan2(y, x)
        angle = radian * 180 / np.pi
        return angle


def runMD(Samname,dt,method,cycnum,if_restart):
        os.system('clear')
        print ('##Running MD simulation process, Method:%s'%(method))
        start = datetime.datetime.now()

        # Setting parameters
        if (if_restart == 0):
                posa,fara,moma,symbol,step = restart(Samname)
                atnum = 4
        else:
                posa,fara,moma,symbol,atnum = initdata(Samname)
        #cout = open('cluster/tmpco','w')
        trajname = Samname + ".zzip"
        cu_path = os.getcwd()+'/'
        dt = dt/constant.fstoau
        If_clev = False
        cycle = 0
        e_tot = 0.0

        #initialzing tensorflow
        if (method == 'nnpes'):
                import lib.softapi.tfcal
                tfcal=lib.softapi.tfcal.nnpes(atnum,symbol)
        gaussian = gaucal.gaussian()
        # Restart setting
        if (if_restart == 0):
                cycle = step
                if (os.path.exists('%s_dft' %(Samname)) == False):
                        tfcal.divcount = 0
                else:
                        tfcal.divcount = get_divnum('%s_dft' %(Samname))

        # Main cycle of MD-simulation
        while (cycle<cycnum):
                cycle +=1
                print ('---------- Cycle %i start ---------'%(cycle))

        #updating coordinate using verlet-algorithm
                pos0 = hop.verletx(posa[-1],fara[-1],moma[-1],symbol,atnum,dt)
                pos0 = np.array(pos0,dtype=float) 
                posa.append(pos0)

        #calculating Force & Energy

                if (method == 'nnpes'):
                        tfcal.getco(pos0)
                        if (tfcal.ifout_z == True or tfcal.ifout_s == True):
                                if (tfcal.divcount == 1):
                                        os.system('mkdir %s_dft' %(Samname))
                                        os.system("cp initial/"+Samname+".gjf "+Samname+"_dft/"+Samname+"_0.gjf")
                                os.chdir(cu_path+Samname+'_dft')
                                newinp = gaussian.reinpotfile(Samname, posa[-1], symbol, tfcal.divcount, atnum, False)
                                oup0 = gaussian.runsoft(newinp, 'r')
                                far0, energy = gaussian.parseout(oup0, atnum, trajname,start)
                                os.chdir(cu_path)
                                print ('>out range, ab-initial count:',tfcal.divcount)
                        else:
                                far0, energy = tfcal.cal(atnum)
                                if (tfcal.ifdiv == True):
                                        if (tfcal.divcount == 1):
                                                os.system('mkdir %s_dft' %(Samname))
                                                os.system("cp initial/"+Samname+".gjf "+Samname+"_dft/"+Samname+"_0.gjf")
                                        os.chdir(cu_path+Samname+'_dft')
                                        newinp = gaussian.reinpotfile(Samname, posa[-1], symbol, tfcal.divcount, atnum, False)
                                        oup0 = gaussian.runsoft(newinp, 'r')
                                        far0, energy = gaussian.parseout(oup0, atnum, trajname,start)
                                        os.chdir(cu_path)
                                        print ('>network div, ab-initial count:',tfcal.divcount)
                                else:
                                        print ('> Using NN-PES...')

                elif (method == 'gaussian'):
                        print ('calculating ab-initio')
                        if (cycle==1):
                                os.system('mkdir %s_dft' %(Samname))
                                os.system("cp initial/"+Samname+".gjf "+Samname+"_dft/"+Samname+"_0.gjf")
                        os.chdir(cu_path+Samname+'_dft')
                        newinp = gaussian.reinpotfile(Samname, posa[-1], symbol, cycle, atnum, False)
                        oup0 = gaussian.runsoft(newinp, 'r')
                        far0, energy, nbores = gaussian.parseout(oup0, atnum, trajname, start,True)
                        os.system('rm %s_%i.*'%(Samname,(cycle - 1)))
                        print ('rm %s_%i.*'%(Samname,(cycle - 1)))
                        os.chdir(cu_path)

                elif (method == 'cross'):
                        tfcal.getco(pos0)
                        far0, energy = tfcal.cal(atnum,True)



        #updating force&momentum with verlet algorithm
                fara.append(far0)
                mom0 = hop.verletp(moma[-1], fara[-1], fara[-2], dt)
                moma.append(mom0)

        #Check for Energy-stability & correction
                #if (cycle != 1):
                #        moma[-1] = hop.corr_etot(energy,moma[-1],symbol,etot)
                etot = hop.cal_etot(moma[-1],symbol,energy)
                print ('Final total energy is %.8f hartree'%(etot))


        #writing log file
                iofile.ncycle(trajname, cycle)
                iofile.outenergy(trajname, energy, 'potential')
                iofile.writecon(trajname, symbol, posa[-1], 'coordinate', atnum)
                iofile.writecon(trajname, symbol, fara[-1], 'Forces', atnum)
                iofile.writecon(trajname, symbol, moma[-1], 'momentum', atnum)
                if (method == 'gaussian'): iofile.outnbo(trajname,nbores)
                iofile.cyclend(trajname, cycle)

        #Ending
                print ('---------- Cycle %i done ----------'%(cycle))
                If_clev = check_bdlen(pos0)
                if (If_clev == True):break

        print('==============MD process normal exit===============')
        iofile.endtime(trajname, 'cycle', cycle, start)

