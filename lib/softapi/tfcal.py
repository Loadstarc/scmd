#/usr/bin/python3
'''
Behler-Parrinello NN-PES calculation API
'''
import sys,os
import math
import numpy as np
import lib.constant as constant
from time import time
from copy import copy
from lib.md import calc_torsion_angle
from lib.softapi.trans import trans
from lib.softapi.gra1 import gra_n1
from lib.softapi.gra2 import gra_n2
from lib.iofile import ptmat,print_info
np.set_printoptions(precision=3, suppress=True)
np.set_printoptions(formatter={'float': '{: 0.9f}'.format})
storge_path = '/home/loadstar/desktop/scmd'

class nnpes:
        def __init__(self,atnum,sbl):
                self = self
                self.ifout_s = False
                self.ifout_z = False
                self.ifdiv = False
                self.cross_stat = 0
                self.divcount = 0
                self.symbol = sbl
                self.atnum = atnum
                self.limit = np.loadtxt(storge_path+"/n1/limit")
                self.slimit = np.loadtxt(storge_path+"/slimit")
                self.zlimit = np.loadtxt(storge_path+'/zlimit',comments='#')
        
        def checkrange(self,symmco):
                co = copy(symmco)
                co = co.reshape(60)
                for i,c in enumerate(co):
                        cdel = (self.slimit[i][0] - self.slimit[i][1])*0.005
                        cmax = self.slimit[i][0] #- cdel
                        cmin = self.slimit[i][1] #+ cdel
                        if (c > cmax or c < cmin):
                                print ("Warning: function number:",str(i+1),"is going out of Confidence range!")
                                self.ifout_s = True
                                return 0
                print ("> Confidence region check was clear.")
                self.ifout_s = False
                return 1

        def check_dihedral(self):
                dihedral_list = [[1,5,8,11],[5,8,11,14],[8,11,14,16],[8,11,14,15]]
                for d in dihedral_list:
                        angle = calc_torsion_angle(self.coord,d)
                        print ('Dihedral of atom: {0[0]}-{0[1]}-{0[2]}-{0[3]} is {1:.2f} '.format(d,angle))

        def getco(self,co):
                #get coordinate & calculate gradient of symmfunction
                self.coord = copy(co) * constant.autoan
                self.symmf,self.gra = trans(self.coord)
                self.checkrange (self.symmf) 
                self.check_zmat ()
                if (self.ifout_s==True or self.ifout_z == True) : self.divcount += 1
        
        def get_zmatrix(self):
                X = copy(self.coord).T
                m,n = X.shape
                G = np.dot(X.T, X)
                H = np.tile(np.diag(G), (n,1))
                self.zmat = np.sqrt(H + H.T - 2*G)

        def check_zmat(self):
#               print ('input coord')
#               print_info(self.coord,self.symbol)
#               print ('check z-matrix:')
                self.get_zmatrix()
                tri_zmat = np.zeros(int((1+self.atnum)*self.atnum/2))
#               ptmat(self.zmat)
                c = 0
#               print ('Check tri:')
                for i in range(self.atnum):
                        for j in range(i+1):
                                tri_zmat[c] = self.zmat[i][j]
#                               print ("%.9f "%(tri_zmat[c]),end='')
                                c += 1
                for i,z in enumerate(tri_zmat):
                        zmax = self.zlimit[i][0] #- zdel
                        zmin = self.zlimit[i][1] #+ zdel
#                       print ('%.9f %.9f %.9f'%(z,zmax,zmin))
                        if ((z > zmax or z < zmin)):
                                print ("Warning: Z-matrix number:",str(i+1),"is going out of Confidence range!")
                                nat1 = int(math.floor(0.5*(-1+np.sqrt(1+4*2*(i+1)))))
                                nat2 = int((i+1) - 0.5*(nat1*(nat1+1)))
                                nat1 += 1
                                print ('>Distance between Atom %i - Atom %i is %.2f'%(nat1,nat2,z))
                                self.ifout_z = True
                                return 0
                print ("> Z-matrix check was clear.")
                self.ifout_z = False
                return 1


        def cal(self,cross=None,out=None,check=None):
                cross = False if cross is None else cross
                out = None if out is None else out
                check = False if check is None else check
                E1,ME1 = gra_n1(self.symmf,self.limit)
                self.cross_stat = self.checkdiv(E1,cross)
                if (self.ifdiv == True and check == False): return 0,0
                for i in range(self.atnum):
                        for d in range(3):
                                self.gra[d][i] = np.multiply(self.gra[d][i],ME1)
                force = np.sum(self.gra,axis = 3)
                force = -np.sum(force,axis = 2).T
                force = (force*constant.autoan)
                return force,E1

        def checkdiv(self,E0,cross,out=None):
                out = None if out is None else out
                E2,ME2 = gra_n2(self.symmf,self.limit)
                Dif = abs(E0-E2)*constant.au2kcal
                #print('In checkdiv, we shutdown double NN check, please check!')
                print ('E1 = %f hartree, E2 = %f hartree\nDifE = %f kcal/mol'%(E0,E2,Dif))
                #print ('Network Diverge Check was shut down.')
                if (Dif > 0.2):
                        print ("W:Network diverge")
                        self.ifdiv = True
                        self.divcount += 1
                        if (cross == True):
                                self.add_randompoint(self.coord,out)
                        return 1
                elif (Dif > 2.0 and cross == True):
                        print ('W: DifE over 2.0kcal/mol, stop sampling')
                        return 2
                else:
                        print ("Network donot diverge here.")
                        self.ifdiv = False
                        return 0
