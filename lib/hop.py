#!/usr/bin/env python
'''
verlet algorithm
'''
import os
import math
import datetime
import numpy as np
from lib.constant import evtoau, mas, au2kcal

def verletx(pos1, far1, mom1, symb, numele, dt):
        pos1 = np.array(pos1, dtype = float)
        far1 = np.array(far1, dtype = float)
        mom1 = np.array(mom1, dtype = float)

        for i in range(numele):
                mass = float(mas[symb[i]])
                x1 = pos1[i] + mom1[i] * dt / mass + far1[i] * dt ** 2 / (2 * mass)
                pos1[i] = x1
        pos0 = tuple(map(tuple, pos1))

        return pos0

#-------------------------------------------------------------------------------------

def verletp(mom1, far1, far2, dt):
        mom1 = np.array(mom1, dtype = float)
        far1 = np.array(far1, dtype = float)
        far2 = np.array(far2, dtype = float)

        mom1 = mom1 + dt * (far1 + far2) / 2
        mom0 = tuple(map(tuple, mom1))

        return mom0

#-------------------------------------------------------------------------------------

def calkine(pp, symb):
        pp = np.array(pp, dtype = float)
        kine = 0
        for i in range(len(pp)):
                k = pp[i] * pp[i] / (mas[symb[i]] * 2)
                kine += np.sum(k)

        return kine

#-------------------------------------------------------------------------------------
def cal_etot(mom,symb,ep):
        return calkine(mom,symb) + ep

#-------------------------------------------------------------------------------------
def corr_etot(ep,mom,symb,pre_etot):
        de_thre = 2.0 # test of correction
        cu_ek = calkine(mom,symb)
        cu_etot = cu_ek + ep
        detot = abs(cu_etot - pre_etot) * au2kcal
        print ('Now, total energy is %.8f hartree, previous is %.8f hartree'%(cu_etot,pre_etot))
        if (detot > de_thre):
            corr_ek = pre_etot - ep
            if (corr_ek < 0.0):
                print ('Potential gap is too large, can not do correction.')
                exit()
            corr_ratio = np.sqrt(corr_ek / cu_ek)
            mom = mom * corr_ratio
            print ('after correction, energy is %.8f hartree'%(cal_etot(mom,symb,ep)))
            return mom
        return mom

