#! /usr/bin/python3
'''
here is a Program to refine dataset
'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sklearn 
import sys,os

def show_cutoff(data,cutoff):
    cut_index = []
    out_num = 0
    for ind,i in enumerate(data):
        if i > cutoff:
            out_num += 1
            cut_index.append(ind)
    return cut_index,out_num


if __name__ == '__main__':
    np.set_printoptions(precision=9)
    filename = sys.argv[1]
    infilel = open(filename,'r').readlines()
    data = []
    for line in infilel:
        data.append(line.split())
    data = np.array(data,dtype=float)
    dif = np.abs(data[:,0]-data[:,1])*627.5094
    index,otcount = show_cutoff(dif,2.0)
    ts = pd.read_csv('data.txt',sep=' ',header=None)
    ts = np.array(ts,dtype=float)
    rts = np.delete(ts,index,axis=0)
    for i in range(len(rts)):      
        for j in range(len(rts[i])):
            print('%.9f'%(rts[i][j]), end=' '),
        print (),


