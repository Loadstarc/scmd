import pandas as pd
import numpy as np
from cluster.clus import cluster

def count_ele(arr):
    arr = np.array(arr)
    key = np.unique(arr)
    result = {}
    for k in key:
        mask = (arr == k)
        arr_new = arr[mask]
        v = arr_new.size
        result[k] = v
    return result

def clus_split(snum,sur_point):
    cluster(snum)
    label = pd.read_csv("cluster/cluster_label.dat",header=None,dtype = int)
    data = pd.read_csv("cluster/tmpco",header=None,sep=" ",dtype = float)
    label = np.array(label)
    label = label.reshape(label.shape[0])
    data = np.array(data,dtype=float)
    cluster_num = np.unique(label).shape[0]
    group_label = [[] for i in range (cluster_num)]
    print (count_ele(label))
    outfile = open("cluster/newco.txt",'w')
    for i in range(data.shape[0]):
        group_label[label[i]].append(i) 
    for cluscount in range (cluster_num):
        for i in group_label[cluscount][:sur_point]:
            for j in range(data.shape[1]):
                outfile.write('%.9f' % data[i][j]),
                if j == (data.shape[1]-1):
                    outfile.write("\n"),
                else:
                    outfile.write(" "),
    outfile.close()

