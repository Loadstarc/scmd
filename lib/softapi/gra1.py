#!/usr/bin/env python
import tensorflow.compat.v1 as tf
import numpy as np
import copy
import os,sys
tf.disable_eager_execution()
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
storge_path = '/home/loadstar/desktop/scmd'
Atnum = 4
Elenum = 3
Gnum = 15
cpu_num = 1
E_bal = -114.562288968
#======================= building tensors =====================
# Class of Layer in NN
# NEED TO SET: 
#1.node size of IN/OUT layer;2. NUM of elements;3. Type of activation function
g1 = tf.Graph()
sess1 = tf.Session(graph=g1,config=tf.ConfigProto(intra_op_parallelism_threads=1,inter_op_parallelism_threads=1))
with g1.as_default():
        class layer:
                empCount = 0
                def __init__(self,in_size,out_size):
                        self.N_in = in_size
                        self.N_out = out_size
                        self.Weight = [[] for i in range(Elenum)]
                        self.bias = [[] for i in range(Elenum)]
                        self.Wx_plus_b = [[] for i in range(Atnum)]
                        self.output = [[] for i in range(Atnum)]
                        layer.empCount += 1

                def buildlayer(self,Elenum):
                        for Elecount in range (Elenum):
                                self.Weight[Elecount] = tf.Variable(tf.random.uniform([self.N_in,self.N_out],minval = -np.sqrt(1/(self.N_in)),maxval=np.sqrt(1/(self.N_in)) , dtype = tf.float32, seed = None, name = None))
                                self.bias[Elecount] = tf.Variable(tf.zeros([1,self.N_out]))

                def cal_layer(self,pre_layer,Atnum):
                        for Atcount in range(Atnum):
                                if (Atcount==0) :
                                        Num = 0
                                elif (Atcount==1 or Atcount==2):
                                        Num = 1
                                else :
                                        Num = 2
                                self.Wx_plus_b[Atcount] = tf.matmul(pre_layer[Atcount],self.Weight[Num]) + self.bias[Num]
                                self.output[Atcount] = tf.nn.sigmoid(self.Wx_plus_b[Atcount])

        infilel = open(storge_path+"/n1/input",'r').readlines()
        paradict = {}
        for line in infilel[1:]:
           paradict[line.split("=")[0].strip()] = line.split("=")[1].strip() 
        layer_size = int(paradict['Layer_size'])
        model_path = storge_path+"/n1/tf_save/check.ckpt"
        #Define space & layer
        gs = [[] for i in range(Atnum)]
        Maximums = [tf.Variable(0.0,dtype=tf.float64, name = None) for i in range(Atnum*Gnum+1)]
        Minimums = [tf.Variable(0.0,dtype=tf.float64, name = None) for i in range(Atnum*Gnum+1)]
        for i in range (Atnum):
                gs[i] = tf.placeholder(tf.float32, [None, Gnum])
        ys = tf.placeholder(tf.float32, [None,1])
        Layer1 = layer(Gnum,layer_size)
        Layer1.buildlayer(Elenum)
        Layer2 = layer(layer_size,layer_size)
        Layer2.buildlayer(Elenum)
        Layer3 = layer(layer_size,layer_size)
        Layer3.buildlayer(Elenum)
        Layer4 = layer(layer_size,1)
        Layer4.buildlayer(Elenum)
        #Building Network
        Layer1.cal_layer(gs,Atnum)
        Layer2.cal_layer(Layer1.output,Atnum)
        Layer3.cal_layer(Layer2.output,Atnum)
        Layer4.cal_layer(Layer3.output,Atnum)
        P = Layer4.Wx_plus_b
        pre = tf.concat([P[0],P[1]],1)
        for i in range (Atnum-2):
                pre = tf.concat([pre,P[i+2]],1)
        trans = tf.ones([Atnum,1], tf.float32)
        prediction = tf.matmul(pre,trans)
        #loss&save
        loss = tf.reduce_mean(tf.square(ys - prediction))
        init = tf.global_variables_initializer()
        saver = tf.train.Saver()
        gradient = tf.gradients(prediction,gs)
        tf.get_default_graph().finalize()
        saver.restore(sess1,model_path)

def gra_n1(co,n1l):
        Rmax = n1l[:,0]
        Rmin = n1l[:,1]
        Emax = Rmax[60] - E_bal
        Emin = Rmin[60] - E_bal
        Edif = Emax - Emin
        #==================== Preconditioning data =====================
        rawdata = copy.copy(co)
        rawdata = rawdata.reshape(1,Atnum*Gnum)
        Spnum = rawdata.shape[0]
        Arynum = rawdata.shape[1]
        #print (rawdata)
        for i in range (Arynum):
                for j in range (Spnum):
                        rawdata[j][i] = (2*(rawdata[j][i] - Rmin[i]) / (Rmax[i] -Rmin[i])) - 1
        #==================== Calculating ========================================
        G = [[]] *(Atnum)
        for i in range (Atnum):
                G[i] = rawdata[: , (Gnum*i) : (Gnum*(i+1))]
        y_data = rawdata[ : ,(Arynum-1) : Arynum]
        d = {}
        for i in range (Atnum):
                d[gs[i]] = G[i]
        d[ys] = y_data
        ### Symmetry functions were NORMALIZED!!
        Dif = (Rmax - Rmin)[:Gnum*Atnum]
        Dif = Dif.reshape(Atnum,Gnum)
        Dif = 2/Dif 
        gra = np.array(sess1.run(gradient,feed_dict=d)).reshape(Atnum,Gnum)
        gra = gra*(0.5*Edif)### Energy was NORMALIZED!!
        gra = np.multiply(gra,Dif)
        pre = np.array(sess1.run(prediction,feed_dict=d),dtype=float)
        pre = (((pre+1)*0.5*Edif)+Emin) + E_bal
        tf.reset_default_graph()
        return pre,gra
