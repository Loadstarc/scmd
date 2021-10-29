#!/usr/bin/python3
'''
Cross sampling method
'''
import copy
import sys,os
import numpy as np
install_path = '/share/home/linsc/bin/scmd'

def rungau(Npoint,turn_num=None):
    turn_num = 8 if turn_num is None else turn_num
    print ('>Starting Gaussian single point calculation.')
    end_symbol = []
    cu_path = os.getcwd()+'/'
    ms_num = Npoint//turn_num
    remain_num =  Npoint%turn_num
    title = ['#!/bin/bash\n','#PBS -l walltime=999:00:00','#PBS -l nodes=1:ppn=16',\
            '#PBS -q opa','#PBS -j oe','#PBS -V\n','cd $PBS_O_WORKDIR\n']

    #for the previous queue
    for turn in range(turn_num-1):
        start_num = 1+turn*ms_num
        end_num = (turn+1)*ms_num+1
        end_symbol.append(cu_path+'tmp/'+str(turn+1)+'/'+str(end_num-1)+'.log')
        print ('In turn %i, sumbit from point %i to point %i' %((turn+1),start_num,(end_num-1)))
        for i in range(start_num,end_num):
            os.system('mv %stmp/%s.gjf %stmp/%s/' %(cu_path,str(i),cu_path,str(turn+1)))
        os.chdir(cu_path+'tmp/'+str(turn+1))
        pbs_out = open(str(turn+1)+'.pbs','w')
        for line in title:
            pbs_out.write(line+'\n')
        pbs_out.write("for ((i=%i;i<%i;i++))\n" %(start_num,end_num))
        pbs_out.write("do\n")
        pbs_out.write("    g16 $i.gjf\n" )
        pbs_out.write("done\n")
        pbs_out.close()
        os.system("qsub %s.pbs" %(str(turn+1)))
        os.chdir (cu_path+'tmp/')

    #the last queue
    start_num = 1+(turn_num-1)*ms_num
    end_num = (turn_num)*ms_num+1+remain_num
    print ('In turn %i, sumbit from point %i to point %i' %(turn_num,start_num,(end_num-1)))
    for i in range(start_num,end_num):
        os.system('mv %stmp/%s.gjf %stmp/%s/' %(cu_path,str(i),cu_path,str(turn_num)))
    os.chdir(cu_path+'tmp/'+str(turn_num))
    pbs_out = open(str(turn_num)+'.pbs','w')
    for line in title:
        pbs_out.write(line+'\n')
    pbs_out.write("for ((i=%i;i<%i;i++))\n" %(start_num,end_num))
    pbs_out.write("do\n")
    pbs_out.write("    g16 $i.gjf\n" )
    pbs_out.write("done\n")
    pbs_out.close()
    os.system("qsub %s.pbs" %(str(turn_num)))

    #sumbit&check
    last_name = str(Npoint)+'.log'
    end_symbol.append(cu_path+'tmp/'+str(turn_num)+'/'+last_name)
    final_msg = 'Normal termination'
    while (True):
        if (Npoint > 7):
            if (os.path.isfile(end_symbol[0])==False or os.path.isfile(end_symbol[1])==False or os.path.isfile(end_symbol[2])==False or os.path.isfile(end_symbol[3])==False or os.path.isfile(end_symbol[4])==False or os.path.isfile(end_symbol[5])==False or os.path.isfile(end_symbol[6])==False or os.path.isfile(end_symbol[7])==False):
                print ('Calculating previous gaussian single points, previous pbs hasnt finished.')
                time.sleep(10)
                continue
        if (len(os.popen("tail -n 1 "+last_name,'r').readlines())==0):
            print ('Calculating previous gaussian single points')
            time.sleep(10)
        else:
            if ((final_msg in os.popen("tail -n 1 "+last_name,'r').readlines()[0]) == True):
                print ('Last point has finished')
                break
            else:
                print ('Last point was calculating')
                time.sleep(10)
	# Ending
    print ('Energy calculation process finished.')
    time.sleep(30)
    os.chdir(cu_path+'tmp/')
    os.system("for i in 1 2 3 4 5 6 7 8 ; do mv $i/*.gjf ./; mv $i/*.log ./; rm $i/*.*; done")
    os.chdir(cu_path)
