#/usr/bin/env python
#Molecular Type: Pentenal
#--------------- Precondition Matrixs  --------------------
import copy
import numpy as np
atnum = 4
Granum = 3 
Gannum = 12
eta = 0.8
lmbda = 1.0
AtElec = np.array([6.0,1.0,1.0,8.0]) #eletron number/ Molecular type:water
F1 = np.zeros((atnum,atnum,atnum))
F2 = np.zeros((atnum,atnum,atnum))
F3 = np.zeros((atnum,atnum,atnum))
# !!!!!! Xr[i][j][k] i is the center atom (was derivativing (FEN MU))
Xr = np.zeros((atnum,atnum,Granum)) # Derivative for radium function (x-axis)
Yr = np.zeros((atnum,atnum,Granum)) # ~ y-axis
Zr = np.zeros((atnum,atnum,Granum)) # ~ Z-azis
Xa = np.zeros((atnum,atnum,Gannum)) # Derivative for angular function (x-axis)
Ya = np.zeros((atnum,atnum,Gannum)) # ~ y-azis
Za = np.zeros((atnum,atnum,Gannum)) # ~ z-azis

#Cutoff Matrix
def Fc(dist):
    fc,mfc = np.zeros((atnum,atnum)),np.zeros((atnum,atnum))
    Rc = 6.0
    for i in range(atnum):
        for j in range(atnum):
            if (dist[i][j]<=Rc):
                fc[i][j] = 0.5*(np.cos(dist[i][j]*np.pi/Rc)+1)
                mfc[i][j] = (-0.5*np.pi/Rc)*np.sin(np.pi*dist[i][j]/Rc)
    return fc,mfc
#Distance cal
def dis(co):
    Dis = np.zeros((atnum,atnum))
    Dis_az = np.zeros((3,atnum,atnum))
    md_x = np.zeros((3,atnum,atnum))
    for i in range(atnum):
        for j in range(i+1,atnum):
                Dis[i][j] = np.linalg.norm(co[i]-co[j])
                Dis[j][i] = Dis[i][j]
                for d in range(3):
                    Dis_az[d][i][j] = co[i][d] - co[j][d]
                    Dis_az[d][j][i] = -Dis_az[d][i][j]
                    md_x[d][i][j] = Dis_az[d][i][j]/Dis[i][j]
                    md_x[d][j][i] = -md_x[d][i][j]
    return Dis,Dis_az,md_x
    '''
    X/Y/Z distance can be repalced by difference of vector
    '''
#Cosine cal& F2/F3 Basis
def cosine(co,Dis):
    Costheta = np.zeros((atnum,atnum,atnum))
    M_cos_ij = np.zeros((3,atnum,atnum,atnum))
    M_cos_ik = np.zeros((3,atnum,atnum,atnum))
    for i in range(atnum):
        for j in range (atnum):
            if (j != i):
                for k in range(atnum):
                    if (k != i and k !=j):
                        Rij = co[i] - co[j]
                        Rik = co[i] - co[k]
                        Costheta[i][j][k] = np.dot(Rij,Rik)/(Dis[i][j]*Dis[i][k])
                        for d in range(3):
                            M_cos_ij[d][i][j][k] = (Rik[d]/(Dis[i][j]*Dis[i][k])) - ((Rij[d]*Costheta[i][j][k])/(Dis[i][j]**2))
                            M_cos_ik[d][i][j][k] = (Rij[d]/(Dis[i][j]*Dis[i][k])) - ((Rik[d]*Costheta[i][j][k])/(Dis[i][k]**2))
    return Costheta,M_cos_ij,M_cos_ik
'''
Costheta(i,j,k) is for the angle of atom i,j,k, atom i is the center atom.
'''

#============= calculating symmfunction&derivative ==========
def trans(cco):
    Atco = copy.copy(cco)
    #Prepare for some Matrixs
    Dis,Dis_az,M_r_c = dis(Atco) #Distance,Dis_az is the distance in x/y/z axis
    M_X_x = np.array([[0,0,-1],[0,0,0],[0,0,0]])#for a == i, deravative from Xij to xi
    Costheta,M_cos_ij,M_cos_ik = cosine(Atco,Dis) #Cosine and Gradients
    fc,M_fc_r = Fc(Dis) #cutoff
    M_G = np.zeros((3,atnum,atnum,Granum+Gannum))
    Ga = np.zeros((atnum,Gannum)) # List of angular function
    # Calculate Radium function || CHECKED!
    Gr = np.zeros((atnum,Granum))
    for s in range (Granum):
        deltaRs = 2.5/(Granum-1)
        Rs = 0.5 + deltaRs * s
        gr = np.exp(-eta*np.power((Dis - Rs),2))
        g = np.multiply(gr,fc)
        # deravative of radium function
        M_gr_r = np.multiply(gr,(-2*eta*(Dis-Rs))) 
        for i in range(atnum):
            for j in range(atnum):
                if (j != i):
                    Gr[i][s] += AtElec[j]*gr[i][j]*fc[i][j]
                    for a in range(atnum):
                        if (a == i):
                            for d in range(3):
                                M_G[d][a][i][s] += AtElec[j]*(fc[i][j]*M_gr_r[i][j] + gr[i][j]*M_fc_r[i][j])*M_r_c[d][i][j]*1
                        elif (a == j):
                            for d in range(3):
                                M_G[d][a][i][s] += AtElec[j]*(fc[i][j]*M_gr_r[i][j] + gr[i][j]*M_fc_r[i][j])*M_r_c[d][i][j]*(-1)
    M_F1_x = np.zeros((3,3,3))
    # Calculate Angular function 
    for s in range(Gannum):
        ksai = 1.0
        deltaRa = 2.5/((Gannum/2)-1)
        if (0<=s<(Gannum/2)):
            lmbda = 1.0
            Ms = 0.5+deltaRa*s
        else:
            lmbda = -1.0
            Ms = 0.5+deltaRa*(s-(Gannum/2))
        F1 = np.power((1+lmbda*Costheta),ksai)
        F2b = np.power((Dis-Ms),2)
        MF1_cos = lmbda*ksai*np.power((1+lmbda*Costheta),(ksai-1))
        Ds = Dis-Ms
        for i in range(atnum):
            for j in range(atnum):
                for k in range(atnum):
                    if (i != j and i != k and j < k):
                        F2[i][j][k] = np.exp(-eta*(F2b[i][j]+F2b[i][k]+F2b[j][k]))
                        F3[i][j][k] = fc[i][j]*fc[i][k]*fc[j][k]
                        Ga[i][s] += AtElec[j]*AtElec[k]*pow(2,(1-ksai))*F1[i][j][k]*F2[i][j][k]*F3[i][j][k]
                        for a in range(atnum):
                            sa = s+Granum
                            for d in range(3):
                                if (a == i):
                                    count = 0
                                    count += F2[i][j][k]*F3[i][j][k]*MF1_cos[i][j][k]*(M_cos_ij[d][i][j][k]+M_cos_ik[d][i][j][k]) 
                                    count += F1[i][j][k]*F3[i][j][k]*F2[i][j][k]*(-2*eta)*(Ds[i][j]*M_r_c[d][i][j]+Ds[i][k]*M_r_c[d][i][k])
                                    count += F1[i][j][k]*F2[i][j][k]*(fc[i][k]*fc[j][k]*M_fc_r[i][j]*M_r_c[d][i][j]+fc[i][j]*fc[j][k]*M_fc_r[i][k]*M_r_c[d][i][k])
                                    count = count*AtElec[j]*AtElec[k]
                                    M_G[d][a][i][sa] += count
                                elif ( a == j ):
                                    count = 0
                                    count += F2[i][j][k]*F3[i][j][k]*MF1_cos[i][j][k]*(M_cos_ij[d][i][j][k])*(-1) 
                                    count += F1[i][j][k]*F3[i][j][k]*F2[i][j][k]*(-2*eta)*(Ds[i][j]*M_r_c[d][i][j]*(-1)+Ds[j][k]*M_r_c[d][j][k])
                                    count += F1[i][j][k]*F2[i][j][k]*(fc[i][k]*fc[j][k]*M_r_c[d][i][j]*M_fc_r[i][j]*(-1)+fc[i][j]*fc[i][k]*M_r_c[d][j][k]*M_fc_r[j][k])
                                    count = count*AtElec[j]*AtElec[k]
                                    M_G[d][a][i][sa] += count
                                elif ( a == k ):
                                    count = 0
                                    count += F2[i][j][k]*F3[i][j][k]*MF1_cos[i][j][k]*(M_cos_ik[d][i][j][k])*(-1) #part 1
                                    count += F1[i][j][k]*F3[i][j][k]*F2[i][j][k]*(-2*eta)*(Ds[i][k]*M_r_c[d][i][k]*(-1)+Ds[j][k]*M_r_c[d][j][k]*(-1))
                                    count += F1[i][j][k]*F2[i][j][k]*(fc[i][j]*fc[j][k]*M_r_c[d][i][k]*M_fc_r[i][k]*(-1)+fc[i][j]*fc[i][k]*M_r_c[d][j][k]*M_fc_r[j][k]*(-1))
                                    count = count*AtElec[j]*AtElec[k]
                                    M_G[d][a][i][sa] += count
    G = np.hstack((Gr,Ga))
    return G,M_G
