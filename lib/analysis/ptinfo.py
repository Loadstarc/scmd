import sys,os
import numpy as np
import matplotlib.pyplot as plt

bdinfo = [[1,2],[1,3],[1,4],[1,5],[5,6],[5,7],[5,8],[8,9],\
        [8,10],[8,11],[11,12],[11,14],[14,15],[14,16]]
dihinfo = [[4,1,2,3],[5,1,2,3],[6,5,1,2],[7,5,1,2],[8,5,1,2],[9,8,5,1]\
        ,[10,8,5,1],[11,8,5,1],[12,11,8,5],[13,11,8,5],[14,11,8,5],[15,14,11,8]\
        ,[16,14,11,8]]


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

def print_bdlen(co,bout):
    for n,bond in enumerate(bdinfo):
        bond_len = np.linalg.norm(co[bond[0]-1]-co[bond[1]-1])
        bout[n].write('%.9f\n'%(bond_len))

def print_dih(co,dout):
    for n,dih in enumerate(dihinfo):
        dih_angle = calc_torsion_angle(co,dih)
        dout[n].write('%.9f\n'%(dih_angle))

def get_pos(filename):
    traj = open(filename+'.zzip','r').readlines()
    cline = []
    for c,line in enumerate(traj):
        if ('The coordinate' in line): cline.append(c+4)
    tco = []
    for i in cline:
        co = []
        for l in traj[i:i+16]:
            co.append(l.split()[1:])
        co = np.array(co,dtype=float)
        tco.append(co)
    return tco

def ptinfo(filename,step):
    tco = get_pos(filename)[:step]
    os.mkdir('info_'+filename)
    os.chdir('info_'+filename)
    os.system('chk_conservation ../%s.zzip > check_energy'%(filename))
    energy = np.loadtxt('check_energy',dtype=float,comments='#')[:step,:]
    for i in range(1,4):
        energy[:,i] = (energy[:,i] - np.min(energy[:,i])) * 627.5
#
    plt.plot(energy[:,1],label='Ek')
    plt.plot(energy[:,2],label='Ep')
    plt.plot(energy[:,3],label='Etot')
    plt.xlabel('Step')
    plt.ylabel('Energy(kcal/mol)')
    plt.legend(loc='upper left')
    plt.savefig('energy.png')
    plt.close()
#
    plt.plot(energy[:,1])
    plt.title('kinectic energy')
    plt.xlabel('Step')
    plt.ylabel('Energy(kcal/mol)')
    plt.savefig('Ek.png')
    plt.close()
#
    plt.plot(energy[:,2])
    plt.title('potential energy')
    plt.xlabel('Step')
    plt.ylabel('Energy(kcal/mol)')
    plt.savefig('Ep.png')
    plt.close()
#
    plt.plot(energy[:,3])
    plt.title('total energy')
    plt.xlabel('Step')
    plt.ylabel('Energy(kcal/mol)')
    plt.savefig('Etot.png')
    plt.close()
#
    os.mkdir('bondlenth')
    os.mkdir('dihedrals')
    bout = []
    dout = []
    for bond in bdinfo:
        bout.append(open('bondlenth/%i-%i.dat'%(bond[0],bond[1]),'w'))
    for dih in dihinfo:
        dout.append(open('dihedrals/%i-%i-%i-%i.dat'%(dih[0],dih[1],dih[2],dih[3]),'w'))
    for co in tco:
        print_bdlen(co,bout)
        print_dih(co,dout)
    for i in range(len(bout)):
        bout[i].close()
        bond = bdinfo[i]
        data = np.loadtxt('bondlenth/%i-%i.dat'%(bond[0],bond[1]))
        plt.plot(data)
        plt.title('%i-%i bond length'%(bond[0],bond[1]))
        plt.xlabel('Step')
        plt.ylabel('Bond length(Angstrom)')
        plt.savefig('bondlenth/%i-%i.png'%(bond[0],bond[1]))
        plt.close()
    for i in range(len(dout)):
        dout[i].close()
        dih = dihinfo[i]
        data = np.loadtxt('dihedrals/%i-%i-%i-%i.dat'%(dih[0],dih[1],dih[2],dih[3]))
        plt.plot(data)
        plt.title('%i-%i-%i-%i dihdrals'%(dih[0],dih[1],dih[2],dih[3]))
        plt.xlabel('Step')
        plt.ylabel('dihedral(angle)')
        plt.savefig('dihedrals/%i-%i-%i-%i.png'%(dih[0],dih[1],dih[2],dih[3]))
        plt.close()



