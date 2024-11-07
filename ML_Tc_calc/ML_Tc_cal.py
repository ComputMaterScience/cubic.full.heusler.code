# Created by Thi Ho ver. 1.5 - 2022/01/13

# import library
import os
import shutil
import numpy as np
import re


# setup OpenMX and Vampire
openmx_link = '/home/user2/thiho/app/openmx3.9'
vampire_link = '/home/user2/thiho/app/vampire'

vps_pao = ['H_PBE19', 'He_PBE19', 'Li_PBE19', 'Be_PBE19', 'B_PBE19', 'C_PBE19', 'N_PBE19', 'O_PBE19',\
      'F_PBE19', 'Ne_PBE19', 'Na_PBE19', 'Mg_PBE19', 'Al_PBE19', 'Si_PBE19', 'P_PBE19', 'S_PBE19',\
      'Cl_PBE19', 'Ar_PBE19', 'K_PBE19', 'Ca_PBE19', 'Sc_PBE19', 'Ti_PBE19', 'V_PBE19', 'Cr_PBE19',\
      'Mn_PBE19', 'Fe_PBE19H', 'Co_PBE19H',	'Ni_PBE19H', 'Cu_PBE19H', 'Zn_PBE19H', 'Ga_PBE19',\
      'Ge_PBE19', 'As_PBE19', 'Se_PBE19', 'Br_PBE19', 'Kr_PBE19', 'Rb_PBE19', 'Sr_PBE19', 'Y_PBE19',\
      'Zr_PBE19', 'Nb_PBE19', 'Mo_PBE19', 'Tc_PBE19', 'Ru_PBE19', 'Rh_PBE19', 'Pd_PBE19', 'Ag_PBE19',\
      'Cd_PBE19', 'In_PBE19', 'Sn_PBE19', 'Sb_PBE19', 'Te_PBE19', 'I_PBE19', 'Xe_PBE19', 'Cs_PBE19',\
      'Ba_PBE19', 'La_PBE19', 'Ce_PBE19', 'Pr_PBE19', 'Nd_PBE19', 'Pm_PBE19', 'Sm_PBE19', 'Dy_PBE19',\
      'Ho_PBE19', 'Lu_PBE19', 'Hf_PBE19', 'Ta_PBE19', 'W_PBE19', 'Re_PBE19', 'Os_PBE19', 'Ir_PBE19',\
      'Pt_PBE19', 'Au_PBE19', 'Hg_PBE19', 'Tl_PBE19', 'Pb_PBE19', 'Bi_PBE19']

std_pao = ['H6.0-s2p1', 'He8.0-s2p1', 'Li8.0-s3p2', 'Be7.0-s2p2', 'B7.0-s2p2d1', 'C6.0-s2p2d1', 'N6.0-s2p2d1',\
       'O6.0-s2p2d1', 'F6.0-s2p2d1', 'Ne9.0-s2p2d1', 'Na9.0-s3p2d1', 'Mg9.0-s3p2d1', 'Al7.0-s2p2d1',\
       'Si7.0-s2p2d1', 'P7.0-s2p2d1f1', 'S7.0-s2p2d1f1', 'Cl7.0-s2p2d1f1', 'Ar9.0-s2p2d1f1', 'K10.0-s3p2d1',\
       'Ca9.0-s3p2d1', 'Sc9.0-s3p2d1', 'Ti7.0-s3p2d1', 'V6.0-s3p2d1', 'Cr6.0-s3p2d1', 'Mn6.0-s3p2d1', 'Fe5.5H-s3p2d1',\
       'Co6.0H-s3p2d1', 'Ni6.0H-s3p2d1', 'Cu6.0H-s3p2d1', 'Zn6.0H-s3p2d1', 'Ga7.0-s3p2d2', 'Ge7.0-s3p2d2',\
       'As7.0-s3p2d2', 'Se7.0-s3p2d2', 'Br7.0-s3p2d2', 'Kr10.0-s3p2d2', 'Rb11.0-s3p2d2', 'Sr10.0-s3p2d2',\
       'Y10.0-s3p2d2', 'Zr7.0-s3p2d2', 'Nb7.0-s3p2d2', 'Mo7.0-s3p2d2', 'Tc7.0-s3p2d2', 'Ru7.0-s3p2d2',\
       'Rh7.0-s3p2d2', 'Pd7.0-s3p2d2', 'Ag7.0-s3p2d2', 'Cd7.0-s3p2d2', 'In7.0-s3p2d2', 'Sn7.0-s3p2d2',\
       'Sb7.0-s3p2d2', 'Te7.0-s3p2d2f1', 'I7.0-s3p2d2f1', 'Xe11.0-s3p2d2', 'Cs12.0-s3p2d2', 'Ba10.0-s3p2d2',\
       'La8.0-s3p2d2f1', 'Ce8.0-s3p2d2f1', 'Pr8.0-s3p2d2f1', 'Nd8.0-s3p2d2f1', 'Pm8.0-s3p2d2f1', 'Sm8.0-s3p2d2f1',\
       'Dy8.0-s3p2d2f1', 'Ho8.0-s3p2d2f1', 'Lu8.0-s3p2d2f1', 'Hf9.0-s3p2d2f1', 'Ta7.0-s3p2d2f1', 'W7.0-s3p2d2f1',\
       'Re7.0-s3p2d2f1', 'Os7.0-s3p2d2f1', 'Ir7.0-s3p2d2f1', 'Pt7.0-s3p2d2f1', 'Au7.0-s3p2d2f1', 'Hg8.0-s3p2d2f1',\
       'Tl8.0-s3p2d2f1', 'Pb8.0-s3p2d2f1', 'Bi8.0-s3p2d2f1']

val_pao = [1, 2, 3, 2, 3, 4, 5, 6, 7, 8, 9, 8, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,\
    	   19, 20, 13, 4, 15, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 14, 15, 16, 17, 12, 13, 14, 15, 16, 7, 8,\
           9, 10, 11, 12, 13, 14, 15, 16, 20, 21, 11, 12, 13, 12, 15, 14, 15, 16, 17, 18, 19, 14, 15]

# constant
meV2J = 1.60217662e-22

# override variables
use_kpoint = True
mag_thres = 0.1 # threshold for magnetic atoms

use_mca = False
mca_const = [9.0e-24, 9.0e-24, 9.0e-24, 9.0e-24]
mca_direct = [[0.0,0.0,1.0], [0.0,0.0,-1.0], [0.0,0.0,1.0], [0.0,0.0,-1.0]]

use_ucell = True
ucell = [1,1,1]

# openMX
cutoff = 400.0
scf_step = 200
scf_conv = 3.6749308136648877e-9
kpoints_openmx = [31, 31, 31]

# jx code
rcut = 8.0
num_poles = 120
kpoints_jx = [31, 31, 31]
use_bunch = True
bunch_pair = 52

# vampire
system_size_x = 3.0
system_size_y = 3.0
system_size_z = 3.0

max_temp = 1500
min_temp = 0

damp_const = 0.5

use_field_cool = True

use_random_mag = True

################################################# IO Functions ###########################################################

def get_folders(folder = './'): # get list of folder names
    return [name for name in os.listdir(folder) if os.path.isdir(os.path.join(folder, name))]

def read_checkpoint(path,value): # checkpoint of calculation
    if os.path.isfile(path+'/checkpoint'):
        with open(path+'/checkpoint', 'r') as f:
            s = int(f.readline())
            if s == value:
                return True
            else:
                return False
    else:
        return True

def write_checkpoint(path,value):
    with open(path+'/checkpoint', 'w') as f:
        f.write(str(value))

def read_poscar(path):
    with open(path+'/step1/POSCAR', 'r') as f:
        f.readline() #ignore comment line
        a0 = float(f.readline())
        a  = [float(x)*a0 for x in f.readline().split()]
        b  = [float(x)*a0 for x in f.readline().split()]
        c  = [float(x)*a0 for x in f.readline().split()]
        lat = [a, b, c]
        atom_label = [x for x in f.readline().split()]
        atom_num = [int(x) for x in f.readline().split()]
        f.readline() # ignore line
        atom_pos = list()
        for i in range(len(atom_label)):
            tmp_pos = list()
            for j in range(atom_num[i]):
                s = f.readline().split()
                tmp_pos.append([float(s[0]), float(s[1]), float(s[2])])
            atom_pos.append(tmp_pos)
        return lat, atom_label, atom_num, atom_pos

def read_outcar(path, atom_label, atom_num):
    with open(path+'/step1/OUTCAR', 'r') as f:
        for line in f:
            if (line.find('magnetization (x)') != -1): # get magnetic moments
                mag = list()
                for i in range(3): f.readline() # ignore lines
                for i in range(len(atom_label)):
                    mag_type = list()
                    for j in range(atom_num[i]):
                        s = f.readline()
                        s = s.split()
                        s = [float(i) for i in s]
                        if np.abs(s[4]) > mag_thres: # magnetic atom
                            mag_type.append([int(s[0]), s[4], 1])
                        else: # non-magnetic atom
                            mag_type.append([int(s[0]), s[4], 0])
                    mag.append(mag_type)
    return mag

def read_kpoints(path):
    with open(path+'/step1/KPOINTS', 'r') as f:
        for i in range(3): f.readline() # ignore comment line
        return [int(x) for x in f.readline().split()]

def get_element_pao(atom_label):
    for i in range(len(vps_pao)):
        if (vps_pao[i].find(atom_label) != -1):
            break
    return ' {0}  {1}   {2}'.format(atom_label,std_pao[i],vps_pao[i])

def get_val_pao(atom_label, mag):
    for i in range(len(vps_pao)):
        if (vps_pao[i].find(atom_label) != -1):
            break
    return '{0}  {1}'.format((val_pao[i]+mag)/2,(val_pao[i]-mag)/2)

######################################################################################################################
# functions for step2

def check_mag(path): # check system magnetic or nonmagnetic
    # read data
    lat, atom_label, atom_num, atom_pos = read_poscar(path)
    mag = read_outcar(path,atom_label,atom_num)
    num_atom_mag = 0
    for i in range(len(atom_label)):
        for j in range(atom_num[i]):
            if mag[i][j][2] == 1:
                return True # run calculations
    return False # stop calculation

def write_step2_datfile(path): # write .dat file for step2
    # read data
    lat, atom_label, atom_num, atom_pos = read_poscar(path)
    mag = read_outcar(path,atom_label,atom_num)
    if (not use_kpoint): 
        kpoint = read_kpoints(path)
    else:
        kpoint = kpoints_openmx
    # open file to write data
    with open(path+'/step2/input.dat', 'w') as f1:
        f1.write('System.CurrrentDirectory         ./    # default=./\n')
        f1.write('DATA.PATH                      {0}/DFT_DATA19     # default=../DFT_DATA19\n'.format(openmx_link))
        f1.write('System.Name                  {0}\n'.format(path))
        f1.write('level.of.stdout                   1    # default=1 (1-3)\n')
        f1.write('level.of.fileout                  1    # default=1 (1-3)\n')
        f1.write('\n')
        f1.write('HS.fileout                        on   # on|off, default=off\n')
        f1.write('\n')
        f1.write('Species.Number  {0}\n'.format(len(atom_label)))
        f1.write('<Definition.of.Atomic.Species\n')
        for i in atom_label:
            f1.write('{0}\n'.format(get_element_pao(i)))
        f1.write('Definition.of.Atomic.Species>\n')
        f1.write('\n')
        f1.write('Atoms.Number   {0}\n'.format(np.sum(atom_num)))
        f1.write('Atoms.SpeciesAndCoordinates.Unit  FRAC\n')
        f1.write('<Atoms.SpeciesAndCoordinates\n')
        for i in range(len(atom_label)):
            for j in range(atom_num[i]):
                f1.write('{0}  {1}    {2}    {3}    {4}  {5}\n'.format(mag[i][j][0], atom_label[i],\
                    atom_pos[i][j][0], atom_pos[i][j][1], atom_pos[i][j][2], get_val_pao(atom_label[i], mag[i][j][1])))
        f1.write('Atoms.SpeciesAndCoordinates>\n')
        f1.write('\n')
        f1.write('Atoms.UnitVectors.Unit  Ang\n')
        f1.write('<Atoms.UnitVectors\n')
        for i in range(3):
            f1.write('{0} {1} {2}\n'.format(lat[i][0], lat[i][1], lat[i][2]))
        f1.write('Atoms.UnitVectors>\n')
        f1.write('\n')
        f1.write('scf.XcType                 GGA-PBE                    # LDA|LSDA-CA|LSDA-PW\n')
        f1.write('scf.SpinPolarization        on                        # On|Off\n')
        f1.write('scf.energycutoff           {0}                      # default=150 (Ry)\n'.format(cutoff))
        f1.write('scf.maxIter                 {0}                       # default=40\n'.format(scf_step))
        f1.write('scf.EigenvalueSolver      band                        # Recursion|Cluster|Band\n')
        f1.write('scf.Kgrid                 {0} {1} {2}\n'.format(kpoint[0],kpoint[1],kpoint[2]))
        f1.write('scf.criterion             {0}       # default=1.0e-6 (Hartree)\n'.format(scf_conv))
        f1.write('\n')
        # Mixing
        f1.write('scf.Mixing.Type            rmm-diish   # Simple|Rmm-Diis|Gr-Pulay\n')
        f1.write('scf.Init.Mixing.Weight     0.0100       # default=0.30\n')
        f1.write('scf.Min.Mixing.Weight      0.0001       # default=0.001\n')
        f1.write('scf.Max.Mixing.Weight      0.3000       # default=0.40\n')
        f1.write('scf.Mixing.History         30          # default=5\n')
        f1.write('scf.Mixing.StartPulay      40          # default=6\n')
        f1.write('\n')

def run_step2(path):
    # check checkpoint
    if read_checkpoint(path+'/step2',0) == True:
        os.makedirs(path+'/step2',exist_ok=True)
        write_checkpoint(path+'/step2',0)
        # make .dat file for openMX
        write_step2_datfile(path)
    
        # --------> run openMX
        os.system('cd ' + path + '/step2 && srun '+ openmx_link +'/work/openmx input.dat > stdout && cd .. && cd ..')
        # ------------------

        write_checkpoint(path+'/step2',1)

######################################################################################################################

# functions for step3

def read_openmx_output(path):
    # read data
    lat, atom_label, atom_num, atom_pos = read_poscar(path)
    with open(path+'/step2/{0}.out'.format(path), 'r') as f:
        for line in f:
            if (line.find('Total spin moment') != -1): # get magnetic moments
                mag = list()
                for i in range(2): f.readline() # ignore lines
                for i in range(len(atom_label)):
                    mag_type = list()
                    for j in range(atom_num[i]):
                        s = f.readline().split()
                        s = float(s[5])
                        if (np.abs(s) > mag_thres): # magnetic element
                            mag_type.append([s, 1])
                        else:
                            mag_type.append([s, 0])
                    mag.append(mag_type)
    return mag

def write_jx_config(path):
    # read data
    lat, atom_label, atom_num, atom_pos = read_poscar(path)
    if (not use_kpoint): 
        kpoint = read_kpoints(path)
    else:
        kpoint = kpoints_jx

    mag = read_openmx_output(path)
    # convert pos to cartesian coordinate
    atom_pos_cart = list()
    atom_pos_mag = list()
    for i in range(len(atom_label)):
        for j in range(atom_num[i]):
            atom_pos_cart.append([atom_pos[i][j][0]*lat[0][0]+atom_pos[i][j][1]*lat[1][0]+atom_pos[i][j][2]*lat[2][0],\
                atom_pos[i][j][0]*lat[0][1]+atom_pos[i][j][1]*lat[1][1]+atom_pos[i][j][2]*lat[2][1],\
                    atom_pos[i][j][0]*lat[0][2]+atom_pos[i][j][1]*lat[1][2]+atom_pos[i][j][2]*lat[2][2]])
            atom_pos_mag.append(mag[i][j][1])

    # calcutate maximum number Rmax for R1, R2 and R3
    def calc_R():
        R_list = list()
        R_dis = list()
        R1_cut = int(np.ceil(rcut/np.linalg.norm(lat[0]))) + 1
        R2_cut = int(np.ceil(rcut/np.linalg.norm(lat[1]))) + 1
        R3_cut = int(np.ceil(rcut/np.linalg.norm(lat[2]))) + 1
        for i in range(len(atom_pos_cart)):
            if (atom_pos_mag[i] == 1):
                for R1 in range(-R1_cut,R1_cut+1):
                    for R2 in range(-R2_cut,R2_cut+1):
                        for R3 in range(-R3_cut,R3_cut+1):     
                            r_dis = R1*np.array(lat[0]) + R2*np.array(lat[1]) + R3*np.array(lat[2]) + np.array(atom_pos_cart[i])
                            for j in range(len(atom_pos_cart)):
                                if (atom_pos_mag[j] == 1):
                                    if (i == j) and (R1 == 0 and R2 == 0 and R3 == 0):
                                        pass
                                    else:
                                        r = np.linalg.norm(r_dis-np.array(atom_pos_cart[j]))
                                        R_list.append([i+1, j+1, R1, R2, R3])
                                        R_dis.append(r)
        R_dis = np.array(R_dis)
        id = np.where(R_dis <= rcut)
        id = id[0]
        Jij = list()
        for i in id:
            Jij.append(R_list[i])
        return Jij

    # calculate Jij matrix
    def get_Jij():
        Jij = list()
        if (not use_ucell):
            Jij = calc_R()
        else:
            Rmax = ucell
            for i in range(1, len(atom_pos_cart)+1):
                for j in range(1, len(atom_pos_cart)+1):
                    for R1 in range(-int(Rmax[0]), int(Rmax[0])+1):
                        for R2 in range(-int(Rmax[1]), int(Rmax[1])+1):
                            for R3 in range(-int(Rmax[2]), int(Rmax[2])+1):
                                if (atom_pos_mag[i-1] == 1):
                                    if (atom_pos_mag[j-1] == 1):
                                        if (i == j) and (R1 == 0 and R2 == 0 and R3 == 0):
                                            pass
                                        else:
                                            Jij.append([i,j,R1,R2,R3])
        return Jij

    Jij = get_Jij()

    num_J = len(Jij)

    with open(path+'/step3/jx.config', 'w') as f1:
        f1.write('Flag.PeriodicSum           off\n')
        f1.write('Num.Poles                  {0}\n'.format(num_poles))
        f1.write('Num.Kgrid               {0} {1} {2}\n'.format(kpoint[0], kpoint[1], kpoint[2]))
        f1.write('Num.ij.pairs	             {0}\n'.format(len(Jij)))
        if (use_bunch):
            f1.write('Bunch.ij.pairs	         {0}\n'.format(bunch_pair))
        else:
            f1.write('Bunch.ij.pairs	         {0}\n'.format(len(Jij)))
        f1.write('\n')

        f1.write('<ijpairs.cellid\n')
        for i in range(len(Jij)):
            f1.write('{0} {1} {2} {3} {4}\n'.format(Jij[i][0], Jij[i][1], Jij[i][2], Jij[i][3], Jij[i][4]))
        f1.write('ijpairs.cellid>\n')


def run_step3(path):
    # check checkpoint
    if read_checkpoint(path+'/step3',0) == True:
        os.makedirs(path+'/step3',exist_ok=True)
        write_checkpoint(path+'/step3',0)
        # copy files to step3
        shutil.copyfile(path+'/step2/{0}.scfout'.format(path),path+'/step3/{0}.scfout'.format(path))
        
        # write jx.config
        write_jx_config(path)

        # --------> run vasp
        os.system('cd ' + path + '/step3 && srun '+ openmx_link +'/work/jx {0}.scfout jx.config >jx.log  && cd .. && cd ..'.format(path))
        # ------------------

        write_checkpoint(path+'/step3',1)

######################################################################################################################

def read_jx_output(path):

    mag = read_openmx_output(path)
    lat, atom_label, atom_num, atom_pos = read_poscar(path)

    mag_list = list()

    for i in range(len(atom_label)):
        for j in range(atom_num[i]):
            mag_list.append(mag[i][j][0])

    with open(path+'/step3/jx.config'.format(path), 'r') as f:
        for i in range(4): s = f.readline() # ignore lines
        s = s.split()
    num_J = int(s[1])

    with open(path+'/step3/jx.log'.format(path), 'r') as f:
        for line in f:
            if (line.find('J [meV]') != -1): # get Jij matrix
                Jij = list()
                f.readline() # ignore lines
                for i in range(num_J):
                    s = f.readline().split()
                    Jij.append([int(s[0]), int(s[1]), int(s[2]), int(s[3]), int(s[4]), \
                    2*np.sign(mag_list[int(s[0])-1])*np.sign(mag_list[int(s[1])-1])*meV2J*float(s[5])])
    return Jij

# functions for step4
def write_vampire_input(path): # write input
    # read data
    lat, atom_label, atom_num, atom_pos = read_poscar(path)
    mag = read_openmx_output(path)
    # write input
    with open(path+'/step4/input', 'w') as f1:
        f1.write('#------------------------------------------\n')
        f1.write('# Creation attributes:\n')
        f1.write('#------------------------------------------\n')
        f1.write('create:full\n')
        f1.write('create:periodic-boundaries-x\n')
        f1.write('create:periodic-boundaries-y\n')
        f1.write('create:periodic-boundaries-z\n')
        f1.write('#------------------------------------------\n')
        f1.write('material:file=vampire.mat\n')
        f1.write('material:unit-cell-file = "vampire.UCF"\n')
        f1.write('#------------------------------------------\n')
        f1.write('# System Dimensions:\n')
        f1.write('#------------------------------------------\n')
        f1.write('dimensions:unit-cell-size-x = {0}\n'.format(np.linalg.norm(lat[0])))
        f1.write('dimensions:unit-cell-size-y = {0}\n'.format(np.linalg.norm(lat[1])))
        f1.write('dimensions:unit-cell-size-z = {0}\n'.format(np.linalg.norm(lat[2])))
        f1.write('\n')
        f1.write('dimensions:system-size-x = {0} !nm\n'.format(system_size_x))
        f1.write('dimensions:system-size-y = {0} !nm\n'.format(system_size_y))
        f1.write('dimensions:system-size-z = {0} !nm\n'.format(system_size_z))
        f1.write('#------------------------------------------\n')
        f1.write('# Simulation attributes: \n')
        f1.write('#------------------------------------------\n')
        # use field-cool or curie-temperature
        mag_sum = 0
        for i in range(len(atom_label)): 
            for j in range(atom_num[i]):
                mag_sum += mag[i][j][0]
        if (use_field_cool):
            f1.write('sim:time-step=1.0E-16\n')
            f1.write('sim:time-steps-increment=1000\n')
            f1.write('sim:cooling-time=0.1 !ns\n')
            f1.write('sim:cooling-function=linear\n')
            f1.write('sim:total-time-steps=1000000\n')
            f1.write('#sim:equilibration-applied-field-strength=1 !T\n')
            f1.write('#sim:applied-field-strength=1 !T\n')
            f1.write('sim:temperature = 0\n')
        else:
            if (mag_sum > mag_thres):
                f1.write('sim:temperature=300 \n')
                f1.write('sim:temperature-increment=25\n')
                f1.write('sim:time-steps-increment=1 \n')
                f1.write('sim:equilibration-time-steps=2500\n')
                f1.write('sim:loop-time-steps=3000\n')
            else:
                f1.write('sim:time-step=1.0E-16\n')
                f1.write('sim:time-steps-increment=1000\n')
                f1.write('sim:cooling-time=0.1 !ns\n')
                f1.write('sim:cooling-function=linear\n')
                f1.write('sim:total-time-steps=1000000\n')
                f1.write('#sim:equilibration-applied-field-strength=1 !T\n')
                f1.write('#sim:applied-field-strength=1 !T\n')
                f1.write('sim:temperature = 0\n')
        f1.write('sim:minimum-temperature={0}\n'.format(min_temp))
        f1.write('sim:maximum-temperature={0}\n'.format(max_temp))
        f1.write('#------------------------------------------\n')
        f1.write('# Program and integrator details\n')
        f1.write('#------------------------------------------\n')
        if (use_field_cool):
            f1.write('sim:program=field-cool\n')
        else:
            if (mag_sum > mag_thres):
                f1.write('sim:program=curie-temperature\n')
            else:
                f1.write('sim:program=field-cool\n')
        f1.write('sim:integrator=llg-heun\n')
        f1.write('#------------------------------------------\n')
        f1.write('# Data output\n')
        f1.write('#------------------------------------------\n')
        f1.write('output:temperature\n')
        f1.write('output:mean-magnetisation-length\n')
        f1.write('output:material-magnetisation\n')
        f1.write('screen:temperature\n')
        f1.write('screen:mean-magnetisation-length\n')
        f1.write('screen:material-magnetisation\n')

    with open(path+'/step4/vampire.mat', 'w') as f1:
        f1.write('material:num-materials = {0}\n'.format(np.sum(atom_num)))
        cnt = 0
        for i in range(len(atom_label)):
            for j in range(atom_num[i]):
                cnt += 1
                if (mag[i][j][1] == 1):
                    f1.write('#---------------------------------------------------\n')
                    f1.write('# Material {0}\n'.format(cnt))
                    f1.write('#---------------------------------------------------\n')
                    f1.write('material[{0}]:material-name={1}\n'.format(cnt, atom_label[i]))
                    f1.write('material[{0}]:damping-constant={1}\n'.format(cnt, damp_const))
                    f1.write('material[{0}]:atomic-spin-moment={1} !muB\n'.format(cnt,np.abs(mag[i][j][0])))
                    if (use_mca):
                        f1.write('material[{0}]:uniaxial-anisotropy-constant={1}\n'.format(cnt, mca_const[cnt-1]))
                    else:
                        f1.write('material[{0}]:uniaxial-anisotropy-constant=0.0\n'.format(cnt))
                    f1.write('material[{0}]:material-element={1}\n'.format(cnt, atom_label[i]))
                    if (use_field_cool):
                        f1.write('material[{0}]:initial-spin-direction = random\n'.format(cnt))
                    else:
                        if (mag_sum > mag_thres):
                            if (use_random_mag):
                                f1.write('material[{0}]:initial-spin-direction = random\n'.format(cnt))
                            else:
                                f1.write('material[{0}]:initial-spin-direction = 0.0, 0.0, {1}\n'.format(cnt, np.sign(mag[i][j][0])))
                        else:
                            f1.write('material[{0}]:initial-spin-direction = random\n'.format(cnt))
                    if (use_mca):
                        if len(mca_direct) > 0:
                            f1.write('material[{0}]:uniaxial-anisotropy-direction = {1} , {2}, {3}\n'.format(cnt, mca_direct[cnt-1][0],\
                                mca_direct[cnt-1][1], mca_direct[cnt-1][2]))
                        else:
                            f1.write('material[{0}]:uniaxial-anisotropy-direction = 0.0 , 0.0, 1.0\n'.format(cnt))
                    else:
                        f1.write('material[{0}]:uniaxial-anisotropy-direction = 0.0 , 0.0, 1.0\n'.format(cnt))
                    f1.write('#---------------------------------------------------\n')
                else:
                    f1.write('#---------------------------------------------------\n')
                    f1.write('# Material {0}\n'.format(cnt))
                    f1.write('#---------------------------------------------------\n')
                    f1.write('material[{0}]:material-name={1}\n'.format(cnt, atom_label[i]))
                    f1.write('material[{0}]:non-magnetic\n'.format(cnt))
                    f1.write('#---------------------------------------------------\n')
        f1.write('#---------------------------------------------------\n')
        f1.write('# Interactions\n')

    Jij = read_jx_output(path)

    with open(path+'/step4/vampire.UCF', 'w') as f1:
        f1.write('# Unit cell size (Angstrom):\n')
        f1.write('1 1 1\n')
        f1.write('# Unit cell lattice vectors:\n')
        for i in range(3):
            f1.write('{0} {1} {2}\n'.format(lat[i][0], lat[i][1], lat[i][2]))
        f1.write('# Atoms\n')
        f1.write('{0} {0}\n'.format(np.sum(atom_num),np.sum(atom_num)))
        cnt = 0
        for i in range(len(atom_label)):
            for j in range(atom_num[i]):
                f1.write('{0} {1} {2} {3} {4}\n'.format(cnt, atom_pos[i][j][0], atom_pos[i][j][1], atom_pos[i][j][2], cnt))
                cnt += 1
        f1.write('# Interactions\n')
        f1.write('{0} isotropic\n'.format(len(Jij)))
        for i in range(len(Jij)):
            f1.write('{0} {1} {2} {3} {4} {5} {6}\n'.format(i, Jij[i][0]-1, Jij[i][1]-1, Jij[i][2], Jij[i][3], Jij[i][4], Jij[i][5]))
        

def run_step4(path):
    # check checkpoint
    if read_checkpoint(path+'/step4',0) == True:
        os.makedirs(path+'/step4',exist_ok=True)
        write_checkpoint(path+'/step4',0)
        
        #write input files for Vampire
        write_vampire_input(path)

        # --------> run vampire
        os.system('cd ' + path + '/step4 && srun '+ vampire_link +'/vampire-parallel > stdout && cd .. && cd ..')
        # ------------------

        write_checkpoint(path+'/step4',1)

#-----------------------------
def save_list(name,data):
    with open(name, 'w') as f:
        for i in data:
            f.write(i+'\n')

def save_log(info):
    with open('output.log', 'a') as f:
        f.write(info)

def load_list(name):
    data = list()
    with open(name, 'r') as f:
        for i in f:
            data.append(i.strip())
    return data

#main program
if __name__ == '__main__':
    # get folders to run
    if os.path.isfile('structures.in') == False:
        list_struct = get_folders()
        list_struct.sort()
        save_list('structures.in',list_struct)
    else:
        list_struct = load_list('structures.in')
    # perform calculations
    for i in range(len(list_struct)):
        save_log(f'Calculating {list_struct[i]}\n')
        if check_mag(list_struct[i]): # magnetic system -> run calculation
            #step2
            run_step2(list_struct[i])
            #step3
            run_step3(list_struct[i])
            #step4
            run_step4(list_struct[i])
        #save checkpoint
        if i + 1 < len(list_struct):
            save_list('structures.in',list_struct[i+1:])
        else:
            save_list('structures.in',[])
        save_log('Done.\n')
        save_log('--------------------------------------\n')
    