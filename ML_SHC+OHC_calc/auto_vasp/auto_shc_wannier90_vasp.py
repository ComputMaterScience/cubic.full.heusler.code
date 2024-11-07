import os
import shutil
import numpy as np
import re


vasp_link = '/home/user2/thiho/app/vasp.5.4.4_wan3'
wannier90_link = '/home/user2/thiho/app/wannier90_OHC/wannier90-3.1.0'
wan4hall_link = '/home/user2/thiho/app/wannier4hall'

# Input Variables

pot_type = 1 # 1: PBE, 2: LDA

# step1
encut   = 500

# step2
kmesh_bcc_fcc_step2 = [12, 12, 12]
kmesh_hcp_step2     = [12, 12, 8]

# step3
num_iter = 0
dis_mix_ratio = 1.0
dis_num_iter = 10000
up_fermi = 5
ws_distance = False

# step4
kmesh_bcc_fcc_step4 = [100, 100, 100]
kmesh_hcp_step4     = [100, 100, 80]

# step4_tb
kmesh_bcc_fcc_step4_tb = [300, 300, 300]
kmesh_hcp_step4_tb     = [300, 300, 185]


########################### IO Functions ##############################

# def get_folders(folder = './'): # get list of folder names
#     return [name for name in os.listdir(folder) if os.path.isdir(os.path.join(folder, name))]

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

def write_checkpoint(path,value): # write checkpoint
    with open(path+'/checkpoint', 'w') as f:
        f.write(str(value))

def read_elements(): # read elements for calculation
    elements_list = list()
    with open('elements.list','r') as f:
        for s in f:
            s = s.split()
            elements_list.append([s[0].capitalize(), s[1].lower(), float(s[2]), int(s[3]), int(s[4])])
    return elements_list

def read_poscar(path): # read poscar
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

def write_step1_KPOINTS(path,kmesh): # write KPOINTS for step1
    with open(path+'/step1/KPOINTS', 'w') as f:
        f.write('Automatic Mesh # Generates Automatically the K-mesh\n')
        f.write('  0\n')
        f.write('Gamma # Automatic Distribution of the mesh\n')
        f.write(f'  {kmesh[0]}  {kmesh[1]}  {kmesh[2]}\n')

def write_potcar(path): # write POTCAR for step1
    lat, atom_label, atom_num, atom_pos = read_poscar(path)
    if (pot_type == 1): # PBE
        link = '/Potential/potpaw_PBE.54/'
    else: # LDA
        link = '/Potential/potpaw_LDA.54/'
    with open(path+'/step1/POTCAR', 'w') as f:
        for element in atom_label:
            if os.path.exists(vasp_link + link + element + '/POTCAR'): 
                pot_link = vasp_link + link + element + '/POTCAR'
            else:
                if os.path.exists(vasp_link + link + element + '_pv/POTCAR'): 
                    pot_link = vasp_link + link + element + '_pv/POTCAR'
                else:
                    pot_link = vasp_link + link + element + '_sv/POTCAR'
            with open(pot_link) as potfile:
                for line in potfile:
                    f.write(line)

def write_poscar(path, type, lat): # write POSCAR for step1
    with open(path+'/step1/POSCAR', 'w') as f:
        f.write('{0}\n'.format(path))
        f.write('   {0}\n'.format(lat))
        if (type == 'bcc'):
            f.write('     0.5000000000000000    0.5000000000000000   -0.5000000000000000\n')
            f.write('    -0.5000000000000000    0.5000000000000000    0.5000000000000000\n')
            f.write('     0.5000000000000000   -0.5000000000000000    0.5000000000000000\n')
            f.write('   {0}\n'.format(path))
            f.write('     1\n')
            f.write('Direct\n')
            f.write('  0.0000000000000000  0.0000000000000000  0.0000000000000000\n')
        elif (type == 'fcc'):
            f.write('     0.0000000000000000    0.5000000000000000    0.5000000000000000\n')
            f.write('     0.5000000000000000    0.0000000000000000    0.5000000000000000\n')
            f.write('     0.5000000000000000    0.5000000000000000    0.0000000000000000\n')
            f.write('   {0}\n'.format(path))
            f.write('     1\n')
            f.write('Direct\n')
            f.write('  0.0000000000000000  0.0000000000000000  0.0000000000000000\n')
        else:
            f.write('     1.0000000000000000    0.0000000000000000    0.0000000000000000\n')
            f.write('    -0.5000000000000000    0.8660254037844390    0.0000000000000000\n')
            f.write('     0.0000000000000000    0.0000000000000000    1.6200000000000000\n')
            f.write('   {0}\n'.format(path))
            f.write('     2\n')
            f.write('Direct\n')
            f.write('  0.6666666666666666  0.3333333333333333  0.7500000000000000\n')
            f.write('  0.3333333333333333  0.6666666666666666  0.2500000000000000\n')

# functions for step1
def write_step1_INCAR(path, opt=False): # write INCAR for step1
    lat, atom_label, atom_num, atom_pos = read_poscar(path)
    with open(path+'/step1/INCAR', 'w') as f:
        f.write('# SCF Job\n')
        f.write('#===========================================\n\n')
        f.write('# General Setup\n')
        f.write('  System = OPT Job       	# Calculation Title\n')
        f.write('  ENCUT  = {0}       		# Kinetic Energy Cutoff in eV\n'.format(encut))
        f.write('  ISTART = 0         	    # Job: 0-new  1-cont  2-samecut\n')
        f.write('  ICHARG = 2               # initial charge density: 1-file 2-atom 10-cons 11-DOS\n')
        f.write('  ISPIN  = 2\n')
        f.write('  ISIF   = 3\n')
        if (opt):
            f.write('  NSW    = 200\n')
            f.write('  IBRION = 2\n')
        else:
            f.write('  NSW    =  0\n')
            f.write('  IBRION = -1\n')
        f.write('  POTIM  = 0.4\n')
        f.write('  MAGMOM = {0}*2\n'.format(np.sum(atom_num)))
        f.write('  LORBIT = 11\n')
        f.write('  LWAVE   = .FALSE.\n')
        f.write('  ADDGRID = .TRUE.\n\n')
        f.write('# Electronic Relaxation (SCF)\n')
        f.write('  NELMDL =     -6      # Number of delayed ELM steps\n')
        f.write('  NELM   =     300     # Number of ELM steps\n')
        f.write('  EDIFF  =     1.0E-7  # Stopping criteria for ESC\n')
        f.write('  EDIFFG =     1.0E-5\n')
        f.write('  ALGO   =     Fast    # Electronic algorithm minimization\n')
        f.write('  SIGMA  =     0.050000     # Broadening in eVb\n')
        f.write('  ISMEAR = 0           # Partial Occupancies for Each Orbital\n')
        f.write('  NEDOS  = 3000\n')
        f.write('# Parallelization Scheme\n')
        f.write('  LPLANE = .TRUE.\n')
        f.write('  NPAR   = 4\n')

def run_step1(path, type, lat):
    # check checkpoint
    if read_checkpoint(path+'/step1',0) == True:
        os.makedirs(path+'/step1',exist_ok=True)
        write_checkpoint(path+'/step1',0)

        # relax structure

        # make POSCAR
        write_poscar(path, type, lat)
        # make INCAR
        write_step1_INCAR(path,opt=True)
        # make KPOINTS
        if (type == 'hcp'):
            write_step1_KPOINTS(path,[15,15,9])
        else:
            write_step1_KPOINTS(path,[15,15,15])
        # make POTCAR
        write_potcar(path)
        # --------> run vasp
        os.system('cd ' + path + '/step1 && srun '+ vasp_link +'/bin/vasp_std > stdout && cd .. && cd ..')
        # ------------------

        write_checkpoint(path+'/step1', 1)

    if read_checkpoint(path+'/step1', 1) == True:

        # do scf

        # make INCAR
        write_step1_INCAR(path,opt=False)
        # make KPOINTS
        if (type == 'hcp'):
            write_step1_KPOINTS(path,[25,25,15])
        else:
            write_step1_KPOINTS(path,[25,25,25])
        # copy CONTCAR to POSCAR
        shutil.copyfile(path+'/step1/CONTCAR',path+'/step1/POSCAR')

        # --------> run vasp
        os.system('cd ' + path + '/step1 && srun '+ vasp_link +'/bin/vasp_std > stdout && cd .. && cd ..')
        # ------------------

        write_checkpoint(path+'/step1', 2)

# functions for step2
def write_step2_INCAR(path, proj, wannier90 = False): # write INCAR for step2
    tag_ISYM = False
    lat, atom_label, atom_num, atom_pos = read_poscar(path)
    # read INCAR from step1
    with open(path+'/step1/INCAR', 'r') as f1:
        # modify INCAR and save to step2
        with open(path+'/step2/INCAR', 'w') as f2:
            for s in f1:
                if re.search('MAGMOM',s): continue
                if re.search('ISTART',s):
                    if (wannier90):
                        f2.write('  ISTART = 1         	    # Job: 0-new  1-cont  2-samecut\n')
                    else:
                        f2.write('  ISTART = 0         	    # Job: 0-new  1-cont  2-samecut\n')
                elif re.search('ICHARG',s):
                    f2.write('  ICHARG = 11               # initial charge density: 1-file 2-atom 10-cons 11-DOS\n')
                elif re.search('ISYM',s):
                    f2.write('  ISYM   = -1\n')
                    tag_ISYM = True
                elif re.search('LWAVE',s):
                    f2.write('  LWAVE   = .TRUE.\n')
                elif re.search('NELM ',s):
                    if (wannier90):
                        f2.write('  NELM   =     0     # Number of ELM steps\n')
                    else:
                        f2.write('  NELM   =     200     # Number of ELM steps\n')
                elif re.search('ALGO',s):
                    if (wannier90):
                        f2.write('  ALGO   =     None    # Electronic algorithm minimization\n')
                    else:
                        f2.write('  ALGO   =     Fast    # Electronic algorithm minimization\n')
                elif re.search('NPAR',s):
                    if (wannier90):
                        f2.write('  NPAR   = 1\n')
                    else:
                        f2.write('  NPAR   = 4\n')
                else:
                    f2.write(s)

    #append tags to INCAR
    with open(path+'/step2/INCAR', 'a') as f2:
        f2.write('\n')
        f2.write('SOC calculation:\n')
        f2.write('  LMAXMIX       = 4\n')
        f2.write('  SAXIS         = 0 0 1\n')
        f2.write('  LNONCOLLINEAR = .TRUE.\n')
        f2.write('  LSORBIT       = .TRUE.\n')
        f2.write('  GGA_COMPAT    = F\n')
        f2.write('  LORBMOM       = T\n')
        if (tag_ISYM == False):
            if (proj == 1): 
                if (np.sum(atom_num) == 1):
                    f2.write('  NBANDS        = 36\n')
                else:
                    f2.write('  NBANDS        = 72\n')
            else:
                if (np.sum(atom_num) == 1):
                    f2.write('  NBANDS        = 16\n')
                else:
                    f2.write('  NBANDS        = 32\n')
        f2.write('  ISYM          = -1\n')

        if (wannier90 == True):
            f2.write('\n')
            f2.write('Wannier90:\n')
            f2.write('  LWANNIER90   =  True\n')
            f2.write('  LWRITE_SPN  =  True\n')

def write_step2_KPOINTS(path,kmesh): # write KPOINTS for step2
    with open(path+'/step2/KPOINTS', 'w') as f:
        f.write('Automatic Mesh # Generates Automatically the K-mesh\n')
        f.write('  0\n')
        f.write('Gamma # Automatic Distribution of the mesh\n')
        f.write(f'  {kmesh[0]}  {kmesh[1]}  {kmesh[2]}\n')

def read_outcar(path):
    with open(path+'/step1/OUTCAR', 'r') as f:
        for line in f:
            if (line.find('E-fermi') != -1): # get fermi energy
                s = line.split(':')
                s = s[1].split()
                efermi = float(s[0])
    return efermi

def read_kpoints(path):
    with open(path+'/step2/KPOINTS', 'r') as f:
        for i in range(3): f.readline() # ignore comment line
        return [int(x) for x in f.readline().split()]

def write_step2_winfile(path, type, proj): # write wannier90.in for step2
    # read data
    lat, atom_label, atom_num, atom_pos = read_poscar(path)
    efermi = read_outcar(path)
    kpoints = read_kpoints(path)
    with open(path+'/step2/wannier90.win', 'w') as f:
        # lattice vectors
        f.write('begin unit_cell_cart\n')
        f.write('ang\n')
        for vec in lat:
            f.write(f'  {vec[0]}  {vec[1]}  {vec[2]}\n')
        f.write('end unit_cell_cart\n')
        f.write('\n')

        # positions
        f.write('begin atoms_frac\n')
        for i in range(len(atom_label)):
            for j in range(atom_num[i]):
                f.write(f'{atom_label[i]}   {atom_pos[i][j][0]}   {atom_pos[i][j][1]}   {atom_pos[i][j][2]}\n')
        f.write('end atoms_frac\n')
        f.write('\n')

        # use spinor tag
        f.write('spinors    =  true\n')
        f.write('\n')

        # projections
        f.write('begin projections\n')
        if proj == 1:
            f.write(f'{atom_label[0]} : s, p, d\n')
        else:
            f.write(f'{atom_label[0]} : s, p\n')
        f.write('end projections\n')
        f.write('\n')

        # num_wann and num_bands tags
        if (proj == 1):
            if (np.sum(atom_num) == 1):
                f.write('num_bands  =  36\n')
                f.write('num_wann   =  18\n')
            else:
                f.write('num_bands  =  72\n')
                f.write('num_wann   =  36\n')
            f.write('\n')
        else:
            if (np.sum(atom_num) == 1):
                f.write('num_bands  =  16\n')
                f.write('num_wann   =  8\n')
            else:
                f.write('num_bands  =  32\n')
                f.write('num_wann   =  16\n')
            f.write('\n')

        # kpath tags
        f.write('#kpath                   = true\n')
        f.write('#kpath_task              = bands\n')
        f.write('#kpath_num_points        = 400\n')
        f.write('\n')

        # fermi energy tag
        f.write(f'fermi_energy = {efermi}\n')
        f.write('\n')

        # kpoint_path
        f.write('begin kpoint_path\n')
        if (type == 'bcc'):
            f.write('G  0.00000  0.00000  0.00000    H -0.50000  0.50000  0.50000\n')
            f.write('H -0.50000  0.50000  0.50000    N  0.00000  0.50000  0.00000\n')
            f.write('N  0.00000  0.50000  0.00000    P  0.25000  0.25000  0.25000\n')
            f.write('P  0.25000  0.25000  0.25000    G  0.00000  0.00000  0.00000\n')
            f.write('G  0.00000  0.00000  0.00000    N  0.00000  0.50000  0.00000\n')
        elif (type == 'fcc'):
            f.write('K  0.37500  0.75000  0.37500    G  0.00000  0.00000  0.00000\n')
            f.write('G  0.00000  0.00000  0.00000    L  0.50000  0.50000  0.50000\n')
            f.write('L  0.50000  0.50000  0.50000    W  0.25000  0.75000  0.50000\n')
            f.write('W  0.25000  0.75000  0.50000    X  0.00000  0.50000  0.50000\n')
            f.write('X  0.00000  0.50000  0.50000    G  0.00000  0.00000  0.00000\n')
        else:
            f.write('G  0.00000  0.00000  0.00000    M  0.00000  0.50000  0.00000\n')
            f.write('M  0.00000  0.50000  0.00000    K  0.33333  0.33333  0.00000\n')
            f.write('K  0.33333  0.33333  0.00000    G  0.00000  0.00000  0.00000\n')
            f.write('G  0.00000  0.00000  0.00000    A  0.00000  0.00000  0.50000\n')
            f.write('A  0.00000  0.00000  0.50000    L  0.00000  0.50000  0.50000\n')
            f.write('L  0.00000  0.50000  0.50000    H  0.33333  0.33333  0.50000\n')
            f.write('H  0.33333  0.33333  0.50000    A  0.00000  0.00000  0.50000\n')
        f.write('end kpoint_path\n')
        f.write('\n')

        # kslice tags
        f.write('#kslice         = true\n')
        f.write('#kslice_task    = shc+fermi_lines\n')
        f.write('#kslice_2dkmesh = 400\n')
        f.write('#kslice_corner  = 0.0  0.0  0.0\n')
        f.write('#kslice_b1      = 1.00  0.00  0.00\n')
        f.write('#kslice_b2      = 0.00  1.00  0.00\n')
        f.write('\n')

        # berry tags
        f.write('#berry                        = true\n')
        f.write('#berry_task                   = eval_ahc\n')
        f.write('#berry_kmesh                  =  100\n')
        f.write('#berry_curv_unit              = ang2\n')
        f.write('#berry_curv_adpt_kmesh        = 5\n')
        f.write('#berry_curv_adpt_kmesh_thresh = 100.0\n')
        f.write('\n')

        # fermi scan tags
        f.write('#fermi_energy_min  = 6\n')
        f.write('#fermi_energy_max  = 26\n')
        f.write('#fermi_energy_step = 0.1\n')
        f.write('\n')

        # wannierise tags
        f.write('search_shells = 200\n')
        f.write('trial_step = 1.0\n')
        f.write('\n')

        f.write('#dis_win_min       =   0.0\n')
        f.write('#dis_win_max       =  60.0\n')
        f.write('#dis_froz_min      =  0.0\n')
        if (proj == 1):
            f.write(f'dis_froz_max      =  {efermi+up_fermi+2}\n')
        else:
            f.write(f'dis_froz_max      =  {efermi+up_fermi}\n')
        f.write('dis_num_iter       =  {0}\n'.format(dis_num_iter))
        f.write('dis_conv_tol       = 1.0e-10\n')
        f.write('conv_tol           = 1.0e-10\n')
        f.write('conv_window        = 10\n')
        f.write('num_iter           = {0}\n'.format(num_iter))
        f.write('dis_mix_ratio = {0}\n'.format(dis_mix_ratio))
        f.write('guiding_centres = T\n')
        f.write('\n')

        f.write(f'mp_grid           = {kpoints[0]} {kpoints[1]} {kpoints[2]}\n')
        f.write('\n')

def run_step2(path, type, proj):
    # check checkpoint
    if read_checkpoint(path+'/step2',0) == True:
        os.makedirs(path+'/step2',exist_ok=True)
        write_checkpoint(path+'/step2',0)
        #copy files to step2
        shutil.copyfile(path+'/step1/CHGCAR',path+'/step2/CHGCAR')
        shutil.copyfile(path+'/step1/POSCAR',path+'/step2/POSCAR')
        shutil.copyfile(path+'/step1/POTCAR',path+'/step2/POTCAR')
        # make INCAR
        write_step2_INCAR(path, proj, wannier90 = False)
        # make KPOINTS
        if (type == 'hcp'):
            write_step2_KPOINTS(path, kmesh_hcp_step2)
        else:
            write_step2_KPOINTS(path, kmesh_bcc_fcc_step2)
        # make wannier90.win
        write_step2_winfile(path, type, proj)
    
        # --------> run vasp
        os.system('cd ' + path + '/step2 && srun '+ vasp_link +'/bin/vasp_ncl > stdout && cd .. && cd ..')
        # ------------------

        write_checkpoint(path+'/step2',1)

    if read_checkpoint(path+'/step2',1):
        # make INPUT
        write_step2_INCAR(path, proj, wannier90 = True)

        # --------> run vasp
        os.system('cd ' + path + '/step2 && srun --ntasks 1 '+ vasp_link +'/bin/vasp_ncl > stdout && cd .. && cd ..')
        # ------------------

        write_checkpoint(path+'/step2',2)

# functions for step3
def get_dis_win_en(path):
    data = list()
    with open(path+'/step2/wannier90.eig', 'r') as f:
        for s in f:
            data.append([float(x) for x in s.split()])
    data = np.matrix(data)
    id = np.where(data[:,2] > -15)
    return(data[min(id[0]),2])

def write_step3_wannier90(path): # write wannier90.win for step3
    dis_en = get_dis_win_en(path)
    # read INCAR from step1
    with open(path+'/step2/wannier90.win', 'r') as f1:
        # modify INCAR and save to step2
        with open(path+'/step3/wannier90.win', 'w') as f2:
            for s in f1:
                if re.search('#dis_win_min',s):
                    f2.write(f'dis_win_min       =   {dis_en}\n')
                    if (ws_distance == False):
                        f2.write('use_ws_distance = false\n')
                    f2.write('write_hr        = T\n')
                elif re.search('guiding_centres',s):
                    f2.write('guiding_centres = T\n')
                    if os.path.isfile(path+'/step3/wannier90.chk'):
                        f2.write('restart = wannierise\n')
                else:
                    f2.write(s)

def run_step3(path):
    # check checkpoint
    if read_checkpoint(path+'/step3',0) == True:
        os.makedirs(path+'/step3',exist_ok=True)
        write_checkpoint(path+'/step3',0)
        #copy files to step3
        if os.path.isfile(path+'/step3/wannier90.chk') == False:
            shutil.copyfile(path+'/step2/wannier90.amn',path+'/step3/wannier90.amn')
            shutil.copyfile(path+'/step2/wannier90.eig',path+'/step3/wannier90.eig')
            shutil.copyfile(path+'/step2/wannier90.mmn',path+'/step3/wannier90.mmn')
        
        #write wannier90.win
        write_step3_wannier90(path)

        # --------> run vasp
        os.system('cd ' + path + '/step3 && srun '+ wannier90_link +'/wannier90.x wannier90 && cd .. && cd ..')
        # ------------------

        write_checkpoint(path+'/step3',1)

# functions for step4
def write_step4_wannier90(path, type, ahc): # write wannier90.win for step4
    dis_en = get_dis_win_en(path)
    # read INCAR from step1
    with open(path+'/step3/wannier90.win', 'r') as f1:
        # modify INCAR and save to step2
        with open(path+'/step4/wannier90.win', 'w') as f2:
            for s in f1:
                if re.search('#berry ',s):
                    f2.write('berry                        = true\n')
                elif re.search('#berry_task ',s):
                    if (ahc == 0):
                        f2.write('berry_task                   = shc\n')
                    else:
                        f2.write('berry_task                   = ahc+shc\n')
                elif re.search('#berry_kmesh ',s):
                    if (type == 'hcp'):
                        f2.write('berry_kmesh                  =  {0}, {1}, {2}\n'.format(\
                            kmesh_hcp_step4[0], kmesh_hcp_step4[1], kmesh_hcp_step4[2]))
                    else:
                        f2.write('berry_kmesh                  =  {0}, {1}, {2}\n'.format(\
                            kmesh_bcc_fcc_step4[0],kmesh_bcc_fcc_step4[1],kmesh_bcc_fcc_step4[2]))
                elif re.search('#berry_curv_unit ',s):
                    f2.write('berry_curv_unit              = ang2\n')
                elif re.search('#berry_curv_adpt_kmesh ',s):
                    f2.write('berry_curv_adpt_kmesh        = 5\n')
                elif re.search('#berry_curv_adpt_kmesh_thresh ',s):
                    f2.write('berry_curv_adpt_kmesh_thresh = 100.0\n')
                else:
                    f2.write(s)

def run_step4(path, type, ahc):
    # check checkpoint
    if read_checkpoint(path+'/step4',0) == True:
        os.makedirs(path+'/step4',exist_ok=True)
        write_checkpoint(path+'/step4',0)
        #copy files to step3
        shutil.copyfile(path+'/step3/wannier90.amn',path+'/step4/wannier90.amn')
        shutil.copyfile(path+'/step3/wannier90.eig',path+'/step4/wannier90.eig')
        shutil.copyfile(path+'/step3/wannier90.mmn',path+'/step4/wannier90.mmn')
        shutil.copyfile(path+'/step3/wannier90.chk',path+'/step4/wannier90.chk')
        shutil.copyfile(path+'/step2/wannier90.spn',path+'/step4/wannier90.spn')
        
        #write wannier90.win
        write_step4_wannier90(path, type, ahc)

        # --------> run vasp
        os.system('cd ' + path + '/step4 && srun '+ wannier90_link +'/postw90.x wannier90 && cd .. && cd ..')
        # ------------------

        write_checkpoint(path+'/step4',1)


# functions for step4_tb
def write_step4_tb_wannier90(path, type, proj, ahc): # write wannier90.win for step4
    # read POSCAR from step1
    lat, atom_label, atom_num, atom_pos = read_poscar(path)
    efermi = read_outcar(path)
    if (type == 'hcp'):
        kmesh = kmesh_hcp_step4_tb
    else:
        kmesh = kmesh_bcc_fcc_step4_tb
    with open(path+'/step3/wannier90.win', 'r') as f1:
        # modify INCAR and save to step2
        with open(path+'/step4_tb/wannier90.win', 'w') as f2:
            for s in f1:
                if re.search('begin unit_cell_cart',s):
                    f2.write('begin unit_cell_cart\n')
                    for i in range(4):
                        f2.write(f1.readline())
                    f2.write('end unit_cell_cart\n\n')

            f2.write('begin atoms_frac\n')
            for i in range(len(atom_label)):
                for j in range(atom_num[i]):
                    f2.write('{0} {1} {2} {3}\n'.format(atom_label[i], atom_pos[i][j][0], atom_pos[i][j][1], atom_pos[i][j][2]))
            f2.write('end atoms_frac\n\n')

            f2.write('begin projections\n')
            for i in range(len(atom_label)):
                if (proj == 1):
                    f2.write('{0} : s;p;d\n'.format(atom_label[i]))
                else:
                    f2.write('{0} : s;p\n'.format(atom_label[i]))
            f2.write('end projections\n\n')

            f2.write('&kpath\n')
            if (type == 'bcc'):
                f2.write('5 200\n')
                f2.write('G  0.00000  0.00000  0.00000\n')
                f2.write('H -0.50000  0.50000  0.50000\n')
                f2.write('N  0.00000  0.50000  0.00000\n')
                f2.write('P  0.25000  0.25000  0.25000\n')
                f2.write('G  0.00000  0.00000  0.00000\n')
            elif (type == 'fcc'):
                f2.write('5 200\n')
                f2.write('K  0.37500  0.75000  0.37500\n')
                f2.write('G  0.00000  0.00000  0.00000\n')
                f2.write('L  0.50000  0.50000  0.50000\n')
                f2.write('W  0.25000  0.75000  0.50000\n')
                f2.write('X  0.00000  0.50000  0.50000\n')
            else:
                f2.write('7 200\n')
                f2.write('G  0.00000  0.00000  0.000000\n')
                f2.write('M  0.00000  0.50000  0.00000\n')
                f2.write('K  0.33333  0.33333  0.00000\n')
                f2.write('G  0.00000  0.00000  0.00000\n')
                f2.write('A  0.00000  0.00000  0.50000\n')
                f2.write('L  0.00000  0.50000  0.50000\n')
                f2.write('H  0.33333  0.33333  0.50000\n')
            f2.write('/\n\n')

            f2.write('&calc_type\n')
            if (ahc == 1):
                f2.write('calc_ahc        =   T\n')
            else:
                f2.write('calc_ahc        =   F\n')
            f2.write('calc_shc        =   T\n')
            f2.write('calc_ohc        =   T\n')
            f2.write('plot_band       =   T\n')
            f2.write('plot_band_ahc   =   F\n')
            f2.write('plot_band_shc   =   F\n')
            f2.write('plot_band_ohc   =   F\n')
            f2.write('plot_fermi	  =	  F\n')
            f2.write('plot_fermi_ahc  =   F\n')
            f2.write('plot_fermi_shc  =   F\n')
            f2.write('plot_fermi_ohc  =   F\n')
            if (ahc == 1):
                f2.write('scan_ahc	      =   T\n')
            else:
                f2.write('scan_ahc	      =   F\n')
            f2.write('scan_shc	      =   F\n')
            f2.write('scan_ohc	      =   F\n')
            f2.write('/\n\n')

            f2.write('&plotpar\n')
            f2.write('band_emin  = 0.0\n')
            f2.write('band_emax  = 0.0\n')
            f2.write('fermi_nkx  = 10\n')
            f2.write('fermi_nky  = 10\n')
            f2.write('fermi_nkz  = 10\n')
            f2.write('berry_type = 4\n')
            f2.write('/\n\n')

            f2.write('&kmesh\n')
            f2.write('nkx  =  {0}\n'.format(kmesh[0]))
            f2.write('nky  =  {0}\n'.format(kmesh[1]))
            f2.write('nkz  =  {0}\n'.format(kmesh[2]))
            f2.write('/\n\n')

            f2.write('&adpt_kmesh\n')
            f2.write('use_adpt_kmesh  = F\n')
            f2.write('adpt_nkx  		=  5 \n')
            f2.write('adpt_nky  		=  5 \n')
            f2.write('adpt_nkz  		=  5\n')
            f2.write('adpt_div  		=  5\n')
            f2.write('/\n\n')

            f2.write('&enpar\n')
            f2.write('en_min        =  {0}\n'.format(efermi-1))
            f2.write('en_max        =  {0}\n'.format(efermi+1))
            f2.write('en_num        =  80\n')
            f2.write('efermi        =  {0}\n'.format(efermi))
            f2.write('scale_efermi  =  F\n')
            f2.write('proj_type     =  1\n')
            f2.write('/\n\n')

def run_step4_tb(path, type, proj, ahc):
    # check checkpoint
    if read_checkpoint(path+'/step4_tb',0) == True:
        os.makedirs(path+'/step4_tb',exist_ok=True)
        write_checkpoint(path+'/step4_tb',0)
        #copy files to step3
        shutil.copyfile(path+'/step3/wannier90_hr.dat',path+'/step4_tb/wannier90_hr.dat')
        
        #write wannier90.win
        write_step4_tb_wannier90(path, type, proj, ahc)

        # --------> run vasp
        os.system('cd ' + path + '/step4_tb && srun '+ wan4hall_link +'/w4h_mpi > stdout && cd .. && cd ..')
        # ------------------

        write_checkpoint(path+'/step4_tb',1)

# functions for band_wannier90
def write_band_wannier90(path): # write wannier90.win for step4
    with open(path+'/step3/wannier90.win', 'r') as f1:
        # modify INCAR and save to step2
        with open(path+'/band_w90/wannier90.win', 'w') as f2:
            for s in f1:
                if re.search('#kpath ',s):
                    f2.write('kpath                   = true\n')
                elif re.search('#kpath_task ',s):
                    f2.write('kpath_task              = bands\n')
                elif re.search('#kpath_num_points ',s):
                    f2.write('kpath_num_points        = 400\n')
                else:
                    f2.write(s)

def run_band_w90(path):
    # check checkpoint
    if read_checkpoint(path+'/band_w90',0) == True:
        os.makedirs(path+'/band_w90',exist_ok=True)
        write_checkpoint(path+'/band_w90',0)
        #copy files to step3
        shutil.copyfile(path+'/step3/wannier90.amn',path+'/band_w90/wannier90.amn')
        shutil.copyfile(path+'/step3/wannier90.eig',path+'/band_w90/wannier90.eig')
        shutil.copyfile(path+'/step3/wannier90.mmn',path+'/band_w90/wannier90.mmn')
        shutil.copyfile(path+'/step3/wannier90.chk',path+'/band_w90/wannier90.chk')
        
        #write wannier90.win
        write_band_wannier90(path)

        # --------> run vasp
        os.system('cd ' + path + '/band_w90 && srun '+ wannier90_link +'/postw90.x wannier90 && cd .. && cd ..')
        # ------------------

        write_checkpoint(path+'/band_w90',1)

        if read_checkpoint(path+'/band_w90',1) == True:
            os.system('cd ' + path + '/band_w90 && python wannier90-bands.py && cd .. && cd ..')

# functions for band_vasp
def write_band_KPOINTS(path, type): # write KPOINTS for band_vasp
    with open(path+'/band_vasp/KPOINTS', 'w') as f:
        f.write('k-points along high symmetry line\n')
        f.write('  150\n')
        f.write('line-mode\n')
        f.write('rec\n')

        if (type == 'bcc'):
            f.write('  0.00000  0.00000  0.00000  !G\n')
            f.write(' -0.50000  0.50000  0.50000  !H\n')
            f.write('\n')

            f.write(' -0.50000  0.50000  0.50000  !H\n')
            f.write('  0.00000  0.50000  0.00000  !N\n')
            f.write('\n')

            f.write('  0.00000  0.50000  0.00000  !N\n')
            f.write('  0.25000  0.25000  0.25000  !P\n')
            f.write('\n')

            f.write('  0.25000  0.25000  0.25000  !P\n')
            f.write('  0.00000  0.00000  0.00000  !G\n')
            f.write('\n')

            f.write('  0.00000  0.00000  0.00000  !G\n')
            f.write('  0.00000  0.50000  0.00000  !N\n')
            f.write('\n')

        elif (type == 'fcc'):
            f.write('  0.37500  0.75000  0.37500  !K\n')
            f.write('  0.00000  0.00000  0.00000  !G\n')
            f.write('\n')

            f.write('  0.00000  0.00000  0.00000  !G\n')
            f.write('  0.50000  0.50000  0.50000  !L\n')
            f.write('\n')

            f.write('  0.50000  0.50000  0.50000  !L\n')
            f.write('  0.25000  0.75000  0.50000  !W\n')
            f.write('\n')

            f.write('  0.25000  0.75000  0.50000  !W\n')
            f.write('  0.00000  0.50000  0.50000  !X\n')
            f.write('\n')

            f.write('  0.00000  0.50000  0.50000  !X\n')
            f.write('  0.00000  0.00000  0.00000  !G\n')
            f.write('\n')
        else:
            f.write('  0.00000  0.00000  0.00000  !G\n')
            f.write('  0.00000  0.50000  0.00000  !M\n')
            f.write('\n')

            f.write('  0.00000  0.50000  0.00000  !M\n')
            f.write('  0.33333  0.33333  0.00000  !K\n')
            f.write('\n')

            f.write('  0.33333  0.33333  0.00000  !K\n')
            f.write('  0.00000  0.00000  0.00000  !G\n')
            f.write('\n')

            f.write('  0.00000  0.00000  0.00000  !G\n')
            f.write('  0.00000  0.00000  0.50000  !A\n')
            f.write('\n')

            f.write('  0.00000  0.00000  0.50000  !A\n')
            f.write('  0.00000  0.50000  0.50000  !L\n')
            f.write('\n')

            f.write('  0.00000  0.50000  0.50000  !L\n')
            f.write('  0.33333  0.33333  0.50000  !H\n')
            f.write('\n')

            f.write('  0.33333  0.33333  0.50000  !H\n')
            f.write('  0.00000  0.00000  0.50000  !A\n')
            f.write('\n')

def write_band_INCAR(path): # write INCAR for band_vasp
    # read INCAR from step2
    with open(path+'/step2/INCAR', 'r') as f1:
        # modify INCAR and save to step2
        with open(path+'/band_vasp/INCAR', 'w') as f2:
            for s in f1:
                if re.search('ISTART',s):
                    f2.write('  ISTART = 0         	    # Job: 0-new  1-cont  2-samecut\n')
                elif re.search('NELM ',s):
                    f2.write('  NELM   =     200     # Number of ELM steps\n')
                elif re.search('ALGO',s):
                    f2.write('  ALGO   =     Fast    # Electronic algorithm minimization\n')
                elif re.search('LWAVE',s):
                    f2.write('  LWAVE   = .FALSE.\n')
                elif re.search('NPAR',s):
                    f2.write('  NPAR   = 4\n')
                elif re.search('LWANNIER90',s):
                    f2.write('  LWANNIER90   =  False\n')
                else:
                    f2.write(s)

def run_band_vasp(path, type):
    # check checkpoint
    if read_checkpoint(path+'/band_vasp',0) == True:
        os.makedirs(path+'/band_vasp',exist_ok=True)
        write_checkpoint(path+'/band_vasp',0)
        #copy files to step2
        shutil.copyfile(path+'/step2/CHGCAR',path+'/band_vasp/CHGCAR')
        shutil.copyfile(path+'/step2/POSCAR',path+'/band_vasp/POSCAR')
        shutil.copyfile(path+'/step2/POTCAR',path+'/band_vasp/POTCAR')
        # make INCAR
        write_band_INCAR(path)
        # make KPOINTS
        write_band_KPOINTS(path, type)
    
        # --------> run vasp
        os.system('cd ' + path + '/band_vasp && srun '+ vasp_link +'/bin/vasp_ncl > stdout && cd .. && cd ..')
        # ------------------

        write_checkpoint(path+'/band_vasp',1)
#-----------------------------

def save_log(info):
    with open('output.log', 'a') as f:
        f.write(info)

#main program
if __name__ == '__main__':
    # get info of elements to run calculation
    list_struct = read_elements()
    # perform calculations
    for i in range(len(list_struct)):
        save_log(f'Calculating {list_struct[i][0]}\n')
        #step1
        run_step1(list_struct[i][0],list_struct[i][1],list_struct[i][2])
        #step2
        run_step2(list_struct[i][0],list_struct[i][1], list_struct[i][3])
        #step3
        run_step3(list_struct[i][0])
        #band_w90
        run_band_w90(list_struct[i][0])
        #band_vasp
        run_band_vasp(list_struct[i][0], list_struct[i][1])
        #step4
        run_step4(list_struct[i][0], list_struct[i][1], list_struct[i][4])
        #step4_tb
        run_step4_tb(list_struct[i][0], list_struct[i][1], list_struct[i][3], list_struct[i][4])
        # print log
        save_log('Done.\n')
        save_log('--------------------------------------\n')
    