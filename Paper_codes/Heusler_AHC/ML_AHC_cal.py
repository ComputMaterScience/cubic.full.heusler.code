import os
import shutil
import numpy as np
import re


vasp_link = '/home/user2/thiho/app/vasp.6.3.1_wan3/bin'
vasp_link2 = '/home/user2/thiho/app/vasp.5.4.4_wan3/bin'
wannier90_link = '/home/user2/thiho/app/wannier90-3.1.0'

########################### IO Functions ##############################

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

# functions for step2
def write_step2_INCAR(path,wannier90 = False): # write INCAR for step2
    tag_ISYM = False
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
                        f2.write('  NPAR   = 1\n')
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
            f2.write('  NBANDS        = 98\n')
        f2.write('  ISYM          = -1\n')

        if (wannier90 == True):
            f2.write('\n')
            f2.write('Wannier90:\n')
            f2.write('  LWANNIER90   =  True\n')
            # make wannier90.win
            write_step2_winfile(path)

def write_step2_spn_INCAR(path,wannier90 = False): # write INCAR for step2
    tag_ISYM = False
    # read INCAR from step1
    with open(path+'/step1/INCAR', 'r') as f1:
        # modify INCAR and save to step2
        with open(path+'/step2_spn/INCAR', 'w') as f2:
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
                        f2.write('  NPAR   = 1\n')
                else:
                    f2.write(s)

    #append tags to INCAR
    with open(path+'/step2_spn/INCAR', 'a') as f2:
        f2.write('\n')
        f2.write('SOC calculation:\n')
        f2.write('  LMAXMIX       = 4\n')
        f2.write('  SAXIS         = 0 0 1\n')
        f2.write('  LNONCOLLINEAR = .TRUE.\n')
        f2.write('  LSORBIT       = .TRUE.\n')
        f2.write('  GGA_COMPAT    = F\n')
        f2.write('  LORBMOM       = T\n')
        if (tag_ISYM == False):
            f2.write('  NBANDS        = 98\n')
        f2.write('  ISYM          = -1\n')

        if (wannier90 == True):
            f2.write('\n')
            f2.write('Wannier90:\n')
            f2.write('  LWANNIER90   =  True\n')
            f2.write('  LCALC_MMN    =  .FALSE.\n')
            f2.write('  LCALC_AMN    =  .FALSE.\n')
            f2.write('  LWRITE_MMN   =  .FALSE.\n')
            f2.write('  LWRITE_AMN   =  .FALSE.\n')
            f2.write('  LWRITE_EIG   =  .FALSE.\n')
            f2.write('  LWRITE_SPN   =  .TRUE.\n')
            # make wannier90.win
            write_step2_winfile(path)

def write_step2_KPOINTS(path,kmesh): # write KPOINTS for step2
    with open(path+'/step2/KPOINTS', 'w') as f:
        f.write('Automatic Mesh # Generates Automatically the K-mesh\n')
        f.write('  0\n')
        f.write('Gamma # Automatic Distribution of the mesh\n')
        f.write(f'  {kmesh[0]}  {kmesh[1]}  {kmesh[2]}\n')

def read_poscar(path):
    with open(path+'/step1/POSCAR', 'r') as f:
        f.readline() # ignore comment line
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
                tmp_pos.append([float(x) for x in f.readline().split()])
            atom_pos.append(tmp_pos)
        return lat, atom_label, atom_num, atom_pos

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

def write_step2_winfile(path): # write wannier90.in for step2
    # read data
    lat, atom_label, atom_num, atom_pos = read_poscar(path)
    efermi = read_outcar(path)
    kpoints = read_kpoints(path)
    with open(path+'/step2/INCAR', 'a') as f:
        f.write('NUM_WANN = 62\n')
        f.write('WANNIER90_WIN = "\n')
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
        f.write(f'{atom_label[0]} : sp3d2;dxy;dxz;dyz\n')
        if atom_label[0] != atom_label[1]:
            f.write(f'{atom_label[1]} : sp3d2;dxy;dxz;dyz\n')
        f.write(f'{atom_label[2]} : sp3\n')
        f.write('end projections\n')
        f.write('\n')

        # num_wann and num_bands tags
        f.write('num_bands  =  98\n')
        f.write('num_wann   =  62\n')
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
        f.write('G  0.00000  0.00000  0.00000    X  0.00000  0.50000  0.50000\n')
        f.write('X  0.00000  0.50000  0.50000    W  0.25000  0.75000  0.50000\n')
        f.write('W  0.25000  0.75000  0.50000    K  0.37500  0.75000  0.37500\n')
        f.write('K  0.37500  0.75000  0.37500    G  0.00000  0.00000  0.00000\n')
        f.write('G  0.00000  0.00000  0.00000    L  0.50000  0.50000  0.50000\n')
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
        f.write(f'dis_froz_max      =  {efermi+5}\n')
        f.write('dis_num_iter       =  800000\n')
        f.write('dis_conv_tol       = 1.0e-10\n')
        f.write('conv_tol           = 1.0e-10\n')
        f.write('conv_window        = 10\n')
        f.write('num_iter           = 800000\n')
        f.write('guiding_centres = T\n')
        f.write('\n')

        f.write(f'mp_grid           = {kpoints[0]} {kpoints[1]} {kpoints[2]}\n')
        f.write('\n')
        f.write('"\n')

def run_step2(path):
    # check checkpoint
    if read_checkpoint(path+'/step2',0) == True:
        os.makedirs(path+'/step2',exist_ok=True)
        write_checkpoint(path+'/step2',0)
        #copy files to step2
        shutil.copyfile(path+'/step1/CHGCAR',path+'/step2/CHGCAR')
        shutil.copyfile(path+'/step1/POSCAR',path+'/step2/POSCAR')
        shutil.copyfile(path+'/step1/POTCAR',path+'/step2/POTCAR')
        # make INCAR
        write_step2_INCAR(path,wannier90 = False)
        # make KPOINTS
        write_step2_KPOINTS(path,[4,4,4])
    
        # --------> run vasp
        os.system('cd ' + path + '/step2 && srun --ntasks 12 '+ vasp_link +'/vasp_ncl > stdout && cd .. && cd ..')
        # ------------------

        write_checkpoint(path+'/step2',1)

    if read_checkpoint(path+'/step2',1):
        # make INPUT
        write_step2_INCAR(path,wannier90 = True)

        # --------> run vasp
        os.system('cd ' + path + '/step2 && srun --ntasks 1 '+ vasp_link +'/vasp_ncl > stdout && cd .. && cd ..')
        # ------------------

        write_checkpoint(path+'/step2',2)
        
def run_step2_spn(path):
    # check checkpoint
    if read_checkpoint(path+'/step2_spn',0) == True:
        os.makedirs(path+'/step2_spn',exist_ok=True)
        write_checkpoint(path+'/step2_spn',0)
        #copy files to step2_spn
        src_files = os.listdir(path+'/step2')
        for file_name in src_files:
            full_file_name = os.path.join(path+'/step2', file_name)
            if os.path.isfile(full_file_name):
                shutil.copy(full_file_name, path+'/step2_spn')

        # make INCAR
        write_step2_spn_INCAR(path,wannier90 = True)
    
        # --------> run vasp
        os.system('cd ' + path + '/step2_spn && srun --ntasks 1 '+ vasp_link2 +'/vasp_ncl > stdout && cd .. && cd ..')
        # ------------------

        write_checkpoint(path+'/step2_spn',1)

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
        os.system('cd ' + path + '/step3 && srun --ntasks 52  '+ wannier90_link +'/wannier90.x wannier90 && cd .. && cd ..')
        # ------------------

        write_checkpoint(path+'/step3',1)

# functions for step4
def write_step4_wannier90(path): # write wannier90.win for step4
    dis_en = get_dis_win_en(path)
    # read INCAR from step1
    with open(path+'/step3/wannier90.win', 'r') as f1:
        # modify INCAR and save to step2
        with open(path+'/step4/wannier90.win', 'w') as f2:
            for s in f1:
                if re.search('#berry ',s):
                    f2.write('berry                        = true\n')
                elif re.search('#berry_task ',s):
                    f2.write('berry_task                   = eval_ahc\n')
                elif re.search('#berry_kmesh ',s):
                    f2.write('berry_kmesh                  =  100\n')
                elif re.search('#berry_curv_unit ',s):
                    f2.write('berry_curv_unit              = ang2\n')
                elif re.search('#berry_curv_adpt_kmesh ',s):
                    f2.write('berry_curv_adpt_kmesh        = 5\n')
                elif re.search('#berry_curv_adpt_kmesh_thresh ',s):
                    f2.write('berry_curv_adpt_kmesh_thresh = 100.0\n')
                else:
                    f2.write(s)

def run_step4(path):
    # check checkpoint
    if read_checkpoint(path+'/step4',0) == True:
        os.makedirs(path+'/step4',exist_ok=True)
        write_checkpoint(path+'/step4',0)
        #copy files to step3
        shutil.copyfile(path+'/step3/wannier90.amn',path+'/step4/wannier90.amn')
        shutil.copyfile(path+'/step3/wannier90.eig',path+'/step4/wannier90.eig')
        shutil.copyfile(path+'/step3/wannier90.mmn',path+'/step4/wannier90.mmn')
        shutil.copyfile(path+'/step3/wannier90.chk',path+'/step4/wannier90.chk')
        
        #write wannier90.win
        write_step4_wannier90(path)

        # --------> run vasp
        os.system('cd ' + path + '/step4 && srun --ntasks 52  '+ wannier90_link +'/postw90.x wannier90 && cd .. && cd ..')
        # ------------------

        write_checkpoint(path+'/step4',1)
        
# functions for step5
def write_step5_wannier90(path): # write wannier90.win for step5
    dis_en = get_dis_win_en(path)
    efermi = read_outcar(path)
    # read INCAR from step1
    with open(path+'/step3/wannier90.win', 'r') as f1:
        # modify INCAR and save to step2
        with open(path+'/step5/wannier90.win', 'w') as f2:
            for s in f1:
                if re.search('#berry ',s):
                    f2.write('berry                        = true\n')
                elif re.search('#berry_task ',s):
                    f2.write('berry_task                   = ahc\n')
                elif re.search('#berry_kmesh ',s):
                    f2.write('berry_kmesh                  =  100\n')
                elif re.search('#berry_curv_unit ',s):
                    f2.write('berry_curv_unit              = ang2\n')
                elif re.search('#berry_curv_adpt_kmesh ',s):
                    f2.write('berry_curv_adpt_kmesh        = 5\n')
                elif re.search('#berry_curv_adpt_kmesh_thresh ',s):
                    f2.write('berry_curv_adpt_kmesh_thresh = 100.0\n')
                elif re.search('#fermi_energy_min ',s):
                    f2.write(f'fermi_energy_min = {efermi-10}\n')
                elif re.search('#fermi_energy_max ',s):
                    f2.write(f'fermi_energy_max = {efermi+4}\n')
                elif re.search('#fermi_energy_step ',s):
                    f2.write('fermi_energy_step = 0.1\n')
                elif re.search('fermi_energy =',s):
                    pass # keep writing EF
                else:
                    f2.write(s)

def run_step5(path):
    # check checkpoint
    if read_checkpoint(path+'/step5',0) == True:
        os.makedirs(path+'/step5',exist_ok=True)
        write_checkpoint(path+'/step5',0)
        #copy files to step3
        shutil.copyfile(path+'/step3/wannier90.amn',path+'/step5/wannier90.amn')
        shutil.copyfile(path+'/step3/wannier90.eig',path+'/step5/wannier90.eig')
        shutil.copyfile(path+'/step3/wannier90.mmn',path+'/step5/wannier90.mmn')
        shutil.copyfile(path+'/step3/wannier90.chk',path+'/step5/wannier90.chk')
        
        #write wannier90.win
        write_step5_wannier90(path)

        # --------> run vasp
        os.system('cd ' + path + '/step5 && srun --ntasks 52  '+ wannier90_link +'/postw90.x wannier90 && cd .. && cd ..')
        # ------------------

        write_checkpoint(path+'/step5',1)
        
# functions for step6
def write_step6_wannier90(path): # write wannier90.win for step6
    dis_en = get_dis_win_en(path)
    efermi = read_outcar(path)
    # read INCAR from step1
    with open(path+'/step3/wannier90.win', 'r') as f1:
        # modify INCAR and save to step2
        with open(path+'/step6/wannier90.win', 'w') as f2:
            for s in f1:
                if re.search('#kpath_num_points ',s):
                    f2.write('\n')
                    f2.write('kpath = true\n')
                    f2.write('kpath_task = bands+curv\n')
                    f2.write('kpath_bands_colour = spin\n')
                    f2.write('kpath_num_points = 400\n')
                    f2.write('kubo_adpt_smr = false\n')
                    f2.write('kubo_smr_fixed_en_width = 0.05\n')
                else:
                    f2.write(s)

def run_step6(path):
    # check checkpoint
    if read_checkpoint(path+'/step6',0) == True:
        os.makedirs(path+'/step6',exist_ok=True)
        write_checkpoint(path+'/step6',0)
        #copy files to step3
        shutil.copyfile(path+'/step3/wannier90.amn',path+'/step6/wannier90.amn')
        shutil.copyfile(path+'/step3/wannier90.eig',path+'/step6/wannier90.eig')
        shutil.copyfile(path+'/step3/wannier90.mmn',path+'/step6/wannier90.mmn')
        shutil.copyfile(path+'/step3/wannier90.chk',path+'/step6/wannier90.chk')
        shutil.copyfile(path+'/step2_spn/wannier90.spn',path+'/step6/wannier90.spn')
        
        #write wannier90.win
        write_step6_wannier90(path)

        # --------> run vasp
        os.system('cd ' + path + '/step6 && srun --ntasks 52  '+ wannier90_link +'/postw90.x wannier90 && cd .. && cd ..')
        # ------------------

        write_checkpoint(path+'/step6',1)

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
        os.system('cd ' + path + '/band_w90 && srun --ntasks 52 '+ wannier90_link +'/postw90.x wannier90 && cd .. && cd ..')
        # ------------------

        write_checkpoint(path+'/band_w90',1)

        if read_checkpoint(path+'/band_w90',1) == True:
            os.system('cd ' + path + '/band_w90 && python wannier90-bands.py && cd .. && cd ..')

# functions for band_vasp
def write_band_KPOINTS(path): # write KPOINTS for band_vasp
    with open(path+'/band_vasp/KPOINTS', 'w') as f:
        f.write('k-points along high symmetry line\n')
        f.write('  150\n')
        f.write('line-mode\n')
        f.write('rec\n')

        f.write('  0.00000  0.00000  0.00000  !G\n')
        f.write('  0.00000  0.50000  0.50000  !X\n')
        f.write('\n')

        f.write('  0.00000  0.50000  0.50000  !X\n')
        f.write('  0.25000  0.75000  0.50000  !W\n')
        f.write('\n')

        f.write('  0.25000  0.75000  0.50000  !W\n')
        f.write('  0.37500  0.75000  0.37500  !K\n')
        f.write('\n')

        f.write('  0.37500  0.75000  0.37500  !K\n')
        f.write('  0.00000  0.00000  0.00000  !G\n')
        f.write('\n')

        f.write('  0.00000  0.00000  0.00000  !G\n')
        f.write('  0.50000  0.50000  0.50000  !L\n')
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

def run_band_vasp(path):
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
        write_band_KPOINTS(path)
    
        # --------> run vasp
        os.system('cd ' + path + '/band_vasp && srun --ntasks 52 '+ vasp_link +'/vasp_ncl > stdout && cd .. && cd ..')
        # ------------------

        write_checkpoint(path+'/band_vasp',1)
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
        #step2
        run_step2(list_struct[i])
        #step2_spn
        run_step2_spn(list_struct[i])
        #step3
        run_step3(list_struct[i])
        #band_w90
        run_band_w90(list_struct[i])
        #band_vasp
        run_band_vasp(list_struct[i])
        #step4
        run_step4(list_struct[i])
        #step5
        run_step5(list_struct[i])
        #step6
        run_step6(list_struct[i])
        #save checkpoint
        if i + 1 < len(list_struct):
            print(list_struct[i+1:])
            save_list('structures.in',list_struct[i+1:])
        else:
            save_list('structures.in',[])
        save_log('Done.\n')
        save_log('--------------------------------------\n')
    