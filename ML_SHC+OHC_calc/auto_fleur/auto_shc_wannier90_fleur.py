import os
import shutil
import numpy as np
import re
import xml.etree.ElementTree as et


vasp_link = '/home/user2/thiho/app/vasp.5.4.4_wan3'
fleur_link = '/home/user2/thiho/app/fleur_d'
wannier90_link = '/home/user2/thiho/app/wannier90_OHC/wannier90-3.1.0_mod'

element_name = ['H', 'He', 'Li', 'Be', 'B',	'C', 'N', 'O', 'F',	'Ne', 'Na',	'Mg', 'Al',	'Si', 'P',\
    	'S', 'Cl', 'Ar', 'K', 'Ca',	'Sc', 'Ti',	'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga',\
        'Ge', 'As',	'Se', 'Br',	'Kr', 'Rb',	'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag',\
        'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu',\
        'Gd', 'Tb',	'Dy', 'Ho',	'Er', 'Tm',	'Yb', 'Lu',	'Hf', 'Ta',	'W', 'Re', 'Os', 'Ir', 'Pt', 'Au',\
        'Hg', 'Tl', 'Pb', 'Bi',	'Po', 'At',	'Rn', 'Fr',	'Ra', 'Ac',	'Th', 'Pa',	'U', 'Np', 'Pu', 'Am',\
        'Cm', 'Bk',	'Cf', 'Es',	'Fm', 'Md',	'No', 'Lr',	'Rf', 'Db', 'Sg', 'Bh',	'Hs', 'Mt',	'Ds', 'Rg',\
        'Cn ', 'Uut', 'Uuq', 'Uup', 'Uuh', 'Uus', 'Uuo']

element_number = [1,2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,	17,	18,	19,	20,	21,	22,\
    	23,	24,	25,	26,	27,	28,	29,	30,	31,	32,	33,	34,	35,	36,	37,	38,	39,	40,	41,	42,	43,	44,\
        45,	46,	47,	48,	49,	50,	51,	52,	53,	54,	55,	56,	57,	58,	59,	60,	61,	62,	63,	64,	65,	66,\
        67,	68,	69,	70,	71,	72,	73,	74,	75,	76,	77,	78,	79,	80,	81,	82,	83,	84,	85,	86,	87,	88,\
        89,	90,	91,	92,	93,	94,	95,	96,	97,	98,	99,	100, 101, 102, 103,	104, 105, 106, 107, 108,\
        109, 110, 111, 112, 113, 114, 115, 116, 117, 118]

# Input Variables

pot_type = 1 # 1: PBE, 2: LDA

# step1
encut   = 500

# step2
kmesh_bcc_fcc_step2 = [8, 8, 8]
kmesh_hcp_step2     = [10, 10, 8]

# step3
num_iter = 0
num_dis = 1000
dis_mix_ratio = 1.0
up_fermi = 5.0

# step4
kmesh_bcc_fcc_step4 = [100, 100, 100]
kmesh_hcp_step4     = [100, 100, 80]

# step4_tb
tb_task = 'shc+ohc'

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
            elements_list.append([s[0].capitalize(), s[1].lower(), float(s[2])])
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

def write_step2_inp(path, type): # make input file for fleur
    # read POSCAR from step1
    lat, atom_label, atom_num, atom_pos = read_poscar(path)
    with open(path+'/step2/inp', 'w') as f:
        f.write('{0}\n\n'.format(path))
        f.write('&input film=f /\n\n')
        if (type == 'hcp'):
            f.write('&lattice latsys={0} a0=1.8897259886 a={1} c={2}/\n\n'.format(type,np.linalg.norm(lat[0]),\
                np.linalg.norm(lat[2])))
        elif (type == 'bcc'):
            f.write('&lattice latsys={0} a0=1.8897259886 a={1}/\n\n'.format(type,np.linalg.norm(lat[0])/\
                np.linalg.norm([0.5,0.5,-0.5])))
        else:
            f.write('&lattice latsys={0} a0=1.8897259886 a={1}/\n\n'.format(type,np.linalg.norm(lat[0])/\
                np.linalg.norm([0.0,0.5,0.5])))
        f.write('{0}\n'.format(np.sum(atom_num)))
        for i in range(len(atom_label)):
            id = element_name.index(atom_label[i].capitalize())
            for j in range(atom_num[i]):
                f.write('{0} {1} {2} {3}\n'.format(element_number[id],\
                    atom_pos[i][j][0], atom_pos[i][j][1], atom_pos[i][j][2]))

        f.write('\n')

        f.write('&soc 0.0 0.0 /\n\n')
        f.write('&comp kmax=4.5 gmaxxc=12.5 gmax=15.0/\n\n')
        if (pot_type == 1):
            f.write('&exco xctyp=\'pbe\' /\n\n')
        else:
            f.write('&exco xctyp=\'vwn\' /\n\n')

        f.write('&comp kmax=4.5 gmaxxc=12.5 gmax=15.0/\n\n')

        if (type == 'hcp'):
            f.write('&kpt div1=25 div2=25 div3=15/\n\n')
        else:
            f.write('&kpt div1=25 div2=25 div3=25/\n\n')
        
def modify_step2_xml(path): # change parameters in inp.xml
    # read inp.xml
    inp = et.parse(path+'/step2/inp.xml')
    root = inp.getroot()
    for i in root.iter('scfLoop'):
        i.set('itmax','99')
        i.set('alpha','.01000000')
    inp.write(path+'/step2/inp.xml')


def run_step2(path, type):
    # check checkpoint
    if read_checkpoint(path+'/step2',0) == True:
        os.makedirs(path+'/step2',exist_ok=True)
        write_checkpoint(path+'/step2',0)
        # make inp file
        write_step2_inp(path, type)

        if (not os.path.exists(path+'/step2/inp.xml')):
            # --------> run inpgen
            os.system('cd ' + path + '/step2 && '+ fleur_link +'/build/inpgen -f inp && cd .. && cd ..')
            # ------------------

        # modify inp.xml
        modify_step2_xml(path)

        # --------> run fleur
        os.system('cd ' + path + '/step2 && srun --ntasks 1 '+ fleur_link +'/build/fleur_MPI -last_extra > stdout \
                                        && cd .. && cd ..')
        # ------------------

        write_checkpoint(path+'/step2',1)

# functions for step3

def write_sym_xml(path): # write sym.xml for step3
    with open(path+'/step3/kpoints/sym.xml', 'w') as f:
        f.write('      <symmetryOperations>\n')
        f.write('         <symOp>\n')
        f.write('            <row-1>1 0 0 .0000000000</row-1>\n')
        f.write('            <row-2>0 1 0 .0000000000</row-2>\n')
        f.write('            <row-3>0 0 1 .0000000000</row-3>\n')
        f.write('         </symOp>\n')
        f.write('      </symmetryOperations>\n')

def write_kpts_xml(path, type): # write kpts.xml for step3
    if (os.path.exists(path+'/step3/kpoints') and os.path.isdir(path+'/step3/kpoints')):
        shutil.rmtree(path+'/step3/kpoints')
    os.makedirs(path+'/step3/kpoints')
    shutil.copyfile(path+'/step2/inp.xml',path+'/step3/kpoints/inp.xml')
    shutil.copyfile(path+'/step2/kpts.xml',path+'/step3/kpoints/kpts.xml')
    write_sym_xml(path)

    if (type == 'hcp'):
        kmesh = kmesh_hcp_step2
    else:
        kmesh = kmesh_bcc_fcc_step2

    # --------> run inpgen
    os.system('cd ' + path + '/step3/kpoints && '+ fleur_link +'/build/inpgen -inp.xml -kpt wan#gamma@grid={0},{1},{2} \
        && cd .. && cd ..'.format(kmesh[0], kmesh[1], kmesh[2]))
    # ------------------

def write_projgen_inp(path): # write projgen_inp file
    # read POSCAR from step1
    lat, atom_label, atom_num, atom_pos = read_poscar(path)
    with open(path+'/step3/projgen_inp', 'w') as f:
        for i in range(len(atom_label)):
            for j in range(3):
                f.write('{0} {1} 0 0\n'.format(atom_label[i], j))

def modify_step3_xml(path, step): # change parameters in inp.xml
    # read POSCAR from step1
    lat, atom_label, atom_num, atom_pos = read_poscar(path)
    # read inp.xml
    inp = et.parse(path+'/step3/inp.xml')
    root = inp.getroot()
    # determine number of bands
    num_lo = 0
    for i in root.iter('lo'): 
        if int(i.attrib['l']) == 0:
            num_lo += 2*np.sum(atom_num)
        elif int(i.attrib['l']) == 1:
            num_lo += 6*np.sum(atom_num)
        elif int(i.attrib['l']) == 3:
            num_lo += 14*np.sum(atom_num)

    # change parameters
    if step == 1:
        for i in root.iter('scfLoop'): i.set('itmax','1')
        for i in root.iter('kPointListSelection'): i.set('listName','wan')
        for i in root.iter('output'): 
            i.set('wannier','T')
            i.append(et.Element('wannier'))
            i.set('eig66','F')
        for i in root.iter('wannier'):
            i.append(et.Element('bandSelection'))
            i.append(et.Element('jobList'))
        for i in root.iter('bandSelection'): 
            i.set('minSpinUp',str(int(num_lo+1)))
            i.set('maxSpinUp',str(int(num_lo+np.sum(atom_num)*36)))
        for i in root.iter('cutoffs'): i.set('numbands',str(int(num_lo+np.sum(atom_num)*36)))
        for i in root.iter('jobList'): i.text = 'projgen prepwan90 stopopt'
    elif step == 2:
        for i in root.iter('jobList'): i.text = 'matrixmmn matrixamn'
    elif step == 3:
        for i in root.iter('output'): i.set('eig66','T')
        for i in root.iter('jobList'): i.text = 'updown mmn0'
    else:
        for i in root.iter('jobList'): i.text = 'mmn0_to_spn_unf'
    
    inp.write(path+'/step3/inp.xml')

def run_step3(path, type):
    # check checkpoint
    if read_checkpoint(path+'/step3',0) == True:
        os.makedirs(path+'/step3',exist_ok=True)
        write_checkpoint(path+'/step3',0)
        # make kpts.xml
        write_kpts_xml(path, type)
        write_checkpoint(path+'/step3',1)

    if read_checkpoint(path+'/step3',1) == True:
        # copy files
        shutil.copyfile(path+'/step3/kpoints/kpts.xml',path+'/step3/kpts.xml')
        shutil.copyfile(path+'/step2/inp.xml',path+'/step3/inp.xml')
        shutil.copyfile(path+'/step2/cdn_last.hdf',path+'/step3/cdn.hdf')
        shutil.copyfile(path+'/step2/sym.xml',path+'/step3/sym.xml')

        # modify inp.xml
        modify_step3_xml(path,1)

        # make projgen_inp
        write_projgen_inp(path)

        # --------> run fleur
        os.system('cd ' + path + '/step3 && srun --ntasks 1 '+ fleur_link +'/build/fleur_MPI -eig hdf > stdout && \
              cd .. && cd ..')
        # ------------------

        write_checkpoint(path+'/step3',2)

    if read_checkpoint(path+'/step3',2) == True:

        # modify inp.xml
        modify_step3_xml(path,2)

        # --------> run fleur
        os.system('cd ' + path + '/step3 && srun --ntasks 1 '+ fleur_link +'/build/fleur_MPI -eig hdf > stdout && \
              cd .. && cd ..')
        # ------------------

        write_checkpoint(path+'/step3',3)

    if read_checkpoint(path+'/step3',3) == True:

        # modify inp.xml
        modify_step3_xml(path,3)

        # --------> run fleur
        os.system('cd ' + path + '/step3 && srun --ntasks 1 '+ fleur_link +'/build/fleur_MPI -eig hdf > stdout && \
              cd .. && cd ..')
        # ------------------

        write_checkpoint(path+'/step3',4)

    if read_checkpoint(path+'/step3',4) == True:

        # modify inp.xml
        modify_step3_xml(path,4)

        # --------> run fleur
        os.system('cd ' + path + '/step3 && srun --ntasks 1 '+ fleur_link +'/build/fleur_MPI -eig hdf > stdout && \
              cd .. && cd ..')
        # ------------------

        write_checkpoint(path+'/step3',5)

# functions for step4

def read_efermi(path): # read Fermi energy from step2
    with open(path+'/step2/out', 'r') as f:
        for line in f:
            if (line.find('-->  new fermi energy') != -1): # get fermi energy
                s = line.split(':')
                s = s[1].split()
                efermi = float(s[0])*27.211396641308 # convert htr to eV
    return efermi

def write_step4_wannier90(path): # write wannier90.win for step4
    # read INCAR from step1
    with open(path+'/step3/WF1.win', 'r') as f1:
        # modify INCAR and save to step2
        with open(path+'/step4/WF1.win', 'w') as f2:
            for s in f1:
                if re.search('num_iter=         100',s):
                    f2.write(' num_iter = {0}\n'.format(num_iter))
                    f2.write(' dis_mix_ratio = 1.0\n')
                elif re.search('dis_num_iter=10000',s): 
                    f2.write(' dis_num_iter={0}\n'.format(num_dis))
                elif re.search('!restart=wannierise',s):
                    if os.path.isfile(path+'/step4/WF1.chk'):
                        f2.write(' restart = wannierise\n')
                elif re.search('!dis_froz_max=',s):
                    f2.write(' dis_froz_max = {0}\n'.format(read_efermi(path)+up_fermi))
                else:
                    f2.write(s)

def run_step4(path):
    # check checkpoint
    if read_checkpoint(path+'/step4',0) == True:
        os.makedirs(path+'/step4',exist_ok=True)
        write_checkpoint(path+'/step4',0)
        #copy files to step4
        if os.path.isfile(path+'/step4/WF1.chk') == False:
            shutil.copyfile(path+'/step3/WF1.amn',path+'/step4/WF1.amn')
            shutil.copyfile(path+'/step3/WF1.eig',path+'/step4/WF1.eig')
            shutil.copyfile(path+'/step3/WF1.mmn',path+'/step4/WF1.mmn')
        
        #modify WF1.win
        write_step4_wannier90(path)

        # --------> run wannier90
        os.system('cd ' + path + '/step4 && srun '+ wannier90_link +'/wannier90.x WF1 && cd .. && cd ..')
        # ------------------

        write_checkpoint(path+'/step4',1)

# functions for step5
def write_step5_wannier90(path, type): # write wannier90.win for step5
    # read INCAR from step1
    with open(path+'/step4/WF1.win', 'r') as f1:
        # modify INCAR and save to step2
        with open(path+'/step5/WF1.win', 'w') as f2:
            for s in f1:
                if re.search('end atoms_frac',s):
                    f2.write(' end atoms_frac\n\n')
                    f2.write(' berry                        = true\n')
                    f2.write(' berry_task                   = eval_shc\n')
                    if (type == 'hcp'):
                        f2.write(' berry_kmesh                  =  {0}, {1}, {2}\n'.format(\
                            kmesh_hcp_step4[0], kmesh_hcp_step4[1], kmesh_hcp_step4[2]))
                    else:
                        f2.write(' berry_kmesh                  =  {0}, {1}, {2}\n'.format(\
                            kmesh_bcc_fcc_step4[0],kmesh_bcc_fcc_step4[1],kmesh_bcc_fcc_step4[2]))
                    f2.write(' berry_curv_unit              = ang2\n')
                    f2.write(' berry_curv_adpt_kmesh        = 5\n')
                    f2.write(' berry_curv_adpt_kmesh_thresh = 100.0\n\n')
                    f2.write(' fermi_energy = {0}\n'.format(read_efermi(path)))
                else:
                    f2.write(s)

def run_step5(path, type):
    # check checkpoint
    if read_checkpoint(path+'/step5',0) == True:
        os.makedirs(path+'/step5',exist_ok=True)
        write_checkpoint(path+'/step5',0)
        #copy files to step5
        shutil.copyfile(path+'/step4/WF1.amn',path+'/step5/WF1.amn')
        shutil.copyfile(path+'/step4/WF1.eig',path+'/step5/WF1.eig')
        shutil.copyfile(path+'/step4/WF1.mmn',path+'/step5/WF1.mmn')
        shutil.copyfile(path+'/step4/WF1.chk',path+'/step5/WF1.chk')
        shutil.copyfile(path+'/step3/WF1.spn',path+'/step5/WF1.spn')
        
        #modify WF1.win
        write_step5_wannier90(path, type)

        # --------> run wannier90
        os.system('cd ' + path + '/step5 && srun '+ wannier90_link +'/postw90.x WF1 && cd .. && cd ..')
        # ------------------

        write_checkpoint(path+'/step5',1)


# functions for step5_tb
def write_step5_tb_wannier90(path, type): # write wannier90.win for step5
    # read INCAR from step1
    with open(path+'/step4/WF1.win', 'r') as f1:
        # modify INCAR and save to step2
        with open(path+'/step5_tb/WF1.win', 'w') as f2:
            for s in f1:
                if re.search('end atoms_frac',s):
                    f2.write(' end atoms_frac\n\n')
                    f2.write(' tb_berry                        = true\n')
                    f2.write(' tb_berry_task                   = {0}\n'.format(tb_task))
                    f2.write(' tb_berry_type                   = 1\n')
                    if (type == 'hcp'):
                        f2.write(' berry_kmesh                  =  {0}, {1}, {2}\n'.format(\
                            kmesh_hcp_step4[0], kmesh_hcp_step4[1], kmesh_hcp_step4[2]))
                    else:
                        f2.write(' berry_kmesh                  =  {0}, {1}, {2}\n'.format(\
                            kmesh_bcc_fcc_step4[0],kmesh_bcc_fcc_step4[1],kmesh_bcc_fcc_step4[2]))
                    f2.write(' berry_curv_unit              = ang2\n')
                    f2.write(' berry_curv_adpt_kmesh        = 5\n')
                    f2.write(' berry_curv_adpt_kmesh_thresh = 100.0\n\n')
                    f2.write(' fermi_energy = {0}\n'.format(read_efermi(path)))
                else:
                    f2.write(s)

def run_step5_tb(path, type):
    # check checkpoint
    if read_checkpoint(path+'/step5_tb',0) == True:
        os.makedirs(path+'/step5_tb',exist_ok=True)
        write_checkpoint(path+'/step5_tb',0)
        #copy files to step5
        shutil.copyfile(path+'/step4/WF1.amn',path+'/step5_tb/WF1.amn')
        shutil.copyfile(path+'/step4/WF1.eig',path+'/step5_tb/WF1.eig')
        shutil.copyfile(path+'/step4/WF1.mmn',path+'/step5_tb/WF1.mmn')
        shutil.copyfile(path+'/step4/WF1.chk',path+'/step5_tb/WF1.chk')
        
        #modify WF1.win
        write_step5_tb_wannier90(path, type)

        # --------> run wannier90
        os.system('cd ' + path + '/step5_tb && srun '+ wannier90_link +'/postw90.x WF1 && cd .. && cd ..')
        # ------------------

        write_checkpoint(path+'/step5_tb',1)

# functions for band_wannier90
def write_band_wannier90(path, type): # write wannier90.win for band_w90
    # read INCAR from step1
    with open(path+'/step4/WF1.win', 'r') as f1:
        # modify INCAR and save to band_w90
        with open(path+'/band_w90/WF1.win', 'w') as f2:
            for s in f1:
                if re.search('end atoms_frac',s):
                    f2.write(' end atoms_frac\n\n')
                    f2.write(' kpath                   = true\n')
                    f2.write(' kpath_task              = bands\n')
                    f2.write(' kpath_num_points        = 400\n')
                    f2.write(' fermi_energy = {0}\n\n'.format(read_efermi(path)))
                    # kpoint_path
                    f2.write('begin kpoint_path\n')
                    if (type == 'bcc'):
                        f2.write('G  0.00000  0.00000  0.00000    H -0.50000  0.50000  0.50000\n')
                        f2.write('H -0.50000  0.50000  0.50000    N  0.00000  0.50000  0.00000\n')
                        f2.write('N  0.00000  0.50000  0.00000    P  0.25000  0.25000  0.25000\n')
                        f2.write('P  0.25000  0.25000  0.25000    G  0.00000  0.00000  0.00000\n')
                        f2.write('G  0.00000  0.00000  0.00000    N  0.00000  0.50000  0.00000\n')
                    elif (type == 'fcc'):
                        f2.write('K  0.37500  0.75000  0.37500    G  0.00000  0.00000  0.00000\n')
                        f2.write('G  0.00000  0.00000  0.00000    L  0.50000  0.50000  0.50000\n')
                        f2.write('L  0.50000  0.50000  0.50000    W  0.25000  0.75000  0.50000\n')
                        f2.write('W  0.25000  0.75000  0.50000    X  0.00000  0.50000  0.50000\n')
                        f2.write('X  0.00000  0.50000  0.50000    G  0.00000  0.00000  0.00000\n')
                    else:
                        f2.write('G  0.00000  0.00000  0.00000    M  0.00000  0.50000  0.00000\n')
                        f2.write('M  0.00000  0.50000  0.00000    K  0.33333  0.33333  0.00000\n')
                        f2.write('K  0.33333  0.33333  0.00000    G  0.00000  0.00000  0.00000\n')
                        f2.write('G  0.00000  0.00000  0.00000    A  0.00000  0.00000  0.50000\n')
                        f2.write('A  0.00000  0.00000  0.50000    L  0.00000  0.50000  0.50000\n')
                        f2.write('L  0.00000  0.50000  0.50000    H  0.33333  0.33333  0.50000\n')
                        f2.write('H  0.33333  0.33333  0.50000    A  0.00000  0.00000  0.50000\n')
                    f2.write('end kpoint_path\n')
                    f2.write('\n')
                else:
                    f2.write(s)

def run_band_w90(path, type):
    # check checkpoint
    if read_checkpoint(path+'/band_w90',0) == True:
        os.makedirs(path+'/band_w90',exist_ok=True)
        write_checkpoint(path+'/band_w90',0)
        #copy files to step3
        shutil.copyfile(path+'/step4/WF1.amn',path+'/band_w90/WF1.amn')
        shutil.copyfile(path+'/step4/WF1.eig',path+'/band_w90/WF1.eig')
        shutil.copyfile(path+'/step4/WF1.mmn',path+'/band_w90/WF1.mmn')
        shutil.copyfile(path+'/step4/WF1.chk',path+'/band_w90/WF1.chk')
        
        #write wannier90.win
        write_band_wannier90(path, type)

        # --------> run vasp
        os.system('cd ' + path + '/band_w90 && srun '+ wannier90_link +'/postw90.x WF1 && cd .. && cd ..')
        # ------------------

        write_checkpoint(path+'/band_w90',1)

        if read_checkpoint(path+'/band_w90',1) == True:
            os.system('cd ' + path + '/band_w90 && python WF1-bands.py && cd .. && cd ..')

# functions for band_fleur
def write_kpts_xml_band(path, type): # write kpts.xml for step3
    if (os.path.exists(path+'/band_fleur/kpoints') and os.path.isdir(path+'/band_fleur/kpoints')):
        shutil.rmtree(path+'/band_fleur/kpoints')
    os.makedirs(path+'/band_fleur/kpoints')
    shutil.copyfile(path+'/step2/inp.xml',path+'/band_fleur/kpoints/inp.xml')
    shutil.copyfile(path+'/step2/kpts.xml',path+'/band_fleur/kpoints/kpts.xml')
    shutil.copyfile(path+'/step2/sym.xml',path+'/band_fleur/kpoints/sym.xml')

    if (type == 'bcc'):
        knum = 900
        kpath = 'gamma=0,0,0;h=-0.5,0.5,0.5;n=0,0.5,0;p=0.25,0.25,0.25;gamma=0,0,0;n=0,0.5,0'
    elif (type == 'fcc'):
        knum = 900
        kpath = 'k=0.375,0.75,0.375;gamma=0,0,0;l=0.5,0.5,0.5;w=0.25,0.75,0.5;x=0,0.5,0.5;gamma=0,0,0'
    else:
        knum = 1200
        kpath = 'gamma=0,0,0;m=0,0.5,0;k=0.33333,0.33333,0;gamma=0,0,0;a=0,0,0.5;l=0,0.5,0.5;h=0.33333,0.33333,0.5;a=0,0,0.5'

    # --------> run inpgen
    os.system('cd ' + path + '/band_fleur/kpoints && '+ fleur_link +'/build/inpgen -inp.xml -kpt band={0} -kptsPath "{1}" \
        && cd .. && cd ..'.format(knum, kpath))

def modify_band_fleur_xml(path): # change parameters in inp.xml
    # read POSCAR from step1
    lat, atom_label, atom_num, atom_pos = read_poscar(path)
    # read inp.xml
    inp = et.parse(path+'/band_fleur/inp.xml')
    root = inp.getroot()
    # determine number of bands
    num_lo = 0
    for i in root.iter('lo'): 
        if int(i.attrib['l']) == 0:
            num_lo += 2
        elif int(i.attrib['l']) == 1:
            num_lo += 6

    # change parameters
    for i in root.iter('kPointListSelection'): i.set('listName','path-3')
    for i in root.iter('output'): 
        i.set('band','T')
    for i in root.iter('cutoffs'): i.set('numbands',str(int(num_lo+np.sum(atom_num)*36)))
    
    inp.write(path+'/band_fleur/inp.xml')

def run_band_fleur(path, type):
    # check checkpoint
    if read_checkpoint(path+'/band_fleur',0) == True:
        os.makedirs(path+'/band_fleur',exist_ok=True)
        write_checkpoint(path+'/band_fleur',0)
        # make kpts.xml
        write_kpts_xml_band(path, type)
        write_checkpoint(path+'/band_fleur',1)

    if read_checkpoint(path+'/band_fleur',1) == True:
        # copy files
        shutil.copyfile(path+'/band_fleur/kpoints/kpts.xml',path+'/band_fleur/kpts.xml')
        shutil.copyfile(path+'/step2/inp.xml',path+'/band_fleur/inp.xml')
        shutil.copyfile(path+'/step2/cdn_last.hdf',path+'/band_fleur/cdn.hdf')
        shutil.copyfile(path+'/step2/sym.xml',path+'/band_fleur/sym.xml')

        # modify inp.xml
        modify_band_fleur_xml(path)

        # --------> run fleur
        os.system('cd ' + path + '/band_fleur && srun --ntasks 1 '+ fleur_link +'/build/fleur_MPI > stdout && \
              cd .. && cd ..')
        # ------------------

        write_checkpoint(path+'/band_fleur',2)

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
        # step1
        run_step1(list_struct[i][0],list_struct[i][1],list_struct[i][2])
        # step2
        run_step2(list_struct[i][0],list_struct[i][1])
        # step3
        run_step3(list_struct[i][0], list_struct[i][1])
        # step4
        run_step4(list_struct[i][0])
        # band_w90
        run_band_w90(list_struct[i][0], list_struct[i][1])
        # band_fleur
        run_band_fleur(list_struct[i][0], list_struct[i][1])
        # step5
        run_step5(list_struct[i][0], list_struct[i][1])
        # step5_tb
        run_step5_tb(list_struct[i][0], list_struct[i][1])
        # print log
        save_log('Done.\n')
        save_log('--------------------------------------\n')
    