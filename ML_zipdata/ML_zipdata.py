import os
import shutil
import numpy as np

struct_type = ['l21','xa']
list_files = ['CHGCAR','INCAR','KPOINTS','OUTCAR','POSCAR','POTCAR']

def get_folders(folder = './'):
    return [name for name in os.listdir(folder) if os.path.isdir(os.path.join(folder, name))]

def read_outcar(path):
    # define variables
    # total energy
    tot_en = 0

    with open(path,'r') as openfileobject:
        for line in openfileobject:
            if (line.find('energy  without entropy=') != -1): # get total energy
                s = line.split('=')
                tot_en = float(s[2])

    return tot_en

#-----------------------------
#main program
if __name__ == '__main__':
    # remove old data
    if os.path.isdir('data'):
        shutil.rmtree('data', ignore_errors=True)
    # read input parameters
    list_struct = get_folders()
    num_struct = len(list_struct)
    # make folder for data storage
    os.makedirs('data')
    # loop for copying files
    for i in range(num_struct):
        en = list()
        for j in struct_type:
            en.append(read_outcar(list_struct[i]+'/'+j+'/'+'scf/OUTCAR'))
        id_min = np.argmin(en)
        # make folder for most stable structure
        os.makedirs('data/'+list_struct[i]+'/step1')
        # copy files
        for j in list_files:
            shutil.copyfile(list_struct[i]+'/'+struct_type[id_min]+'/'+'scf/'+j, 'data/'+list_struct[i]+'/step1/'+j)
    # zip files
    shutil.make_archive('data', 'zip', 'data')
    # remove old data
    if os.path.isdir('data'):
        shutil.rmtree('data', ignore_errors=True)

    