import os
import shutil
import numpy as np

list_folders = ['band_vasp','band_w90','step4']

list_files = [['DOSCAR','EIGENVAL','INCAR','KPOINTS','OUTCAR','POSCAR','PROCAR'],\
    ['wannier90.win','wannier90-bands.dat','wannier90-bands.gnu','wannier90-bands.pdf','wannier90-bands.py','wannier90-path.kpt'],\
        ['wannier90.amn','wannier90.chk','wannier90.eig','wannier90.mmn','wannier90.win','wannier90.wpout']]

def get_folders(folder = './'):
    return [name for name in os.listdir(folder) if os.path.isdir(os.path.join(folder, name))]

#-----------------------------
#main program
if __name__ == '__main__':
    # remove old data
    if os.path.isdir('data_ahc'):
        shutil.rmtree('data_ahc', ignore_errors=True)
    # read input parameters
    list_struct = get_folders()
    num_struct = len(list_struct)
    # make folder for data storage
    os.makedirs('data_ahc')
    # loop for copying files
    for i in range(num_struct):
        # make folder for data storage
        for j in range(len(list_folders)):
            os.makedirs('data_ahc/'+list_struct[i]+'/'+list_folders[j])
            for k in list_files[j]:
            # copy files
                shutil.copyfile(list_struct[i]+'/'+list_folders[j]+'/'+k, 'data_ahc/'+list_struct[i]+'/'+list_folders[j]+'/'+k)
    # zip files
    shutil.make_archive('data_ahc', 'zip', 'data_ahc')
    # remove old data
    if os.path.isdir('data_ahc'):
        shutil.rmtree('data_ahc', ignore_errors=True)

    