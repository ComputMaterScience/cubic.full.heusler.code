# import libraries
import os
import json
import requests
import pandas as pd
from string import digits
import time

path1 = 'https://www.materialsproject.org/rest/v2/materials/'

path2 = '/vasp?API_KEY=TpNpzhEvC1i2wBwOS'

remove_digits = str.maketrans('', '', digits)

def get_data(X,Y,Z):
    print(f'Get data for ',X+Y+Z)
    cnt = 0
    while 1:
        try:
            response = requests.get(path1+X+Y+Z+path2)
            break
        except:
            print("Connection refused by the server... Retrying ...")
            time.sleep(5)
    data = response.json()
    list_data = list()
    if (data["response"]):
        data = data['response']
        for id in data:
            # check xa or l21 structures
            if (id['spacegroup']['number'] == 216 or id['spacegroup']['number'] == 225):
                cnt += 1
                list_data.append({'pretty_formula':id['pretty_formula'],
                                  'X':X.translate(remove_digits),'Y':Y.translate(remove_digits),'Z':Z.translate(remove_digits),
                                  'energy':id['energy'],'energy_per_atom':id['energy_per_atom'],
                                  'volume':id['volume'],'formation_energy_per_atom':id['formation_energy_per_atom'],
                                  'spacegroup':id['spacegroup']['number'],'band_gap':id['band_gap'],'density':id['density'],
                                  'total_magnetization':id['total_magnetization']})
    print(f'Found: ',cnt)
    print(f'-----------------------------------------')
    return list_data
    

######################## Main Program #################################

X_list = ['Mn','Fe','Co','Ni','Zn','Sc','Ti','V','Cr','Y','Zr','Nb','Mo',
          'Hf','W','Li','Be']

Y_list = ['Mn','Fe','Co','Ni','Cu','Ru','Rh','Pd','Ag','Cd','Ir','Pt','Au','Li','Mg']

Z_list = ['B','Al','Si','Ga','Ge','As','In','Sn','Sb','Pb','Bi','Zn','Mg']

#-----------------------------
#main program
if __name__ == '__main__':
    list_data = list()
    for x in X_list:
        for y in Y_list:
            for z in Z_list:
                for i in get_data(x+'2',y,z):
                    time.sleep(5)
                    list_data.append(i)
    df = pd.DataFrame(list_data)
    df.to_csv('data.csv', index=False, encoding='utf-8-sig')
