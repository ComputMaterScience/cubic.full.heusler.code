# import libraries
import os
import json
import requests
import pandas as pd
from string import digits
import time
from unicodedata import normalize
from bs4 import BeautifulSoup

remove_digits = str.maketrans('', '', digits)

####################### String Functions ##############################
def clean_normalize_whitespace(x):
    """ Normalize unicode characters and strip trailing spaces
    """
    if isinstance(x, str):
        return normalize('NFKC', x).strip()
    else:
        return x

####################### Get Data Functions ############################
def get_calculation_data(X, Y, Z, link):
    cnt_try = 0
    while 1:
        try:
            table_df = pd.read_html(link)
            break
        except:
            print("Connection refused by the server... Retrying ...")
            time.sleep(5)
            cnt_try += 1
            if (cnt_try > 5):
                return False, 0, 0, 0, 0, 0, 0, 0
    
    # get summary info
    df_info = table_df[0]
    # Clean up the DataFrame and Columns
    df_info = df_info.applymap(clean_normalize_whitespace)
    df_info.columns = df_info.columns.to_series().apply(clean_normalize_whitespace)

    energy   = df_info.iloc[0,1]
    volume   = df_info.iloc[1,1]
    tot_mag  = df_info.iloc[2,1]
    band_gap = df_info.iloc[3,1]

    # get magmom info
    df_mag = table_df[4]
    # Clean up the DataFrame and Columns
    df_mag = df_mag.applymap(clean_normalize_whitespace)
    df_mag.columns = df_mag.columns.to_series().apply(clean_normalize_whitespace)

    mag_X = list()
    cnt = -1
    for i in df_mag.iloc[:,0]:
        cnt += 1
        if (i.find(X) != -1):
            mag_X.append(df_mag.iloc[cnt,8])
        if (i.find(Y) != -1):
            mag_Y = df_mag.iloc[cnt,8]
        if (i.find(Z) != -1):
            mag_Z = df_mag.iloc[cnt,8]

    return True, energy, volume, tot_mag, band_gap, mag_X, mag_Y, mag_Z

def get_entry_link(id):
    cnt_try = 0
    while 1:
        try:
            table_df = pd.read_html('http://oqmd.org/materials/entry/'+str(id))
            break
        except:
            print("Connection refused by the server... Retrying ...")
            time.sleep(5)
            cnt_try += 1
            if (cnt_try > 5):
                return ''
    # get summary info
    df_info = table_df[0]
    # Clean up the DataFrame and Columns
    df_info = df_info.applymap(clean_normalize_whitespace)
    df_info.columns = df_info.columns.to_series().apply(clean_normalize_whitespace)

    cnt_try = 0
    while 1:
        try:
            response = requests.get('http://oqmd.org/materials/entry/'+str(id))
            break
        except:
            print("Connection refused by the server... Retrying ...")
            time.sleep(5)
            cnt_try += 1
            if (cnt_try > 5):
                return ''

    soup = BeautifulSoup(response.text, 'html.parser')
    table = soup.find('table')

    links = []
    for tr in table.findAll("tr"):
        trs = tr.findAll("td")
        for each in trs:
            try:
                link = each.find('a')['href']
                links.append(link)
            except:
                pass

    return str('http://oqmd.org'+links[-1])

def search_data(X, Y, Z):
    cnt_try = 0
    while 1:
        try:
            table_df = pd.read_html('http://oqmd.org/materials/composition/'+X+'2'+Y+Z)
            break
        except:
            print("Connection refused by the server... Retrying ...")
            time.sleep(5)
            cnt_try += 1
            if (cnt_try > 5):
                return False, 0, '', 0, 0, '', 0, 0, 0, 0, 0, 0, 0
    # get summary info
    df_info = table_df[0]
    # Clean up the DataFrame and Columns
    df_info = df_info.applymap(clean_normalize_whitespace)
    df_info.columns = df_info.columns.to_series().apply(clean_normalize_whitespace)
    
    # get formation energy
    if (int(df_info.iloc[0,6]) == 4):
        if ((df_info.iloc[0,2].find('Fm-3m') != -1) or (df_info.iloc[0,2].find('F-43m') != -1)):
            form_en = df_info.iloc[0,3]
            spacegroup = df_info.iloc[0,2]
            id = df_info.iloc[0,0]
            name = df_info.iloc[0,1]
            link_id = get_entry_link(id)
            if (len(link_id) == 0):
                return False, 0, '', 0, 0, '', 0, 0, 0, 0, 0, 0, 0
            else:
                status, energy, volume, tot_mag, band_gap, mag_X, mag_Y, mag_Z = get_calculation_data(X, Y, Z, link_id)
                if (status):
                    return True, id, name, energy, form_en, spacegroup, volume, tot_mag, band_gap, mag_X[0], mag_X[1], mag_Y, mag_Z
                else:
                    return False, 0, '', 0, 0, '', 0, 0, 0, 0, 0, 0, 0
        else:
            return False, 0, '', 0, 0, '', 0, 0, 0, 0, 0, 0, 0
    else:
        return False, 0, '', 0, 0, '', 0, 0, 0, 0, 0, 0, 0

def add_data(X,Y,Z):
    print(f'Get data for ',X+'2'+Y+Z)
    status, id, name, energy, form_en, spacegroup, volume, tot_mag, band_gap, mag_X1, mag_X2, mag_Y, mag_Z = search_data(X,Y,Z)
    if (status):
        print(f'Found!')
        print(f'-----------------------------------------')
        data_dict = {'id':id, 'formula':name,
                    'X':X.translate(remove_digits),'Y':Y.translate(remove_digits),'Z':Z.translate(remove_digits),
                    'energy':energy,'formation':form_en,
                    'volume':volume,'spacegroup':spacegroup,
                    'band_gap':band_gap,'tot_mag':tot_mag,
                    'X1_mag':mag_X1, 'X2_mag':mag_X2,
                    'Y_mag':mag_Y, 'Z_mag':mag_Z}
        return status, data_dict
    else:
        print(f'Not Found!')
        print(f'-----------------------------------------')
        return status, {}

######################## Main Program #################################

X_list = ['Mn','Fe','Co','Ni','Zn','Sc','Ti','V','Cr','Y','Zr','Nb','Mo', 'Hf','W','Li','Be']

Y_list = ['Mn','Fe','Co','Ni','Cu','Ru','Rh','Pd','Ag','Cd','Ir','Pt','Au','Li','Mg']

Z_list = ['B','Al','Si','Ga','Ge','As','In','Sn','Sb','Pb','Bi','Zn','Mg']

#-----------------------------
#main program
if __name__ == '__main__':
    for x in X_list:
        for y in Y_list:
            for z in Z_list:
                status, data_dict = add_data(x,y,z)
                if (status):
                    df = pd.DataFrame([data_dict])
                    df.to_csv(f'data_oqmd_{x}2{y}{z}.csv', index=False, encoding='utf-8-sig')
