# import libraries
import pandas as pd
import os
import numpy as np

data_csv = list([['atom_num_X','symbol_X','mass_X','period_X','group_X','stable_phase_X','atomic_radius_X','negativity_X','density_X','tot_val_X','s_val_X','p_val_X','d_val_X','f_val_X',\
    'atom_num_Y','symbol_Y','mass_Y','period_Y','group_Y','stable_phase_Y','atomic_radius_Y','negativity_Y','density_Y','tot_val_Y','s_val_Y','p_val_Y','d_val_Y','f_val_Y',\
    'atom_num_Z','symbol_Z','mass_Z','period_Z','group_Z','stable_phase_Z','atomic_radius_Z','negativity_Z','density_Z','tot_val_Z','s_val_Z','p_val_Z','d_val_Z','f_val_Z',\
    'ratio_XY', 'ratio_XZ', 'ratio_YZ',\
    'lat_type','lat', 'x1_mag', 'x2_mag', 'y_mag', 'z_mag', 'tot_mag', 'tot_en','form_en','formula', 'num_valence']])

column_names = data_csv.pop(0)

######################## IO functions #################################

def get_filenames(folder):
    # get filenames in database folder
    list_files = list()
    for path in os.listdir(folder):
        full_path = os.path.join(folder, path)
        if os.path.isfile(full_path):
            list_files.append(full_path)
    return list_files

def combineData(input,output):
    # combine csv files into one csv file
    list_files = get_filenames(input)
    #combine all files in the list
    combined_csv = pd.concat([pd.read_csv(f) for f in list_files ])
    #export to csv
    combined_csv.to_csv(output, index=False, encoding='utf-8-sig')

def processing_data(input,output):
    # read csv file
    data = pd.read_csv(input)
    data.pop('band_gap')

    # convert volume to lattice constant
    data['volume'] = np.power(data['volume']*4*4,1/3)
    data.rename({'volume': 'lat'}, axis=1, inplace=True)

    # set na value to zero
    data.loc[pd.isnull(data['tot_mag']),'tot_mag'] = 0
    data['tot_mag'] = data['tot_mag']*4

    data.to_csv(output, index=False, encoding='utf-8-sig')

def make_input(input,output):
    # read OMDB file
    data_omdb = pd.read_csv(input)
    # read Elements file
    df = pd.read_csv('Elements.csv',encoding='ANSI')
    list_struct = list()
    num_struct = len(data_omdb.index)
    for i in range(num_struct):
        list_struct.append([data_omdb.iloc[i,2],data_omdb.iloc[i,3],data_omdb.iloc[i,4]])

    # add data to list
    for i in range(num_struct):
        info = list()
        for k in range(3):
            info.append(df[df['Symbol'] == list_struct[i][k]])
        info_X = pd.DataFrame(info[0])
        info_Y = pd.DataFrame(info[1])
        info_Z = pd.DataFrame(info[2])

        lat      = data_omdb.loc[i,'lat']
        x1_mag   = data_omdb.loc[i,'X1_mag']
        x2_mag   = data_omdb.loc[i,'X2_mag']
        y_mag    = data_omdb.loc[i,'Y_mag']
        z_mag    = data_omdb.loc[i,'Z_mag']
        tot_mag  = data_omdb.loc[i,'tot_mag']
        tot_en   = data_omdb.loc[i,'energy']
        form_en  = data_omdb.loc[i,'formation']

        if (data_omdb.loc[i,'spacegroup'].find('Fm-3m') != -1):
            lat_type = 'l21'
        else:
            lat_type = 'xa'
        
        # get info atom X
        atom_num_X = info_X['Atomic Number'].values[0]
        symbol_X = info_X['Symbol'].values[0]
        mass_X = info_X['Atomic Weight'].values[0]
        period_X = info_X['Period'].values[0]
        group_X = info_X['Group'].values[0]
        stable_phase_X = info_X['Most Stable Crystal'].values[0]
        atomic_radius_X = info_X['Atomic Radius'].values[0]
        negativity_X = info_X['Electronegativity'].values[0]
        density_X = info_X['Density'].values[0]
        num_val_s = info_X['s Valence'].values[0]
        num_val_p = info_X['p Valence'].values[0]
        num_val_d = info_X['d Valence'].values[0]
        num_val_f = info_X['f Valence'].values[0]
        if ((num_val_d == 10) & (num_val_p > 0)):
            num_val_tot = num_val_s + num_val_p
        else:
            num_val_tot = num_val_s + num_val_p + num_val_d
        tot_val_X = num_val_tot
        s_val_X = num_val_s
        p_val_X = num_val_p
        d_val_X = num_val_d
        f_val_X = num_val_f

        # get info atom Y
        atom_num_Y = info_Y['Atomic Number'].values[0]
        symbol_Y = info_Y['Symbol'].values[0]
        mass_Y = info_Y['Atomic Weight'].values[0]
        period_Y = info_Y['Period'].values[0]
        group_Y = info_Y['Group'].values[0]
        stable_phase_Y = info_Y['Most Stable Crystal'].values[0]
        atomic_radius_Y = info_Y['Atomic Radius'].values[0]
        negativity_Y = info_Y['Electronegativity'].values[0]
        density_Y = info_Y['Density'].values[0]
        num_val_s = info_Y['s Valence'].values[0]
        num_val_p = info_Y['p Valence'].values[0]
        num_val_d = info_Y['d Valence'].values[0]
        num_val_f = info_Y['f Valence'].values[0]
        if ((num_val_d == 10) & (num_val_p > 0)):
            num_val_tot = num_val_s + num_val_p
        else:
            num_val_tot = num_val_s + num_val_p + num_val_d
        tot_val_Y = num_val_tot
        s_val_Y = num_val_s
        p_val_Y = num_val_p
        d_val_Y = num_val_d
        f_val_Y = num_val_f

        # get info atom Z
        atom_num_Z = info_Z['Atomic Number'].values[0]
        symbol_Z = info_Z['Symbol'].values[0]
        mass_Z = info_Z['Atomic Weight'].values[0]
        period_Z = info_Z['Period'].values[0]
        group_Z = info_Z['Group'].values[0]
        stable_phase_Z = info_Z['Most Stable Crystal'].values[0]
        atomic_radius_Z = info_Z['Atomic Radius'].values[0]
        negativity_Z = info_Z['Electronegativity'].values[0]
        density_Z = info_Z['Density'].values[0]
        num_val_s = info_Z['s Valence'].values[0]
        num_val_p = info_Z['p Valence'].values[0]
        num_val_d = info_Z['d Valence'].values[0]
        num_val_f = info_Z['f Valence'].values[0]
        if ((num_val_d == 10) & (num_val_p > 0)):
            num_val_tot = num_val_s + num_val_p
        else:
            num_val_tot = num_val_s + num_val_p + num_val_d
        tot_val_Z = num_val_tot
        s_val_Z = num_val_s
        p_val_Z = num_val_p
        d_val_Z = num_val_d
        f_val_Z = num_val_f

        # get ratio of atomic radius
        ratio_XY = atomic_radius_X/atomic_radius_Y
        ratio_XZ = atomic_radius_X/atomic_radius_Z
        ratio_YZ = atomic_radius_Y/atomic_radius_Z

        # convert magnetic moments total valence < 24
        total_val = tot_val_X*2+tot_val_Y+tot_val_Z
        if (total_val < 24):
            if (tot_mag > 0):
                tot_mag *= -1
                x1_mag *= -1
                x2_mag *= -1
                y_mag  *= -1
                z_mag  *= -1
        else:
            if (tot_mag < 0):
                tot_mag *= -1
                x1_mag *= -1
                x2_mag *= -1
                y_mag  *= -1
                z_mag  *= -1

        data_csv.append([atom_num_X, symbol_X, mass_X, period_X, group_X, stable_phase_X, atomic_radius_X, negativity_X, density_X, tot_val_X, s_val_X, p_val_X, d_val_X, f_val_X,\
            atom_num_Y, symbol_Y, mass_Y, period_Y, group_Y, stable_phase_Y, atomic_radius_Y, negativity_Y, density_Y, tot_val_Y, s_val_Y, p_val_Y, d_val_Y, f_val_Y,\
            atom_num_Z, symbol_Z, mass_Z, period_Z, group_Z, stable_phase_Z, atomic_radius_Z, negativity_Z, density_Z, tot_val_Z, s_val_Z, p_val_Z, d_val_Z, f_val_Z,\
            ratio_XY, ratio_XZ, ratio_YZ,\
            lat_type, lat, x1_mag, x2_mag, y_mag, z_mag, tot_mag, tot_en, form_en, list_struct[i][0]+'2'+list_struct[i][1]+list_struct[i][2],total_val])

    dw = pd.DataFrame(data_csv,columns = column_names)
    # save data to csv
    dw.to_csv(output, index=False, encoding='utf-8-sig')

#-----------------------------
#main program
if __name__ == '__main__':
    combineData('databases','oqmd_data.csv')
    processing_data('oqmd_data.csv','processed_oqmd_data.csv')
    make_input('processed_oqmd_data.csv','oqmd_data_plus.csv')