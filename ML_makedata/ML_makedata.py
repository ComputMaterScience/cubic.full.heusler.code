
import pandas as pd

type_struct = ['l21','xa']
num_type_struct = len(type_struct)

filename_save = 'data.csv'

dash = '/'

data = list([['atom_num_X','symbol_X','mass_X','period_X','group_X','stable_phase_X','atomic_radius_X','negativity_X','density_X','tot_val_X','s_val_X','p_val_X','d_val_X','f_val_X',\
    'atom_num_Y','symbol_Y','mass_Y','period_Y','group_Y','stable_phase_Y','atomic_radius_Y','negativity_Y','density_Y','tot_val_Y','s_val_Y','p_val_Y','d_val_Y','f_val_Y',\
    'atom_num_Z','symbol_Z','mass_Z','period_Z','group_Z','stable_phase_Z','atomic_radius_Z','negativity_Z','density_Z','tot_val_Z','s_val_Z','p_val_Z','d_val_Z','f_val_Z',\
    'ratio_XY', 'ratio_XZ', 'ratio_YZ',\
    'lat_type','lat', 'x1_mag', 'x2_mag', 'y_mag', 'z_mag', 'tot_mag', 'tot_en','formula', 'num_valence']])

column_names = data.pop(0)

def read_outcar(path):
    # define variables
    # lattice constant
    lat = 0
    # Magnetization
    x1_mag = 0
    x2_mag = 0
    y_mag = 0
    z_mag = 0
    tot_mag = 0
    # total energy
    tot_en = 0

    with open(path,'r') as openfileobject:
        for line in openfileobject:
            if (line.find('ALAT') != -1): # get lattice constant
                s = line.split('='); lat = float(s[1])
            elif (line.find('energy  without entropy=') != -1): # get total
                s = line.split('=')
                tot_en = float(s[2])
            elif (line.find('magnetization (x)') != -1): # get magnetization
                openfileobject.readline(); openfileobject.readline(); openfileobject.readline()
                s = openfileobject.readline(); s = s.split()
                x1_mag = float(s[4])
                s = openfileobject.readline(); s = s.split()
                x2_mag = float(s[4])
                s = openfileobject.readline(); s = s.split()
                y_mag = float(s[4])
                s = openfileobject.readline(); s = s.split()
                z_mag = float(s[4])
                openfileobject.readline()
                s = openfileobject.readline(); s = s.split()
                tot_mag = float(s[4])
    return lat, x1_mag, x2_mag, y_mag, z_mag, tot_mag, tot_en

# read INPUT file
f = open('INPUT','r')
num_struct = int(f.readline())
list_struct = list()
for i in range(num_struct):
    s = f.readline().split()
    list_struct.append(s[0:3])
f.close()

# read Elements file
df = pd.read_csv('Elements.csv')

for i in range(num_struct):
    info = list()
    for k in range(3):
        info.append(df[df['Symbol'] == list_struct[i][k]])
    info_X = pd.DataFrame(info[0])
    info_Y = pd.DataFrame(info[1])
    info_Z = pd.DataFrame(info[2])
    for j in range(num_type_struct):
        path = list_struct[i][0] + list_struct[i][1] + list_struct[i][2] + dash + type_struct[j] + dash + 'scf' + dash
        lat, x1_mag, x2_mag, y_mag, z_mag, tot_mag, tot_en = read_outcar(path+'OUTCAR')
        
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

        data.append([atom_num_X, symbol_X, mass_X, period_X, group_X, stable_phase_X, atomic_radius_X, negativity_X, density_X, tot_val_X, s_val_X, p_val_X, d_val_X, f_val_X,\
            atom_num_Y, symbol_Y, mass_Y, period_Y, group_Y, stable_phase_Y, atomic_radius_Y, negativity_Y, density_Y, tot_val_Y, s_val_Y, p_val_Y, d_val_Y, f_val_Y,\
            atom_num_Z, symbol_Z, mass_Z, period_Z, group_Z, stable_phase_Z, atomic_radius_Z, negativity_Z, density_Z, tot_val_Z, s_val_Z, p_val_Z, d_val_Z, f_val_Z,\
            ratio_XY, ratio_XZ, ratio_YZ,\
            type_struct[j], lat, x1_mag, x2_mag, y_mag, z_mag, tot_mag, tot_en, list_struct[i][0]+'2'+list_struct[i][1]+list_struct[i][2],total_val])

    dw = pd.DataFrame(data,columns = column_names)
    # save data to csv
    dw.to_csv(filename_save, index=False, encoding='utf-8-sig')
