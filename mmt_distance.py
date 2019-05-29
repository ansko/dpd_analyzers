# Originates from:
# https://github.com/ansko/dpd_analyzers
# /commit/de492c7789873720169df59a1bb4f4f717fae52f
# /mmt_distance.py

'''
A tool to calculate the distance between clay platelets. Distance is 
calculated as minimum from the distances between atom from platelet 1
and atom from platelet 2. This helps to study exfoliation.
- Proper averaaging is not ready, so works only for 2-stacked MMT.
'''


import os
import sys

from datafile_content import DatafileContent


def platelets_distance(dfc, mmt_atom_type, platelets_count):
    # extract mmt atoms
    atoms = [atom for atom in dfc.atoms if atom['atom_type_id'] == mmt_atom_type]

    # map atoms onto platelets
    platelet_atoms = {idx: [] for idx in range(platelets_count)}
    for atom in atoms:
        platelet_idx = (atom['atom_id'] - 1) // (len(atoms) // platelets_count)
        platelet_atoms[platelet_idx].append({
            'x': atom['x'], 'y': atom['y'], 'z': atom['z']})
    # get data
    data = {idx_1: {idx_2: None} for idx_1 in range(platelets_count)
        for idx_2 in range(idx_1 + 1, platelets_count)}
    lx = dfc.xhi - dfc.xlo
    ly = dfc.yhi - dfc.ylo
    lz = dfc.zhi - dfc.zlo
    big_distance = max(lx, ly, lz)
    for plat_idx_1 in range(platelets_count):
        n1 = len(platelet_atoms[plat_idx_1])
        for plat_idx_2 in range(plat_idx_1 + 1, platelets_count):
            n2 = len(platelet_atoms[plat_idx_2])
            data[plat_idx_1][plat_idx_2] = 0
            for atom1 in platelet_atoms[plat_idx_1]:
                closest = big_distance
                for atom2 in platelet_atoms[plat_idx_2]:
                    dx = abs(atom1['x'] - atom2['x']); dx = min(dx, lx - dx)
                    dy = abs(atom1['y'] - atom2['y']); dy = min(dy, ly - dy)
                    dz = abs(atom1['z'] - atom2['z']); dz = min(dz, lz - dz)
                    dr = (dx**2 + dy**2 + dz**2)**0.5
                    closest = min(closest, dr)
                data[plat_idx_1][plat_idx_2] += closest / n1
    return data


if __name__ == '__main__':
    mmt_atom_type = 1
    platelets_count = 2

    #fnames = ['dfs/{0}'.format(fname) for fname in sorted(os.listdir('dfs'))]
    hd_dir = ('/media/anton/Seagate Expansion Drive/dpd_calculations/'
              'q_10/coul_cut_25/')

    if len(sys.argv) < 2:
        print('usage: ./{0} [tail_length] count'.format(__file__))
        sys.exit()
    elif len(sys.argv) > 3:
        r = int(sys.argv[1])
        tail = int(sys.argv[2])
        count = int(sys.argv[3])
    elif len(sys.argv) == 3:
        r = 10
        tail = int(sys.argv[1])
        count = int(sys.argv[2])
    else:
        r = 10
        tail = 5
        count = int(sys.argv[1])

    if tail == 5 and r == 15:
        hd_dir += {
            160: 'isolated_mmt_r15_n2_mod_n160_tail5_poly_p10_n50767',
            240: 'isolated_mmt_r15_n2_mod_n240_tail5_poly_p10_n50650',
            360: 'isolated_mmt_r15_n2_mod_n360_tail5_poly_p10_n50474'
        }[count]
    if tail == 2 and r == 10:
        hd_dir += {
            4:   'clay_r10_n2_int_5_mod_tail2_n4_poly_p10_n15043',
            10:  'clay_r10_n2_int_5_mod_tail2_n10_poly_p10_n15039',
            20:  'clay_r10_n2_int_5_mod_tail2_n20_poly_p10_n15031',
            50:  'clay_r10_n2_int_5_mod_tail2_n50_poly_p10_n15009',
            100: 'isolated_mmt_r10_n2_mod_n100_tail2_poly_p10_n14973',
            150: 'isolated_mmt_r10_n2_mod_n150_tail2_poly_p10_n14936'
        }[count]
    elif tail == 5 and r == 10:
        hd_dir += {
            4:  'isolated_mmt_r10_n2_mod_n4_tail5_poly_p10_n15040',
            10: 'isolated_mmt_r10_n2_mod_n10_tail5_poly_p10_n15031',
            20: 'isolated_mmt_r10_n2_mod_n20_tail5_poly_p10_n15017',
            40: 'isolated_mmt_r10_n2_mod_n40_tail5_poly_p10_n14987',
            50: 'isolated_mmt_r10_n2_mod_n50_tail5_poly_p10_n14973',
            60: 'isolated_mmt_r10_n2_mod_n60_tail5_poly_p10_n14958',
            70: 'isolated_mmt_r10_n2_mod_n70_tail5_poly_p10_n14944',
            80: 'isolated_mmt_r10_n2_mod_n80_tail5_poly_p10_n14900'
       }[count]
    hd_dir += '/datafiles'
    print('Tail = {0};  Count = {1}'.format(tail, count))

    fnames = ['{0}/{1}'.format(hd_dir, fname)
        for fname in sorted(os.listdir(hd_dir))]

    distance_vs_time = {}  # distances between platelets
    for idx, fname in enumerate(fnames):
        try:
            assert(len(fname.split('/')[-1].split('.')) == 3)
        except AssertionError:
            continue
        fs = int(fname.split('.')[1]) / 1000
        dfc = DatafileContent(fname)
        data = platelets_distance(dfc, mmt_atom_type, platelets_count)
        distance_vs_time[fs] = data[0][1]
        print(idx, '/', len(fnames), data[0][1], fname.split('/')[-1])

    with open('outputs/mmt_distance_tail{0}_n{1}'.format(tail,count), 'w') as f:
        for k in sorted(distance_vs_time.keys()):
            print(k, distance_vs_time[k], file=f)

    print('\a')  # Beep!
