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
#    regime = 'attractive'
#    regime = 'usual'
#    regime = 'notail'
#    regime = 'longtail_small_coeff'
#    regime = '10tail'
    regime = '7tail'

    mmt_atom_type = 1
    hd_dir = ('/media/anton/Seagate Expansion Drive/dpd_calculations/')

    if regime == 'attractive':
        # Weak repulsion between modifier tail and polymer
        hd_dirs = [
            hd_dir + 'a34=0.5_mod_tail5_coeff3/mod_n400_poly_n50156',
            hd_dir + 'a34=0.5_mod_tail5_coeff3/mod_n500_poly_n50096',
            hd_dir + 'a34=0.5_mod_tail5_coeff3/mod_n600_poly_n50036',
            hd_dir + 'a34=0.5_mod_tail5_coeff3/mod_n700_poly_n49976',
            hd_dir + 'a34=0.5_mod_tail5_coeff3/mod_n800_poly_n49916'
        ]
        for idx in range(4):
            tmp_dir = hd_dirs[idx] + '/datafiles'
            fnames = ['{0}/{1}'.format(tmp_dir, fname)
                for fname in sorted(os.listdir(tmp_dir))]
            out_fname = 'attractive_' + str((idx + 4) * 100)

            distance_vs_time = {}  # distances between platelets
            for idx, fname in enumerate(fnames):
                try:
                    assert(len(fname.split('/')[-1].split('.')) == 3)
                except AssertionError:
                    continue
                fs = int(fname.split('.')[-2]) / 1000
                dfc = DatafileContent(fname)
                data = platelets_distance(dfc, mmt_atom_type, 2)
                distance_vs_time[fs] = data[0][1]
                print(idx, '/', len(fnames), data[0][1], fname.split('/')[-1])

            with open('outputs/{0}'.format(out_fname), 'w') as f:
                for k in sorted(distance_vs_time.keys()):
                    print(k, distance_vs_time[k], file=f)
    elif regime == 'usual':
        # The most usual systems with big expansion coefficients
        hd_dirs = [
            hd_dir + 'coeff3_n1/mod_n400_tail5_poly_n50156',
            hd_dir + 'coeff3_n1/mod_n500_tail5_poly_n50096',
            hd_dir + 'coeff3_n1/mod_n600_tail5_poly_n50036',
            hd_dir + 'coeff3_n1/mod_n700_tail5_poly_n49976',
            hd_dir + 'coeff3_n1/mod_n800_tail5_poly_n49916'
        ]
        for idx in range(4):
            tmp_dir = hd_dirs[idx] + '/datafiles'
            fnames = ['{0}/{1}'.format(tmp_dir, fname)
                for fname in sorted(os.listdir(tmp_dir))]
            out_fname = 'usual_' + str((idx + 4) * 100)

            distance_vs_time = {}  # distances between platelets
            for idx, fname in enumerate(fnames):
                try:
                    assert(len(fname.split('/')[-1].split('.')) == 3)
                except AssertionError:
                    continue
                fs = int(fname.split('.')[-2]) / 1000
                dfc = DatafileContent(fname)
                data = platelets_distance(dfc, mmt_atom_type, 2)
                distance_vs_time[fs] = data[0][1]
                print(idx, '/', len(fnames), data[0][1], fname.split('/')[-1])

            with open('outputs/{0}'.format(out_fname), 'w') as f:
                for k in sorted(distance_vs_time.keys()):
                    print(k, distance_vs_time[k], file=f)
    elif regime == 'notail':
        # 'Modifier' has no tail, i.e. it is cation (like Na+)
        hd_dirs = [
            hd_dir + 'coeff3_n1_tail0/mod_n400_tail0_poly_n50356',
            hd_dir + 'coeff3_n1_tail0/mod_n500_tail0_poly_n50346',
            hd_dir + 'coeff3_n1_tail0/mod_n600_tail0_poly_n50336',
            hd_dir + 'coeff3_n1_tail0/mod_n700_tail0_poly_n50326',
            hd_dir + 'coeff3_n1_tail0/mod_n800_tail0_poly_n50316'
        ]
        for idx in range(4):
            tmp_dir = hd_dirs[idx] + '/datafiles'
            fnames = ['{0}/{1}'.format(tmp_dir, fname)
                for fname in sorted(os.listdir(tmp_dir))]
            out_fname = 'notail_' + str((idx + 4) * 100)

            distance_vs_time = {}  # distances between platelets
            for idx, fname in enumerate(fnames):
                try:
                    assert(len(fname.split('/')[-1].split('.')) == 3)
                except AssertionError:
                    continue
                fs = int(fname.split('.')[-2]) / 1000
                dfc = DatafileContent(fname)
                data = platelets_distance(dfc, mmt_atom_type, 2)
                distance_vs_time[fs] = data[0][1]
                print(idx, '/', len(fnames), data[0][1], fname.split('/')[-1])

            with open('outputs/{0}'.format(out_fname), 'w') as f:
                for k in sorted(distance_vs_time.keys()):
                    print(k, distance_vs_time[k], file=f)
    elif regime == 'longtail_small_coeff':
        # 'Modifier' has long tail (10 beads)
        hd_dirs = [
            hd_dir + 'mod_tail_10_small_coeff/mod_n400_tail10_poly_p10_n5804',
            hd_dir + 'mod_tail_10_small_coeff/mod_n500_tail10_poly_p10_n5694',
            hd_dir + 'mod_tail_10_small_coeff/mod_n600_tail10_poly_p10_n5584',
            hd_dir + 'mod_tail_10_small_coeff/mod_n700_tail10_poly_p10_n5474'
        ]
        for idx in range(4):
            tmp_dir = hd_dirs[idx] + '/datafiles'
            fnames = ['{0}/{1}'.format(tmp_dir, fname)
                for fname in sorted(os.listdir(tmp_dir))]
            out_fname = 'longtail_small_coeff_' + str((idx + 4) * 100)

            distance_vs_time = {}  # distances between platelets
            for idx, fname in enumerate(fnames):
                try:
                    assert(len(fname.split('/')[-1].split('.')) == 3)
                except AssertionError:
                    continue
                fs = int(fname.split('.')[-2]) / 1000
                dfc = DatafileContent(fname)
                data = platelets_distance(dfc, mmt_atom_type, 2)
                distance_vs_time[fs] = data[0][1]
                print(idx, '/', len(fnames), data[0][1], fname.split('/')[-1])

            with open('outputs/{0}'.format(out_fname), 'w') as f:
                for k in sorted(distance_vs_time.keys()):
                    print(k, distance_vs_time[k], file=f)
    elif regime == '10tail':
        counts = [200, 400, 600, 800]
        # 'Modifier' has long tail (10 beads)
        hd_dirs = [
            hd_dir + 'mod_tail_10/mod_n200_tail10_poly_n50176',
            hd_dir + 'mod_tail_10/mod_n400_tail10_poly_n49956',
            hd_dir + 'mod_tail_10/mod_n600_tail10_poly_n49736',
            hd_dir + 'mod_tail_10/mod_n800_tail10_poly_n49516'
        ]
        for idx in range(4):
            tmp_dir = hd_dirs[idx] + '/datafiles'
            fnames = ['{0}/{1}'.format(tmp_dir, fname)
                for fname in sorted(os.listdir(tmp_dir))]
            out_fname = '10tail_' + str(counts[idx])

            distance_vs_time = {}  # distances between platelets
            for idx, fname in enumerate(fnames):
                try:
                    assert(len(fname.split('/')[-1].split('.')) == 3)
                except AssertionError:
                    continue
                fs = int(fname.split('.')[-2]) / 1000
                dfc = DatafileContent(fname)
                data = platelets_distance(dfc, mmt_atom_type, 2)
                distance_vs_time[fs] = data[0][1]
                print(idx, '/', len(fnames), data[0][1], fname.split('/')[-1])

            with open('outputs/{0}'.format(out_fname), 'w') as f:
                for k in sorted(distance_vs_time.keys()):
                    print(k, distance_vs_time[k], file=f)
    elif regime == '7tail':
        counts = [200, 400, 600, 800]
        # 'Modifier' has long tail (10 beads)
        hd_dirs = [
            hd_dir + 'mod_tail_7/mod_n200_tail7_poly_n50236',
            hd_dir + 'mod_tail_7/mod_n400_tail7_poly_n50076',
            hd_dir + 'mod_tail_7/mod_n600_tail7_poly_n49916',
            hd_dir + 'mod_tail_7/mod_n800_tail7_poly_n49756'
        ]
        for idx in range(4):
            tmp_dir = hd_dirs[idx] + '/datafiles'
            fnames = ['{0}/{1}'.format(tmp_dir, fname)
                for fname in sorted(os.listdir(tmp_dir))]
            out_fname = '7tail_' + str(counts[idx])

            distance_vs_time = {}  # distances between platelets
            for idx, fname in enumerate(fnames):
                try:
                    assert(len(fname.split('/')[-1].split('.')) == 3)
                except AssertionError:
                    continue
                fs = int(fname.split('.')[-2]) / 1000
                dfc = DatafileContent(fname)
                data = platelets_distance(dfc, mmt_atom_type, 2)
                distance_vs_time[fs] = data[0][1]
                print(idx, '/', len(fnames), data[0][1], fname.split('/')[-1])

            with open('outputs/{0}'.format(out_fname), 'w') as f:
                for k in sorted(distance_vs_time.keys()):
                    print(k, distance_vs_time[k], file=f)
    else:
        print('Unknow regime')
        print(regime)
