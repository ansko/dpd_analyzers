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
    regime = 'mod_tail_5_uniform'
    mmt_atom_type = 1
    hd_dir = ('/media/anton/Seagate Expansion Drive/dpd_calculations/')

    hd_dirs = [
            hd_dir + 'mod_tail_5_uniform/mod_n100_tail5_poly_n50336',
            hd_dir + 'mod_tail_5_uniform/mod_n200_tail5_poly_n50276',
            hd_dir + 'mod_tail_5_uniform/mod_n300_tail5_poly_n50216',
            hd_dir + 'mod_tail_5_uniform/mod_n400_tail5_poly_n50156',
    ]
    out_fnames = []
    for idx in range(4):
            tmp_dir = hd_dirs[idx] + '/datafiles'
            fnames = ['{0}/{1}'.format(tmp_dir, fname)
                for fname in sorted(os.listdir(tmp_dir))]
            out_fname = regime + '_' + str((idx + 4) * 100)

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
            out_fnames += 'outputs/{0}'.format(out_fname)
            with open('outputs/{0}'.format(out_fname), 'w') as f:
                for k in sorted(distance_vs_time.keys()):
                    print(k, distance_vs_time[k], file=f)

    print(out_fnames)
