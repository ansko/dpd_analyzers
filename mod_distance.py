# Originates from:
# https://github.com/ansko/dpd_analyzers
# /commit/de492c7789873720169df59a1bb4f4f717fae52f
# /mod_distance.py

'''
The idea was to study whether there are some modifiers far from the
clay (i.e. CEC is too high).
***
It seems, this is fail.
'''

from datafile_content import DatafileContent


def modifier_distance(dfc, mmt_atom_type, modifier_head_atom_type):
    '''
    Get the distances from the modifier molecules to the
    closest clay platelet.
    '''

    # extract mmt atoms
    mmt_atoms = [atom
        for atom in dfc.atoms if atom['atom_type_id'] == mmt_atom_type]

    # extract modifier head atoms
    mod_atoms = [atom
        for atom in dfc.atoms if atom['atom_type_id'] == modifier_head_atom_type]

    # get data
    data = []  # distances

    lx = dfc.xhi - dfc.xlo
    ly = dfc.yhi - dfc.ylo
    lz = dfc.zhi - dfc.zlo
    big_distance = max(lx, ly, lz)
    for mod_atom in mod_atoms:
        min_distance = big_distance
        for mmt_atom in mmt_atoms:
            dx = abs(mod_atom['x'] - mmt_atom['x']); dx = min(dx, lx - dx)
            dy = abs(mod_atom['y'] - mmt_atom['y']); dy = min(dy, ly - dy)
            dz = abs(mod_atom['z'] - mmt_atom['z']); dz = min(dz, lz - dz)
            dr = (dx**2 + dy**2 + dz**2)**0.5
            min_distance = min(min_distance, dr)
        data.append(min_distance)
    return data


if __name__ == '__main__':
    mmt_atom_type = 1
    modifier_head_atom_type = 2

    import os
    #fnames = ['dfs/{0}'.format(fname) for fname in sorted(os.listdir('dfs'))]
    hd_dir = ('/media/anton/Seagate Expansion Drive/dpd_calculations/'
              'q_10/coul_cut_25/')


#    hd_dir += 'clay_r10_n2_int_5_mod_tail2_n4_poly_p10_n15043'
#    hd_dir += 'clay_r10_n2_int_5_mod_tail2_n10_poly_p10_n15039'
#    hd_dir += 'clay_r10_n2_int_5_mod_tail2_n20_poly_p10_n15031'
#    hd_dir += 'clay_r10_n2_int_5_mod_tail2_n50_poly_p10_n15009'

#    hd_dir += 'isolated_mmt_r10_n2_mod_n4_tail5_poly_p10_n15040'
#    hd_dir += 'isolated_mmt_r10_n2_mod_n10_tail5_poly_p10_n15031'
#    hd_dir += 'isolated_mmt_r10_n2_mod_n20_tail5_poly_p10_n15017'
    hd_dir += 'isolated_mmt_r10_n2_mod_n80_tail5_poly_p10_n14900'

    hd_dir += '/datafiles'

    fnames = ['{0}/{1}'.format(hd_dir, fname)
        for fname in sorted(os.listdir(hd_dir))]

    distance_vs_time = {}

    for idx, fname in enumerate(fnames):
        try:
            assert(len(fname.split('/')[-1].split('.')) == 3)
        except AssertionError:
            continue
        fs = int(fname.split('.')[1]) / 1000
        dfc = DatafileContent(fname)
        data = modifier_distance(dfc, mmt_atom_type, modifier_head_atom_type)
        if idx == 0:
            print(len(data), 'modifiers in system')

        distance_vs_time[fs] = sum(data)/len(data)
        print(idx, '/', len(fnames), sum(data)/len(data))


    with open('mod_distance_tail5_n80', 'w') as f:
        for k in sorted(distance_vs_time.keys()):
            print(k, distance_vs_time[k], file=f)
