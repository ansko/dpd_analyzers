# Originates from:
# https://github.com/ansko/dpd_analyzers
# /commit/de492c7789873720169df59a1bb4f4f717fae52f
# /mod_distance.py


'''
The idea was to study whether there are some modifiers far from the
clay (i.e. CEC is too high). This is an improvement if mod_distance
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


def modifiers_further(dfc, mmt_atom_type, modifier_head_atom_type, threshold):
    '''
    Get count of modifiers that are further from the clay than threshold.
    '''

    # extract mmt atoms
    mmt_atoms = [atom
        for atom in dfc.atoms if atom['atom_type_id'] == mmt_atom_type]

    # extract modifier head atoms
    mod_atoms = [atom
        for atom in dfc.atoms if atom['atom_type_id'] == modifier_head_atom_type]

    # get data
    count = 0

    lx = dfc.xhi - dfc.xlo
    ly = dfc.yhi - dfc.ylo
    lz = dfc.zhi - dfc.zlo
    #big_distance = max(lx, ly, lz)
    for mod_atom in mod_atoms:
        #min_distance = big_distance
        for mmt_atom in mmt_atoms:
            dx = abs(mod_atom['x'] - mmt_atom['x']); dx = min(dx, lx - dx)
            dy = abs(mod_atom['y'] - mmt_atom['y']); dy = min(dy, ly - dy)
            dz = abs(mod_atom['z'] - mmt_atom['z']); dz = min(dz, lz - dz)
            dr = (dx**2 + dy**2 + dz**2)**0.5
            if dr > threshold:
                count += 1
    return count


if __name__ == '__main__':
    mmt_atom_type = 1
    modifier_head_atom_type = 2

    import os
    #fnames = ['dfs/{0}'.format(fname) for fname in sorted(os.listdir('dfs'))]
    hd_dir = ('/media/anton/Seagate Expansion Drive/dpd_calculations/'
              'q_10/coul_cut_25/')

#    hd_dir += 'isolated_mmt_r10_n2_mod_n4_tail5_poly_p10_n15040'
    hd_dir += 'isolated_mmt_r10_n2_mod_n10_tail5_poly_p10_n15031'
#    hd_dir += 'isolated_mmt_r10_n2_mod_n20_tail5_poly_p10_n15017'
#    hd_dir += 'isolated_mmt_r10_n2_mod_n80_tail5_poly_p10_n14900'

    fnames = ['{0}/{1}'.format(hd_dir, fname)
        for fname in sorted(os.listdir(hd_dir))]

    count_vs_time = {}

    count = []
    for idx, fname in enumerate(fnames):
        fs = int(fname.split('.')[1])
        try:
            assert(len(fname.split('/')[-1].split('.')) == 3)
        except AssertionError:
            continue
        dfc = DatafileContent(fname)
        data = modifier_distance(dfc, mmt_atom_type, modifier_head_atom_type)
        if idx == 0:
            print(len(data), 'modifiers in system')

        threshold = sum(data)/len(data) * 3  # why so? why not so?
        count = modifiers_further(dfc, mmt_atom_type, modifier_head_atom_type,
            threshold)

        count_vs_time[idx] = count
        print(idx, '/', len(fnames), fs, data[0][1], count)
        counts.append(count)

    print('Average:', sum(counts) / len(count))

    with open('mod_count_tail2_n4', 'w') as f:
        for k in sorted(count_vs_time.keys()):
            print(k, count_vs_time[k], file=f)
