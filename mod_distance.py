# Originates from:
# https://github.com/ansko/dpd_analyzers
# /commit/de492c7789873720169df59a1bb4f4f717fae52f
# /mod_distance.py


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
    fnames = ['dfs/{0}'.format(fname) for fname in sorted(os.listdir('dfs'))]

    for idx, fname in enumerate(fnames):
        dfc = DatafileContent(fname)
        data = modifier_distance(dfc, mmt_atom_type, modifier_head_atom_type)
        if idx == 0:
            print(len(data), 'modifiers in system')
        print(idx, sum(data)/len(data))
