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
    fnames = ['clay_r10_n2_mod_tail2_poly_p10_n14078.data']
    fnames = ['dfs/dpd_d.{0}.data'.format(idx*1000) for idx in range(4)]
    mmt_atom_type = 1
    platelets_count = 2

    for idx, fname in enumerate(fnames):
#    for idx in range(21):
#        fname = 'dfs/dpd.{0}.data'.format((idx + 1)*5000)
        dfc = DatafileContent(fname)
        data = platelets_distance(dfc, mmt_atom_type, platelets_count)
        print(idx, data[0][1])