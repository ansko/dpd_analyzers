'''
This is done to analyze the density in the MMT interlayer since there 
was a guess that the exfoliation may occur due to the very high value 
of the interlayer density. Interlayer is the space between parallel 
platelets with dx**2 + dy**2 < r**2 where dx, dy are distances from 
the platelets' centers; z = Const for platelets.
- This is done only for 2-stacked platelet now
- Both platelets are parallel to xy plane
***
This quess was not confirmed for now, but this functionality may be
useful later.
'''


import math
import os
import sys

from datafile_content import DatafileContent


if __name__ == '__main__':
    mmt_atom_type = 1
    modifier_head_type = 2
    platelets_count = 2
    r_c = 1  # in lj units!
    R = 10
    R2 = R**2 * r_c**2 * 4

    hd_dir = ('/media/anton/Seagate Expansion Drive/dpd_calculations/'
              'q_10/coul_cut_25/')

    if len(sys.argv) < 2:
        print('usage: ./{0} [platelet_radius] [tail_length] count'.format(
            __file__))
        sys.exit()
    elif len(sys.argv) == 4:
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

    if tail == 2:
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
    elif tail == 5 and r == 15:
        hd_dir += {
            100: 'isolated_mmt_r15_n2_mod_n100_tail5_poly_p10_n50855',
            160: 'isolated_mmt_r15_n2_mod_n160_tail5_poly_p10_n50767',
            240: 'isolated_mmt_r15_n2_mod_n240_tail5_poly_p10_n50650',
            360: 'isolated_mmt_r15_n2_mod_n360_tail5_poly_p10_n50474'
        }[count]
    hd_dir += '/datafiles'
    print('Tail = {0};  Count = {1}'.format(tail, count))

    #fname = '{0}/dpd_d.0.data'.format(hd_dir)
    #fname = '{0}/dpd_d.1000.data'.format(hd_dir)
    fname = '{0}/dpd_d.5000.data'.format(hd_dir)
    #fname = '{0}/dpd_d.10000.data'.format(hd_dir)

    dfc = DatafileContent(fname)

    plat_thickness = 5 * r_c  # really 4 not 5, 1 is for margin
    mmt_atoms = [atom for atom in dfc.atoms
        if atom['atom_type_id'] == mmt_atom_type]
    mod_head_atoms = [atom for atom in dfc.atoms
        if atom['atom_type_id'] == modifier_head_type]

    min_plats_x = dfc.xhi
    max_plats_x = dfc.xlo
    min_plats_y = dfc.yhi
    max_plats_y = dfc.ylo
    # TODO: not only stacking == 2
    ave_z_lower = 0
    ave_z_higher = 0
    for atom in mmt_atoms:
        if atom['z'] < 0:
            ave_z_lower += atom['z']
        else:
            ave_z_higher += atom['z']
        min_plats_x = min(min_plats_x, atom['x'])
        max_plats_x = max(max_plats_x, atom['x'])
        min_plats_y = min(min_plats_y, atom['y'])
        max_plats_y = max(max_plats_y, atom['y'])
    ave_z_lower /= len(mmt_atoms) / 2
    ave_z_higher /= len(mmt_atoms) / 2
    plats_xc = (min_plats_x + max_plats_x) / 2
    plats_yc = (min_plats_y + max_plats_y) / 2

    interlayer_atoms = 0
    for atom in dfc.atoms:
        if atom['atom_type_id'] == mmt_atom_type:
            continue
        if (atom['x'] - plats_xc)**2 + (atom['y'] - plats_yc)**2 > R2:
            continue
        if atom['atom_type_id'] != modifier_head_type:
            continue
        if ave_z_lower <= atom['z'] <= ave_z_higher:
            interlayer_atoms += 1

    interlayer_volume = math.pi * R2 * (ave_z_higher - ave_z_lower - 4*r_c)
    dpd_rho = interlayer_atoms / interlayer_volume

    print('interlayer_atoms:', interlayer_atoms)
    print('interlayer_volume:', interlayer_volume)
    print('dpd_rho:', dpd_rho)

    lx = dfc.xhi - dfc.xlo
    ly = dfc.yhi - dfc.ylo
    lz = dfc.zhi - dfc.zlo
    print('overall rho:', len(dfc.atoms) / lx / ly / lz)
    print(lx, ly, lz)
