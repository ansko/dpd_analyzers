# Originates from:
# https://github.com/ansko/dpd_analyzers
# /commit/???
# /plot_mmt_distance.py

'''
This is a toll to compare exfoliation tempos.
'''


import math
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt


def cec(r, modifiers, platelets=2, bead_radius=1.35):
    mmt_square = math.pi * (bead_radius*r)**2 * platelets
    modifiers_should_be = mmt_square * 0.015
    return int(92.6 * modifiers / modifiers_should_be)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Tail length not specified!')
        sys.exit()
    elif len(sys.argv) == 3 and int(sys.argv[1]) == 15 and int(sys.argv[2]) == 5:
        fnames = [
            'outputs/mmt_distance_tail5_n360',
            'outputs/mmt_distance_tail5_n240',
            'outputs/mmt_distance_tail5_n160',
            'outputs/mmt_distance_tail5_n100'
        ]
        legends = [
            'N = 360 (26->31->51), CEC = {0}'.format(cec(r=15, modifiers=360)),
            'N = 240 (16->21->25), CEC = {0}'.format(cec(r=15, modifiers=240)),
            'N = 160 (10->15->17), CEC = {0}'.format(cec(r=15, modifiers=160)),
            'N = 100 (3->13->15), CEC = {0}'.format(cec(r=15, modifiers=100))

        ]
        ylims = [0, 25]
    elif int(sys.argv[1]) == 2:
        fnames = [
            'outputs/mmt_distance_tail2_n150',
            'outputs/mmt_distance_tail2_n100',
            'outputs/mmt_distance_tail2_n50',
            'outputs/mmt_distance_tail2_n20',
            'outputs/mmt_distance_tail2_n10',
            'outputs/mmt_distance_tail2_n4'
        ]
        legends = [
            'N = 150 (32->75->46), CEC = {0}'.format(cec(r=10, modifiers=150)),
            'N = 100 (19->48->41), CEC = {0}'.format(cec(r=10, modifiers=100)),
            'N = 50 (12->?->25), CEC = {0}'.format(cec(r=10, modifiers=50)),
            'N = 20 (2->?->5), CEC = {0}'.format(cec(r=10, modifiers=20)),
            'N = 10 (3-?->2), CEC = {0}'.format(cec(r=10, modifiers=10)),
            'N = 4 (0-?->2), CEC = {0}'.format(cec(r=10, modifiers=4))
        ]
        ylims = [0, 25]
    elif int(sys.argv[1]) == 5:
        fnames = [
            'outputs/mmt_distance_tail5_n80',
            'outputs/mmt_distance_tail5_n70',
            'outputs/mmt_distance_tail5_n60',
            'outputs/mmt_distance_tail5_n50',
            'outputs/mmt_distance_tail5_n40',
            'outputs/mmt_distance_tail5_n20',
            'outputs/mmt_distance_tail5_n10',
            'outputs/mmt_distance_tail5_n4']
        legends = [
            'N = 80 (16->39->39), CEC = {0}'.format(cec(r=10, modifiers=80)),
            'N = 70 (7->23->31), CEC = {0}'.format(cec(r=10, modifiers=70)),
            'N = 60 (8->25->30), CEC = {0}'.format(cec(r=10, modifiers=60)),
            'N = 50 (11->26->27), CEC = {0}'.format(cec(r=10, modifiers=50)),
            'N = 40 (6->13->16), CEC = {0}'.format(cec(r=10, modifiers=40)),
            'N = 20 (3->5->8), CEC = {0}'.format(cec(r=10, modifiers=20)),
            'N = 10 (1->2->5), CEC = {0}'.format(cec(r=10, modifiers=10)),
            'N = 4 (1->1->1), CEC = {0}'.format(cec(r=10, modifiers=4))
        ]
        ylims = [0, 35]
    else:
        print('Unknown value for tail length:', sys.argv[1])
        sys.exit()

    plotted_lines = []
    for idx, fname in enumerate(fnames):
        xs = []
        ys = []
        for line in open(fname).readlines():
            xs.append(float(line.split()[0]) * 1000)
            ys.append(float(line.split()[1]))
        new_line, = plt.plot(xs, ys, label=legends[idx])
        plotted_lines.append(new_line)

    if len(sys.argv) == 3:
        plt.title('R = {0}, tail = {0}'.format(sys.argv[1], sys.argv[2]))
    else:
        plt.title('Tail = {0}'.format(sys.argv[1]))

    plt.gca().set_ylim(ylims)

    plt.xlabel('Steps')
    plt.ylabel(r'Distance in $r_c$ ($r_c$ = 1.35$\AA$)')

    plt.legend(plotted_lines, legends, fontsize=8)

    if len(sys.argv) == 3:
        plt.savefig('figures/mmt_distances_tail{0}_r{1}.eps'.format(
            sys.argv[2], sys.argv[1]))
    else:
        plt.savefig('figures/mmt_distances_tail{0}.eps'.format(sys.argv[1]))
