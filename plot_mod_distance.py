# Originates from:
# https://github.com/ansko/dpd_analyzers
# /commit/???
# /plot_mmt_distance.py

'''
This is a toll to compare modifier distances.
***
Seems to be ueless for now.
'''

import sys

import matplotlib as mpl
import matplotlib.pyplot as plt


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Tail length not specified!')
        sys.exit()
    elif int(sys.argv[1]) == 2:
        fnames = [
            'mmt_distance_tail2_n50',
            'mmt_distance_tail2_n20',
            'mod_distance_tail2_n10',
            'mod_distance_tail2_n4']
        legends = ['N = 50', 'N = 20', 'N = 10', 'N = 4']
        ylims = [0, 10]
    elif int(sys.argv[1]) == 5:
        fnames = [
            'mod_distance_tail5_n80',
            'mod_distance_tail5_n20',
            'mod_distance_tail5_n10',
            'mod_distance_tail5_n4']
        legends = ['N = 80', 'N = 20', 'N = 10', 'N = 4']
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

    plt.title('Tail = {0}'.format(sys.argv[1]))
    plt.gca().set_ylim(ylims)

    plt.xlabel('Steps')
    plt.ylabel(r'Distance in $r_c$ ($r_c$ = 1.35$\AA$)')

    plt.legend(plotted_lines, legends)

    plt.savefig('mod_distances_tail{0}.eps'.format(sys.argv[1]))
