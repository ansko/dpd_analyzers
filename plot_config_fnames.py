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


if __name__ == '__main__':
    regime = 'mod_tail_5_uniform'
    fnames = [
            'mod_tail_5_uniform_400',
            'mod_tail_5_uniform_500',
            'mod_tail_5_uniform_600',
            'mod_tail_5_uniform_700',
    ]
    legends = [
            '100 modifiers',
            '200 modifiers',
            '300 modifiers',
            '400 modifiers',
    ]

    '''
    regime = 'mod_tail_10'
    fnames = [
            'longtail_200',
            'longtail_400',
            'longtail_600',
            'longtail_800',
    ]
    legends = [
            '200 modifiers',
            '400 modifiers',
            '600 modifiers',
            '800 modifiers',
    ]
    '''

    ylims = [0, 6]
    title = regime
    xlabel = r'Время, $\tau$'
    ylabel = r'Расстояние в $r_c$ ($r_c$ = 5.23$\AA$)'
    out_fname = regime + '.pdf'

    plotted_lines = []
    for idx, fname in enumerate(fnames):
        xs = []
        ys = []
        for line in open('outputs/' + fname).readlines():
            xs.append(float(line.split()[0]) * 1000)
            ys.append(float(line.split()[1]))
        new_line, = plt.plot(xs, ys, label=legends[idx])
        plotted_lines.append(new_line)

    plt.title(title)
    plt.gca().set_ylim(ylims)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(plotted_lines, legends, fontsize=8)
    plt.savefig(out_fname)

    print('Finished for', regime)
