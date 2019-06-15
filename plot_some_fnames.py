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
#    regime = 'attractive'
#    regime = 'usual'
#    regime = 'notail'
#    regime = 'longtail_small_coeff'
#    regime = '10tail'
    regime = '7tail'

    if regime == 'attractive':
        # Weak repulsion between modifier tail and polymer
        fnames = [
            'a34=0.5_mod_tail5_coeff3_mod_n400_poly_n50156',
            'a34=0.5_mod_tail5_coeff3_mod_n500_poly_n50096',
            'a34=0.5_mod_tail5_coeff3_mod_n600_poly_n50036',
            'a34=0.5_mod_tail5_coeff3_mod_n700_poly_n49976',
            'a34=0.5_mod_tail5_coeff3_mod_n800_poly_n49916'
        ]
        legends = [
            '400 modifiers',
            '500 modifiers',
            '600 modifiers',
            '700 modifiers',
            '800 modifiers'
        ]
        ylims = [0, 6]
        title = r'Attractive: $a_{ij}=0.5$ (for modifier tail and polymer)'
        xlabel = 'Time, ktaus'
        ylabel = r'Distance in $r_c$ ($r_c$ = 5.23$\AA$)'
        out_fname = 'attractive.pdf'
    elif regime == 'usual':
        # The most usual systems with big expansion coefficients
        fnames = [
            'usual_mod_tail5_coeff3_mod_n400_poly_n50156',
            'usual_mod_tail5_coeff3_mod_n500_poly_n50096',
            'usual_mod_tail5_coeff3_mod_n600_poly_n50036',
            'usual_mod_tail5_coeff3_mod_n700_poly_n49976',
            'usual_mod_tail5_coeff3_mod_n800_poly_n49916',
        ]
        legends = [
            '400 modifiers',
            '500 modifiers',
            '600 modifiers',
            '700 modifiers',
            '800 modifiers'
        ]
        ylims = [0, 6]
        title = 'Usual (expansion coeff = 3)'
        xlabel = 'Time, ktaus'
        ylabel = r'Distance in $r_c$ ($r_c$ = 5.23$\AA$)'
        out_fname = 'usual.pdf'
    elif regime == 'notail':
        # The most usual systems with big expansion coefficients
        fnames = [
            'notail_mod_n400_poly_n50356',
            'notail_mod_n500_poly_n50346',
            'notail_mod_n600_poly_n50336',
            'notail_mod_n700_poly_n50326',
            'notail_mod_n800_poly_n50316',
        ]
        legends = [
            '400 modifiers',
            '500 modifiers',
            '600 modifiers',
            '700 modifiers',
            '800 modifiers'
        ]
        ylims = [0, 6]
        title = 'No tail (expansion coeff = 3)'
        xlabel = 'Time, ktaus'
        ylabel = r'Distance in $r_c$ ($r_c$ = 5.23$\AA$)'
        out_fname = 'notail.pdf'
    elif regime == 'longtail_small_coeff':
        # The most usual systems with big expansion coefficients
        fnames = [
            'longtail_small_coeff_400',
            'longtail_small_coeff_500',
            'longtail_small_coeff_600',
            'longtail_small_coeff_700',
        ]
        legends = [
            '400 modifiers',
            '500 modifiers',
            '600 modifiers',
            '700 modifiers'
        ]
        ylims = [0, 10]
        title = 'Long (10 beads) tail (expansion coeff = 1.5)'
        xlabel = 'Time, ktaus'
        ylabel = r'Distance in $r_c$ ($r_c$ = 5.23$\AA$)'
        out_fname = 'longtail.pdf'
    elif regime == '10tail':
        fnames = [
            '10tail_200',
            '10tail_400',
            '10tail_600',
            '10tail_800',
        ]
        legends = [
            '200 modifiers',
            '400 modifiers',
            '600 modifiers',
            '800 modifiers'
        ]
        ylims = [0, 15]
        title = '10-beads tail (expansion coeff = 1.5)'
        xlabel = 'Time, ktaus'
        ylabel = r'Distance in $r_c$ ($r_c$ = 5.23$\AA$)'
        out_fname = '10tail.pdf'
    elif regime == '7tail':
        fnames = [
            '7tail_200',
            '7tail_400',
            '7tail_600',
            '7tail_800',
        ]
        legends = [
            '200 modifiers',
            '400 modifiers',
            '600 modifiers',
            '800 modifiers'
        ]
        ylims = [0, 15]
        title = '10-beads tail (expansion coeff = 1.5)'
        xlabel = 'Time, ktaus'
        ylabel = r'Distance in $r_c$ ($r_c$ = 5.23$\AA$)'
        out_fname = '7tail.pdf'
    else:
        print('Unknown regime!')
        print(regime)

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
