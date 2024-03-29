#!/usr/bin/env python3

"""A script to create the Kohomoto butterfly plot.

usage: Butterfly.py [-h] [-l LINEWIDTH] [--driver DRIVER] [-d DPI] [-s SAVE_EVERY] [-v V]
                    [--merge] [--merge-prefix MERGE_PREFIX] [--merge-per MERGE_PER] n

Plot the Kohomoto butterfly and save as a PNG

positional arguments:
  n                     Use the nth farey sequence as the alphas

optional arguments:
  -h, --help            show this help message and exit
  -l LINEWIDTH, --linewidth LINEWIDTH
                        line width for each spectrum line
  --driver DRIVER       eigenvalue driver (Read more in scipy.linalg.eigvalsh docs)
  -d DPI, --dpi DPI     output image dpi
  -s SAVE_EVERY, --save-every SAVE_EVERY
                        save every x alphas, -1 for a single output
  -v V                  Butterfly potential amplitude
  --merge               Prints the command to merge the split plot using ImageMagick CLI
  --merge-prefix MERGE_PREFIX
                        Add a prefix to the png images
  --merge-per MERGE_PER
                        ImageMagick takes lots of RAM so we do multiple merges
"""

__author__ = "Barak Biber"
__license__ = "MIT"
__copyright__ = "Copyright (c) 2020 Barak Biber"


import os
import tqdm
import contfrac
import argparse
import sortednp
import fractions
import numpy as np
import scipy.sparse
import scipy.linalg
import scipy.special
import scipy.optimize
import matplotlib.pyplot as plt
from timeit import default_timer as timer


def get_matrix_for_alpha(p, q, v=1):
    alpha = fractions.Fraction(p,q)
    p = alpha.numerator
    q = alpha.denominator
    na = np.arange(q) * alpha 
    na -= np.floor(na)
    v_n = np.where(na >= 1-alpha, v, 0)
    return scipy.sparse.diags([1, v_n, 1], [-1, 0, 1], shape=(q, q), dtype=np.int32)


def get_intervals_from_matrix(_mat, driver='ev'):
    _q = _mat.shape[0]
    
    diag_bl = np.diag([1], -_q+1)
    diag_tr = np.diag([1], _q-1)

    _mat += diag_bl
    _mat += diag_tr
    eig_zero = scipy.linalg.eigvalsh(_mat, driver=driver)
    _mat -= 2*diag_bl
    _mat -= 2*diag_tr
    eig_pi = scipy.linalg.eigvalsh(_mat, driver=driver)
    _mat += diag_bl
    _mat += diag_tr

    _intervals = sortednp.merge(eig_zero, eig_pi)
    _intervals.shape = (_q, 2)

    return _intervals

def get_intervals_from_alpha(alpha, v=1, driver='ev'):
    mat = get_matrix_for_alpha(alpha.numerator, alpha.denominator, v)
    return get_intervals_from_matrix(mat, driver=driver)

def get_spectrum_min_max(_mat):
    _q = _mat.shape[0]
    
    diag_bl = np.diag([1], -_q+1)
    diag_tr = np.diag([1], _q-1)

    _mat += diag_bl
    _mat += diag_tr
    min_zero = scipy.linalg.eigvalsh(_mat, subset_by_index=[0, 0], driver='evr')[0]
    max_zero = scipy.linalg.eigvalsh(_mat, subset_by_index=[_q-1, _q-1], driver='evr')[0]
    _mat -= 2*diag_bl
    _mat -= 2*diag_tr
    min_pi = scipy.linalg.eigvalsh(_mat, subset_by_index=[0, 0], driver='evr')[0]
    max_pi = scipy.linalg.eigvalsh(_mat, subset_by_index=[_q-1, _q-1], driver='evr')[0]
    _mat += diag_bl
    _mat += diag_tr
    
    _max = np.maximum(max_zero, max_pi)
    _min = np.maximum(min_zero, min_pi)

    return _min, _max


def set_butterfly_axis(ax, v):
    ax.set_yticks(np.linspace(0, 1, 21))
    if v > 0:
        ax.set_xlim([-2.05, 2.05+v])
    else:
        ax.set_xlim([-2.05+v, 2.05])

    ax.set_ylim([-0.05, 1.05])
    ax.set_ylabel('Frequencies \u03B1')
    ax.set_xlabel('Spectrum')


def is_interval_contained(interval, intervals):
    for a, b in intervals:
        if interval[0] >= a and interval[1] <= b:
            return True

    return False


def get_alpha_1_2(alpha, coefficients=None):
    if coefficients is None:
        coefficients = list(contfrac.continued_fraction(alpha))
    alpha_1 = contfrac_to_fraction(coefficients[:-1])
    coefficients[-1] -= 1
    alpha_2 = contfrac_to_fraction(coefficients)

    return alpha_1, alpha_2


def plot_alpha(alpha, ax, linewidth=0.1, v=1, driver='ev', coloring_scheme=None, overide_y=None, coeff=None, verbose=False):
    should_color = (coloring_scheme is not None)
    color_by_up_down = (coloring_scheme == 'up_down')

    color = 'k'
    if should_color:
        alpha_1, alpha_2 = get_alpha_1_2(alpha, coefficients=coeff)
        alpha_1_intervals = get_intervals_from_alpha(alpha_1, v=v, driver=driver)
        alpha_2_intervals = get_intervals_from_alpha(alpha_2, v=v, driver=driver)
    intervals = get_intervals_from_alpha(alpha, v=v, driver=driver)
    if verbose:
        print(f"{alpha} has {len(intervals)} bands")
    for x1, x2 in intervals:
        interval = (x1, x2)
        if color_by_up_down: 
            if is_interval_contained(interval, alpha_1_intervals):
                if alpha < alpha_1:
                    color = 'b'
                else:
                    color = 'r'
            elif is_interval_contained(interval, alpha_2_intervals):
                if alpha < alpha_2:
                    color = 'b'
                else:
                    color = 'r'
            else:
                color = 'k'
        elif should_color:
            if is_interval_contained(interval, alpha_1_intervals):
                color = 'r'
                if verbose:
                    print('A', end=' ')
            elif is_interval_contained(interval, alpha_2_intervals):
                color = 'b'
                if verbose:
                    print('B', end=' ')
            else:
                color = 'k'
                if verbose:
                    print('?', end=' ')

        if overide_y is None:
            y = alpha
        else:
            y = overide_y


        ax.plot(interval, [y, y], color + '-', linewidth=linewidth)

    if verbose:
        print()

def plot_butterfly(alphas, name, *, linewidth=0.1, dpi=400, v=1, driver='evr', save_every=500, coloring_scheme=None, ax=None):
    if coloring_scheme is not None:
        name = name + '_colored_' + coloring_scheme
    if not os.path.exists(f'output/{name}'):
        os.mkdir(f'output/{name}')
    saved = False
    do_not_save = ax is not None

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 10))
        set_butterfly_axis(ax, v)
    else:
        save_every = -1        

    ax.plot([-2, 2], [0, 0], 'k-', linewidth=linewidth)
    ax.plot([-2 + v, 2 + v], [1, 1], 'k-', linewidth=linewidth)

    if save_every != -1:
        fig.savefig(f'output/{name}/background.png', dpi=dpi)

        ax.cla()
        set_butterfly_axis(ax, v)

    counter = 1
    for i, alpha in enumerate(tqdm.tqdm(alphas)):
        saved = False
        plot_alpha(alpha, ax=ax, linewidth=linewidth, v=v, driver=driver, coloring_scheme=coloring_scheme)
        if save_every != -1 and i > 0 and i % save_every == 0:
            fig.savefig(f'output/{name}/{counter}.png', dpi=dpi, transparent=True)
            ax.cla()
            set_butterfly_axis(ax, v)
            counter += 1
            saved = True

    if not do_not_save:
        if save_every == -1:
            fig.savefig(f'output/{name}/{name}.png', dpi=dpi)
        elif not saved:
            fig.savefig(f'output/{name}/{counter}.png', dpi=dpi, transparent=True)


def plot_spectral_radius(alphas, *, markersize=2):
    fig, ax = plt.subplots(2, figsize=(10,10))
    counter = 0
    max_radius = []
    min_radius = []
    for alpha in tqdm.tqdm(alphas):
        mat = get_matrix_for_alpha(alpha.numerator, alpha.denominator)
        _min, _max = get_spectrum_min_max(mat)
        max_radius.append(_max)
        min_radius.append(_min)
    ax[0].plot(alphas, max_radius, 'o-', markersize=markersize)
    ax[0].set_title('max radius')
    ax[1].plot(alphas, min_radius, 'o-', markersize=markersize)
    ax[1].set_title('min radius')


def farey_sequence(n: int, descending: bool = False, print_ends: bool = False):
    """Generator for the n'th Farey sequence. with or without 0 and 1, Allow for either ascending or descending."""
    (a, b, c, d) = (0, 1, 1, n)
    c_end = n if print_ends else n-1
    if descending:
        (a, c) = (1, n - 1)
    if print_ends:
        yield fractions.Fraction(a, b)
    while (c <= c_end and not descending) or (a > 0 and descending):
        k = (n + b) // d
        (a, b, c, d) = (c, d, k * c - a, k * d - b)
        if a != 0 and a != b:
            yield fractions.Fraction(a, b)


def contfrac_to_fraction(coefficients):
    numerator_2_ago = 0
    numerator_1_ago = 1
    denominator_2_ago = 1
    denominator_1_ago = 0
    for coefficient in coefficients:
        numerator = coefficient * numerator_1_ago + numerator_2_ago
        numerator_2_ago = numerator_1_ago
        numerator_1_ago = numerator
        denominator = coefficient * denominator_1_ago + denominator_2_ago
        denominator_2_ago = denominator_1_ago
        denominator_1_ago = denominator
    return fractions.Fraction(numerator, denominator)
    

def merge(count, prefix='', per=30):
    cmd_format = "convert {{{}}}.png -set page '+%[fx:0]+%[fx:0]' -background none -layers merge +repage {}.png"

    for i, start in enumerate(range(0, count, per)):
        n = min(count - start, per)
        cmd = cmd_format.format(",".join([prefix+str(v+1) for v in np.arange(start, start+n)]), f'merged_{i}')
        print(cmd)


def format_time(dt):
    precision = 3
    units = {"nsec": 1e-9, "usec": 1e-6, "msec": 1e-3, "sec": 1.0}

    scales = [(scale, unit) for unit, scale in units.items()]
    scales.sort(reverse=True)
    for scale, unit in scales:
        if dt >= scale:
            break

    return "%.*g %s" % (precision, dt / scale, unit)


def main():
    parser = argparse.ArgumentParser("Butterfly", description="Plot the Kohomoto butterfly and save as a PNG")
    parser.add_argument('n', type=int, help='Use the nth farey sequence as the alphas')
    parser.add_argument('-l', '--linewidth', type=float, help='line width for each spectrum line', default=0.1)
    parser.add_argument('--driver', help='eigenvalue driver (Read more in scipy.linalg.eigvalsh docs)', default='evr')
    parser.add_argument('-d', '--dpi', type=int, help='output image dpi', default=400)
    parser.add_argument('-s', '--save-every', type=int, help='save every x alphas, -1 for a single output', default=-1)
    parser.add_argument('-v', type=float, help='Butterfly potential amplitude', default=1)
    parser.add_argument('-c', '--color', choices=['up_down', 'a_b_type'], help='Color the butterfly according to type A/B intervals')
    parser.add_argument('--merge', action='store_true', help='Prints the command to merge the split plot using '
                                                             'ImageMagick CLI')
    parser.add_argument('--merge-prefix', default='', help="Add a prefix to the png images")
    parser.add_argument('--merge-per', type=int, default=30, help="ImageMagick takes lots of RAM so we do multiple "
                                                                  "merges")
    args = parser.parse_args()
    sequence = list(farey_sequence(args.n))

    if args.merge:
        num_files = len(sequence)//args.save_every
        if len(sequence)%args.save_every > 0:
            num_files += 1
        merge(num_files, args.merge_prefix, per=args.merge_per)
        return

    start = timer()
    plot_butterfly(sequence, f"farey_{args.n}", linewidth=args.linewidth, dpi=args.dpi, driver=args.driver,
                   save_every=args.save_every, v=args.v, coloring_scheme=args.color)
    print(f'Plotting took {format_time(timer()-start)}')


if __name__ == '__main__':
    main()
