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
import argparse
import sortednp
import fractions
import numpy as np
import scipy.sparse
import scipy.linalg
import scipy.special
import scipy.optimize
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


def _set_butterfly_axis(ax, v):
    ax.set_yticks(np.linspace(0, 1, 21))
    if v > 0:
        ax.set_xlim([-2.05, 2.05+v])
    else:
        ax.set_xlim([-2.05+v, 2.05])

    ax.set_ylim([-0.05, 1.05])
    ax.set_ylabel('Frequencies \u03B1')
    ax.set_xlabel('Spectrum')


def plot_butterfly(alphas, name, *, linewidth=0.1, dpi=400, v=1, driver='evr', save_every=500):
    if not os.path.exists(f'output/{name}'):
        os.mkdir(f'output/{name}')
    fig, ax = plt.subplots(figsize=(10, 10))
    _set_butterfly_axis(ax, v)

    ax.plot([-2, 2], [0, 0], 'k-', linewidth=linewidth)
    ax.plot([-2 + v, 2 + v], [1, 1], 'k-', linewidth=linewidth)

    if save_every != -1:
        fig.savefig(f'output/{name}/background.png', dpi=dpi)

        ax.cla()
        _set_butterfly_axis(ax, v)

    saved = False
    counter = 1
    for i, alpha in enumerate(tqdm.tqdm(alphas)):
        saved = False
        mat = get_matrix_for_alpha(alpha.numerator, alpha.denominator, v)
        for x1, x2 in get_intervals_from_matrix(mat, driver=driver):
            ax.plot([x1, x2], [alpha, alpha], 'k-', linewidth=linewidth)
        if save_every != -1 and i > 0 and i % save_every == 0:
            fig.savefig(f'output/{name}/{counter}.png', dpi=dpi, transparent=True)
            ax.cla()
            _set_butterfly_axis(ax, v)
            counter += 1
            saved = True

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


def merge(count, prefix='', per=30):
    cmd_format = "convert {{{}}}.png -set page '+%[fx:0]+%[fx:0]' -background none -layers merge +repage {}.png"

    for i, start in enumerate(range(0, count, per)):
        n = min(count - start, per)
        cmd = cmd_format.format(",".join([prefix+str(v+1) for v in np.arange(start, start+n)]), f'merged_{i}')
        print(cmd)


def main():
    parser = argparse.ArgumentParser("Butterfly", description="Plot the Kohomoto butterfly and save as a PNG")
    parser.add_argument('n', type=int, help='Use the nth farey sequence as the alphas')
    parser.add_argument('-l', '--linewidth', type=float, help='line width for each spectrum line', default=0.1)
    parser.add_argument('--driver', help='eigenvalue driver (Read more in scipy.linalg.eigvalsh docs)', default='evr')
    parser.add_argument('-d', '--dpi', type=int, help='output image dpi', default=400)
    parser.add_argument('-s', '--save-every', type=int, help='save every x alphas, -1 for a single output', default=-1)
    parser.add_argument('-v', type=int, help='Butterfly potential amplitude', default=1)
    parser.add_argument('--merge', action='store_true', help='Prints the command to merge the split plot using '
                                                             'ImageMagick CLI')
    parser.add_argument('--merge-prefix', default='', help="Add a prefix to the png images")
    parser.add_argument('--merge-per', type=int, default=30, help="ImageMagick takes lots of RAM so we do multiple "
                                                                  "merges")
    args = parser.parse_args()
    sequence = list(farey_sequence(args.n))

    if args.merge:
        merge(args.n, args.merge_prefix, per=args.merge_per)
        return

    start = timer()
    plot_butterfly(sequence, f"farey_{args.n}", linewidth=args.linewidth, dpi=args.dpi, driver=args.driver,
                   save_every=args.save_every, v=args.v)
    print(f'Plotting took {timer()-start} sec')


if __name__ == '__main__':
    main()
