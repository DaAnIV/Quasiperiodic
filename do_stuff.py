#!/usr/bin/env python3
from Butterfly import *

# call plot_butterfly(alphas, name, *, linewidth=0.1, dpi=400, v=1, driver='evr', save_every=500, coloring_scheme=None)
# coloring_scheme = "up_down"/"a_b_type"


def main():
    n = 50

    # Plot farey sequnece butterfly
    sequence = list(farey_sequence(n))
    plot_butterfly(sequence, "test", save_every=-1)



if __name__ == '__main__':
    main()
