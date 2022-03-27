#!/usr/bin/env python3
import tqdm
from Butterfly import *

# call plot_butterfly(alphas, name, *, linewidth=0.1, dpi=400, v=1, driver='evr', save_every=500, coloring_scheme=None)
# coloring_scheme = "up_down"/"a_b_type"


def set_axis(ax, x_min, x_max):
    ax.set_xlim([x_min, x_max])
    # ax.get_yaxis().set_visible(False)

    # ax.set_ylabel('Frequencies \u03B1')
    ax.set_xlabel('Spectrum')

def c_to_fraction(c):
    return contfrac_to_fraction(c[1:])

def main():
    ################### Parameters ######################
    n = 100
    v = 4.2
    name = "test"
    dpi = 1000
    butterfly_line_width = 0.1
    special_line_width = 0.3
    coloring="a_b_type" # only special_case

    should_plot_buttefly = True
    should_plot_zoom = True

    x_min = -2.05
    # x_max = 2.05+v
    x_max = 1.1
    special_alpha= [(1,     [0, 0, 3]), 
                    (1.1,   [0,0,3,2]),
                    (1.2,   [0,0,3,2,1]),
                    (1.3,   [0,0,3,2,2]),
                    (1.4,   [0,0,3,2,3]),
                    (1.5,   [0,0,3,2,3,1])
                ]

    ################### Code ######################

    for y, c in special_alpha:  
        a = c_to_fraction(c)  
        print(f'{y}: {c}={a}')

    print()

    sequence = list(farey_sequence(n))

    ################### Zoomed image ######################
    if should_plot_zoom:
        s_fig, s_ax = plt.subplots(figsize=(10, 5))
        set_axis(s_ax, x_min, x_max)

        ticks = []
        tick_labels = []

        for y, c in tqdm.tqdm(special_alpha):  
            a = c_to_fraction(c)  
            plot_alpha(a, s_ax, v=v, linewidth=0.5, coloring_scheme=coloring, overide_y=y, coeff=c[1:], verbose=True)

            ticks.append(y)
            tick_labels.append(str(c))
        
        s_ax.set_yticks(ticks)
        s_ax.set_yticklabels(tick_labels)
            

        s_fig.savefig(f'output/{name}/{name}_zoom.png', dpi=dpi)

    ################### Butterfly image ######################
    if should_plot_buttefly:
        b_fig, b_ax = plt.subplots(figsize=(10, 10))
        set_butterfly_axis(b_ax, v)
    
        # Plot farey sequnece butterfly
        plot_butterfly(sequence, name, linewidth=butterfly_line_width, save_every=-1, v=v, ax=b_ax)

        for y, c in special_alpha:    
            a = c_to_fraction(c)
            plot_alpha(a, b_ax, v=v, linewidth=special_line_width, coloring_scheme=coloring)

        b_fig.savefig(f'output/{name}/{name}_butterfly.png', dpi=dpi)




if __name__ == '__main__':
    main()
