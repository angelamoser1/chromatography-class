# -*- coding: utf-8 -*-
"""
This script contains the functions used for plotting
"""
import matplotlib.pyplot as plt


def set_x_axis(s, ms): 
    x_exp = []
    if ms['xaxis_plot_set'].lower() == 'time':
        x_sim = s['time_sec'] 
        label = 'Time [s]'  
    elif ms['xaxis_plot_set'].lower() == 'volume':
        x_sim = s['volume']
        label = 'Volume [mL]'   
    else:
        x_sim = s['CV']
        label = 'Volume [CV]'
        
    return x_sim, x_exp, label


def set_plot_format(fig, ax, ax2, x_label):
    # label axes
    ax.set_xlabel(x_label)
    ax.set_ylabel('Concentration [g/L]')
    ax2.set_ylabel('salt [mM]')
    
    # make it pretty
    fig.tight_layout()
    
    ax.tick_params(direction="in")
    ax2.tick_params(direction="in")
    
    ax.tick_params(axis='x', colors='black')
    ax.tick_params(axis='y', colors='black')
    ax2.tick_params(axis='y', colors='red')
    
    ax2.yaxis.label.set_color('red')
    
    return fig, ax, ax2


def generate_plot(s, ms):
    plt.rcParams.update({'font.size': 14})
    fig, ax = plt.subplots()
    x_sim, x_exp, x_label = set_x_axis(s, ms)
    cmap = plt.colormaps['viridis']
    for i, comp in enumerate(s['prot_c_mg']):
        color = cmap(i / len(s['prot_c_mg']))
        y_sim = comp
        # plot sim data for each component 
        ax.plot(x_sim, y_sim, color=color)        
            
    # plot salt
    ax2 = ax.twinx()
    ax2.plot(x_sim, s['salt'], color='red')
    
    # plot fraction cutoff line
    ax.plot((0, x_sim[-1]), (ms['cutoff'], ms['cutoff']), color='black', linestyle='dashed')

    fig, ax, ax2 = set_plot_format(fig, ax, ax2, x_label)
    
    return fig, ax, ax2