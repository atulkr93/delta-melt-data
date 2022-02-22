#!/usr/bin/python
import matplotlib as mpl
import pandas as pd
import numpy as np
import sys
from nmrglue import proc_base
from nmrglue import analysis
import math
import matplotlib.pyplot as plt
import sys
lf = 70.94
fs = 32
mpl.rcParams['axes.linewidth'] = 2.0

def plot_profiles(filename, op_filename):
    ''' Given the parameters, intensity from simulation, plot the CEST profile '''
    bf = pd.read_csv(filename)
    number_slps = np.unique(bf['slp(hz)'])
    colors_plot = ['b', 'cyan', 'b', 'g', 'cyan', 'magenta', 'brown', 'yellow', 'teal', 'lightgreen']
    # Plot profile

    fig, ax = plt.subplots(1)
    counter = 0
    for dummy_slp in number_slps:
        bf_new = bf.loc[bf['slp(hz)'] == dummy_slp]
        plt.errorbar(bf_new['offset(hz)'], bf_new['norm_intensity'], yerr=bf_new['norm_intensity_error'], color = colors_plot[counter], label='%d'%math.ceil(float(dummy_slp)), linewidth=3, fmt='o')
        counter = counter + 1
    #plt.xlim([((np.amin(params['offsets'])-(2*math.pi*0.8*params['lf']))/(2*math.pi*params['lf']))+138.4, ((np.amax(params['offsets'])+(2*math.pi*0.8*params['lf']))/(2*math.pi*params['lf']))+138.4])
    #plt.xlim([((np.amin(params['offsets'])-(2*math.pi*0.8*params['lf']))/(2*math.pi*params['lf'])), ((np.amax(params['offsets'])+(2*math.pi*0.8*params['lf']))/(2*math.pi*params['lf']))])
    #plt.ylim([-0.05, 0.6])
    plt.xlabel('$\Omega (kHz)$', fontsize=fs) 
    plt.ylabel('Norm. Intensity', fontsize=fs)
    plt.xticks([-4000., -2000., 0.0, 2000., 4000.], [-4, -2, 0, 2, 4], fontsize=fs)
    plt.ylim([-0.05, 1.0])
    plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], [0.0, 0.2, 0.4, 0.6, 0.8, 1.0], fontsize=fs)
    ax.tick_params(axis='x', direction='out', width=2, top=False) 
    ax.tick_params(axis='y', direction='out', width=2, right=False) 
    
    # Legend
    handles, labels = ax.get_legend_handles_labels()
    handles = [h[0] if type(h) is not mpl.lines.Line2D else h for h in handles]
    leg = ax.legend(handles, labels, title='$\omega_{1}/2\pi$' + ' ' + '$(Hz)$', fancybox=False, ncol=2, handlelength=0, frameon=False, columnspacing=0.0, markerscale=0.0, handletextpad=0.5, borderpad=0.0, handleheight=0.0, labelspacing=0.2, loc=1, fontsize=fs)
    counter = 0
    for text in leg.get_texts():
        text.set_color(colors_plot[counter])
        counter = counter + 1
    for l in leg.get_lines():
        l.set_linestyle('None')
    leg.set_title('$\omega_{1}/2\pi$' + ' ' + '$(Hz)$', prop={'size':fs})
    plt.tight_layout()
    plt.show()

    #plt.savefig(op_filename + "-intensity.pdf")

    fig, ax = plt.subplots(1)
    counter = 0
    for dummy_slp in number_slps:
        bf_new = bf.loc[bf['slp(hz)'] == dummy_slp]
        plt.errorbar(bf_new['offset(hz)']/lf, bf_new['norm_volume'], yerr=bf_new['norm_volume_error'], color = colors_plot[counter], label='%d'%float(dummy_slp), linewidth=3, fmt='o')
        counter = counter + 1
    #plt.xlim([((np.amin(params['offsets'])-(2*math.pi*0.8*params['lf']))/(2*math.pi*params['lf']))+138.4, ((np.amax(params['offsets'])+(2*math.pi*0.8*params['lf']))/(2*math.pi*params['lf']))+138.4])
    #plt.xlim([((np.amin(params['offsets'])-(2*math.pi*0.8*params['lf']))/(2*math.pi*params['lf'])), ((np.amax(params['offsets'])+(2*math.pi*0.8*params['lf']))/(2*math.pi*params['lf']))])
    #plt.ylim([-0.05, 0.6])
    plt.xlabel('$\Omega/2\pi$' + '(kHz)', fontsize=fs) 
    plt.ylabel('Norm. Volume', fontsize=fs)
    plt.xticks([-12, -8, -4, 0, 4, 8, 12])
    ax.tick_params(axis='x', direction='out', width=2, top=False) 
    ax.tick_params(axis='y', direction='out', width=2, right=False) 
 
    # Legend
    handles, labels = ax.get_legend_handles_labels()
    handles = [h[0] if type(h) is not mpl.lines.Line2D else h for h in handles]
    leg = ax.legend(handles, labels, title='$\omega_{1}/2\pi$' + ' ' + '$(Hz)$', fancybox=False, ncol=2, handlelength=0, frameon=False, columnspacing=0.0, markerscale=0.0, handletextpad=0.5, borderpad=0.0, handleheight=0.0, labelspacing=0.2, loc=1, fontsize=fs)
    counter = 0
    for text in leg.get_texts():
        text.set_color(colors_plot[counter])
        counter = counter + 1
    for l in leg.get_lines():
        l.set_linestyle('None')
    leg.set_title('$\omega_{1}/2\pi$' + ' ' + '$(Hz)$', prop={'size':fs})

    #plt.ylim([0.0, 0.35])
    #plt.savefig(op_filename + "-volume.pdf")
    plt.tight_layout()
    plt.show()
 
filename = str(sys.argv[1])
op_filename = filename[:filename.index('.')]
plot_profiles(filename, op_filename) 
