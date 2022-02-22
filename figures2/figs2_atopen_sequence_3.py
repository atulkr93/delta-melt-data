import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import pandas as pd
import math
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

mpl.rcParams['axes.linewidth'] = 4.0
fs=36
#fs=12
R = 8.314/4184
Ct = 3 * math.pow(10, -6)

def remove_spaces(line):
    line_temp = []
    for ele in line:
        if ele != "":
            line_temp.append(ele)
    return line_temp

def read_data(filename):
    #print filename
    # Read the UV data file
    f = open(filename, "r")
    temperature = []
    absorbance = []
    #line = f.readline()
    #line = f.readline()
    # Now read the rest of the lines
    for line in f.readlines():
        if line[0] == 'D' or len(line) <= 2:
            break
        line = line.strip("\n")
        line = line.split("\t")
        line = remove_spaces(line)
        line[0] = float(line[0])
        line[1] = float(line[1])
        temperature.append(line[0])
        absorbance.append(line[1]) 
    return [np.array(temperature), np.array(absorbance)]

def uv_hairpin(t, mds, bds, mss, bss, delH, Tm):
    # Gas constant in kcal/mol
    R = 8.314 / 4184
    # Get the exponential factor
    exp_factor = np.exp(((1/Tm) - (1/(t+273.16))) * delH / R)
    # Get the fraction of hairpin that is intact
    f = exp_factor / (1 + exp_factor)
    # Use f to get absorbance of hairpin/unfolded species
    absorbance = (((mds*t)+bds)*f) + (((mss*t)+bss)*(1-f))
    return absorbance 

def uv_duplex(t, mds, bds, mss, bss, delH, Tm):
    # Gas constant in kcal/mol
    R = 8.314 / 4184
    # Get the exponential factor
    exp_factor = np.exp(((1/Tm) - (1/(t+273.16))) * delH / R)
    # Get the fraction of double stranded species for a duplex
    f = ((1 + (4*exp_factor)) - np.sqrt(1 + (8*exp_factor)))/(4*exp_factor)
    # Use f to get absorbance of duplex/ss mixture
    absorbance = (((mds*t)+bds)*f) + (((mss*t)+bss)*(1-f))
    return absorbance

def plot_uv(ax, dirname, filename, yaxislabel, ylim, yinterval, xlim, xinterval, xaxisticks_status, fit_type, mode):
    ''' plot a 2state data + fit curve '''
    bf = pd.read_csv(dirname + filename)
    colors_plot = ['k', 'b', 'r']
    counter = 0

    for dummy in range(len(bf['filename'])-2):
        t_dummy = np.linspace(xlim[dummy][0], xlim[dummy][1], 500)   
        # Read the data
        dataname = bf['filename'].iloc[dummy]
        [t, a] = read_data(dirname + dataname)
        # If the fit is 2-state
        if fit_type == '2state':
            #def uv_hairpin(t, mds, bds, mss, bss, delH, Tm):
            if mode == 'duplex':
                sim_a = uv_duplex(t, bf['mds'].iloc[dummy], bf['bds'].iloc[dummy], bf['mss'].iloc[dummy], bf['bss'].iloc[dummy], bf['delH'].iloc[dummy], bf['Tm(K)'].iloc[dummy])
            elif mode == 'hairpin':
                sim_a = uv_hairpin(t, bf['mds'].iloc[dummy], bf['bds'].iloc[dummy], bf['mss'].iloc[dummy], bf['bss'].iloc[dummy], bf['delH'].iloc[dummy], bf['Tm(K)'].iloc[dummy])
            else:
                print "ERROR: Sample is neither duplex not hairpin" 
                sys.exit(0)

            # compute chi2 for comparison
            error_absorbance = np.std(a[:20])
            chi2 = np.sum(np.square((a - sim_a)/error_absorbance))
            redchi2 = chi2 / (len(a)-6.0)
            print redchi2, bf['redchi2'].iloc[dummy]

            if mode == 'duplex':
                sim_a = uv_duplex(t_dummy, bf['mds'].iloc[dummy], bf['bds'].iloc[dummy], bf['mss'].iloc[dummy], bf['bss'].iloc[dummy], bf['delH'].iloc[dummy], bf['Tm(K)'].iloc[dummy])
            elif mode == 'hairpin':
                sim_a = uv_hairpin(t_dummy, bf['mds'].iloc[dummy], bf['bds'].iloc[dummy], bf['mss'].iloc[dummy], bf['bss'].iloc[dummy], bf['delH'].iloc[dummy], bf['Tm(K)'].iloc[dummy])
            else:
                print "ERROR: Sample is neither duplex not hairpin" 
                sys.exit(0)

            ax[counter].plot(t_dummy, sim_a, color='b', linewidth=4.0) 

        # data points
        ax[counter].plot(t, a, color='k', linewidth=0.0, marker='o', markersize=8) 
        #if counter == 0:
        title = ax[counter].set_title(yaxislabel, fontsize=fs)
        title.set_position([0.5, 1.02])

        ax[counter].tick_params(axis='x', direction='out', width=4, length=8, top=False)
        ax[counter].tick_params(axis='y', direction='out', width=4, length=8, right=False)
        #ax[counter].set_title('r' + '$\chi^{2}$' + ' ' + '%.2f'%redchi2, fontsize=fs)
        ax[counter].set_ylim(ylim[dummy])
        ax[counter].set_yticks(np.arange(ylim[dummy][0], ylim[dummy][1], yinterval))
        ax[counter].set_yticklabels(['%3.2f'%ele for ele in np.arange(ylim[dummy][0], ylim[dummy][1], yinterval)], fontsize=fs)
        #if counter == 0:
        #else:
        #    ax[counter].set_yticklabels([])
        ax[counter].set_xlim(xlim[dummy])
        ax[counter].set_xticks(np.arange(xlim[dummy][0]+xinterval, xlim[dummy][1]+1, xinterval))
        if xaxisticks_status == 1:
            ax[counter].set_xticklabels(['%d'%ele for ele in np.arange(xlim[dummy][0]+xinterval, xlim[dummy][1]+1, xinterval)], fontsize=fs)
        else:
            ax[counter].set_xticklabels([])
 
        counter = counter + 1
    print
   
fig, axarr = plt.subplots(6, 6, figsize=(45, 42))
#fig, axarr = plt.subplots(2, 6)
fig.suptitle('pH 8.8HS', fontsize=fs)
plot_uv(list(axarr[0, 2:]) + list(axarr[1, :2]), '../figure9/cary_uv/scaf2_gat_m3t_done/final_all/', 'test.csv', 'scaf2_GAT' + '$^{Am3T}$', [[0.57, 0.74], [0.51, 0.70], [0.51, 0.68], [0.57, 0.74], [0.68, 0.87], [0.49, 0.68]], 0.05, [[2., 100.], [2., 100.], [2., 100.], [2., 100.], [2., 100.], [2., 100.]], 20.0, 1, '2state', 'duplex')
plot_uv(axarr[1, 2:5], '../figure9/cary_uv/scaf2_taa_at_done/final_all/', 'test.csv', 'scaf2_TAA' + '$^{AT}$', [[0.54, 0.74], [0.50, 0.70], [0.63, 0.83]], 0.05, [[2., 100.], [2., 100.], [2., 100.]], 20.0, 1, '2state', 'duplex')
plot_uv([axarr[1, 5]] + list(axarr[2, :2]), '../figure9/cary_uv/scaf2_taa_m3t_done/final_all/', 'test.csv', 'scaf2_TAA' + '$^{Am3T}$', [[0.61, 0.78], [0.61, 0.78], [0.61, 0.78]], 0.05, [[2., 100.], [2., 100.], [2., 100.]], 20.0, 1, '2state', 'duplex')
plot_uv(list(axarr[2, 2:5]), '../figure9/cary_uv/scaf2_tac_at_done/final_all/', 'test.csv', 'scaf2_TAC' + '$^{AT}$', [[0.62, 0.82], [0.53, 0.71], [0.62, 0.81]], 0.05, [[2., 100.], [2., 100.], [2., 100.]], 20.0, 1, '2state', 'duplex')
plot_uv([axarr[2, 5]] + list(axarr[3, :3]), '../figure9/cary_uv/scaf2_tac_m3t_done/final_all/', 'test.csv', 'scaf2_TAC' + '$^{Am3T}$', [[0.62, 0.81], [0.62, 0.81], [0.65, 0.82], [0.62, 0.81]], 0.05, [[2., 100.], [2., 100.], [2., 100.], [2., 100.]], 20.0, 1, '2state', 'duplex')
plot_uv(list(axarr[3, 3:]) + [axarr[4, 0]], '../figure9/cary_uv/scaf2_tag_m3t_done/final_all/', 'test.csv', 'scaf2_TAG' + '$^{Am3T}$', [[0.58, 0.76], [0.63, 0.80], [0.58, 0.76], [0.58, 0.76]], 0.05, [[2., 100.], [2., 100.], [2., 100.], [2., 100.]], 20.0, 1, '2state', 'duplex')
plot_uv(axarr[4, 1:4], '../figure9/cary_uv/scaf2_tat_at_done/final_all/', 'test.csv', 'scaf2_TAT' + '$^{AT}$', [[0.62, 0.82], [0.62, 0.82], [0.62, 0.82]], 0.05, [[2., 100.], [2., 100.], [2., 100.]], 20.0, 1, '2state', 'duplex')
plot_uv(list(axarr[4, 4:]) + list(axarr[5, :3]), '../figure9/cary_uv/scaf2_tat_m3t_done/final_all/', 'test.csv', 'scaf2_TAT' + '$^{Am3T}$', [[0.66, 0.84], [0.63, 0.80], [0.65, 0.78], [0.71, 0.88], [0.66, 0.85]], 0.05, [[2., 100.], [2., 100.], [2., 100.], [2., 100.], [2., 100.]], 20.0, 1, '2state', 'duplex')



for dummy in range(6):
    axarr[dummy, 0].set_ylabel('A' + '$_{260}$', fontsize=fs)

for dummy in range(6):
    axarr[0, dummy].remove()
    if dummy <= 2:
        axarr[5, dummy].set_xlabel('Temperature (' + '$^{\circ}$' + 'C)', fontsize=fs)
    else:
        axarr[5, dummy].remove()
        axarr[4, dummy].set_xlabel('Temperature (' + '$^{\circ}$' + 'C)', fontsize=fs)

    #axarr[1, dummy].remove()

#plt.show()
plt.tight_layout()
plt.savefig('figs2_atopen_sequence_3.pdf')

