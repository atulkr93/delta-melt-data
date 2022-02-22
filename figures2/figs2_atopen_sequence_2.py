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
   
fig, axarr = plt.subplots(9, 6, figsize=(45, 63))
#fig, axarr = plt.subplots(2, 6)
fig.suptitle('pH 8.8HS', fontsize=fs)
plot_uv(list(axarr[0, :3]), '../figure9/cary_uv/scaf2_aat_m3t_done/final_all/', 'test.csv', 'scaf2_AAT' + '$^{Am3T}$', [[-0.16, 0.02], [0.46, 0.65], [0.69, 0.88]], 0.05, [[15, 100.], [15., 100.], [15., 100.]], 20.0, 1, '2state', 'duplex')
plot_uv(axarr[0, 3:], '../figure9/cary_uv/scaf2_caa_at_done/final_all/', 'test.csv', 'scaf2_CAA' + '$^{AT}$', [[0.61, 0.82], [0.61, 0.82], [0.61, 0.82]], 0.05, [[10., 100.], [10., 100.], [10., 100.]], 20.0, 1, '2state', 'duplex')
plot_uv(axarr[1, :4], '../figure9/cary_uv/scaf2_caa_m3t_done/final_all/', 'test.csv', 'scaf2_CAA' + '$^{Am3T}$', [[0.64, 0.82], [0.63, 0.81], [0.66, 0.85], [0.63, 0.81]], 0.05, [[10., 100.], [10., 100.], [10., 100.], [10., 100.]], 20.0, 1, '2state', 'duplex')
plot_uv(list(axarr[1, 4:]) + [axarr[2, 0]], '../figure9/cary_uv/scaf2_cac_m3t_done/final_all/', 'test.csv', 'scaf2_CAC' + '$^{Am3T}$', [[0.40, 0.58], [0.50, 0.67], [0.40, 0.58]], 0.05, [[10., 100.], [10., 100.], [10., 100.]], 20.0, 1, '2state', 'duplex')
plot_uv(axarr[2, 1:4], '../figure9/cary_uv/scaf2_cag_at_done/final_all/', 'test.csv', 'scaf2_CAG' + '$^{AT}$', [[-0.14, 0.04], [0.72, 0.89], [0.72, 0.90]], 0.05, [[10., 100.], [10., 100.], [10., 100.]], 20.0, 1, '2state', 'duplex')
plot_uv(list(axarr[2, 4:]) + [axarr[3, 0]], '../figure9/cary_uv/scaf2_cag_m3t_done/final_all/', 'test.csv', 'scaf2_CAG' + '$^{Am3T}$', [[0.62, 0.79], [0.62, 0.79], [0.62, 0.79]], 0.05, [[10., 100.], [10., 100.], [10., 100.]], 20.0, 1, '2state', 'duplex')
plot_uv(axarr[3, 1:4], '../figure9/cary_uv/scaf2_cat_at_done/final_all/', 'test.csv', 'scaf2_CAT' + '$^{AT}$', [[0.42, 0.62], [0.42, 0.62], [0.42, 0.62]], 0.05, [[10., 100.], [10., 100.], [10., 100.]], 20.0, 1, '2state', 'duplex')
plot_uv(list(axarr[3, 4:]) + list(axarr[4, :2]), '../figure9/cary_uv/scaf2_cat_m3t_done/final_all/', 'test.csv', 'scaf2_CAT' + '$^{Am3T}$', [[0.62, 0.82], [0.56, 0.72], [0.59, 0.76], [0.62, 0.82]], 0.05, [[10., 100.], [10., 100.], [10., 100.], [10., 100.]], 20.0, 1, '2state', 'duplex')
plot_uv(axarr[4, 2:5], '../figure9/cary_uv/scaf2_gaa_at_done/final_all/', 'test.csv', 'scaf2_GAA' + '$^{AT}$', [[0.60, 0.81], [0.63, 0.82], [0.60, 0.81]], 0.05, [[10., 100.], [10., 100.], [10., 100.]], 20.0, 1, '2state', 'duplex')
plot_uv([axarr[4, 5]] + list(axarr[5, :2]), '../figure9/cary_uv/scaf2_gaa_m3t_done/final_all/', 'test.csv', 'scaf2_GAA' + '$^{Am3T}$', [[0.62, 0.80], [0.62, 0.80], [0.5, 0.67]], 0.05, [[10., 100.], [10., 100.], [10., 100.]], 20.0, 1, '2state', 'duplex')
plot_uv(list(axarr[5, 2:]), '../figure9/cary_uv/scaf2_gac_m3t_done/final_all/', 'test.csv', 'scaf2_GAC' + '$^{Am3T}$', [[0.65, 0.84], [0.65, 0.84], [0.67, 0.84], [0.65, 0.84]], 0.05, [[10., 100.], [10., 100.], [10., 100.], [10., 100.]], 20.0, 1, '2state', 'duplex')
plot_uv(axarr[6, :3], '../figure9/cary_uv/scaf2_gag_at_done/final_all/', 'test.csv', 'scaf2_GAG' + '$^{AT}$', [[0.68, 0.87], [0.68, 0.87], [0.70, 0.90]], 0.05, [[2., 100.], [2., 100.], [2., 100.]], 20.0, 1, '2state', 'duplex')
plot_uv(list(axarr[6, 3:]) + list(axarr[7, :2]), '../figure9/cary_uv/scaf2_gat_at_done/final_all/', 'test.csv', 'scaf2_GAT' + '$^{AT}$', [[0.46, 0.67], [0.46, 0.67], [0.47, 0.63], [0.47, 0.65], [0.60, 0.77]], 0.05, [[2., 100.], [2., 100.], [2., 100.], [2., 100.], [2., 100.]], 20.0, 1, '2state', 'duplex')
plot_uv(list(axarr[7, 2:]) + list(axarr[8, :2]), '../figure9/cary_uv/scaf2_gat_m3t_done/final_all/', 'test.csv', 'scaf2_GAT' + '$^{Am3T}$', [[0.57, 0.74], [0.51, 0.70], [0.51, 0.68], [0.57, 0.74], [0.68, 0.87], [0.49, 0.68]], 0.05, [[2., 100.], [2., 100.], [2., 100.], [2., 100.], [2., 100.], [2., 100.]], 20.0, 1, '2state', 'duplex')



for dummy in range(9):
    axarr[dummy, 0].set_ylabel('A' + '$_{260}$', fontsize=fs)

for dummy in range(6):
    axarr[8, dummy].remove()


#plt.show()
plt.tight_layout()
plt.savefig('figs2_atopen_sequence_2.pdf')
