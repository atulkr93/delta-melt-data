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
    #line = line.strip('\n')
    #print line

    #line = f.readline()
    #line = line.strip('\n')
    #print line
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

def plot_uv(ax, dirname, filename, yaxislabel, ylim, yinterval, xlim, xinterval, xaxisticks_status, fit_type, mode, xlims, ylims, xtickins, ytickins):
    ''' plot a 2state data + fit curve '''
    bf = pd.read_csv(dirname + filename)
    colors_plot = ['k', 'b', 'r']
    t_dummy = np.linspace(0.0, 120.0, 500) 
    counter = 0

    for dummy in range(len(bf['filename'])-2):
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

            # inset
            #axins = inset_axes(ax[counter], width="82%", height="82%", bbox_to_anchor=[0.44, -0.01, 0.55, 0.55], bbox_transform=ax[counter].transAxes, borderpad=0.2)
            axins = inset_axes(ax[counter], width="82%", height="82%", bbox_to_anchor=[0.20, 0.49, 0.5, 0.5], bbox_transform=ax[counter].transAxes, borderpad=0.2)
            axins.plot(t_dummy, sim_a, color='b', linewidth=4.0)
            axins.plot(t, a, color='k', linewidth=0.0, marker='o', markersize=8) 
            axins.set_xlim([xlims[dummy][0], xlims[dummy][1]])
            axins.set_ylim([ylims[dummy][0], ylims[dummy][1]])
            axins.tick_params(axis='x', direction='out', width=4, length=4, pad=0.0, top=False)
            axins.tick_params(axis='y', direction='out', width=4, length=4, pad=0.0, right=False)
            axins.set_yticks(np.arange(ytickins[dummy][0], ytickins[dummy][1], ytickins[dummy][2]))
            axins.set_yticklabels(['%3.3f'%ele for ele in np.arange(ytickins[dummy][0], ytickins[dummy][1], ytickins[dummy][2])], fontsize=fs)
            axins.set_xticks([int(ele) for ele in np.arange(xtickins[dummy][0], xtickins[dummy][1], xtickins[dummy][2])])
            axins.set_xticklabels(['%d'%ele for ele in np.arange(xtickins[dummy][0], xtickins[dummy][1], xtickins[dummy][2])], fontsize=fs)

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
        ax[counter].set_xticks(np.arange(xlim[dummy][0], xlim[dummy][1]+1, xinterval))
        if xaxisticks_status == 1:
            ax[counter].set_xticklabels(['%d'%ele for ele in np.arange(xlim[dummy][0], xlim[dummy][1]+1, xinterval)], fontsize=fs)
        else:
            ax[counter].set_xticklabels([])
 
        counter = counter + 1
    print
   
fig, axarr = plt.subplots(5, 6, figsize=(46, 35))
#fig, axarr = plt.subplots(2, 6)
fig.suptitle('pH 8.8HS', fontsize=fs)
plot_uv(axarr[0, :3], '../figure9/cary_uv/scaf2_aaa_m3t_done/final_all/', 'test.csv', "scaf2_AAA" + "$^{Am3T}$", [[0.67, 0.91], [0.66, 0.945], [0.66, 0.95]], 0.05, [[5., 105.0], [5., 105.], [5., 105.]], 20.0, 1, '2state', 'duplex', [[8., 22.], [7., 22.], [6., 22.]], [[0.6835, 0.692], [0.672, 0.685], [0.674, 0.685]], [[10., 21., 10.], [10., 21., 10.], [10., 21., 10.]], [[0.685, 0.692, 0.005], [0.675, 0.685, 0.007], [0.675, 0.685, 0.007]])
plot_uv(axarr[0, 3:], '../figure9/cary_uv/scaf2_aag_at_done/final_all/', 'test.csv', "scaf2_AAG" + "$^{AT}$", [[0.60, 0.84], [0.60, 0.84], [0.60, 0.84]], 0.05, [[0, 105.0], [0., 105.], [0., 105.]], 20.0, 1, '2state', 'duplex', [[3., 40.], [3., 40.], [3., 40.]], [[0.615, 0.645], [0.617, 0.645], [0.616, 0.645]], [[10., 41., 20.], [10., 41., 20.], [10., 41., 20.]], [[0.62, 0.65, 0.02], [0.62, 0.65, 0.02], [0.62, 0.65, 0.02]])
plot_uv(axarr[1, :3], '../figure9/cary_uv/scaf2_aag_m1a_done/final_all/', 'test.csv', "scaf2_AAG" + "$^{m1AT}$", [[0.59, 0.83], [0.59, 0.83], [0.59, 0.84]], 0.05, [[0, 105.0], [0., 105.], [0., 105.]], 20.0, 1, '2state', 'duplex', [[70., 100.], [70., 100.], [70., 100.]], [[0.72, 0.775], [0.70, 0.76], [0.71, 0.77]], [[70., 101., 20.], [70., 101., 20.], [70., 101., 20.]], [[0.72, 0.79, 0.05], [0.70, 0.78, 0.05], [0.71, 0.77, 0.05]])
plot_uv(axarr[1, 3:], '../figure9/cary_uv/scaf2_aag_m3t_done/final_all/', 'test.csv', "scaf2_AAG" + "$^{Am3T}$", [[0.61, 0.79], [0.60, 0.81], [0.62, 0.82]], 0.05, [[3., 105.0], [3., 105.], [3., 105.]], 20.0, 1, '2state', 'duplex', [[5., 25.], [6., 25.], [5., 25.]], [[0.615, 0.628], [0.615, 0.628], [0.627, 0.640]], [[10., 21., 10.], [10., 21., 10.], [10., 21., 10.]], [[0.615, 0.628, 0.007], [0.615, 0.628, 0.008], [0.63, 0.641, 0.007]])
plot_uv(axarr[2, :3], '../figure9/cary_uv/scaf2_cac_at_done/final_all/', 'test.csv', "scaf2_CAC" + "$^{AT}$", [[0.60, 0.89], [0.63, 0.91], [0.61, 0.9]], 0.05, [[2., 105.0], [2., 105.], [2., 105.]], 20.0, 1, '2state', 'duplex', [[2., 35.], [4., 35.], [3., 35.]], [[0.615, 0.630], [0.640, 0.657], [0.620, 0.638]], [[10., 31., 10.], [10., 31., 10.], [10., 31., 10.]], [[0.615, 0.630, 0.008], [0.640, 0.660, 0.015], [0.620, 0.640, 0.015]])
plot_uv(axarr[2, 3:], '../figure9/cary_uv/scaf2_gac_at_done/final_all/', 'test.csv', "scaf2_GAC" + "$^{AT}$", [[0.62, 0.86], [0.62, 0.86], [0.62, 0.86]], 0.05, [[5., 105.0], [5., 105.], [5., 105.]], 20.0, 1, '2state', 'duplex', [[8., 40.], [5., 40.], [5., 40.]], [[0.631, 0.645], [0.627, 0.643], [0.627, 0.644]], [[10., 31., 10.], [10., 31., 10.], [10., 31., 10.]], [[0.632, 0.645, 0.01], [0.628, 0.644, 0.01], [0.627, 0.645, 0.01]])
plot_uv(axarr[3, :3], '../figure9/cary_uv/scaf2_gag_m3t_done/final_all/', 'test.csv', "scaf2_GAG" + "$^{Am3T}$", [[0.67, 0.94], [0.60, 0.84], [0.66, 0.92]], 0.05, [[5., 105.0], [5., 105.], [5., 105.]], 20.0, 1, '2state', 'duplex', [[20., 45.], [20., 45.], [20., 45.]], [[0.695, 0.74], [0.61, 0.65], [0.685, 0.725]], [[20., 51., 20.], [20., 51., 20.], [20., 51., 20.]], [[0.70, 0.74, 0.03], [0.61, 0.65, 0.03], [0.69, 0.72, 0.02]])
plot_uv(list(axarr[3, 3:]) + [axarr[4, 0]], '../figure9/cary_uv/scaf2_tag_at_done/final_all/', 'test.csv', "scaf2_TAG" + "$^{AT}$", [[0.61, 0.88], [0.61, 0.89], [0.61, 0.87], [0.60, 0.88]], 0.05, [[2., 105.0], [2., 105.], [2., 105.], [2., 105.]], 20.0, 1, '2state', 'duplex', [[4., 35.], [4., 35.], [4., 35.], [4., 35.]], [[0.622, 0.640], [0.623, 0.643], [0.617, 0.638], [0.612, 0.634]], [[10., 35., 20.], [10., 35., 20.], [10., 35., 20.], [10., 35., 20.]], [[0.625, 0.640, 0.01], [0.625, 0.643, 0.015], [0.620, 0.636, 0.015], [0.615, 0.634, 0.015]])
plot_uv(axarr[4, 1:4], '../figure9/cary_uv/scaf2_tag_m1a_done/final_all/', 'test.csv', "scaf2_TAG" + "$^{m1AT}$", [[0.64, 0.91], [0.66, 0.93], [0.74, 1.02]], 0.05, [[2., 105.0], [2., 105.], [2., 105.]], 20.0, 1, '2state', 'duplex', [[4., 35.], [4., 35.], [4., 35.]], [[0.652, 0.678], [0.666, 0.700], [0.752, 0.777]], [[10., 35., 20.], [10., 35., 20.], [10., 35., 20.]], [[0.655, 0.678, 0.02], [0.67, 0.70, 0.02], [0.755, 0.775, 0.020]])





for dummy in range(2):
    axarr[dummy, 0].set_ylabel('A' + '$_{260}$', fontsize=fs)

for dummy in range(6):
    if dummy <= 3:
        axarr[4, dummy].set_xlabel('Temperature (' + '$^\circ$' + 'C)', fontsize=fs)
    else:
        axarr[3, dummy].set_xlabel('Temperature (' + '$^\circ$' + 'C)', fontsize=fs)
        axarr[4, dummy].remove()
#plt.show()
plt.tight_layout()
plt.savefig('figs3_b_2.pdf')
