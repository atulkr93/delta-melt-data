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
    line = f.readline()
    line = line.strip('\n')
    #print line
    line = f.readline()
    line = line.strip('\n')
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

def plot_uv(ax, dirname, filename, yaxislabel, ylim, yinterval, xlim, xinterval, xaxisticks_status, fit_type, mode):
    ''' plot a 2state data + fit curve '''
    bf = pd.read_csv(dirname + filename)
    colors_plot = ['k', 'b', 'r']
    t_dummy = np.linspace(0.0, 105.0, 500)   
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
   
fig, axarr = plt.subplots(4, 6, figsize=(44, 28.3))
fig.suptitle('pH 5.4HS', fontsize=fs)
plot_uv(axarr[0, :3], '../figure8/a6dna_ph5.4nbhs/', 'test.csv', 'A' + '$_{6}$' + '-DNA', [[0.56, 0.78], [0.48, 0.68], [0.48, 0.68]], 0.05, [[10, 90.0], [10., 90.], [10., 90.]], 15.0, 1, '2state', 'duplex')
plot_uv(axarr[0, 3:], '../figure8/a6dna_m1a16_ph5.4hsnb/', 'test.csv', 'A' + '$_{6}$' + '-DNA' + '$^{m1A16}$', [[0.57, 0.78], [0.57, 0.78], [0.55, 0.76]], 0.05, [[10, 90.0], [10., 90.], [10., 90.]], 15.0, 1, '2state', 'duplex')
plot_uv(axarr[1, :3], '../figure8/a6dna_m1g10_keck_ph5.4nbhs/', 'test.csv', 'A' + '$_{6}$' + '-DNA' + '$^{m1G10}$', [[0.55, 0.75], [0.56, 0.76], [0.55, 0.75]], 0.05, [[10, 90.0], [10., 90.], [10., 90.]], 15.0, 1, '2state', 'duplex')
plot_uv(axarr[1, 3:], '../figure8/a6dna_m1a16m1g10_ph5.4hsnb/', 'test.csv', 'A' + '$_{6}$' + '-DNA' + '$^{m1G10,m1A16}$', [[0.55, 0.76], [0.52, 0.72], [0.52, 0.75]], 0.05, [[10, 90.0], [10., 90.], [10., 90.], [10., 90.]], 15.0, 1, '2state', 'duplex')

plot_uv(axarr[2, :3], '../figure8/a2dna_ph5.4hsnb/', 'test.csv', 'A' + "$_{2}$" + "-DNA", [[0.59, 0.77], [0.59, 0.76], [0.60, 0.80]], 0.05, [[0, 95.0], [0., 95.], [0., 95.]], 20.0, 1, '2state', 'duplex')
plot_uv(axarr[2, 3:], '../figure8/a2dna_m1g6_ph5.4hsnb/', 'test.csv', 'A' + "$_{2}$" + "-DNA" + "$^{m1G6}$", [[0.80, 0.97], [0.57, 0.74], [0.58, 0.74]], 0.05, [[0, 95.0], [0., 95.], [0., 95.]], 20.0, 1, '2state', 'duplex')
plot_uv(axarr[3, :3], '../figure8/a2dna_m1g20_ph5.4hsnb/', 'test.csv', 'A' + "$_{2}$" + "-DNA" + "$^{m1G20}$", [[0.62, 0.82], [0.62, 0.82], [0.62, 0.82]], 0.05, [[0, 95.0], [0., 95.], [0., 95.]], 20.0, 1, '2state', 'duplex')
plot_uv(axarr[3, 3:], '../figure8/a2dna_m1g6m1g20_ph5.4hsnb/', 'test.csv', 'A' + "$_{2}$" + "-DNA" + "$^{m1G6,20}$", [[0.54, 0.72], [0.54, 0.72], [0.54, 0.72]], 0.05, [[0, 95.0], [0., 95.], [0., 95.]], 20.0, 1, '2state', 'duplex')

for dummy in range(2):
    axarr[dummy, 0].set_ylabel('A' + '$_{260}$', fontsize=fs)

for dummy in range(6):
    axarr[3, dummy].set_xlabel('Temperature (' + '$^{\circ}$' + 'C)', fontsize=fs)
#plt.show()
plt.tight_layout()
plt.savefig('figs1_hoogsteen_coop_ph5p4hs.pdf')
