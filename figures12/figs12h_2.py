import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from os import path
from uncertainties import unumpy
import sys
from matplotlib import rcParams
from uncertainties import ufloat

rcParams['axes.linewidth'] = 2.0
rcParams['font.family']='Sans-serif'
rcParams['font.sans-serif']=['Arial']

cpsize=5
linew=2
fs = 18

nontwo_state = ['aaa', 'aag', 'gag']
sequence_contexts = []
delG_wc_values = []
delG_m1g_values = []
delG_wc_error_values = []
delG_m1g_error_values = []
delH_wc_values = []
delH_m1g_values = []
delH_wc_error_values = []
delH_m1g_error_values = []
delS_wc_values = []
delS_m1g_values = []
delS_wc_error_values = []
delS_m1g_error_values = []

comp_nts = {'a':'t', 't':'a', 'g':'c', 'c':'g'}

for dummy_outer in ['a', 'c', 'g', 't']:
    for dummy_inner in ['a', 'c', 'g', 't']:
        # get sequence context
        sequence_context = dummy_outer + 'a' + dummy_inner

        #if sequence_context not in ["AGG", "TGG"]:
        # Get wc delG
        wc_file = '../figure9/cary_uv/scaf2_' + sequence_context + '_at_done/final_all/' + 'test.csv'
        bf_wc = pd.read_csv(wc_file)
        # Get values
        #delG_wc = bf_wc['delG_37C'].values
        delS_wc = bf_wc['delS'].values
        delH_wc = bf_wc['delH'].values
        delG_wc = delH_wc - ((273.16+25)*delS_wc)
        delG_wc[-2] = np.average(delG_wc[:-2])
        delG_wc[-1] = np.std(delG_wc[:-2])

        #print sequence_context, delG_wc[-2], delG_wc[-1]
    
        # Get m1GC delG
        m1g_file = '../figure9/cary_uv/scaf2_' + sequence_context + '_m3t_done/final_all/' + 'test.csv'
        bf_m1g = pd.read_csv(m1g_file)
        # Get values
        #delG_m1g = bf_m1g['delG_37C'].values
        delS_m1g = bf_m1g['delS'].values
        delH_m1g = bf_m1g['delH'].values
        delG_m1g = delH_m1g - ((273.16+25)*delS_m1g)
        delG_m1g[-2] = np.average(delG_m1g[:-2])
        delG_m1g[-1] = np.std(delG_m1g[:-2])
        #print sequence_context, delG_m1g[-2], delG_m1g[-1]
        #print
    
        # Shift seq context to signature notation
        seq_context_new = comp_nts[dummy_inner] + 't' + comp_nts[dummy_outer]
        #sequence_contexts.append("5'-" + seq_context_new + "-3'")

        #sequence_contexts.append("5'-" + sequence_context + "-3'")
        sequence_contexts.append(sequence_context)

        if sequence_context in nontwo_state:
            delG_wc_values.append(5.0)
            delG_wc_error_values.append(0.0)
            delG_m1g_values.append(1.0)
            delG_m1g_error_values.append(0.0)

            delH_wc_values.append(0.0)
            delH_wc_error_values.append(0.0)
            delH_m1g_values.append(0.0)
            delH_m1g_error_values.append(0.0)

            delS_wc_values.append(0.0)
            delS_wc_error_values.append(0.0)
            delS_m1g_values.append(0.0)
            delS_m1g_error_values.append(0.0)

        else:
            delG_wc_values.append(delG_wc[-2])
            delG_wc_error_values.append(delG_wc[-1])
            delG_m1g_values.append(delG_m1g[-2])
            delG_m1g_error_values.append(delG_m1g[-1])

            delH_wc_values.append(delH_wc[-2])
            delH_wc_error_values.append(delH_wc[-1])
            delH_m1g_values.append(delH_m1g[-2])
            delH_m1g_error_values.append(delH_m1g[-1])

            delS_wc_values.append(delS_wc[-2])
            delS_wc_error_values.append(delS_wc[-1])
            delS_m1g_values.append(delS_m1g[-2])
            delS_m1g_error_values.append(delS_m1g[-1])

print "wc stats", min(delG_wc_values), max(delG_wc_values)
print "m1g stats", min(delG_m1g_values), max(delG_m1g_values)

# Convert lists to arrays
delG_wc_values = unumpy.uarray(delG_wc_values, delG_wc_error_values)
delG_m1g_values = unumpy.uarray(delG_m1g_values, delG_m1g_error_values)
delH_wc_values = unumpy.uarray(delH_wc_values, delH_wc_error_values)
delH_m1g_values = unumpy.uarray(delH_m1g_values, delH_m1g_error_values)
delS_wc_values = unumpy.uarray(delS_wc_values, delS_wc_error_values)
delS_m1g_values = unumpy.uarray(delS_m1g_values, delS_m1g_error_values)
deldelG = delG_m1g_values - delG_wc_values
deldelH = delH_m1g_values - delH_wc_values
deldelS = delS_m1g_values - delS_wc_values
deldelG_deltamelt = deldelG + ufloat(2.3917732231117785, 0.13772230621098872)
delG_m1g_values = delG_m1g_values * -1.0

print "deldelG stats", np.amin(deldelG), np.amax(deldelG)
print len(delG_m1g_values)
#f = open('m3t_37c.csv', 'w')
#f.write('context,deldelG_25C,deldelG_25C_error,delG_m1g_25C,delG_m1g_25C_error,delG_wc_25C,delG_wc_25C_error\n')
#for dummy in range(len(deldelG)):
#    print sequence_contexts[dummy], deldelG[dummy] 
#    f.write(sequence_contexts[dummy][3:6] + "," + str(unumpy.nominal_values(deldelG_deltamelt[dummy])) + "," + str(unumpy.std_devs(deldelG_deltamelt[dummy])) + "," + str(-1.0*unumpy.nominal_values(delG_m1g_values[dummy])) + "," + str(unumpy.std_devs(delG_m1g_values[dummy])) + "," + str(-1.0*unumpy.nominal_values(delG_wc_values[dummy])) + "," + str(unumpy.std_devs(delG_wc_values[dummy])) + "\n")
#f.close()

def plot_subplot(ax, xlabels, yvals, ylabel, ylims, fmt, xaxislabel):
    ax.bar(range(1, len(yvals)+1), unumpy.nominal_values(yvals), yerr=unumpy.std_devs(yvals), width=0.5, align='center', color='k', edgecolor='k', linewidth=linew, ecolor='k', error_kw=dict(lw=linew, capsize=cpsize, capthick=linew)) 
    ax.set_xlim([0, len(delG_wc_values)+1])
    ax.set_ylabel(ylabel, fontsize=fs)
    ax.set_yticks(np.arange(ylims[0], ylims[1]+0.02, ylims[2]))
    ax.set_yticklabels([fmt%ele for ele in np.arange(ylims[0], ylims[1]+0.02, ylims[2])], fontsize=fs)
    ax.set_xticks(range(1, len(yvals)+1))
    ax.set_xticklabels([ele.upper() for ele in xlabels], rotation='vertical', fontsize=fs)
    ax.tick_params(axis='x', direction='out', width=2, top=False)
    ax.tick_params(axis='y', direction='out', width=2, right=False)
    ax.set_ylim([ylims[0], ylims[1]])
    ax.set_xlabel(xaxislabel, fontsize=fs) 

    for dummy in range(len(xlabels)):
        if xlabels[dummy] in nontwo_state:
            ax.text(dummy+1, 0.2, '*')
 
def plot_corr(ax, xvals, yvals, xlims, xticklims, ylims, yticklims, xlabel, ylabel, xfmt, yfmt):
    ax.errorbar(1.0 * unumpy.nominal_values(xvals), 1.0 * unumpy.nominal_values(yvals), xerr=unumpy.std_devs(xvals), yerr=unumpy.std_devs(yvals), fmt='o', color='k', mec='k', markersize=1, elinewidth=linew, mew=0, capthick=linew, capsize=cpsize)
    ax.errorbar(1.0 * unumpy.nominal_values(xvals), 1.0 * unumpy.nominal_values(yvals), fmt='o', color='r', mec='k', markersize=15, elinewidth=linew, mew=linew, capthick=linew, capsize=cpsize)
    #ax.set_xlim([xlims[0], xlims[1]])
    #ax.set_ylim([ylims[0], ylims[1]])
    ax.set_xlabel(xlabel, fontsize=fs)
    ax.set_ylabel(ylabel, fontsize=fs)
    #ax.set_yticks(np.arange(yticklims[0], yticklims[1]+0.001, yticklims[2]))
    #ax.set_yticklabels([yfmt%ele for ele in np.arange(yticklims[0], yticklims[1]+0.001, yticklims[2])], fontsize=fs)
    #ax.set_xticks(np.arange(xticklims[0], xticklims[1]+0.001, xticklims[2]))
    #ax.set_xticklabels([xfmt%ele for ele in np.arange(xticklims[0], xticklims[1]+0.001, xticklims[2])], fontsize=fs)
    ax.tick_params(axis='x', direction='out', width=2, top=False)
    ax.tick_params(axis='y', direction='out', width=2, right=False)

# Plotting
fig, ax = plt.subplots()
plot_subplot(ax, sequence_contexts, delG_m1g_values, '$\Delta G^{\circ}_{melt,25 ^{\circ} C}$' + ' (kcal/mol)', [0.0, 14.0, 3.0], '%d', 'Sequence Contexts')
title = ax.set_title("A-m3T base pairs", fontsize=fs)
title.set_position([0.5, 1.02])


plt.tight_layout()
#plt.show()
plt.savefig('figs12h_2.pdf')
