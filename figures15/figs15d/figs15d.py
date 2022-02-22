import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from uncertainties import ufloat
from uncertainties import unumpy
from matplotlib import rcParams
from scipy.optimize import curve_fit
from scipy.stats import pearsonr
import matplotlib.gridspec as gridspec
import os

rcParams['axes.linewidth'] = 4.0
rcParams['font.family']='Sans-serif'
rcParams['font.sans-serif']=['Arial']

elw = 4.0
fs = 27
tick_length = 8.0
msize = 16
cpthick=4
cpsize=8
mewal = 1.0
mc_iterations = 10000

# Get the delH and delS values from Chen and Russu 2004
os.system('python get_tbp_opening_thermodynamics.py')

# Get the delH and delS values for introducing A-m3T base pairs
# in TBP DNA
os.system('python get_m3t_thermodynamics.py')

# delG = delH - (T*delS)
def linear_fit(x, m, b):
    return ((-1.0*m*x)+b)

def compute_deldelG(wt_filename, mt_filename, label):
    ''' Compute deldelG at temperature t '''
    wt_data = pd.read_csv(wt_filename)
    mt_data = pd.read_csv(mt_filename)
    wt_delg = ufloat(wt_data[label].iloc[-2], wt_data[label].iloc[-1])
    mt_delg = ufloat(mt_data[label].iloc[-2], mt_data[label].iloc[-1])
    deldelg = mt_delg - wt_delg
    return deldelg

def compute_delG(p_arr, t_arr):
    ''' Compute delGs '''
    delg = []
    for dummy in range(len(p_arr)):
        print p_arr[dummy], t_arr[dummy]
        delg.append((-8.314*t_arr[dummy]/4184)*unumpy.log(p_arr[dummy]/(1-p_arr[dummy])))
    delg = np.array(delg)
    return delg

   
fig = plt.figure(figsize=(18, 6.2))
gs = gridspec.GridSpec(ncols=7, nrows=2)
ax_one = fig.add_subplot(gs[:, :3])
#ax_two = fig.add_subplot(gs[:, 3:5])
#ax_three = fig.add_subplot(gs[:, 5:])

rd = pd.read_csv('tbp_opening_delh_dels.csv')
rd_delh = unumpy.uarray(rd['delh(kcal/mol)'], rd['delh_err(kcal/mol)'])
rd_dels = unumpy.uarray(rd['dels(kcal/mol)'], rd['dels_err(kcal/mol)']) * 1000.
rd_resi = list(rd['resi'].values)

chartt = pd.read_csv('m3t_delh_dels.csv')
chartt_delh = unumpy.uarray(chartt['ddh(kcal/mol)'], chartt['ddh_error(kcal/mol)'])
chartt_dels = unumpy.uarray(chartt['dds(cal/mol/k)'], chartt['dds_error(cal/mol/k)']) * 1000.
chartt_resi = list(chartt['resi'].values)
print rd_resi, chartt_resi

delh_lim = [np.amin(np.array([np.amin(rd_delh), np.amin(chartt_delh)])), np.amax(np.array([np.amax(rd_delh), np.amax(chartt_delh)]))]
delh_lim = [unumpy.nominal_values(delh_lim[0]) - unumpy.std_devs(delh_lim[0]), unumpy.nominal_values(delh_lim[1]) + unumpy.std_devs(delh_lim[1])]
dels_lim = [np.amin(np.array([np.amin(rd_dels), np.amin(chartt_dels)])), np.amax(np.array([np.amax(rd_dels), np.amax(chartt_dels)]))]
dels_lim = [unumpy.nominal_values(dels_lim[0]) - unumpy.std_devs(dels_lim[0]), unumpy.nominal_values(dels_lim[1]) + unumpy.std_devs(dels_lim[1])]
delh_lim[0] = delh_lim[0] - 5.
delh_lim[1] = delh_lim[1] + 2.
dels_lim[0] = dels_lim[0] - 5.
dels_lim[1] = dels_lim[1] + 2.

def draw_corr(ax, x, y, xlabel, ylabel, lim, ticks, fmt_string, text_labels):
    ax.errorbar(unumpy.nominal_values(x), unumpy.nominal_values(y), markersize=msize, capthick=cpthick, capsize=cpsize, mec='k', mfc='k', mew=mewal, fmt='o', c='k', yerr=unumpy.std_devs(y), xerr=unumpy.std_devs(x), linewidth=elw)
    for dummy in range(len(text_labels)):
        ax.text(unumpy.nominal_values(x[dummy]), unumpy.nominal_values(y[dummy]), 'T' + str(text_labels[dummy]), fontsize=fs)
    ax.set_xlabel(xlabel, fontsize=fs)
    ax.set_ylabel(ylabel, fontsize=fs)
    ax.tick_params(axis='x', which='minor', direction='out', pad=6, width=elw, length=0, top=False)
    ax.tick_params(axis='x', which='major', direction='out', pad=6, width=elw, length=tick_length, top=False)
    ax.tick_params(axis='y', which='minor', direction='out', pad=6, width=elw, length=0, right=False)
    ax.tick_params(axis='y', which='major', direction='out', pad=6, width=elw, length=tick_length, right=False)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xticks(ticks)
    ax.set_xticklabels([fmt_string%ele for ele in ticks], fontsize=fs) 
    ax.set_yticks(ticks)
    ax.set_yticklabels([fmt_string%ele for ele in ticks], fontsize=fs) 
    ax.plot([lim[0], lim[1]], [lim[0], lim[1]], color='k', linewidth=elw)
    title = ax.set_title('Base Opening', fontsize=fs)
    title.set_position([0.5, 1.02])

draw_corr(ax_one, chartt_delh, rd_delh, '$\Delta\Delta H^{\circ}_{melt}$' + '(i) (kcal/mol)', '$\Delta H^{\circ}_{conf}$' + '(i) (kcal/mol)', delh_lim, np.arange(-10, 31, 10), '%d', rd_resi)

#plt.tight_layout()
plt.savefig('figs15d.pdf')
#plt.show()
