import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib

matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = 'Arial'
matplotlib.rcParams['axes.linewidth'] = 2.0
fs = 24
def linear_fit(x, b):
    return x + b

m3t = pd.read_csv('../figure9/cary_uv/m1a_25C_new.csv')
#m3t = pd.read_csv('m1a_37c.csv')

context_exclude_m3t = ['aag', 'gac', 'tag', 'cac']

m3t = m3t[~m3t['context'].isin(context_exclude_m3t)]
#for dummy in range(len(m1a['context'].values)):
#    print m1a['context'].iloc[dummy], m3t['context'].iloc[dummy]

popt1, pcov1 = curve_fit(linear_fit, m3t['delG_m1g_25C'].values,  m3t['deldelG_25C'].values, p0=[-1])
print popt1
print np.corrcoef(m3t['deldelG_25C'].values, m3t['delG_m1g_25C'].values)[0,1]
print np.corrcoef(m3t['delG_m1g_25C'].values, m3t['deldelG_25C'].values)[0,1]
eqn1 = 'y = x + ' + '%2.1f'%popt1 + ', r = ' + '%3.2f'%np.corrcoef(m3t['deldelG_25C'].values, m3t['delG_m1g_25C'].values)[0,1]

popt2, pcov2 = curve_fit(linear_fit, m3t['delG_wc_25C'].values, m3t['deldelG_25C'].values, p0=[1])
print popt2
print np.corrcoef(m3t['deldelG_25C'].values, m3t['delG_wc_25C'].values)[0,1]
print np.corrcoef(m3t['delG_wc_25C'].values, m3t['deldelG_25C'].values)[0,1]
eqn2 = 'y = x + ' + '%2.1f'%popt2 + ', r = ' + '%3.2f'%np.corrcoef(m3t['deldelG_25C'].values, m3t['delG_wc_25C'].values)[0,1]


fig, ax = plt.subplots(1, 2, figsize=(14, 6))

ms = 12

# Correlate m3AT deldelG with AT melting
#ax[0].set_title(eqn2)
#ax[0].plot([12, 17], [12 + popt2[0], 17 + popt2[0]], color='k')
ax[0].errorbar(m3t['delG_wc_25C'].values, m3t['deldelG_25C'].values, xerr=m3t['delG_wc_25C_error'].values, yerr=m3t['deldelG_25C_error'].values, fmt='o', color='k', markersize=ms, elinewidth=2, capthick=2)
#for dummy in range(len(m3t['context'].values)):
#    ax[0].text(m3t['delG_wc_25C'].iloc[dummy]+0.05, m3t['deldelG_25C'].iloc[dummy]+0.05, m3t['context'].iloc[dummy], fontsize=16)
ax[0].text(15.5, 2.1, '$r$' + '=%3.2f'%(np.corrcoef(m3t['deldelG_25C'].values, m3t['delG_wc_25C'].values)[0,1]), fontsize=fs)
ax[0].set_xlabel('$\Delta G_{melting, 25^{\\circ}C}^{\\circ}$' + ' A-T\n(kcal/mol)', fontsize=fs)
#ax[0].set_ylabel('A-T Hoogsteen\n' + '$\Delta \Delta G_{delta-melt, 25 ^{\\circ}C}^{\\circ}$' + ' (kcal/mol)', fontsize=fs)
ax[0].set_ylabel('A-T Hoogsteen\n' + '$\Delta G_{conf, 25 ^{\\circ}C}^{\\circ}$' + ' (kcal/mol)', fontsize=fs)
ax[0].set_xlim([12, 17.2])
ax[0].set_ylim([2, 6.5])
ax[0].tick_params(axis='x', which='minor', direction='out', length=0.0)
ax[0].tick_params(axis='x', which='major', direction='out', length=6.0, width=2.0, top=False)
ax[0].tick_params(axis='y', which='minor', direction='out', length=0.0)
ax[0].tick_params(axis='y', which='major', direction='out', length=6.0, width=2.0, right=False)
ax[0].set_xticks([12., 14., 16.])
ax[0].set_xticklabels(['12', '14', '16'], fontsize=fs)
ax[0].set_yticks([2.0, 4.0, 6.0])
ax[0].set_yticklabels(['2', '4', '6'], fontsize=fs)


# Correlate m3AT deldelG with m3AT melting
#ax[1].set_title(eqn1)
#ax[1].plot([9.5, 12.0], [9.5 + popt1[0], 12.0 + popt1[0]], color='k')
#ax[1].plot(m3t['delG_m1g_25C'].values, m3t['deldelG_25C'].values, linewidth=0, marker='o', color='k', markersize=ms)
ax[1].errorbar(m3t['delG_m1g_25C'].values, m3t['deldelG_25C'].values, xerr=m3t['delG_m1g_25C_error'].values, yerr=m3t['deldelG_25C_error'].values, fmt='o', color='k', markersize=ms, elinewidth=2, capthick=2)
#for dummy in range(len(m3t['context'].values)):
#    ax[1].text(m3t['delG_m1g_25C'].iloc[dummy]+0.05, m3t['deldelG_25C'].iloc[dummy]+0.05, m3t['context'].iloc[dummy], fontsize=16)

ax[1].text(11.3, 2.1, '$r$' + '=%3.2f'%(np.corrcoef(m3t['deldelG_25C'].values, m3t['delG_m1g_25C'].values)[0,1]), fontsize=fs)
ax[1].set_xlabel('$\Delta G_{melting, 25^{\\circ}C}^{\\circ}$' + ' m' + '$^{1}$' + 'A-T\n(kcal/mol)', fontsize=fs)
ax[1].set_xlim([9.5, 12.0])
ax[1].set_ylim([2, 6.5])
ax[1].set_xticks([10., 11., 12.])
ax[1].set_xticklabels(['10', '11', '12'], fontsize=fs)
ax[1].set_yticks([2.0, 4.0, 6.0])
ax[1].set_yticklabels(['2', '4', '6'], fontsize=fs)
ax[1].tick_params(axis='x', which='minor', direction='out', length=0.0)
ax[1].tick_params(axis='x', which='major', direction='out', length=6.0, width=2.0, top=False)
ax[1].tick_params(axis='y', which='minor', direction='out', length=0.0)
ax[1].tick_params(axis='y', which='major', direction='out', length=6.0, width=2.0, right=False)

plt.savefig('figs14a.pdf')
#plt.show()
