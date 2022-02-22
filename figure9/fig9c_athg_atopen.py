import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib

matplotlib.rcParams['axes.linewidth'] = 2.0

def linear_fit(x, b):
    return x + b

m1a = pd.read_csv('cary_uv/m1a_25c_new.csv')
m3t = pd.read_csv('cary_uv/m3t_25c_new.csv')

context_exclude_m1a = ['aag', 'gac', 'tag', 'cac']
context_exclude_m3t = ['aaa', 'aag', 'gac', 'tag', 'cac', 'gag']
context_exclude = context_exclude_m1a + context_exclude_m3t
context_exclude = list(set(context_exclude))

m1a = m1a[~m1a['context'].isin(context_exclude)]
m3t = m3t[~m3t['context'].isin(context_exclude)]
#for dummy in range(len(m1a['context'].values)):
#    print m1a['context'].iloc[dummy], m3t['context'].iloc[dummy]

popt, pcov = curve_fit(linear_fit, m1a['deldelG_25C'].values, m3t['deldelG_25C'].values, p0=[1])
#print popt
print 'r = ', np.corrcoef(m1a['deldelG_25C'].values, m3t['deldelG_25C'].values)[0,1]
#print np.corrcoef(m3t['deldelG_25C'].values, m1a['deldelG_25C'].values)[0,1]
eqn = 'y = x + ' + '%2.1f'%popt + ', r = ' + '%3.2f'%np.corrcoef(m1a['deldelG_25C'].values, m3t['deldelG_25C'].values)[0,1]

fs = 24
plt.figure(1)
plt.tick_params(axis='x', direction='out', which='major', length=4, width=2, top=False, pad=2)
plt.tick_params(axis='y', direction='out', which='major', length=4, width=2, right=False, pad=2)
#plt.title(eqn)
plt.text(5.8, 6.1, '$r$' + '=' + '%3.2f'%np.corrcoef(m1a['deldelG_25C'].values, m3t['deldelG_25C'].values)[0,1], fontsize=fs)
#plt.plot([2.0, 6.5], [2.0 + popt[0], 6.5 + popt[0]], color='k', linewidth=2)
#plt.plot(m1a['deldelG_25C'].values, m3t['deldelG_25C'].values, linewidth=0, marker='o', color='r', mec='k', mew=2, markersize=16)
plt.errorbar(m1a['deldelG_25C'].values, m3t['deldelG_25C'].values, xerr=m1a['deldelG_25C_error'].values, yerr=m3t['deldelG_25C_error'].values, fmt='o', color='k', mec='k', mew=2, markersize=16, elinewidth=4.0, capthick=4.0, capsize=8)
#for dummy in range(len(m1a['context'].values)):
#    plt.text(m1a['deldelG_25C'].iloc[dummy]+0.05, m3t['deldelG_25C'].iloc[dummy]+0.05, m1a['context'].iloc[dummy], fontsize=16)
plt.xlabel('$\Delta G^{\\circ}_{conf,25^{\\circ} C}$' + ' (kcal/mol)\n' + ' A-T Hoogsteen', fontsize=fs)
plt.ylabel('$\Delta G^{\\circ}_{conf,25^{\\circ} C}$' + ' (kcal/mol)\n' + ' A-T Opening', fontsize=fs)
#plt.xlim([2, 5.])
plt.xlim([3, 6.5])
plt.ylim([6, 9.])
#plt.xticks([2.0, 3.0, 4.0, 5.0], fontsize=fs)
plt.xticks([3.0, 4.0, 5.0, 6.0], fontsize=fs)
plt.yticks([6.0, 7.0, 8.0, 9.0], fontsize=fs)
plt.savefig('fig9c_athg_atopen.pdf')
#plt.show()
