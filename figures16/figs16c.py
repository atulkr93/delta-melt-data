import math
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

matplotlib.rcParams['axes.linewidth'] = 2.0

delH = -84.17
delS = -0.231473
conc = np.logspace(0.0, math.log10(150.), 10.0)
print conc
tm_array = []
conc_log = []

for dummy_conc in conc:
    conc_log.append(math.log(dummy_conc * math.pow(10, -6) * 0.25))
    temp = ((8.314 / 4184 / delH) * math.log(dummy_conc * math.pow(10, -6) * 0.25)) + (delS/delH)
    tm = 1. / temp
    print dummy_conc, tm
    tm_array.append(tm)

fs=16
conc_log = np.array(conc_log)
tm_array = np.array(tm_array)
plt.figure(1)
plt.plot(-1.0*conc_log, np.reciprocal(tm_array) * math.pow(10, 4), color='b', linewidth=2.0, label='2-state Fit')
plt.plot(-1.0*conc_log, np.reciprocal(tm_array) * math.pow(10, 4), marker='o', color='k', linewidth=0.0, markersize=8, label='Data')
plt.tick_params(axis='x', which='both', width=2.0, length=6.0, direction='out', top=False)
plt.tick_params(axis='y', which='both', width=2.0, length=6.0, direction='out', right=False)
plt.xticks([10., 12., 14., 16.], ['%d'%ele for ele in [10., 12., 14., 16.]], fontsize=fs)
plt.yticks([29.8, 30.2, 30.6, 31.0], ['%.1f'%ele for ele in [29.8, 30.2, 30.6, 31.0]], fontsize=fs)
plt.xlabel('-ln(' + '$C_{t}$' + '/4)', fontsize=fs)
plt.ylabel('$10^{4}/T_{m}$' + ' (1/K)', fontsize=fs)
plt.legend(fontsize=fs, loc=2)
plt.savefig('figs16c.pdf')
#plt.show()
