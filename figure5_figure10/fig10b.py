import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from uncertainties import ufloat
from uncertainties import unumpy

linew = 1
cpsize = 8
fs = 12

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial']
mpl.rcParams['axes.linewidth'] = 1.0

fig, ax = plt.subplots()

# declare rd and uv arrays
# CGT_5.4LS_30c
cgt_5p4ls_30c_rd = ufloat(2.954, 0.012)
cgt_5p4ls_30c_pred = ufloat(2.53, 0.21)

# CGT_5.4HS_30c
cgt_5p4hs_30c_rd = ufloat(3.734, 0.006)
cgt_5p4hs_30c_pred = ufloat(3.24, 0.22)

# TGC_5.4LS_30c
tgc_5p4ls_30c_rd = ufloat(3.755, 0.015)
tgc_5p4ls_30c_pred = ufloat(4.40, 0.12)

# TGC_5.4HS_30c
tgc_5p4hs_30c_rd = ufloat(4.38, 0.04)
tgc_5p4hs_30c_pred = ufloat(5.19, 0.21)

# TGC_5.4HS_40c
tgc_5p4hs_40c_rd = ufloat(4.51, 0.17)
tgc_5p4hs_40c_pred = ufloat(4.46, 0.15)

# TGC_5.4HS_15C
tgc_5p4hs_15c_rd = ufloat(0.0, 0.0)
tgc_5p4hs_15c_pred = ufloat(6.30, 0.33)

rd_net = np.array([cgt_5p4ls_30c_rd, cgt_5p4hs_30c_rd, tgc_5p4ls_30c_rd, tgc_5p4hs_30c_rd, tgc_5p4hs_40c_rd, tgc_5p4hs_15c_rd])
pred_net = np.array([cgt_5p4ls_30c_pred, cgt_5p4hs_30c_pred, tgc_5p4ls_30c_pred, tgc_5p4hs_30c_pred, tgc_5p4hs_40c_pred, tgc_5p4hs_15c_pred])

rmse = np.sqrt(unumpy.nominal_values(np.sum(np.square(rd_net[:-1] - pred_net[:-1]))/len(list(rd_net[:-1]))))
#print 'rmse = ', rmse

bar_width = 1.0
sample_spacing = 3.0
plt.tick_params(axis='x', top=False, length=12.0, width=1, direction='out', pad=5)
plt.tick_params(axis='y', right=False, length=8.0, width=1, direction='out', pad=10)

ax.bar(np.array(range(len(rd_net)))*sample_spacing, unumpy.nominal_values(rd_net), yerr=unumpy.std_devs(rd_net), color='#FFA500', width=bar_width, align='center', edgecolor='#FFA500', linewidth=0, ecolor='#FFA500', error_kw=dict(lw=linew, capsize=cpsize, capthick=linew), label='$\Delta G^{\circ}_{conf}(i)$')
ax.bar((np.array(range(len(rd_net)))*sample_spacing)+bar_width, unumpy.nominal_values(pred_net), yerr=unumpy.std_devs(pred_net), color='r', width=bar_width, align='center', edgecolor='r', linewidth=0, ecolor='r', error_kw=dict(lw=linew, capsize=cpsize, capthick=linew), label='$\Delta\Delta G^{\circ}_{delta-Melt}(i)$')
ax.text(3*sample_spacing, unumpy.nominal_values(rd_net)[3] + 0.2, '*')
ax.text(4*sample_spacing, unumpy.nominal_values(rd_net)[4] + 0.2, '*')

ax.set_xticks((np.array(range(len(rd_net)))*sample_spacing) + (bar_width*0.5))
ax.set_xticklabels([str(ele+1) for ele in range(len(rd_net))], fontsize=fs)
plt.xlim([-1.0, 17.0])
plt.legend(fontsize=fs, loc=2)
plt.ylim([0., 7.])
ax.set_yticks([0., 2., 4., 6.])
ax.set_yticklabels([str(int(ele)) for ele in [0., 2., 4., 6.]], fontsize=fs)
plt.ylabel('Free Energy\n' + ' (kcal/mol)', fontsize=fs)
plt.savefig('fig10b.pdf')
#plt.show()
