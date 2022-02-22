# Plot delta-melt correlation for Hoogsteen cooperativity
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from uncertainties import ufloat
from uncertainties import unumpy
import pandas as pd

def compute_deldelG(wt_filename, mt_filename, label, add_status):
    ''' Compute deldelG at temperature t '''
    global uv_values, uv_error
    wt_data = pd.read_csv(wt_filename)
    mt_data = pd.read_csv(mt_filename)
    wt_delg = ufloat(wt_data[label].iloc[-2], wt_data[label].iloc[-1])
    mt_delg = ufloat(mt_data[label].iloc[-2], mt_data[label].iloc[-1])
    deldelg = mt_delg - wt_delg
    return deldelg

linew = 1
cpsize = 8
fs = 12

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial']
mpl.rcParams['axes.linewidth'] = 1.0

fig, ax = plt.subplots()

offset_athg = ufloat(0.33, 0.06)
offset_gchg = ufloat(-0.12, 0.07)

# Data
# a6dna m1a16m1g10 
a6dna_m1a16m1g10_ph5p4hs_double = compute_deldelG('a6dna_ph5.4nbhs/test.csv', 'a6dna_m1a16m1g10_ph5.4hsnb/test.csv', 'delG_25C', 0)
a6dna_m1a16m1g10_ph5p4hs_single1 = compute_deldelG('a6dna_ph5.4nbhs/test.csv', 'a6dna_m1g10_keck_ph5.4nbhs/test.csv', 'delG_25C', 0)
a6dna_m1a16m1g10_ph5p4hs_single2 = compute_deldelG('a6dna_ph5.4nbhs/test.csv', 'a6dna_m1a16_ph5.4hsnb/test.csv', 'delG_25C', 0)
a6dna_m1a16m1g10_ph5p4hs_coop = a6dna_m1a16m1g10_ph5p4hs_single1 + a6dna_m1a16m1g10_ph5p4hs_single2 - a6dna_m1a16m1g10_ph5p4hs_double
#print a6dna_m1a16m1g10_ph5p4hs_single1, a6dna_m1a16m1g10_ph5p4hs_single2, a6dna_m1a16m1g10_ph5p4hs_double
print "A6DNA m1A16-m1G10 coop = ", a6dna_m1a16m1g10_ph5p4hs_coop, ' kcal/mol'

# a2dna m1g6-m1g20
a2dna_m1g6m1g20_ph5p4ns_double = compute_deldelG('a2dna_ph5.4hsnb/test.csv', 'a2dna_m1g6m1g20_ph5.4hsnb/test.csv', 'delG_25C', 0)
a2dna_m1g6m1g20_ph5p4ns_single1 = compute_deldelG('a2dna_ph5.4hsnb/test.csv', 'a2dna_m1g6_ph5.4hsnb/test.csv', 'delG_25C', 0)
a2dna_m1g6m1g20_ph5p4ns_single2 = compute_deldelG('a2dna_ph5.4hsnb/test.csv', 'a2dna_m1g20_ph5.4hsnb/test.csv', 'delG_25C', 0)
a2dna_m1g6m1g20_ph5p4ns_coop = a2dna_m1g6m1g20_ph5p4ns_single1 + a2dna_m1g6m1g20_ph5p4ns_single2 - a2dna_m1g6m1g20_ph5p4ns_double
#print a2dna_m1g6m1g20_ph5p4ns_single1, a2dna_m1g6m1g20_ph5p4ns_single2, a2dna_m1g6m1g20_ph5p4ns_double
print "A2DNA m1G6-m1G20 coop = ", a2dna_m1g6m1g20_ph5p4ns_coop, ' kcal/mol'

# a2dna m1a16m1a17
a2dna_m1a16m1a17_ph6p8hs_double = compute_deldelG('a2dna_ph6.8hsnb/test.csv', 'a2dna_m1a16m1a17_pH6.8hsnb/test.csv', 'delG_25C', 0)
a2dna_m1a16m1a17_ph6p8hs_single1 = compute_deldelG('a2dna_ph6.8hsnb/test.csv', 'a2dna_m1a16_pH6.8hsnb/test.csv', 'delG_25C', 0)
a2dna_m1a16m1a17_ph6p8hs_single2 = compute_deldelG('a2dna_ph6.8hsnb/test.csv', 'a2dna_m1a17_ph6.8hsnb/test.csv', 'delG_25C', 0)
a2dna_m1a16m1a17_ph6p8hs_coop = a2dna_m1a16m1a17_ph6p8hs_single1 + a2dna_m1a16m1a17_ph6p8hs_single2 - a2dna_m1a16m1a17_ph6p8hs_double
#print a2dna_m1a16m1a17_ph6p8hs_single1, a2dna_m1a16m1a17_ph6p8hs_single2, a2dna_m1a16m1a17_ph6p8hs_double
print 'a2dna m1a16-m1a17 coop = ', a2dna_m1a16m1a17_ph6p8hs_coop, ' kcal/mol'

# a2dna m1a17, 7
a2dna_m1a7m1a17_ph6p8hs_double = compute_deldelG('a2dna_ph6.8hsnb/test.csv', 'a2dna_m1a7m1a17_ph6.8hsnb/test.csv', 'delG_25C', 0)
a2dna_m1a7m1a17_ph6p8hs_single1 = compute_deldelG('a2dna_ph6.8hsnb/test.csv', 'a2dna_m1a7_ph6.8hsnb/test.csv', 'delG_25C', 0)
a2dna_m1a7m1a17_ph6p8hs_single2 = compute_deldelG('a2dna_ph6.8hsnb/test.csv', 'a2dna_m1a17_ph6.8hsnb/test.csv', 'delG_25C', 0)
a2dna_m1a7m1a17_ph6p8hs_coop = a2dna_m1a7m1a17_ph6p8hs_single1 + a2dna_m1a7m1a17_ph6p8hs_single2 - a2dna_m1a7m1a17_ph6p8hs_double
#print a2dna_m1a7m1a17_ph6p8hs_single1, a2dna_m1a7m1a17_ph6p8hs_single2, a2dna_m1a7m1a17_ph6p8hs_double
print 'a2dna m1a7, 17 coop = ', a2dna_m1a7m1a17_ph6p8hs_coop, ' kcal/mol'

uv_values = np.array([a6dna_m1a16m1g10_ph5p4hs_coop, a2dna_m1g6m1g20_ph5p4ns_coop, a2dna_m1a16m1a17_ph6p8hs_coop, a2dna_m1a7m1a17_ph6p8hs_coop])

bar_width = 1.0
sample_spacing = 2.0
plt.tick_params(axis='x', top=False, length=12.0, width=1, direction='out', pad=5)
plt.tick_params(axis='y', right=False, length=8.0, width=1, direction='out', pad=10)

ax.bar([0], unumpy.nominal_values(uv_values)[0], yerr=unumpy.std_devs(uv_values)[0], color='r', width=bar_width, align='center', edgecolor='r', linewidth=0, ecolor='r', error_kw=dict(lw=linew, capsize=cpsize, capthick=linew))
ax.bar(sample_spacing, unumpy.nominal_values(uv_values)[1], yerr=unumpy.std_devs(uv_values)[1], color='magenta', width=bar_width, align='center', edgecolor='magenta', linewidth=0, ecolor='magenta', error_kw=dict(lw=linew, capsize=cpsize, capthick=linew))
ax.bar(sample_spacing*2, unumpy.nominal_values(uv_values)[2], yerr=unumpy.std_devs(uv_values)[2], color='blue', width=bar_width, align='center', edgecolor='blue', linewidth=0, ecolor='blue', error_kw=dict(lw=linew, capsize=cpsize, capthick=linew))
ax.bar(sample_spacing*3, unumpy.nominal_values(uv_values)[3], yerr=unumpy.std_devs(uv_values)[3], color='brown', width=bar_width, align='center', edgecolor='brown', linewidth=0, ecolor='brown', error_kw=dict(lw=linew, capsize=cpsize, capthick=linew))

#ax.set_xticks((np.array(range(len(uv_values)))*sample_spacing))
#ax.set_xticklabels([str(ele+1) for ele in range(len(uv_values))], fontsize=fs)
ax.set_xticks([])
plt.xlim([-1.0, 7.0])
#plt.legend(fontsize=fs)
plt.ylim([0., 4.])
ax.set_yticks([1., 2., 3., 4.])
ax.set_yticklabels([str(int(ele)) for ele in [1., 2., 3., 4.]], fontsize=fs)
plt.ylabel('$\Delta G^{\circ}_{Coop,HG}$' + '\n(kcal/mol)', fontsize=fs)
plt.savefig('fig8b.pdf')
##plt.show()
