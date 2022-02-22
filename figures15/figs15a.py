import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from uncertainties import ufloat
from uncertainties import unumpy
from matplotlib import rcParams
from scipy.optimize import curve_fit
from scipy.stats import pearsonr
import matplotlib.gridspec as gridspec

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

def draw_bar(ax, val1, val2, label1, label2, yticks, fmt_string, ylabel):
    ax.bar([0], unumpy.nominal_values(val2), yerr=unumpy.std_devs(val2), width=1, align='center', color='orange', label=label2, edgecolor='orange', linewidth=0, ecolor='orange', error_kw=dict(lw=cpthick, capsize=cpsize, capthick=cpthick))
    ax.bar([1], unumpy.nominal_values(val1), yerr=unumpy.std_devs(val1), width=1, align='center', color='red', label=label1, edgecolor='red', linewidth=0, ecolor='red', error_kw=dict(lw=cpthick, capsize=cpsize, capthick=cpthick))
    ax.set_xlim([-1.0, 2.0])
    ax.tick_params(axis='x', which='both', direction='out', pad=6, width=elw, length=0, top=False)
    ax.tick_params(axis='y', which='minor', direction='out', pad=6, width=elw, length=0, right=False)
    ax.tick_params(axis='y', which='major', direction='out', pad=6, width=elw, length=tick_length, right=False)
    ax.set_xticks([])
    ax.set_ylim([yticks[0], yticks[1]])
    ax.set_yticks(np.arange(yticks[0], yticks[1], yticks[2]))
    ax.set_yticklabels([fmt_string%ele for ele in np.arange(yticks[0], yticks[1], yticks[2])], fontsize=fs)
    leg=ax.legend(fontsize=fs)
    leg.get_frame().set_linewidth(elw)
    ax.set_ylabel(ylabel, fontsize=fs)

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

def get_therm(ax, delg_array, t_array, xticks, yticks, title):
    m_vals = []
    b_vals = []
    delg_vals = unumpy.nominal_values(delg_array)
    delg_errors = unumpy.std_devs(delg_array)
    delg_error_matrix = np.zeros((mc_iterations, np.size(delg_array)))
    for dummy in range(np.size(delg_array)):
        op = (np.ones(mc_iterations) * delg_vals[dummy]) + np.random.normal(0.0, delg_errors[dummy], mc_iterations)
        #op = (np.ones(mc_iterations) * delg_vals[dummy])
        delg_error_matrix[:, dummy] = op.transpose()
    for dummy in range(mc_iterations):
        popt, pcov = curve_fit(linear_fit, t_array, delg_error_matrix[dummy, :], p0=[0.1, 1])
        m_vals.append(popt[0])
        b_vals.append(popt[1])
    m_vals = np.array(m_vals)
    b_vals = np.array(b_vals)
    m_avg = np.mean(m_vals)
    m_std = np.std(m_vals)
    b_avg = np.mean(b_vals)
    b_std = np.std(b_vals)
    line_coeff = [m_avg, b_avg]
    m_net = ufloat(m_avg, m_std)
    b_net = ufloat(b_avg, b_std)
    print "M = ", m_avg, "+/-", m_std
    print "b = ", b_avg, "+/-", b_std

    #fig2, ax2 = plt.subplots(2)
    #ax2[0].set_title('m')
    #ax2[0].hist(m_vals, 50, range=(np.amin(m_vals), np.amax(m_vals)))
    #ax2[1].set_title('b')
    #ax2[1].hist(b_vals, 50, range=(np.amin(b_vals), np.amax(b_vals)))
 
    #fig3, ax3 = plt.subplots()
    #H, xedges, yedges = np.histogram2d(m_vals, b_vals, 20)
    #X, Y = np.meshgrid(xedges, yedges)
    #ax3.pcolormesh(X, Y, H)
    #ax3.set_xlabel('m')
    #ax3.set_ylabel('b')

    predicted_delg = np.array([((line_coeff[0] * ele * -1.0) + line_coeff[1]) for ele in t_array])
    print '*', pearsonr(delg_vals, predicted_delg)
    print pearsonr(predicted_delg, delg_vals)

    ax.errorbar(t_array, delg_vals, markersize=msize, capthick=cpthick, capsize=cpsize, mec='k', mfc='k', mew=mewal, fmt='o', c='k', yerr=delg_errors, linewidth=elw)
    ax.plot([t_array[0]-2, t_array[-1]+2], [(line_coeff[0]*(t_array[0]-2)*-1.0)+line_coeff[1], (-1.0*line_coeff[0]*(t_array[-1]+2))+line_coeff[1]], color='k', linewidth=elw)
    m_vals = np.array([list(m_vals)])
    b_vals = np.array([list(b_vals)])
    x_dummy = np.linspace(t_array[0]-2.0, t_array[-1] + 2.0, 100)
    x_dummy_mat = np.tile(x_dummy, (mc_iterations, 1))
    m_matrix = np.tile(m_vals.transpose(), (1, 100))
    b_matrix = np.tile(b_vals.transpose(), (1, 100))
    line_matrix = np.multiply(-1.0*m_matrix, x_dummy_mat) + b_matrix
    ax.fill_between(x_dummy, np.amin(line_matrix, axis=0), np.amax(line_matrix, axis=0), color='#B0E2FF')

    # tick params
    ax.tick_params(axis='x', which='minor', direction='out', pad=6, width=elw, length=0, top=False)
    ax.tick_params(axis='x', which='major', direction='out', pad=6, width=elw, length=tick_length, top=False)
    ax.tick_params(axis='y', which='minor', direction='out', pad=6, width=elw, length=0, right=False)
    ax.tick_params(axis='y', which='major', direction='out', pad=6, width=elw, length=tick_length, right=False)
    
    ax.set_xlim([t_array[0]-2.0, t_array[-1]+2.0])
    ax.set_xticks(np.arange(xticks[0], xticks[1]+1, xticks[2]))
    ax.set_xticklabels([str(int(ele)) for ele in np.arange(xticks[0], xticks[1]+1, xticks[2])], fontsize=fs)
    ax.set_ylim([yticks[0], yticks[1]])
    ax.set_yticks(np.arange(yticks[0], yticks[1], yticks[2]))
    ax.set_yticklabels(['%2.1f'%ele for ele in np.arange(yticks[0], yticks[1], yticks[2])], fontsize=fs)
    ax.set_xlabel('Temperature (K)', fontsize=fs)
    ax.set_ylabel('$\Delta G^{\circ}_{conf}$' + '(i) ' + '(kcal/mol)', fontsize=fs)
    title = ax.set_title(title, fontsize=fs)
    title.set_position([0.5, 1.02])
    return [m_net, b_net]
   
fig = plt.figure(figsize=(18, 6.5))
gs = gridspec.GridSpec(ncols=7, nrows=2)
ax_one = fig.add_subplot(gs[:, :3])
ax_two = fig.add_subplot(gs[:, 3:5])
ax_three = fig.add_subplot(gs[:, 5:])

## A6DNA-A16 (pH 5.4, 25mM NaCl)
pb_17p5c = ufloat(0.463, 0.019)
pb_20c = ufloat(0.459, 0.030)
pb_22p5c = ufloat(0.559, 0.047)
pb_25c = ufloat(0.629, 0.081)
p_array = np.array([pb_17p5c, pb_20c, pb_22p5c, pb_25c])/100.
t_array = np.array([17.5, 20.0, 22.5, 25.0]) + 273.16
delg_array = compute_delG(p_array, t_array)
a6dna_h = compute_deldelG('../figure3/a6dna_ph5.4/test.csv', '../figure3/a6dna_m1a16_ph5.4/test.csv', 'delH')
a6dna_s = compute_deldelG('../figure3/a6dna_ph5.4/test.csv', '../figure3/a6dna_m1a16_ph5.4/test.csv', 'delS')
[a6dna_m, a6dna_b] = get_therm(ax_one, delg_array, t_array, [290., 300., 4.], [2.6, 3.4, 0.3], 'A-T Hoogsteen')
#print "CHARTT", a6dna_h, a6dna_s*1000.
#print "NMR", a6dna_b, a6dna_m*1000.
draw_bar(ax_two, a6dna_h, a6dna_b, '$\Delta\Delta H^{\circ}_{melt}$' + '(i)', '$\Delta H^{\circ}_{conf}$' + "(i)", [0., 35., 10], '%d', 'Enthalpy (kcal/mol)')
draw_bar(ax_three, a6dna_s*1000, a6dna_m*1000, '$\Delta\Delta S^{\circ}_{melt}$'+"(i)", '$\Delta S^{\circ}_{conf}$'+"(i)", [0., 95., 20], '%d', 'Entropy (cal/mol/K)')


## A6DNA-A16 (pH 5.4, 25mM NaCl)
#pb_10c = ufloat(0.180, 0.008)
#pb_12p5c = ufloat(0.286, 0.046)
#pb_15c = ufloat(0.311, 0.074)
#p_array = np.array([pb_10c, pb_12p5c, pb_15c])/100.
#t_array = np.array([10., 12.5, 15.]) + 273.16
#delg_array = compute_delG(p_array, t_array)
#fig, ax = plt.subplots()
#get_therm(ax, delg_array, t_array)
#plt.show()

 
## hpGGACU m6a6 
#pb_37c = ufloat(0.643, 0.006)
#pb_55c = ufloat(1.205, 0.042)
#pb_65c = ufloat(1.526, 0.013)
#p_array = np.array([pb_37c, pb_55c, pb_65c])/100.
#t_array = np.array([37., 55., 65.]) + 273.16
#delg_array = compute_delG(p_array, t_array)
#hpggacu_h = compute_deldelG('dsgbc_m6a_ph6.8bb/test.csv', 'dsgbc_dm6a_ph6.8bb/test.csv', 'delH')
#hpggacu_s = compute_deldelG('dsgbc_m6a_ph6.8bb/test.csv', 'dsgbc_dm6a_ph6.8bb/test.csv', 'delS')
#print "CHARTT", hpggacu_h, hpggacu_s*1000.
#[hpggacu_m, hpggacu_b] = get_therm(ax_one, delg_array, t_array, [310., 340., 10.], [2.7, 3.15, 0.1], '${N}^{6}$' + '-methylamino rotation')
#print "NMR", hpggacu_b, hpggacu_m*1000.
#draw_bar(ax_two, hpggacu_h, hpggacu_b, '$\Delta\Delta H^{\circ}_{unfold}$', '$\Delta H^{\circ}_{conf,NMR}$', [0., 17., 5], '%d', 'Enthalpy (kcal/mol)')
#draw_bar(ax_three, hpggacu_s*1000, hpggacu_m*1000, '$\Delta\Delta S^{\circ}_{unfold}$', '$\Delta S^{\circ}_{conf,NMR}$', [0., 40., 10], '%d', 'Entropy (cal/mol/K)')
plt.tight_layout()
plt.savefig('figs15a.pdf')
#plt.show()
