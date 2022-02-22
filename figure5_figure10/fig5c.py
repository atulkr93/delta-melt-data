# plot delta-melt correlation for G-C+ Hoogsteen
from scipy.optimize import curve_fit
import math
import pandas as pd
from scipy.odr import *
import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib import rcParams
import matplotlib as mpl
from uncertainties import ufloat
from uncertainties import unumpy
from scipy.stats import pearsonr
from sklearn.metrics import r2_score

# Plot params
fs = 16
fmt_point = 'o'
msize = 16
elw = 1
cpthick =1
cpsize = 8
mewval = 1.0
at_color='cyan'
gc_color='orange'
coop_color = 'magenta'
rntp_color = 'green'
test_color = 'violet'

# Intiialize plots
# Initialize plot
# Plot params
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
rcParams['axes.linewidth'] = 1.

### Global variables
uv_values = []
uv_error = []
rd_values = []
rd_error = []


# Define fit function
def fit_func(B, x):
    return (B[0]*x) + B[1]

def linear_fit(x, m, b):
    return ((m*x)+b)

def linear_fit2(x, b):
    return ((x)+b)

def compute_deldelG(wt_filename, mt_filename, label, add_status):
    ''' Compute deldelG at temperature t '''
    # add status - controls whether data is added to plot correlation curve
    global uv_values, uv_error
    wt_data = pd.read_csv(wt_filename)
    mt_data = pd.read_csv(mt_filename)
    wt_delg = ufloat(wt_data[label].iloc[-2], wt_data[label].iloc[-1])
    mt_delg = ufloat(mt_data[label].iloc[-2], mt_data[label].iloc[-1])
    #print wt_delg, mt_delg
    deldelg = mt_delg - wt_delg
    if add_status == 1:
        uv_values.append(float(unumpy.nominal_values(deldelg)))
        uv_error.append(float(unumpy.std_devs(deldelg)))
    return deldelg

def compute_delG(pop, pop_error, t, add_status):
    ''' given population, compute delG '''
    global rd_values, rd_error
    pop_mag = ufloat(pop, pop_error) / 100.0
    delg = (-8.314 * (273.16 + t) / 4184) * unumpy.log(pop_mag/(1-pop_mag))
    if add_status == 1:
        rd_values.append(float(unumpy.nominal_values(delg)))
        rd_error.append(float(unumpy.std_devs(delg)))
    return delg

gc_color_status = 0
def plot_point(uv, rd, color_pt):
    ''' Plot point on correlation plot '''
    global gc_color_status
    #if color_pt == gc_color:
    #    gc_color_status = 1
    if color_pt == 'unknown':
        plt.errorbar(unumpy.nominal_values(uv), unumpy.nominal_values(rd), yerr=unumpy.std_devs(rd), xerr=unumpy.std_devs(uv), color='r', fmt=fmt_point, markersize=msize, elinewidth=elw, capthick=cpthick, capsize=cpsize, mec='r', mfc='white', mew=mewval, zorder=100)
    elif color_pt == 'blue2':
        plt.errorbar(unumpy.nominal_values(uv), unumpy.nominal_values(rd), yerr=unumpy.std_devs(rd), xerr=unumpy.std_devs(uv), color=test_color, fmt=fmt_point, markersize=msize, elinewidth=elw, capthick=cpthick, capsize=cpsize, mec=test_color, mfc='white', mew=mewval, zorder=50)
    else:
        plt.errorbar(unumpy.nominal_values(uv), unumpy.nominal_values(rd), yerr=unumpy.std_devs(rd), xerr=unumpy.std_devs(uv), color=color_pt, fmt=fmt_point, markersize=msize, elinewidth=elw, capthick=cpthick, capsize=cpsize, mec=color_pt, mfc=color_pt, mew=mewval, zorder=50)

    #pass

### Read data for getting correlation curve ###
a2dna_g10_ph4p4hs_ddg = compute_deldelG('a2dna_ph4.4hs/test.csv', 'a2dna_m1g10_ph4.4hsnb/test.csv', 'delG_25C', 1)
a2dna_g10_ph4p4hs_rd = compute_delG(2.586, 0.052, 25, 1)

actd_g7_ph5p3hs_ddg = compute_deldelG('actd_ph5.3hs/test.csv', 'actd_m1g_ph5.3hs/test.csv', 'delG_25C', 1)
actd_g7_ph5p3hs_rd = compute_delG(0.137, 0.036, 25, 1)

a2dna_g10_ph5p4_ddg = compute_deldelG('a2dna_ph5.4nb/test.csv', 'a2dna_m1g10_ph5.4nb/test.csv', 'delG_25C', 1)
a2dna_g10_ph5p4_rd = compute_delG(1.014, 0.024, 25, 1)

a6dna_g10_ph5p4_ddg = compute_deldelG('a6dna_ph5.4/test.csv', 'a6dna_m1g10_keck_ph5.4nb/test.csv', 'delG_25C', 1)
a6dna_g10_ph5p4_rd = compute_delG(0.717, 0.010, 25, 1)

a6dna_rg10_ph5p4_ddg = compute_deldelG('a6dna_rg10_ph5.4/test.csv', 'a6dna_m1rg10_ph5.4/test.csv', 'delG_25C', 1)
a6dna_rg10_ph5p4_rd = compute_delG(0.577, 0.105, 25, 1)

a6dna_g10_ph5p4hsnb_ddg = compute_deldelG('a6dna_ph5.4nbhs/test.csv', 'a6dna_m1g10_keck_ph5.4nbhs/test.csv', 'delG_25C', 1)
a6dna_g10_ph5p4hsnb_rd = compute_delG(0.251, 0.012, 25, 1)

# weak rd
a6dna_g10_ph6p8_ddg = compute_deldelG('a6dna_ph6.8/test.csv', 'a6dna_m1g10_keck_ph6.8nb/test.csv', 'delG_25C', 0)
a6dna_g10_ph6p8_rd = compute_delG(0.063, 0.015, 25, 0)

print 'Number of UV data points = ', len(uv_values), len(uv_error)
print 'Number of NMR RD data points = ', len(rd_values), len(rd_error)


### Read test data ###
scaf2_tgc_ph5p4ls_30c_ddg = compute_deldelG('scaf2_tgc_gc_ph5.4lsnb/test.csv', 'scaf2_tgc_m1gc_ph5.4lsnb/test.csv', 'delG_30C', 0)
scaf2_tgc_ph5p4ls_30c_rd = compute_delG(0.196, 0.005, 30., 0)

scaf2_tgc_ph5p4hs_30c_ddg = compute_deldelG('scaf2_tgc_gc_ph5.4hsnb/test.csv', 'scaf2_tgc_m1gc_ph5.4hsnb/test.csv', 'delG_30C', 0)
scaf2_tgc_ph5p4hs_30c_rd = compute_delG(0.069, 0.005, 30., 0)

scaf2_tgc_ph5p4hs_40c_ddg = compute_deldelG('scaf2_tgc_gc_ph5.4hsnb/test.csv', 'scaf2_tgc_m1gc_ph5.4hsnb/test.csv', 'delG_40C', 0)
scaf2_tgc_ph5p4hs_40c_rd = compute_delG(0.071, 0.019, 40., 0)

scaf2_cgt_ph5p4ls_30c_ddg = compute_deldelG('scaf2_cgt_gc_ph5.4lsnb/test.csv', 'scaf2_cgt_m1gc_ph5.4lsnb/test.csv', 'delG_30C', 0)
scaf2_cgt_ph5p4ls_30c_rd = compute_delG(0.736, 0.015, 30., 0)

scaf2_cgt_ph5p4hs_30c_ddg = compute_deldelG('scaf2_cgt_gc_ph5.4hsnb/test.csv', 'scaf2_cgt_m1gc_ph5.4hsnb/test.csv', 'delG_30C', 0)
scaf2_cgt_ph5p4hs_30c_rd = compute_delG(0.203, 0.002, 30., 0)

scaf2_tgc_ph5p4hs_15c_ddg = compute_deldelG('scaf2_tgc_gc_ph5.4hsnb/test.csv', 'scaf2_tgc_m1gc_ph5.4hsnb/test.csv', 'delG_15C', 0)

#print len(uv_values), len(uv_error)
#print len(rd_values), len(rd_error)

#print "scaf2_cgt_ph5p4hs_30c_ddg", scaf2_cgt_ph5p4hs_30c_ddg
#print "scaf2_cgt_ph5p4hs_30c_pred", (0.9*scaf2_cgt_ph5p4hs_30c_ddg)+0.3
#print "scaf2_cgt_ph5p4hs_30c_rd", scaf2_cgt_ph5p4hs_30c_rd

# convert to arrays
rd_values = np.array(rd_values)
rd_error = np.array(rd_error)
uv_values = np.array(uv_values)
uv_error = np.array(uv_error)

# Initialize MC iterations
mc_iterations = 10000

# initialize error matrices
uv_error_matrix = np.zeros((mc_iterations, len(uv_values)))
rd_error_matrix = np.zeros((mc_iterations, len(rd_values)))
for dummy in range(len(uv_values)):
    op = (np.ones(mc_iterations) * uv_values[dummy]) + np.random.normal(0.0, uv_error[dummy], mc_iterations)
    #op = (np.ones(mc_iterations) * uv_values[dummy])  
    uv_error_matrix[:, dummy] = op.transpose()
    op2 = (np.ones(mc_iterations) * rd_values[dummy]) + np.random.normal(0.0, rd_error[dummy], mc_iterations)
    #op2 = (np.ones(mc_iterations) * rd_values[dummy]) 
    rd_error_matrix[:, dummy] = op2.transpose()
#print uv_values
#print np.mean(uv_error_matrix, axis=0)
#print rd_values
#print np.mean(rd_error_matrix, axis=0)
#print
#print uv_error
#print np.std(uv_error_matrix, axis=0)
#print rd_error
#print np.std(rd_error_matrix, axis=0)
filename = 'results.csv'
f = open(filename, 'w')
f.write('m,b\n')

for dummy in range(mc_iterations):
    #popt, pcov = curve_fit(linear_fit, uv_error_matrix[dummy, :], rd_error_matrix[dummy, :], p0=[1, 1])
    #line_coeff = popt 
    popt, pcov = curve_fit(linear_fit2, uv_error_matrix[dummy, :], rd_error_matrix[dummy, :], p0=[1])
    line_coeff = list([1]) + list(popt) 

    #eqn = 'y = ' + '%3.2f'%(line_coeff[0]) + 'x + ' + '%3.2f'%(line_coeff[1])
    #predicted_rd = np.array([((line_coeff[0]*ele)+line_coeff[1]) for ele in uv_values])
    #rmse = np.sqrt(np.sum(np.square(rd_values - predicted_rd))/len(list(rd_values)))
    #r2 = 1 - (np.sum(np.square(predicted_rd - rd_values)) / np.sum(np.square(rd_values - np.mean(rd_values))))

    #print "iteration", dummy + 1
    #print "***", eqn, "***"
    #print "see rmse", rmse
    #print "see r2", r2
    f.write(str(line_coeff[0]) + ',' + str(line_coeff[1]) + "\n")
    

f.close()

# Read the output file and plot
bf = pd.read_csv(filename)

m_avg = np.mean(np.array(bf['m'].values))
m_std = np.std(np.array(bf['m'].values))
b_avg = np.mean(np.array(bf['b'].values))
b_std = np.std(np.array(bf['b'].values))
m_net = ufloat(m_avg, m_std)
b_net = ufloat(b_avg, b_std)
#rmse_avg = np.mean(np.array(bf['rmse'].values))
#rmse_std = np.std(np.array(bf['rmse'].values))
#r2_avg = np.mean(np.array(bf['r2'].values))
#r2_std = np.std(np.array(bf['r2'].values))
#print "M = ", m_avg, "+/-", m_std
print "c(i) = ", b_avg, "+/-", b_std
#print "rmse = ", rmse_avg, "+/-", rmse_std
#print "r2 = ", r2_avg, "+/-", r2_std

#fig, axarr = plt.subplots(2)
#axarr[0].set_title('m')
#axarr[0].hist(np.array(bf['m'].values), 50, range=(np.amin(np.array(bf['m'].values)), np.amax(np.array(bf['m'].values))))
#axarr[1].set_title('b')
#axarr[1].hist(np.array(bf['b'].values), 50, range=(np.amin(np.array(bf['b'].values)), np.amax(np.array(bf['b'].values))))
##axarr[1, 0].set_title('rmse')
##axarr[1, 0].hist(np.array(bf['rmse'].values), 50, range=(np.amin(np.array(bf['rmse'].values)), np.amax(np.array(bf['rmse'].values))))
##axarr[1, 1].set_title('r2')
##axarr[1, 1].hist(np.array(bf['r2'].values), 50, range=(np.amin(np.array(bf['r2'].values)), np.amax(np.array(bf['r2'].values))))
#plt.show()

#fig, axarr = plt.subplots()
##axarr.plot(np.array(bf['m'].values), np.array(bf['b'].values), color='k', linewidth=0, marker='o')
#H, xedges, yedges = np.histogram2d(np.array(bf['m'].values), np.array(bf['b'].values), 20)
#X, Y = np.meshgrid(xedges, yedges)
#axarr.pcolormesh(X, Y, H)
#axarr.set_xlabel('m')
#axarr.set_ylabel('b')
#plt.show()

line_coeff = [m_avg, b_avg] 
eqn = 'y = mx + b\nm=' + '%3.2f'%(m_avg) + '+/-' + '%3.2f'%(m_std) + '\nb=' + '%3.2f'%(b_avg) + '+/-' + '%3.2f'%(b_std)
#eqn = 'y = ' + '%3.2f'%(line_coeff[0]) + 'x + ' + '%3.2f'%(line_coeff[1])
predicted_rd = np.array([((line_coeff[0]*ele)+line_coeff[1]) for ele in uv_values])
rmse = np.sqrt(np.sum(np.square(rd_values - predicted_rd))/len(list(rd_values)))
r2 = 1 - (np.sum(np.square(predicted_rd - rd_values)) / np.sum(np.square(rd_values - np.mean(rd_values))))
#for dummy in range(len(predicted_rd)):
#    print predicted_rd[dummy], rd_values[dummy], -1.0*(predicted_rd[dummy] - rd_values[dummy])

# shaded lines
x_dummy = np.linspace(0.0, 6.0, 100)
x_dummy_mat = np.tile(x_dummy, (mc_iterations, 1))
m_vals = np.array([list(bf['m'].values)])
b_vals = np.array([list(bf['b'].values)])
m_matrix = np.tile(m_vals.transpose(), (1, 100))
b_matrix = np.tile(b_vals.transpose(), (1, 100))
line_matrix = np.multiply(m_matrix, x_dummy_mat) + b_matrix

fig, ax = plt.subplots()
plt.tick_params(axis='x', top=False, length=8.0, width=1, direction='out', pad=5)
plt.tick_params(axis='y', right=False, length=8.0, width=1, direction='out', pad=5)
#plt.errorbar(uv_values, rd_values, xerr=uv_error, yerr=rd_error, color='r', fmt='o')
#gc_color = '#FFA500'
gc_color = 'r'

plot_point(a6dna_g10_ph5p4hsnb_ddg, a6dna_g10_ph5p4hsnb_rd, 'r')
plot_point(a6dna_g10_ph5p4_ddg, a6dna_g10_ph5p4_rd, 'r')
plot_point(a6dna_g10_ph6p8_ddg, a6dna_g10_ph6p8_rd, 'unknown')

plot_point(a6dna_rg10_ph5p4_ddg, a6dna_rg10_ph5p4_rd, 'magenta')

plot_point(a2dna_g10_ph5p4_ddg, a2dna_g10_ph5p4_rd, 'blue')
plot_point(a2dna_g10_ph4p4hs_ddg, a2dna_g10_ph4p4hs_rd, 'blue')

plot_point(actd_g7_ph5p3hs_ddg, actd_g7_ph5p3hs_rd, 'brown')



test_color= 'k'
test_color2 = 'orange'
#plot_point(scaf2_tgc_ph5p4ls_30c_ddg, scaf2_tgc_ph5p4ls_30c_rd, test_color2)
#plot_point(scaf2_tgc_ph5p4hs_30c_ddg, scaf2_tgc_ph5p4hs_30c_rd, test_color2)
#plot_point(scaf2_cgt_ph5p4ls_30c_ddg, scaf2_cgt_ph5p4ls_30c_rd, test_color)
#plot_point(scaf2_cgt_ph5p4hs_30c_ddg, scaf2_cgt_ph5p4hs_30c_rd, test_color)
#plot_point(scaf2_tgc_ph5p4hs_40c_ddg, scaf2_tgc_ph5p4hs_40c_rd, test_color2)

#print "SEE", m_avg, b_avg
#print "CGT_5.4LS, RD = ", scaf2_cgt_ph5p4ls_30c_rd, " UV = ", scaf2_cgt_ph5p4ls_30c_ddg, ", Pred = ", (m_net*scaf2_cgt_ph5p4ls_30c_ddg)+b_net 
#print "CGT_5.4HS, RD = ", scaf2_cgt_ph5p4hs_30c_rd, " UV = ", scaf2_cgt_ph5p4hs_30c_ddg, ", Pred = ", (m_net*scaf2_cgt_ph5p4hs_30c_ddg)+b_net
#print "TGC_5.4LS, RD = ", scaf2_tgc_ph5p4ls_30c_rd, " UV = ", scaf2_tgc_ph5p4ls_30c_ddg, ", Pred = ", (m_net*scaf2_tgc_ph5p4ls_30c_ddg)+b_net
#print "TGC_5.4HS 30c, RD = ", scaf2_tgc_ph5p4hs_30c_rd, " UV = ", scaf2_tgc_ph5p4hs_30c_ddg, ", Pred = ", (m_net*scaf2_tgc_ph5p4hs_30c_ddg)+b_net
#print "TGC_5.4HS 40c, RD = ", scaf2_tgc_ph5p4hs_40c_rd, " UV = ", scaf2_tgc_ph5p4hs_40c_ddg, ", Pred = ", (m_net*scaf2_tgc_ph5p4hs_40c_ddg)+b_net
#print "TGC_5.4HS 15c, RD = ", 0.0, " UV = ", scaf2_tgc_ph5p4hs_15c_ddg, ", Pred = ", (m_net*scaf2_tgc_ph5p4hs_15c_ddg)+b_net

plt.ylabel("$\Delta G^{\circ}_{conf}(i)$ (kcal/mol)", fontsize=fs)
plt.xlabel("$\Delta\Delta G^{\circ}_{melt}(i)$ (kcal/mol)", fontsize=fs)
plt.plot([1.65, 6.0], [(1.65*line_coeff[0]) + line_coeff[1], (6*line_coeff[0]) + line_coeff[1]], color='k', linewidth=2, zorder=10, label=eqn)
plt.fill_between(x_dummy, np.amin(line_matrix, axis=0), np.amax(line_matrix, 0), color='#B0E2FF', zorder=0)

plt.xticks([2.0, 3.0, 4.0, 5.0], fontsize=fs)
plt.yticks([2.0, 3.0, 4.0, 5.0], fontsize=fs)
plt.xlim([1.65, 5])
plt.ylim([1.65, 5])
#print "r2", r2
#print 'r2 check = ', r2_score(rd_values, predicted_rd, multioutput='raw_values')
#print 'r2 check = ', r2_score(predicted_rd, rd_values, multioutput='raw_values')
#print 'r2 check2 = ', r2_score(rd_values, uv_values, multioutput='raw_values')
#print 'r2 check2 = ', r2_score(uv_values, rd_values, multioutput='raw_values')

print "rmse", rmse, ' kcal/mol'
print 'r2 = ', pearsonr(rd_values, predicted_rd)
#print pearsonr(predicted_rd, rd_values)
#print np.corrcoef(rd_values, predicted_rd)[0,1]
#print np.corrcoef(predicted_rd, rd_values)[0,1]
#print pearsonr(rd_values, uv_values)
#print pearsonr(uv_values, rd_values)
#print np.corrcoef(rd_values, uv_values)[0,1]
#print np.corrcoef(uv_values, rd_values)[0,1]
r2 = np.corrcoef(uv_values, rd_values)[0,1]
#ax.text(2.0, 4.50, "$R^{2}$ = " + "%3.2f"%r2, fontsize=24)
#ax.text(2.0, 4.25, "$RMSE$ = " + "%3.2f"%rmse, fontsize=24)
#plt.legend(fontsize=18, loc=4)
#plt.show()
plt.savefig('fig5c.pdf')

#fig, ax = plt.subplots()
#plt.errorbar(predicted_rd, rd_values, fmt='o')
#plt.xlabel('predicted')
#plt.ylabel('rd meas')
##plt.xlim([2, 4])
##plt.ylim([2, 4])
#plt.show()
