# delta-melt correlation plot for A-T Opening
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
fs = 12
fmt_point = 'o'
msize = 16
elw = 1
cpthick = 1
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
    deldelg = mt_delg - wt_delg
    if add_status == 1:
        uv_values.append(float(unumpy.nominal_values(deldelg)))
        uv_error.append(float(unumpy.std_devs(deldelg)))
    return deldelg

def compute_delG(pop, pop_error, t, add_status):
    ''' given population, compute delG '''
    global rd_values, rd_error
    pop_mag = ufloat(pop, pop_error) 
    delg = (-8.314 * (273.16 + t) / 4184) * unumpy.log(pop_mag/(1-pop_mag))
    if add_status == 1:
        rd_values.append(float(unumpy.nominal_values(delg)))
        rd_error.append(float(unumpy.std_devs(delg)))
    return delg

g_label_status = 0
k_label_status = 0

def plot_point(uv, rd, color_pt):
    ''' Plot point on correlation plot '''
    global g_label_status 
    global k_label_status 
    if color_pt == 'magenta':
        plt.errorbar(unumpy.nominal_values(uv), unumpy.nominal_values(rd), yerr=unumpy.std_devs(rd), xerr=unumpy.std_devs(uv), color=color_pt, fmt=fmt_point, markersize=msize, elinewidth=elw, capthick=cpthick, capsize=cpsize, mec=color_pt, mew=mewval, zorder=50)
    else:
        plt.errorbar(unumpy.nominal_values(uv), unumpy.nominal_values(rd), yerr=unumpy.std_devs(rd), xerr=unumpy.std_devs(uv), color=color_pt, fmt=fmt_point, markersize=msize, elinewidth=elw, capthick=cpthick, capsize=cpsize, mec=color_pt, mew=mewval, zorder=100)
    #if color_pt == 'g' and g_label_status == 0:
    #    plt.errorbar(unumpy.nominal_values(uv), unumpy.nominal_values(rd), yerr=unumpy.std_devs(rd), xerr=unumpy.std_devs(uv), color=color_pt, fmt=fmt_point, markersize=msize, elinewidth=elw, capthick=cpthick, capsize=cpsize, mec=color_pt, mew=mewval, label='TBP')
    #    g_label_status = 1
    #else:
    #    plt.errorbar(unumpy.nominal_values(uv), unumpy.nominal_values(rd), yerr=unumpy.std_devs(rd), xerr=unumpy.std_devs(uv), color=color_pt, fmt=fmt_point, markersize=msize, elinewidth=elw, capthick=cpthick, capsize=cpsize, mec=color_pt, mew=mewval)

    #if color_pt == 'k' and k_label_status == 0:
    #    plt.errorbar(unumpy.nominal_values(uv), unumpy.nominal_values(rd), yerr=unumpy.std_devs(rd), xerr=unumpy.std_devs(uv), color=color_pt, fmt=fmt_point, markersize=msize, elinewidth=elw, capthick=cpthick, capsize=cpsize, mec=color_pt, mew=mewval, label='$A_{6}$')
    #    k_label_status = 1
    #else:
    #    plt.errorbar(unumpy.nominal_values(uv), unumpy.nominal_values(rd), yerr=unumpy.std_devs(rd), xerr=unumpy.std_devs(uv), color=color_pt, fmt=fmt_point, markersize=msize, elinewidth=elw, capthick=cpthick, capsize=cpsize, mec=color_pt, mew=mewval)
    #pass

### Read data for getting correlation curve ###
# New Error estimation for A6DNA based on MC sampling
a6dna_t5_ph8p8_ddg = compute_deldelG('opening/A6DNA/wt/test.csv', 'opening/A6DNA/m3T5/test.csv', 'delG_25C', 1)
#a6dna_t5_ph8p8_rd = compute_delG(8.02005333753e-07, 4.61041853596e-08, 25, 1)
#print a6dna_t5_ph8p8_rd, "SEEE"
a6dna_t5_ph8p8_rd = ufloat(8.31, 0.04)
rd_values.append(float(unumpy.nominal_values(a6dna_t5_ph8p8_rd)))
rd_error.append(float(unumpy.std_devs(a6dna_t5_ph8p8_rd)))

a6dna_t7_ph8p8_ddg = compute_deldelG('opening/A6DNA/wt/test.csv', 'opening/A6DNA/m3T7/test.csv', 'delG_25C', 1)
#a6dna_t7_ph8p8_rd = compute_delG(1.0618557324e-06, 5.63652186104e-08, 25, 1)
#print a6dna_t7_ph8p8_rd, "SEE"
a6dna_t7_ph8p8_rd = ufloat(8.10, 0.04) 
rd_values.append(float(unumpy.nominal_values(a6dna_t7_ph8p8_rd)))
rd_error.append(float(unumpy.std_devs(a6dna_t7_ph8p8_rd)))

a6dna_t6_ph8p8_ddg = compute_deldelG('cary_uv/a6dna_baseopen/final_all/test.csv', 'cary_uv/a6dna_m3t6_open/final_all/test.csv', 'delG_25C', 1)
a6dna_t6_ph8p8_rd = ufloat(8.34, 0.044) 
rd_values.append(float(unumpy.nominal_values(a6dna_t6_ph8p8_rd)))
rd_error.append(float(unumpy.std_devs(a6dna_t6_ph8p8_rd)))
#print uv_values, uv_error
#print rd_values, rd_error

a6dna_t22_ph8p8_ddg = compute_deldelG('opening/A6DNA/wt/test.csv', 'opening/A6DNA/m3T22/test.csv', 'delG_25C', 1)
#a6dna_t22_ph8p8_ddg = ufloat(4.02, 0.375)
#uv_values.append(4.02)
#uv_error.append(0.375)
#a6dna_t22_ph8p8_rd = compute_delG(7.56305619483e-06, 1.97076341976e-07, 25, 1)
#print a6dna_t22_ph8p8_rd, "SEE"
a6dna_t22_ph8p8_rd = ufloat(7.00, 0.02) 
rd_values.append(float(unumpy.nominal_values(a6dna_t22_ph8p8_rd)))
rd_error.append(float(unumpy.std_devs(a6dna_t22_ph8p8_rd)))

tbp_t3_ph8p0_ddg = compute_deldelG('opening/CDNA/wt/test.csv', 'opening/CDNA/m3T3/test.csv', 'delG_15C', 1)
tbp_t3_ph8p0_rd = compute_delG(2.42089017297e-05, 1.52035399771e-06, 15, 1)
#print tbp_t3_ph8p0_ddg, tbp_t3_ph8p0_rd

tbp_t5_ph8p0_ddg = compute_deldelG('opening/CDNA/wt/test.csv', 'opening/CDNA/m3T5/test.csv', 'delG_15C', 1)
tbp_t5_ph8p0_rd = compute_delG(4.41447119379e-06, 3.68567868016e-07, 15, 1)
#print tbp_t5_ph8p0_ddg, tbp_t5_ph8p0_rd

tbp_t16_ph8p0_ddg = compute_deldelG('opening/CDNA/wt/test.csv', 'opening/CDNA/m3T16/test.csv', 'delG_15C', 1)
tbp_t16_ph8p0_rd = compute_delG(2.89415323649e-06, 9.77986402355e-07, 15, 1)
#print tbp_t16_ph8p0_ddg, tbp_t16_ph8p0_rd

tbp_t17_ph8p0_ddg = compute_deldelG('opening/CDNA/wt/test.csv', 'opening/CDNA/m3T17/test.csv', 'delG_15C', 1)
tbp_t17_ph8p0_rd = compute_delG(1.72597687727e-06, 4.08798431407e-07, 15, 1)
#print tbp_t17_ph8p0_ddg, tbp_t17_ph8p0_rd

tbp_t18_ph8p0_ddg = compute_deldelG('opening/CDNA/wt/test.csv', 'opening/CDNA/m3T18/test.csv', 'delG_15C', 1)
tbp_t18_ph8p0_rd = compute_delG(2.14318271985e-06, 4.71379307794e-07, 15, 1)
#print tbp_t18_ph8p0_ddg, tbp_t18_ph8p0_rd

tbp_t19_ph8p0_ddg = compute_deldelG('opening/CDNA/wt/test.csv', 'opening/CDNA/m3T19/test.csv', 'delG_15C', 1)
tbp_t19_ph8p0_rd = compute_delG(5.00203608876e-06,1.16473568427e-06, 15, 1)
#print tbp_t19_ph8p0_ddg, tbp_t19_ph8p0_rd

tbp_t21_ph8p0_ddg = compute_deldelG('opening/CDNA/wt/test.csv', 'opening/CDNA/m3T21/test.csv', 'delG_15C', 1)
tbp_t21_ph8p0_rd = compute_delG(6.35050497281e-06, 2.05179477338e-06, 15, 1)
#print tbp_t21_ph8p0_ddg, tbp_t21_ph8p0_rd


# convert to arrays
rd_values = np.array(rd_values)
rd_error = np.array(rd_error)
uv_values = np.array(uv_values)
uv_error = np.array(uv_error)
print 'Number of UV data points = ', len(uv_values), len(uv_error)
print 'Number of NMR RD data points = ', len(rd_values), len(rd_error)

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

filename = 'results.csv'
f = open(filename, 'w')
f.write('m,b\n')

for dummy in range(mc_iterations):
    #popt, pcov = curve_fit(linear_fit, uv_error_matrix[dummy, :], rd_error_matrix[dummy, :], p0=[1, 1])
    #line_coeff = popt 
    popt, pcov = curve_fit(linear_fit2, uv_error_matrix[dummy, :], rd_error_matrix[dummy, :], p0=[1])
    line_coeff = list([1]) + list(popt)
    #plt.subplots()
    #plt.plot(uv_error_matrix[dummy, :], rd_error_matrix[dummy, :], color='k', linewidth=0.0, marker='o')
    #plt.errorbar(uv_values, rd_values, xerr=uv_error, yerr=rd_error, color='r', fmt='o')
    #plt.show()

    #eqn = 'y = ' + '%3.2f'%(line_coeff[0]) + 'x + ' + '%3.2f'%(line_coeff[1])
    #predicted_rd = np.array([((line_coeff[0]*ele)+line_coeff[1]) for ele in uv_values])
    #rmse = np.sqrt(np.sum(np.square(rd_values - predicted_rd))/len(list(rd_values)))
    #r2 = 1 - (np.sum(np.square(predicted_rd - rd_values)) / np.sum(np.square(rd_values - np.mean(rd_values))))

    #print "iteration", dummy + 1
    #print "***", eqn, "***"
    #print "see rmse", rmse
    #print "see r2", r2
    #f.write(str(line_coeff[0]) + ',' + str(line_coeff[1]) + "," + str(r2) + "," + str(rmse) + "\n")
    f.write(str(line_coeff[0]) + ',' + str(line_coeff[1]) + "\n")
    
f.close()

# Read the output file and plot
bf = pd.read_csv(filename)

m_avg = np.mean(np.array(bf['m'].values))
m_std = np.std(np.array(bf['m'].values))
b_avg = np.mean(np.array(bf['b'].values))
b_std = np.std(np.array(bf['b'].values))
#rmse_avg = np.mean(np.array(bf['rmse'].values))
#rmse_std = np.std(np.array(bf['rmse'].values))
#r2_avg = np.mean(np.array(bf['r2'].values))
#r2_std = np.std(np.array(bf['r2'].values))
#print "M = ", m_avg, "+/-", m_std
print "c(i) = ", b_avg, "+/-", b_std
#print "rmse = ", rmse_avg, "+/-", rmse_std
#print "r2 = ", r2_avg, "+/-", r2_std

fig, axarr = plt.subplots(2)
axarr[0].set_title('m')
axarr[0].hist(np.array(bf['m'].values), 50, range=(np.amin(np.array(bf['m'].values)), np.amax(np.array(bf['m'].values))))
axarr[1].set_title('b')
axarr[1].hist(np.array(bf['b'].values), 50, range=(np.amin(np.array(bf['b'].values)), np.amax(np.array(bf['b'].values))))
#axarr[1, 0].set_title('rmse')
#axarr[1, 0].hist(np.array(bf['rmse'].values), 50, range=(np.amin(np.array(bf['rmse'].values)), np.amax(np.array(bf['rmse'].values))))
#axarr[1, 1].set_title('r2')
#axarr[1, 1].hist(np.array(bf['r2'].values), 50, range=(np.amin(np.array(bf['r2'].values)), np.amax(np.array(bf['r2'].values))))
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

# shaded lines
x_dummy = np.linspace(2.5, 6.5, 100)
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

plt.ylabel("$\Delta G^{\circ}_{conf}(i)$ (kcal/mol)", fontsize=fs)
plt.xlabel("$\Delta\Delta G^{\circ}_{melt}(i)$ (kcal/mol)", fontsize=fs)
plt.fill_between(x_dummy, np.amin(line_matrix, axis=0), np.amax(line_matrix, 0), color='#B0E2FF', zorder=0)
plt.plot([2.5, 6.5], [(2.5*line_coeff[0]) + line_coeff[1], (6.5*line_coeff[0]) + line_coeff[1]], color='k', linewidth=2, zorder=10)

#a6_color = '#FFA500'
#tbp_color = '#FFA500'
a6_color = 'blue'
tbp_color = 'magenta'

plot_point(a6dna_t6_ph8p8_ddg, a6dna_t6_ph8p8_rd, a6_color)
plot_point(a6dna_t7_ph8p8_ddg, a6dna_t7_ph8p8_rd, a6_color)
plot_point(a6dna_t5_ph8p8_ddg, a6dna_t5_ph8p8_rd, a6_color)
plot_point(a6dna_t22_ph8p8_ddg, a6dna_t22_ph8p8_rd, a6_color)

plot_point(tbp_t3_ph8p0_ddg, tbp_t3_ph8p0_rd, tbp_color)
plot_point(tbp_t5_ph8p0_ddg, tbp_t5_ph8p0_rd, tbp_color)
plot_point(tbp_t16_ph8p0_ddg, tbp_t16_ph8p0_rd, tbp_color)
plot_point(tbp_t17_ph8p0_ddg, tbp_t17_ph8p0_rd, tbp_color)
plot_point(tbp_t18_ph8p0_ddg, tbp_t18_ph8p0_rd, tbp_color)
plot_point(tbp_t19_ph8p0_ddg, tbp_t19_ph8p0_rd, tbp_color)
plot_point(tbp_t21_ph8p0_ddg, tbp_t21_ph8p0_rd, tbp_color)

plt.yticks([5.0, 6.0, 7.0, 8.0], fontsize=fs)
plt.xticks([3.0, 4.0, 5.0, 6.0], fontsize=fs)
#print "R2 = ", r2
#print 'r2 check = ', r2_score(rd_values, predicted_rd, multioutput='raw_values')
#print 'r2 check = ', r2_score(predicted_rd, rd_values, multioutput='raw_values')
#print 'r2 check2 = ', r2_score(rd_values, uv_values, multioutput='raw_values')
#print 'r2 check2 = ', r2_score(uv_values, rd_values, multioutput='raw_values')

print "rmse =", rmse
print 'r2 = ', pearsonr(rd_values, predicted_rd)[0]
#print pearsonr(predicted_rd, rd_values)
#print np.corrcoef(rd_values, predicted_rd)[0,1]
#print np.corrcoef(predicted_rd, rd_values)[0,1]
#print pearsonr(rd_values, uv_values)
#print pearsonr(uv_values, rd_values)
#print np.corrcoef(rd_values, uv_values)[0,1]
#print np.corrcoef(uv_values, rd_values)[0,1]

plt.xlim([2.5, 6.5])
plt.ylim([5.0, 8.5])

#ax.text(4.9, 6.25, "$R^{2}$ = " + "%3.2f"%r2, fontsize=fs)
#ax.text(4.9, 6.0, "$RMSE$ = " + "%3.2f"%rmse, fontsize=fs)
#plt.legend(fontsize=18, loc=4)
plt.savefig('fig6c.pdf')
#plt.show()

#fig, ax = plt.subplots()
#plt.errorbar(predicted_rd, rd_values, fmt='o')
#plt.xlabel('predicted')
#plt.ylabel('rd meas')
##plt.xlim([2, 4])
##plt.ylim([2, 4])
#plt.show()
