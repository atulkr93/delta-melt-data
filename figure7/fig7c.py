# delta-melt correlation plot for m6A-U bps
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
fs=12
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
    pop_mag = ufloat(pop, pop_error) / 100.0
    delg = (-8.314 * (273.16 + t) / 4184) * unumpy.log(pop_mag/(1-pop_mag))
    if add_status == 1:
        rd_values.append(float(unumpy.nominal_values(delg)))
        rd_error.append(float(unumpy.std_devs(delg)))
    return delg

def plot_point(uv, rd, color_pt):
    ''' Plot point on correlation plot '''
    plt.errorbar(unumpy.nominal_values(uv), unumpy.nominal_values(rd), yerr=unumpy.std_devs(rd), xerr=unumpy.std_devs(uv), color=color_pt, fmt=fmt_point, markersize=msize, elinewidth=elw, capthick=cpthick, capsize=cpsize, mec=color_pt, mew=mewval, zorder=50)
    #pass

### Read data for getting correlation curve ###
a6rna_m6a16_ph6p8_ddg = compute_deldelG('a6rna_m6a_ph6.8bb/test.csv', 'a6rna_dm6a_ph6.8bb/test.csv', 'delG_37C', 1)
a6rna_m6a16_ph6p8_rd = compute_delG(1.203, 0.580, 37, 1)
#print 'a6 37c ddg/rd =', a6rna_m6a16_ph6p8_ddg, a6rna_m6a16_ph6p8_rd
 
### NOTE THAT THIS IS MC ERROR (N=100)
gbc_m6a6_ph6p8_37c_ddg = compute_deldelG('dsgbc_m6a_ph6.8bb/test.csv', 'dsgbc_dm6a_ph6.8bb/test.csv', 'delG_37C', 1)
gbc_m6a6_ph6p8_37c_rd = compute_delG(0.643, 0.006, 37, 1)
#print 'gbc 37c ddg/rd = ', gbc_m6a6_ph6p8_37c_ddg, gbc_m6a6_ph6p8_37c_rd

gbc_m6a6_ph6p8_55c_ddg = compute_deldelG('dsgbc_m6a_ph6.8bb/test.csv', 'dsgbc_dm6a_ph6.8bb/test.csv', 'delG_55C', 1)
gbc_m6a6_ph6p8_55c_rd = compute_delG(1.205, 0.042, 55, 1)
#print 'gbc 55c ddg/rd = ', gbc_m6a6_ph6p8_55c_ddg, gbc_m6a6_ph6p8_55c_rd

gbc_m6a6_ph6p8_65c_ddg = compute_deldelG('dsgbc_m6a_ph6.8bb/test.csv', 'dsgbc_dm6a_ph6.8bb/test.csv', 'delG_65C', 1)
gbc_m6a6_ph6p8_65c_rd = compute_delG(1.526, 0.013, 65, 1)
#print 'gbc 65c ddg/rd = ', gbc_m6a6_ph6p8_65c_ddg, gbc_m6a6_ph6p8_65c_rd



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
#gbc_color = '#FFA500'
#a6_color = '#FFA500'
gbc_color = 'blue'
a6_color = 'magenta'

plot_point(gbc_m6a6_ph6p8_65c_ddg, gbc_m6a6_ph6p8_65c_rd, gbc_color)
plot_point(gbc_m6a6_ph6p8_55c_ddg, gbc_m6a6_ph6p8_55c_rd, gbc_color)
plot_point(gbc_m6a6_ph6p8_37c_ddg, gbc_m6a6_ph6p8_37c_rd, gbc_color)
plot_point(a6rna_m6a16_ph6p8_ddg, a6rna_m6a16_ph6p8_rd, a6_color)
plt.ylabel("$\Delta G^{\circ}_{conf}(i)$ (kcal/mol)", fontsize=fs)
plt.xlabel("$\Delta\Delta G^{\circ}_{melt}(i)$ (kcal/mol)", fontsize=fs)
plt.fill_between(x_dummy, np.amin(line_matrix, axis=0), np.amax(line_matrix, 0), color='#B0E2FF', zorder=0)
plt.plot([1.5, 4.0], [(1.5*line_coeff[0]) + line_coeff[1], (4*line_coeff[0]) + line_coeff[1]], color='k', label=eqn, linewidth=2, zorder=10)

plt.xticks([2.0, 3.0], fontsize=fs)
plt.yticks([2.0, 3.0], fontsize=fs)
plt.xlim([1.5, 3.0])
plt.ylim([2., 3.5])
#print "r2 = ", r2
#print 'r2 check = ', r2_score(rd_values, predicted_rd, multioutput='raw_values')
#print 'r2 check = ', r2_score(predicted_rd, rd_values, multioutput='raw_values')
#print 'r2 check2 = ', r2_score(rd_values, uv_values, multioutput='raw_values')
#print 'r2 check2 = ', r2_score(uv_values, rd_values, multioutput='raw_values')

print "rmse = ", rmse, ' kcal/mol'
print 'r2 = ', pearsonr(rd_values, predicted_rd)[0]
#print pearsonr(predicted_rd, rd_values)
#print np.corrcoef(rd_values, predicted_rd)[0,1]
#print np.corrcoef(predicted_rd, rd_values)[0,1]
#print pearsonr(rd_values, uv_values)
#print pearsonr(uv_values, rd_values)
#print np.corrcoef(rd_values, uv_values)[0,1]
#print np.corrcoef(uv_values, rd_values)[0,1]

#ax.text(2.6, 2.50, "$R^{2}$ = " + "%3.2f"%r2, fontsize=24)
#ax.text(2.6, 2.25, "$RMSE$ = " + "%3.2f"%rmse, fontsize=24)
#plt.legend(fontsize=fs, loc=4)
plt.savefig('fig7c.pdf')
#plt.show()


#fig, ax = plt.subplots()
#plt.errorbar(predicted_rd, rd_values, fmt='o')
#plt.xlabel('predicted')
#plt.ylabel('rd meas')
##plt.xlim([2, 4])
##plt.ylim([2, 4])
#plt.show()
