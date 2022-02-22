# Plot delta-melt correlation for TT mismatches
from sklearn.metrics import r2_score
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

# Plot params
fs=12
fmt_point = 'o'
msize = 16
elw = 1
cpthick = 1
cpsize = 8
mewval = 1.0
at_color='r'
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
    #print "WT = ", wt_delg, ", Mutant = ", mt_delg

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
    if color_pt == 'unknown':
        plt.errorbar(unumpy.nominal_values(uv), unumpy.nominal_values(rd), yerr=unumpy.std_devs(rd), xerr=unumpy.std_devs(uv), color='green', fmt=fmt_point, markersize=msize, elinewidth=elw, capthick=cpthick, capsize=cpsize, mec='green', mew=mewval, zorder=100., mfc='white')
    else:
        plt.errorbar(unumpy.nominal_values(uv), unumpy.nominal_values(rd), yerr=unumpy.std_devs(rd), xerr=unumpy.std_devs(uv), color=color_pt, fmt=fmt_point, markersize=msize, elinewidth=elw, capthick=cpthick, capsize=cpsize, mec=color_pt, mew=mewval, zorder=50., mfc=color_pt)
    #pass

### Read data for getting correlation curve ###
a6dna_a16_ph4p4hs_ddg = compute_deldelG('a6dna_ph4.4hsnb/final_all/test.csv', 'a6dna_t16t9_ph4.4hsnb/final_all/test.csv', 'delG_25C', 1)
a6dna_a16_ph4p4hs_rd = compute_delG(0.888, 0.043, 25, 1)

a6dna_a16_ph5p4lsnb_17p5c_ddg = compute_deldelG('a6dna_ph5.4nb/final_all/test.csv', 'a6dna_t16t9_ph5.4nb/final_all/test.csv', 'delG_17.5C', 1)
a6dna_a16_ph5p4lsnb_17p5c_rd = compute_delG(0.463, 0.019, 17.5, 1)

a6dna_a16_ph5p4lsnb_20c_ddg = compute_deldelG('a6dna_ph5.4nb/final_all/test.csv', 'a6dna_t16t9_ph5.4nb/final_all/test.csv', 'delG_20.0C', 1)
a6dna_a16_ph5p4lsnb_20c_rd = compute_delG(0.459, 0.030, 20.0, 1)

a6dna_a16_ph5p4lsnb_22p5c_ddg = compute_deldelG('a6dna_ph5.4nb/final_all/test.csv', 'a6dna_t16t9_ph5.4nb/final_all/test.csv', 'delG_22.5C', 1)
a6dna_a16_ph5p4lsnb_22p5c_rd = compute_delG(0.559, 0.047, 22.5, 1)

a6dna_a16_ph5p4lsnb_25c_ddg = compute_deldelG('a6dna_ph5.4nb/final_all/test.csv', 'a6dna_t16t9_ph5.4nb/final_all/test.csv', 'delG_25C', 1)
a6dna_a16_ph5p4lsnb_25c_rd = compute_delG(0.629, 0.081, 25.0, 1)

a6dna_ra16_ph5p4_10c_ddg = compute_deldelG('a6dna_ra16_ph5.4nb/final_all/test.csv', 'a6dna_rt16_ph5.4nb/final_all/test.csv', 'delG_10C', 1)
a6dna_ra16_ph5p4_10c_rd = compute_delG(0.180, 0.008, 10, 1)

a6dna_ra16_ph5p4_12p5c_ddg = compute_deldelG('a6dna_ra16_ph5.4nb/final_all/test.csv', 'a6dna_rt16_ph5.4nb/final_all/test.csv', 'delG_12.5C', 1)
a6dna_ra16_ph5p4_12p5c_rd = compute_delG(0.286, 0.046, 12.5, 1)

a6dna_ra16_ph5p4_15c_ddg = compute_deldelG('a6dna_ra16_ph5.4nb/final_all/test.csv', 'a6dna_rt16_ph5.4nb/final_all/test.csv', 'delG_15C', 1)
a6dna_ra16_ph5p4_15c_rd = compute_delG(0.311, 0.074, 15, 1)
#print a6dna_ra16_ph5p4_10c_ddg, a6dna_ra16_ph5p4_12p5c_ddg, a6dna_ra16_ph5p4_15c_ddg

a6dna_a16_ph5p4hsnb_25c_ddg = compute_deldelG('a6dna_ph5.4nbhs/final_all/test.csv', 'a6dna_t16t9_ph5.4hsnb/final_all/test.csv', 'delG_25C', 1)
a6dna_a16_ph5p4hsnb_25c_rd = compute_delG(0.498, 0.038, 25., 1)

a6dnam1g_a16_ph5p4hsnb_ddg = compute_deldelG('a6dna_m1g10_ph5.4hsnb/final_all/test.csv', 'a6dna_m1g10_t16t9_ph5.4hsnb/final_all/test.csv', 'delG_25C', 1)
a6dnam1g_a16_ph5p4hsnb_rd = compute_delG(2.460, 0.041, 25, 1)

# weak rd
a2dna_a16_ph6p9hb_25c_ddg = compute_deldelG('a2dna_ph6.8nb/final_all/test.csv', 'a2dna_t16t9_ph6.8nb/final_all/test.csv', 'delG_25C',0)
a2dna_a16_ph6p9hb_25c_rd = compute_delG(0.252, 0.086, 25., 0)

a6dna_a16_ph6p8_ddg = compute_deldelG('a6dna_ph6.8nb/final_all/test.csv', 'a6dna_t16t9_ph6.8nb/final_all/test.csv', 'delG_25C', 1)
a6dna_a16_ph6p8_rd = compute_delG(0.441, 0.058, 25, 1)

a6dna_a21_ph6p8_25c_ddg = compute_deldelG('a6dna_ph6.8nb/final_all/test.csv', 'a6dna_t21t4_ph6.8nb/final_all/test.csv', 'delG_25C', 1)
a6dna_a21_ph6p8_25c_rd = compute_delG(0.204, 0.029, 25, 1)

a6dna_a21_ph6p8_30c_ddg = compute_deldelG('a6dna_ph6.8nb/final_all/test.csv', 'a6dna_t21t4_ph6.8nb/final_all/test.csv', 'delG_30C', 1)
a6dna_a21_ph6p8_30c_rd = compute_delG(0.205, 0.019, 30, 1)

a6dna_a21_ph4p4hs_30c_ddg = compute_deldelG('a6dna_ph4.4hsnb/final_all/test.csv', 'a6dna_t21t4_ph4.4hsnb/final_all/test.csv', 'delG_30C', 1)
a6dna_a21_ph4p4hs_30c_rd = compute_delG(0.369, 0.026, 30, 1)
#print a6dna_a21_ph4p4hs_30c_ddg, a6dna_a21_ph4p4hs_30c_rd

print 'Number of UV data points = ', len(uv_values), len(uv_error)
print 'Number of NMR RD data points = ', len(rd_values), len(rd_error)

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

#fig, axarr = plt.subplots(1, 2)
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
#print len(uv_values), len(uv_error)
#print len(rd_values), len(rd_error)
#print ''.join([str(ele) + "," for ele in rd_values])
#print ''.join([str(ele) + "," for ele in predicted_rd])

#print predicted_rd
#for dummy in range(len(rd_values)):
#    print rd_values[dummy]
#print
#print
#for dummy in range(len(rd_values)):
#    print predicted_rd[dummy]


#print ''.join([str(ele) + "," for ele in rd_values])
#print ''.join([str(ele) + "," for ele in predicted_rd])

rmse = np.sqrt(np.sum(np.square(rd_values - predicted_rd))/len(list(rd_values)))
r2 = 1 - (np.sum(np.square(predicted_rd - rd_values)) / np.sum(np.square(rd_values - np.mean(rd_values))))
#for dummy in range(len(uv_values)):
#    print predicted_rd[dummy], rd_values[dummy], uv_values[dummy], -1.0*(predicted_rd[dummy] - rd_values[dummy])

# shaded lines
x_dummy = np.linspace(0.0, 4.0, 100)
x_dummy_mat = np.tile(x_dummy, (mc_iterations, 1))
m_vals = np.array([list(bf['m'].values)])
b_vals = np.array([list(bf['b'].values)])
m_matrix = np.tile(m_vals.transpose(), (1, 100))
b_matrix = np.tile(b_vals.transpose(), (1, 100))
line_matrix = np.multiply(m_matrix, x_dummy_mat) + b_matrix

fig, ax = plt.subplots()
plt.tick_params(axis='x', top=False, length=8.0, width=elw, direction='out', pad=5)
plt.tick_params(axis='y', right=False, length=8.0, width=elw, direction='out', pad=5)
#plt.errorbar(uv_values, rd_values, xerr=uv_error, yerr=rd_error, color='r', fmt='o')

plot_point(a6dna_ra16_ph5p4_15c_ddg, a6dna_ra16_ph5p4_15c_rd, 'magenta')
plot_point(a6dna_ra16_ph5p4_12p5c_ddg, a6dna_ra16_ph5p4_12p5c_rd, 'magenta')
plot_point(a6dna_ra16_ph5p4_10c_ddg, a6dna_ra16_ph5p4_10c_rd, 'magenta')

plot_point(a6dna_a16_ph4p4hs_ddg, a6dna_a16_ph4p4hs_rd, at_color)
plot_point(a6dna_a16_ph5p4lsnb_17p5c_ddg, a6dna_a16_ph5p4lsnb_17p5c_rd, at_color)
plot_point(a6dna_a16_ph5p4lsnb_20c_ddg, a6dna_a16_ph5p4lsnb_20c_rd, at_color)
plot_point(a6dna_a16_ph5p4lsnb_22p5c_ddg, a6dna_a16_ph5p4lsnb_22p5c_rd, at_color)
plot_point(a6dna_a16_ph5p4lsnb_25c_ddg, a6dna_a16_ph5p4lsnb_25c_rd, at_color)
plot_point(a6dna_a16_ph5p4hsnb_25c_ddg, a6dna_a16_ph5p4hsnb_25c_rd, at_color)
plot_point(a6dna_a16_ph6p8_ddg, a6dna_a16_ph6p8_rd, at_color)

plot_point(a6dna_a21_ph6p8_30c_ddg, a6dna_a21_ph6p8_30c_rd, 'blue')
plot_point(a6dna_a21_ph6p8_25c_ddg, a6dna_a21_ph6p8_25c_rd, 'blue')
plot_point(a6dna_a21_ph4p4hs_30c_ddg, a6dna_a21_ph4p4hs_30c_rd, 'blue')

plot_point(a6dnam1g_a16_ph5p4hsnb_ddg, a6dnam1g_a16_ph5p4hsnb_rd, 'brown')

plot_point(a2dna_a16_ph6p9hb_25c_ddg, a2dna_a16_ph6p9hb_25c_rd, 'unknown')

plt.ylabel("$\Delta G^{\circ}_{conf}(i)$ (kcal/mol)", fontsize=fs)
plt.xlabel("$\Delta\Delta G^{\circ}_{melt}(i)$ (kcal/mol)", fontsize=fs)
plt.plot([1.65, 4.0], [(1.65*line_coeff[0]) + line_coeff[1], (4*line_coeff[0]) + line_coeff[1]], color='k', linewidth=2, zorder=10)
plt.fill_between(x_dummy, np.amin(line_matrix, axis=0), np.amax(line_matrix, 0), color='#B0E2FF', zorder=0)

plt.title('TT mismatches')
plt.xticks([2.0, 3.0, 4.0], fontsize=fs)
plt.yticks([2.0, 3.0, 4.0], fontsize=fs)
plt.xlim([1.65, 4])
plt.ylim([1.65, 4])
#print "r2 = ", r2
eqn2 = ""
print "rmse = ", rmse, ' kcal/mol'
#eqn2 = "rmse = " + "%2.1f"%rmse
#print 'r2 check = ', r2_score(rd_values, predicted_rd, multioutput='raw_values')
print 'r2 = ', pearsonr(rd_values, predicted_rd)
#print pearsonr(predicted_rd, rd_values)
#print np.corrcoef(rd_values, predicted_rd)[0,1]
#print np.corrcoef(predicted_rd, rd_values)[0,1]
#print pearsonr(rd_values, uv_values)
#print pearsonr(uv_values, rd_values)
#print np.corrcoef(rd_values, uv_values)[0,1]
#print np.corrcoef(uv_values, rd_values)[0,1]
eqn2 = eqn + "\n" + eqn2 + "\n" + "r=" + "%2.1f"%np.corrcoef(uv_values, rd_values)[0,1]
#plt.text(3.0, 2.0, eqn2, fontsize=18)

#ax.text(3.0, 2.25, "$R^{2}$ = " + "%3.2f"%r2, fontsize=24)
#ax.text(3.0, 2.0, "$RMSE$ = " + "%3.2f"%rmse, fontsize=24)
plt.legend(fontsize=18, loc=4)
#plt.show()
plt.savefig('fig4d.pdf')

#fig, ax = plt.subplots()
#plt.errorbar(predicted_rd, rd_values, fmt='o')
#plt.xlabel('predicted')
#plt.ylabel('rd meas')
#plt.xlim([2, 4])
#plt.ylim([2, 4])
#plt.show()
