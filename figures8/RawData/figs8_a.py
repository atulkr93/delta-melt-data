import math
import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
from uncertainties import unumpy
import pandas as pd
import matplotlib
from scipy.optimize import curve_fit

matplotlib.rcParams['axes.linewidth']=4.0
matplotlib.rcParams['font.family']='Sans-serif'
matplotlib.rcParams['font.sans-serif']=['Arial']

elw = 4.0
fs = 30
tick_length=8.0


def fun_ex_ir(x,kex,R1n):
    [x_time, E, R1w] = x
    return (1 - (E*kex/(R1w-R1n)*(np.exp(-R1n*x_time)-np.exp(-R1w*x_time))))


def fun_ex_ir2(x,R1n):
    [x_time, E, R1w, kex_val] = x
    return (1 - (E*kex_val/(R1w-R1n)*(np.exp(-R1n*x_time)-np.exp(-R1w*x_time))))


def plot_curve(ax, dirname, data, R1w, E, xlim, xticks, xfmt, ylim, yticks, yfmt, plot_title, add_status, base_conc):
    global xpts, ypts
    global kex_array

    # Read input data
    data = pd.read_csv(dirname + data)
    delays = np.array(data['delay(s)'].values)
    intensity = np.array(data['vol_rel'])

    # Do the fitting
    popt, pcov = curve_fit(fun_ex_ir,[delays, E, R1w],intensity, p0 = [1., 1.], absolute_sigma=False)
    kex = ufloat(popt[0],np.sqrt(np.diag(pcov))[0])
    R1n = ufloat(popt[1],np.sqrt(np.diag(pcov))[1])
    print plot_title
    print 'kex = ', kex
    print 'R1n = ', R1n
    print

    # get optimal trendline
    time = np.linspace(0., xlim[1], 500)
    opt_trendline =  (1 - (E*popt[0]/(R1w-popt[1])*(np.exp(-popt[1]*time)-np.exp(-R1w*time))))

    # Plot
    ax.plot(delays, intensity, linewidth=0.0, markersize=8., color='k', marker='o')
    ax.plot(time, opt_trendline, linewidth=elw, color='r')

    # tick params
    ax.tick_params(axis='x', which='minor', direction='out', pad=6, width=elw, length=0, top=False)
    ax.tick_params(axis='x', which='major', direction='out', pad=6, width=elw, length=tick_length, top=False)
    ax.tick_params(axis='y', which='minor', direction='out', pad=6, width=elw, length=0, right=False)
    ax.tick_params(axis='y', which='major', direction='out', pad=6, width=elw, length=tick_length, right=False)
     
    # Labels
    ax.set_xlim(xlim)
    ax.set_xticks([ele for ele in np.arange(xticks[0], xticks[1]+0.00001, xticks[2])]) 
    ax.set_xticklabels([xfmt%ele for ele in np.arange(xticks[0], xticks[1]+0.00001, xticks[2])], fontsize=fs) 
    ax.set_ylim(ylim)
    ax.set_yticks([ele for ele in np.arange(yticks[0], yticks[1], yticks[2])]) 
    ax.set_yticklabels([yfmt%ele for ele in np.arange(yticks[0], yticks[1], yticks[2])], fontsize=fs) 
    title = ax.set_title(plot_title, fontsize=fs) 
    title.set_position([0.5, 1.02])

    # Perform degeneracy calculations
    #fig2, axarr2 = plt.subplots(4, 5)
    #horiz_counter = 0
    #vert_counter = 0
    kex_array.sort()
    residuals_net = []
    for dummy in kex_array:
        popt, pcov = curve_fit(fun_ex_ir2, [delays, E, R1w, dummy],intensity, p0 = [10.], absolute_sigma=False, bounds=[0.0, math.pow(10, 8)])
        predicted_vals =  (1 - (E*dummy/(R1w-popt[0])*(np.exp(-popt[0]*delays)-np.exp(-R1w*delays))))
        residual = np.sum(np.square(predicted_vals - intensity))
        residuals_net.append(residual)
        #title = 'kex=' + '%.2f'%dummy + ',R1=' + '%.2f'%popt[0] + ', rss=', '%.2f'%residual
        #time2 = np.linspace(0., xlim[1], 500)
        #fix_trendline =  (1 - (E*dummy/(R1w-popt[0])*(np.exp(-popt[0]*time2)-np.exp(-R1w*time2))))
        #axarr2[horiz_counter, vert_counter].set_title(title, fontsize=10) 
        #axarr2[horiz_counter, vert_counter].set_xticks([])
        #axarr2[horiz_counter, vert_counter].plot(delays, intensity, linewidth=0.0, markersize=8., color='k', marker='o')
        #axarr2[horiz_counter, vert_counter].plot(time2, fix_trendline, linewidth=elw, color='r')

        #vert_counter = vert_counter + 1
        #if vert_counter == 5:
        #    vert_counter = 0
        #    horiz_counter = horiz_counter + 1
    #plt.show()
    residuals_net = np.array(residuals_net)
    residuals_net = residuals_net / np.amin(residuals_net)
    #ax2.plot(kex_array, residuals_net, linewidth=elw, color='k')
    #ax2.plot([unumpy.nominal_values(kex), unumpy.nominal_values(kex)], [0., math.pow(10, 4)], linewidth=elw, linestyle="--", color='k')
    #ax2.set_yscale('log')
    #ax2.set_xscale('log')
    #ax2.set_ylim([0.0, math.pow(10, 4)])
    #ax2.set_yticks([math.pow(10, ele) for ele in range(5)])
    #ax2.set_yticklabels(['$10^{%d}$'%int(ele) for ele in range(5)], fontsize=fs)
    #ax2.set_xticks([math.pow(10, ele) for ele in range(int(math.ceil(math.log10(np.amin(kex_array)))), int(math.floor(math.log10(np.amax(kex_array))))+1)])
    #ax2.set_xticklabels(['$10^{%d}$'%int(ele) for ele in range(int(math.ceil(math.log10(np.amin(kex_array)))), int(math.floor(math.log10(np.amax(kex_array))))+1)], fontsize=fs)
    #ax2.set_xlim([np.amin(kex_array), np.amax(kex_array)])

    ## tick params
    #ax2.tick_params(axis='x', which='minor', direction='out', pad=6, width=elw, length=0, top=False)
    #ax2.tick_params(axis='x', which='major', direction='out', pad=6, width=elw, length=tick_length, top=False)
    #ax2.tick_params(axis='y', which='minor', direction='out', pad=6, width=elw, length=0, right=False)
    #ax2.tick_params(axis='y', which='major', direction='out', pad=6, width=elw, length=tick_length, right=False)
    #title = ax2.set_title(plot_title, fontsize=fs) 
    #title.set_position([0.5, 1.01])

    if add_status == 1:
        xpts.append(1/(0.001*base_conc))
        ypts.append(1./kex)

def linear_fit(x, m, b):
    return ((m*x)+b)

def plot_fit(ax, x, y):
    # MC fitting
    mc_iterations = 1000

    # Initialize error matrices
    y_error_matrix = np.zeros((mc_iterations, len(y)))
    yvals = unumpy.nominal_values(y)
    yerrors = unumpy.std_devs(y)

    for dummy in range(len(yvals)):
        op = (np.ones(mc_iterations)*yvals[dummy]) + np.random.normal(0.0, yerrors[dummy], mc_iterations) 
        y_error_matrix[:, dummy] = op.transpose()
    
    m_vals = np.zeros(mc_iterations)
    b_vals = np.zeros(mc_iterations)
    rmse_vals = np.zeros(mc_iterations)
    r2_vals = np.zeros(mc_iterations)

    for dummy in range(mc_iterations):
        popt, pcov = curve_fit(linear_fit, x, y_error_matrix[dummy, :], p0=[1., 1.])
        m_vals[dummy] = popt[0]
        b_vals[dummy] = popt[1]
 
        predicted_vals = np.array([((popt[0]*ele)+popt[1]) for ele in x]) 
        rmse = np.sqrt(np.sum(np.square(yvals - predicted_vals))/len(list(yvals))) 
        r2 = 1 - (np.sum(np.square(predicted_vals - yvals)) / np.sum(np.square(yvals - np.mean(yvals))))
        rmse_vals[dummy] = rmse
        r2_vals[dummy] = r2
     
    fig, axarr = plt.subplots(2, 2)

    axarr[0, 0].set_title('m')
    axarr[0, 0].hist(m_vals, 50, range=(np.amin(m_vals), np.amax(m_vals)))
    axarr[0, 1].set_title('b')
    axarr[0, 1].hist(b_vals, 50, range=(np.amin(b_vals), np.amax(b_vals)))
    axarr[1, 0].set_title('r2')
    axarr[1, 0].hist(r2_vals, 50, range=(np.amin(r2_vals), np.amax(r2_vals)))
    axarr[1, 1].set_title('rmse')
    axarr[1, 1].hist(rmse_vals, 50, range=(np.amin(rmse_vals), np.amax(rmse_vals)))

    m_avg = np.mean(m_vals)
    m_std = np.std(m_vals)
    b_avg = np.mean(b_vals)
    b_std = np.std(b_vals)
    kcoll = 1/(4*math.pow(10, -9))
    m_net = ufloat(m_avg, m_std)
    K = (1/m_net) * kcoll
    delG = -(8.314 * 298.16 / 4184) * unumpy.log(K)
    tau0 = ufloat(b_avg, b_std)
    print "Eqbm constant K = ", K
    print "delG = ", delG
    print "Tau0 = ", tau0

fig, ax = plt.subplots(4, 5, figsize=(32, 24))
#fig, ax = plt.subplots(2, 5)

xpts = []
ypts = []


kex_array = list(np.logspace(-2., 1, 10.)) + [0.26] 
plot_curve(ax[0, 0], 'A6DNA_8.8_25C_0mM_Resi_5/', 'data_A6DNA_8.8_25C_0mM_Resi_5.csv', 0.353183667, 1.900162, [-0.2, 3.7], [0, 3., 1], '%d', [0.73, 1.02], [0.8, 1.01, 0.1], '%2.1f', 'A' + '$_{6}$' + '-DNA T5-H3 0mM NH' + '$_{3}$', 0, 0)
kex_array = list(np.logspace(-1., 2, 10.)) + [3.] 
plot_curve(ax[0, 1], 'A6DNA_8.8_25C_20mM_Resi_5/', 'data_A6DNA_8.8_25C_20mM_Resi_5.csv', 0.334388287, 1.853923, [-0.2, 3.2], [0, 3., 1], '%d', [-0.1, 1.1], [0.0, 1.01, 0.2], '%2.1f', 'A' + '$_{6}$' + '-DNA T5-H3 20mM NH' + '$_{3}$', 1, 20)
kex_array = list(np.logspace(-1., 2, 10.)) + [4.6] 
plot_curve(ax[0, 2], 'A6DNA_8.8_25C_40mM_Resi_5/', 'data_A6DNA_8.8_25C_40mM_Resi_5.csv', 0.329595522, 1.782285, [-0.2, 3.2], [0, 3., 1], '%d', [-0.2, 1.1], [0.0, 1.01, 0.2], '%2.1f', 'A' + '$_{6}$' + '-DNA T5-H3 40mM NH' + '$_{3}$', 1, 40)
kex_array = list(np.logspace(-1., 3, 10.)) + [7.2] 
plot_curve(ax[0, 3], 'A6DNA_8.8_25C_100mM_Resi_5/', 'data_A6DNA_8.8_25C_100mM_Resi_5.csv', 0.331505686, 1.951047, [-0.05, 1.1], [0, 1., 0.5], '%2.1f', [-0.4, 1.1], [0.0, 1.1, 0.5], '%2.1f', 'A' + '$_{6}$' + '-DNA T5-H3 100mM NH' + '$_{3}$', 1, 100)
kex_array = list(np.logspace(-1., 3, 10.)) + [8.2] 
plot_curve(ax[0, 4], 'A6DNA_8.8_25C_150mM_Resi_5/', 'data_A6DNA_8.8_25C_150mM_Resi_5.csv', 0.320788387, 1.929775, [-0.05, 1.1], [0, 1., 0.5], '%2.1f', [-0.4, 1.1], [0.0, 1.1, 0.5], '%2.1f', 'A' + '$_{6}$' + '-DNA T5-H3 150mM NH' + '$_{3}$', 1, 150)


kex_array = list(np.logspace(-2., 1.7, 10.)) + [0.47] 
plot_curve(ax[2, 0], 'A6DNA_8.8_25C_0mM_Resi_7/', 'data_A6DNA_8.8_25C_0mM_Resi_7.csv', 0.353183667, 1.900162, [-0.2, 3.7], [0, 3., 1], '%d', [0.65, 1.02], [0.7, 1.01, 0.1], '%2.1f', 'A' + '$_{6}$' + '-DNA T7-H3 0mM NH' + '$_{3}$', 0, 0)
kex_array = list(np.logspace(-1., 2.7, 10.)) + [4.53] 
plot_curve(ax[2, 1], 'A6DNA_8.8_25C_20mM_Resi_7/', 'data_A6DNA_8.8_25C_20mM_Resi_7.csv', 0.334388287, 1.853923, [-0.2, 3.2], [0, 3., 1], '%d', [-0.2, 1.05], [0.0, 1.01, 0.5], '%2.1f', 'A' + '$_{6}$' + '-DNA T7-H3 20mM NH' + '$_{3}$', 1, 20)
kex_array = list(np.logspace(-1., 2.85, 10.)) + [6.8] 
plot_curve(ax[2, 2], 'A6DNA_8.8_25C_40mM_Resi_7/', 'data_A6DNA_8.8_25C_40mM_Resi_7.csv', 0.329595522, 1.782285, [-0.2, 3.2], [0, 3., 1], '%d', [-0.3, 1.05], [0.0, 1.01, 0.5], '%2.1f', 'A' + '$_{6}$' + '-DNA T7-H3 40mM NH' + '$_{3}$', 1, 40)
kex_array = list(np.logspace(-1., 3., 10.)) + [11.3] 
plot_curve(ax[2, 3], 'A6DNA_8.8_25C_100mM_Resi_7/', 'data_A6DNA_8.8_25C_100mM_Resi_7.csv', 0.331505686, 1.951047, [-0.05, 1.05], [0, 1., 0.5], '%2.1f', [-0.5, 1.05], [-0.5, 1.01, 0.5], '%2.1f', 'A' + '$_{6}$' + '-DNA T7-H3 100mM NH' + '$_{3}$', 1, 100)
kex_array = list(np.logspace(-1., 3., 10.)) + [14.6] 
plot_curve(ax[2, 4], 'A6DNA_8.8_25C_150mM_Resi_7/', 'data_A6DNA_8.8_25C_150mM_Resi_7.csv', 0.320788387, 1.929775, [-0.05, 1.05], [0, 1., 0.5], '%2.1f', [-0.6, 1.1], [-0.5, 1.01, 0.5], '%2.1f', 'A' + '$_{6}$' + '-DNA T7-H3 150mM NH' + '$_{3}$', 1, 150)

kex_array = list(np.logspace(-1., 2., 10.)) + [1.8] 
plot_curve(ax[3, 0], 'A6DNA_8.8_25C_0mM_Resi_22/', 'data_A6DNA_8.8_25C_0mM_Resi_22.csv', 0.353183667, 1.900162, [-0.2, 3.7], [0, 3., 1], '%d', [0.1, 1.05], [0.2, 1.01, 0.2], '%2.1f', 'A' + '$_{6}$' + '-DNA T22-H3 0mM NH' + '$_{3}$', 0, 0)
kex_array = list(np.logspace(0., 3., 10.)) + [31.7] 
plot_curve(ax[3, 1], 'A6DNA_8.8_25C_20mM_Resi_22/', 'data_A6DNA_8.8_25C_20mM_Resi_22.csv', 0.334388287, 1.853923, [-0.2, 3.2], [0, 3., 1], '%d', [-0.8, 1.1], [-0.5, 1.01, 0.5], '%2.1f', 'A' + '$_{6}$' + '-DNA T22-H3 20mM NH' + '$_{3}$', 1, 20)
kex_array = list(np.logspace(0., 3., 10.)) + [50.9] 
plot_curve(ax[3, 2], 'A6DNA_8.8_25C_40mM_Resi_22/', 'data_A6DNA_8.8_25C_40mM_Resi_22.csv', 0.329595522, 1.782285, [-0.2, 3.2], [0, 3., 1], '%d', [-0.9, 1.1], [-0.5, 1.01, 0.5], '%2.1f', 'A' + '$_{6}$' + '-DNA T22-H3 40mM NH' + '$_{3}$', 1, 40)
kex_array = list(np.logspace(1., 3.7, 10.)) + [98.78] 
plot_curve(ax[3, 3], 'A6DNA_8.8_25C_100mM_Resi_22/', 'data_A6DNA_8.8_25C_100mM_Resi_22.csv', 0.331505686, 1.951047, [-0.05, 1.05], [0, 1., 0.5], '%d', [-1.3, 1.1], [-1.0, 1.01, 0.5], '%2.1f', 'A' + '$_{6}$' + '-DNA T22-H3 100mM NH' + '$_{3}$', 1, 100)
kex_array = list(np.logspace(1., 4., 10.)) + [112.38] 
plot_curve(ax[3, 4], 'A6DNA_8.8_25C_150mM_Resi_22/', 'data_A6DNA_8.8_25C_150mM_Resi_22.csv', 0.320788387, 1.929775, [-0.05, 1.05], [0, 1., 0.5], '%d', [-1.3, 1.1], [-1.0, 1.01, 0.5], '%2.1f', 'A' + '$_{6}$' + '-DNA T22-H3 150mM NH' + '$_{3}$', 1, 150)

plot_curve(ax[1, 0], 'A6DNA_8.8_25C_0mM_Resi_6/', 'data_A6DNA_8.8_25C_0mM_Resi_6.csv', 0.353183667, 1.900162, [-0.2, 3.7], [0, 3., 1], '%d', [0.77, 1.02], [0.8, 1.01, 0.1], '%2.1f', 'A' + '$_{6}$' + '-DNA T6-H3 0mM NH' + '$_{3}$', 0, 0)
plot_curve(ax[1, 1], 'A6DNA_8.8_25C_20mM_Resi_6/', 'data_A6DNA_8.8_25C_20mM_Resi_6.csv', 0.334388287, 1.853923, [-0.2, 3.2], [0, 3., 1], '%d', [0., 1.05], [0., 1.01, 0.2], '%2.1f', 'A' + '$_{6}$' + '-DNA T6-H3 20mM NH' + '$_{3}$', 0, 20)
plot_curve(ax[1, 2], 'A6DNA_8.8_25C_40mM_Resi_6/', 'data_A6DNA_8.8_25C_40mM_Resi_6.csv', 0.329595522, 1.782285, [-0.2, 3.2], [0, 3., 1], '%d', [-0.1, 1.05], [0., 1.01, 0.2], '%2.1f', 'A' + '$_{6}$' + '-DNA T6-H3 40mM NH' + '$_{3}$', 0, 40)
plot_curve(ax[1, 3], 'A6DNA_8.8_25C_100mM_Resi_6/', 'data_A6DNA_8.8_25C_100mM_Resi_6.csv', 0.331505686, 1.951047, [-0.1, 1.1], [0, 1.0, 0.5], '%2.1f', [-0.3, 1.1], [0., 1.01, 0.5], '%2.1f', 'A' + '$_{6}$' + '-DNA T6-H3 100mM NH' + '$_{3}$', 0, 100)
plot_curve(ax[1, 4], 'A6DNA_8.8_25C_150mM_Resi_6/', 'data_A6DNA_8.8_25C_150mM_Resi_6.csv', 0.320788387, 1.929775, [-0.1, 1.1], [0, 1.0, 0.5], '%2.1f', [-0.3, 1.1], [0., 1.01, 0.5], '%2.1f', 'A' + '$_{6}$' + '-DNA T6-H3 150mM NH' + '$_{3}$', 0, 150)

for dummy in range(5):
    ax[3, dummy].set_xlabel('Time (s)', fontsize=fs)
#for dummy in range(5):
#    ax[7, dummy].set_xlabel('k' + "$_{ex}$" + ' (s' + '$^{-1}$' + ')', fontsize=fs)

for dummy in range(4):
    ax[dummy, 0].set_ylabel('Norm. Volume', fontsize=fs)
#for dummy in range(4, 8):
#    ax[dummy, 0].set_ylabel('Norm. RSS', fontsize=fs)
#plt.show()

plt.tight_layout()
plt.savefig('figs8_a.pdf')
