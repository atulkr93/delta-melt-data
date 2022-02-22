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


def plot_curve(ax, dirname, data, R1w, E, xlim, xticks, xfmt, ylim, yticks, yfmt, plot_title, add_status, base_conc):
    global xpts, ypts

    # Read input data
    data = pd.read_csv(dirname + data)
    delays = np.array(data['delay(s)'].values)
    intensity = np.array(data['vol_rel'])

    # Do the fitting
    popt, pcov = curve_fit(fun_ex_ir,[delays, E, R1w],intensity, p0 = [1., 1.], absolute_sigma=False)
    kex = ufloat(popt[0],np.sqrt(np.diag(pcov))[0])
    R1n = ufloat(popt[1],np.sqrt(np.diag(pcov))[1])
    print plot_title
    #print 'kex = ', kex
    #print 'R1n = ', R1n
    #print

    # get optimal trendline
    time = np.linspace(0., xlim[1], 500)
    opt_trendline =  (1 - (E*popt[0]/(R1w-popt[1])*(np.exp(-popt[1]*time)-np.exp(-R1w*time))))

    ## Plot
    #ax.plot(delays, intensity, linewidth=0.0, markersize=8., color='k', marker='o')
    #ax.plot(time, opt_trendline, linewidth=elw, color='r')

    ## tick params
    #ax.tick_params(axis='x', which='minor', direction='out', pad=6, width=elw, length=0, top=False)
    #ax.tick_params(axis='x', which='major', direction='out', pad=6, width=elw, length=tick_length, top=False)
    #ax.tick_params(axis='y', which='minor', direction='out', pad=6, width=elw, length=0, right=False)
    #ax.tick_params(axis='y', which='major', direction='out', pad=6, width=elw, length=tick_length, right=False)
    # 
    ## Labels
    #ax.set_xlim(xlim)
    #ax.set_xticks([ele for ele in np.arange(xticks[0], xticks[1]+0.00001, xticks[2])]) 
    #ax.set_xticklabels([xfmt%ele for ele in np.arange(xticks[0], xticks[1]+0.00001, xticks[2])], fontsize=fs) 
    #ax.set_ylim(ylim)
    #ax.set_yticks([ele for ele in np.arange(yticks[0], yticks[1], yticks[2])]) 
    #ax.set_yticklabels([yfmt%ele for ele in np.arange(yticks[0], yticks[1], yticks[2])], fontsize=fs) 
    #title = ax.set_title(plot_title, fontsize=fs) 
    #title.set_position([0.5, 1.01])

    if add_status == 1:
        xpts.append(1/(0.001*base_conc))
        ypts.append(1./kex)

def linear_fit(x, m, b):
    return ((m*x)+b)

def plot_fit(ax, x, y, xmaster, xticks, xfmt, ylim, yticks, yfactor, yfmt, plot_title):
    # MC fitting
    mc_iterations = 100

    # Initialize error matrices
    y_error_matrix = np.zeros((mc_iterations, len(y)))
    yvals = unumpy.nominal_values(y)
    yerrors = unumpy.std_devs(y)
    #print y

    for dummy in range(len(yvals)):
        op = (np.ones(mc_iterations)*yvals[dummy]) + np.random.normal(0.0, yerrors[dummy], mc_iterations) 
        #op = (np.ones(mc_iterations)*yvals[dummy]) 
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
        rmse = np.sqrt(np.sum(np.square(y_error_matrix[dummy, :] - predicted_vals))/len(list(y_error_matrix[dummy, :]))) 
        r2 = 1 - (np.sum(np.square(predicted_vals - y_error_matrix[dummy, :])) / np.sum(np.square(y_error_matrix[dummy, :] - np.mean(y_error_matrix[dummy, :]))))
        rmse_vals[dummy] = rmse
        r2_vals[dummy] = r2
     
    #fig, axarr = plt.subplots(2, 2)

    #axarr[0, 0].set_title('m')
    #axarr[0, 0].hist(m_vals, 50, range=(np.amin(m_vals), np.amax(m_vals)))
    #axarr[0, 1].set_title('b')
    #axarr[0, 1].hist(b_vals, 50, range=(np.amin(b_vals), np.amax(b_vals)))
    #axarr[1, 0].set_title('r2')
    #axarr[1, 0].hist(r2_vals, 50, range=(np.amin(r2_vals), np.amax(r2_vals)))
    #axarr[1, 1].set_title('rmse')
    #axarr[1, 1].hist(rmse_vals, 50, range=(np.amin(rmse_vals), np.amax(rmse_vals)))

    #fig, axarr = plt.subplots()
    #H, xedges, yedges = np.histogram2d(m_vals, b_vals, 20)
    #X, Y = np.meshgrid(xedges, yedges)
    #axarr.pcolormesh(X, Y, H)
    #axarr.set_xlabel('m')
    #axarr.set_ylabel('b')

    m_avg = np.mean(m_vals)
    m_std = np.std(m_vals)
    b_avg = np.mean(b_vals)
    b_std = np.std(b_vals)
    kcoll = 1/(4*math.pow(10, -9))
    m_net = ufloat(m_avg, m_std)
    tau0 = ufloat(b_avg, b_std)
    K = (1/(m_net * kcoll))
    pes = K/(1+K)
    delG = -(8.314 * 298.16 / 4184) * unumpy.log(K)
    delGn = -(8.314 * 298.16 / 4184) * unumpy.log(pes/(1-pes))
    print m_avg, b_avg, "***"
    print "Eqbm constant K = ", unumpy.nominal_values(K), unumpy.std_devs(K)
    print "delG = ", unumpy.nominal_values(delG), unumpy.std_devs(delG)
    print "Tau0 = ", unumpy.nominal_values(tau0), unumpy.std_devs(tau0)
    print "Pes = ", unumpy.nominal_values(pes), unumpy.std_devs(pes)
    print "delGn = ", delGn
    print 1/tau0  
    ax.errorbar(x, yvals, yerr=yerrors, color='k', fmt='o', markersize=5, elinewidth=elw, capsize=elw, capthick=elw, mec='k', mfc='k', mew=2)
    ax.plot([xmaster[0], xmaster[-1]], [(xmaster[0]*m_avg)+b_avg, (xmaster[-1]*m_avg)+b_avg], color='r', linewidth=elw)
    # tick params
    ax.tick_params(axis='x', which='minor', direction='out', pad=6, width=elw, length=0, top=False)
    ax.tick_params(axis='x', which='major', direction='out', pad=6, width=elw, length=tick_length, top=False)
    ax.tick_params(axis='y', which='minor', direction='out', pad=6, width=elw, length=0, right=False)
    ax.tick_params(axis='y', which='major', direction='out', pad=6, width=elw, length=tick_length, right=False)

    predicted_vals = np.array([((m_avg*ele)+b_avg) for ele in x]) 
    rmse = np.sqrt(np.sum(np.square(yvals - predicted_vals))/len(list(yvals))) 
    r2 = 1 - (np.sum(np.square(predicted_vals - yvals)) / np.sum(np.square(yvals - np.mean(yvals))))
   
    m_vals = np.array([list(m_vals)])
    b_vals = np.array([list(b_vals)]) 
    x_dummy = np.linspace(xmaster[0], xmaster[1], 1000)
    x_dummy_mat = np.tile(x_dummy, (mc_iterations, 1))
    m_matrix = np.tile(m_vals.transpose(), (1, 1000))
    b_matrix = np.tile(b_vals.transpose(), (1, 1000))
    line_matrix = np.multiply(m_matrix, x_dummy_mat) + b_matrix
    ax.fill_between(x_dummy, np.amin(line_matrix, axis=0), np.amax(line_matrix, 0), color='cyan')
    ax.set_xlim(xmaster)
    ax.set_xticks([ele for ele in np.arange(xticks[0], xticks[1], xticks[2])])
    ax.set_xticklabels([xfmt%ele for ele in np.arange(xticks[0], xticks[1], xticks[2])], fontsize=fs)
    ax.set_xlabel('1/[NH' + '$_{3}$' + '] (M' + '$^{-1}$' + ')', fontsize=fs)
    ax.set_ylim(ylim)
    ax.set_yticks([ele for ele in np.arange(yticks[0], yticks[1], yticks[2])]) 
    ax.set_yticklabels([yfmt%(ele*yfactor) for ele in np.arange(yticks[0], yticks[1], yticks[2])], fontsize=fs) 
    title = ax.set_title(plot_title, fontsize=fs) 
    title.set_position([0.5, 1.02])
 
    return delG

def draw_bar(val1, val2, ax):
    # tick params
    ax.tick_params(axis='x', which='minor', direction='out', pad=6, width=elw, length=0, top=False)
    ax.tick_params(axis='x', which='major', direction='out', pad=6, width=elw, length=tick_length, top=False)
    ax.tick_params(axis='y', which='minor', direction='out', pad=6, width=elw, length=0, right=False)
    ax.tick_params(axis='y', which='major', direction='out', pad=6, width=elw, length=tick_length, right=False)
    ax.bar([0.0], [unumpy.nominal_values(val1)], yerr=unumpy.std_devs(val1), align='center', width=1.0, label='This study', linewidth=elw, error_kw={'elinewidth':elw, 'capthick':4}, ecolor='k', capsize=elw*6, color='k')
    ax.bar([1.0], [unumpy.nominal_values(val1)], yerr=unumpy.std_devs(val1), align='center', width=1.0, label='Gueron et. al.', linewidth=elw, error_kw={'elinewidth':elw, 'capthick':4}, ecolor='k', capsize=elw*6, color='r')
    ax.set_xlim([-1.0, 2.0])
    ax.set_xticks([])
    ax.set_ylabel('$\Delta$' + 'G' + '$^{\circ}$' + '$_{conf,25 ^{\circ} C}$' + '\n(kcal/mol)', fontsize=fs)
    ax.set_ylim([0., 11.])
    ax.set_yticks([ele for ele in np.arange(0., 10.1, 4.)])
    ax.set_yticklabels([int(ele) for ele in np.arange(0., 10.1, 4.)], fontsize=fs)
    leg=ax.legend(fontsize=fs)
    leg.get_frame().set_linewidth(elw)
    title = ax.set_title('GDNA T6-H3', fontsize=fs)
    title.set_position([0.5, 1.01])

fig, ax = plt.subplots(1, 5, figsize=(32, 6.6))
#fig, ax = plt.subplots(2, 5)

xpts = []
ypts = []
plot_curve(ax[0], 'A6DNA_8.8_25C_0mM_Resi_5/', 'data_A6DNA_8.8_25C_0mM_Resi_5.csv', 0.353183667, 1.900162, [-0.2, 3.7], [0, 3., 1], '%d', [0.73, 1.02], [0.8, 1.01, 0.1], '%2.1f', 'A' + '$_{6}$' + '-DNA T5-H3 0mM NH' + '$_{3}$', 0, 0)
plot_curve(ax[0], 'A6DNA_8.8_25C_20mM_Resi_5/', 'data_A6DNA_8.8_25C_20mM_Resi_5.csv', 0.334388287, 1.853923, [-0.2, 3.2], [0, 3., 1], '%d', [-0.1, 1.1], [0.0, 1.01, 0.2], '%2.1f', 'A' + '$_{6}$' + '-DNA T5-H3 20mM NH' + '$_{3}$', 1, 20)
plot_curve(ax[0], 'A6DNA_8.8_25C_40mM_Resi_5/', 'data_A6DNA_8.8_25C_40mM_Resi_5.csv', 0.329595522, 1.782285, [-0.2, 3.2], [0, 3., 1], '%d', [-0.2, 1.1], [0.0, 1.01, 0.2], '%2.1f', 'A' + '$_{6}$' + '-DNA T5-H3 40mM NH' + '$_{3}$', 1, 40)
plot_curve(ax[0], 'A6DNA_8.8_25C_100mM_Resi_5/', 'data_A6DNA_8.8_25C_100mM_Resi_5.csv', 0.331505686, 1.951047, [-0.05, 1.1], [0, 1., 0.5], '%2.1f', [-0.4, 1.1], [0.0, 1.1, 0.5], '%2.1f', 'A' + '$_{6}$' + '-DNA T5-H3 100mM NH' + '$_{3}$', 1, 100)
plot_curve(ax[0], 'A6DNA_8.8_25C_150mM_Resi_5/', 'data_A6DNA_8.8_25C_150mM_Resi_5.csv', 0.320788387, 1.929775, [-0.05, 1.1], [0, 1., 0.5], '%2.1f', [-0.4, 1.1], [0.0, 1.1, 0.5], '%2.1f', 'A' + '$_{6}$' + '-DNA T5-H3 150mM NH' + '$_{3}$', 1, 150)
a6dna_t5_dg = plot_fit(ax[0], xpts, ypts, [0., 55.], [0., 50., 20.], '%d', [0.0, 0.45], [0.0, 0.45, 0.10], 1000., '%d', 'A' + '$_{6}$' + '-DNA T5-H3')
ax[0].set_ylabel('(1/k' + '$_{ex}$' + ')*10' + '$^{3}$' + ' (s' + '$^{-1}$' + ')', fontsize=fs)


xpts = []
ypts = []
plot_curve(ax[2], 'A6DNA_8.8_25C_0mM_Resi_7/', 'data_A6DNA_8.8_25C_0mM_Resi_7.csv', 0.353183667, 1.900162, [-0.2, 3.7], [0, 3., 1], '%d', [0.65, 1.02], [0.7, 1.01, 0.1], '%2.1f', 'A' + '$_{6}$' + '-DNA T7-H3 0mM NH' + '$_{3}$', 0, 0)
plot_curve(ax[2], 'A6DNA_8.8_25C_20mM_Resi_7/', 'data_A6DNA_8.8_25C_20mM_Resi_7.csv', 0.334388287, 1.853923, [-0.2, 3.2], [0, 3., 1], '%d', [-0.2, 1.05], [0.0, 1.01, 0.5], '%2.1f', 'A' + '$_{6}$' + '-DNA T7-H3 20mM NH' + '$_{3}$', 1, 20)
plot_curve(ax[2], 'A6DNA_8.8_25C_40mM_Resi_7/', 'data_A6DNA_8.8_25C_40mM_Resi_7.csv', 0.329595522, 1.782285, [-0.2, 3.2], [0, 3., 1], '%d', [-0.3, 1.05], [0.0, 1.01, 0.5], '%2.1f', 'A' + '$_{6}$' + '-DNA T7-H3 40mM NH' + '$_{3}$', 1, 40)
plot_curve(ax[2], 'A6DNA_8.8_25C_100mM_Resi_7/', 'data_A6DNA_8.8_25C_100mM_Resi_7.csv', 0.331505686, 1.951047, [-0.05, 1.05], [0, 1., 0.5], '%2.1f', [-0.5, 1.05], [-0.5, 1.01, 0.5], '%2.1f', 'A' + '$_{6}$' + '-DNA T7-H3 100mM NH' + '$_{3}$', 1, 100)
plot_curve(ax[2], 'A6DNA_8.8_25C_150mM_Resi_7/', 'data_A6DNA_8.8_25C_150mM_Resi_7.csv', 0.320788387, 1.929775, [-0.05, 1.05], [0, 1., 0.5], '%2.1f', [-0.6, 1.1], [-0.5, 1.01, 0.5], '%2.1f', 'A' + '$_{6}$' + '-DNA T7-H3 150mM NH' + '$_{3}$', 1, 150)
a6dna_t7_dg = plot_fit(ax[2], xpts, ypts, [0., 55.], [0., 50., 20.], '%d', [0.0, 0.3], [0.0, 0.301, 0.10], 1000., '%d', 'A' + '$_{6}$' + '-DNA T7-H3')

xpts = []
ypts = []
plot_curve(ax[3], 'A6DNA_8.8_25C_0mM_Resi_22/', 'data_A6DNA_8.8_25C_0mM_Resi_22.csv', 0.353183667, 1.900162, [-0.2, 3.7], [0, 3., 1], '%d', [0.1, 1.05], [0.2, 1.01, 0.2], '%2.1f', 'A' + '$_{6}$' + '-DNA T22-H3 0mM NH' + '$_{3}$', 0, 0)
plot_curve(ax[3], 'A6DNA_8.8_25C_20mM_Resi_22/', 'data_A6DNA_8.8_25C_20mM_Resi_22.csv', 0.334388287, 1.853923, [-0.2, 3.2], [0, 3., 1], '%d', [-0.8, 1.1], [-0.5, 1.01, 0.5], '%2.1f', 'A' + '$_{6}$' + '-DNA T22-H3 20mM NH' + '$_{3}$', 1, 20)
plot_curve(ax[3], 'A6DNA_8.8_25C_40mM_Resi_22/', 'data_A6DNA_8.8_25C_40mM_Resi_22.csv', 0.329595522, 1.782285, [-0.2, 3.2], [0, 3., 1], '%d', [-0.9, 1.1], [-0.5, 1.01, 0.5], '%2.1f', 'A' + '$_{6}$' + '-DNA T22-H3 40mM NH' + '$_{3}$', 1, 40)
plot_curve(ax[3], 'A6DNA_8.8_25C_100mM_Resi_22/', 'data_A6DNA_8.8_25C_100mM_Resi_22.csv', 0.331505686, 1.951047, [-0.05, 1.05], [0, 1., 0.5], '%d', [-1.3, 1.1], [-1.0, 1.01, 0.5], '%2.1f', 'A' + '$_{6}$' + '-DNA T22-H3 100mM NH' + '$_{3}$', 1, 100)
plot_curve(ax[3], 'A6DNA_8.8_25C_150mM_Resi_22/', 'data_A6DNA_8.8_25C_150mM_Resi_22.csv', 0.320788387, 1.929775, [-0.05, 1.05], [0, 1., 0.5], '%d', [-1.3, 1.1], [-1.0, 1.01, 0.5], '%2.1f', 'A' + '$_{6}$' + '-DNA T22-H3 150mM NH' + '$_{3}$', 1, 150)
a6dna_t22_dg = plot_fit(ax[3], xpts, ypts, [0., 55.], [0., 50., 20.], '%d', [0.0, 0.04], [0.0, 0.041, 0.01], 1000., '%d', 'A' + '$_{6}$' + '-DNA T22-H3')

xpts = []
ypts = []
plot_curve(ax[1], 'A6DNA_8.8_25C_0mM_Resi_6/', 'data_A6DNA_8.8_25C_0mM_Resi_6.csv', 0.353183667, 1.900162, [-0.2, 3.7], [0, 3., 1], '%d', [0.77, 1.02], [0.8, 1.01, 0.1], '%2.1f', 'A' + '$_{6}$' + '-DNA T6-H3 0mM NH' + '$_{3}$', 0, 0)
plot_curve(ax[1], 'A6DNA_8.8_25C_20mM_Resi_6/', 'data_A6DNA_8.8_25C_20mM_Resi_6.csv', 0.334388287, 1.853923, [-0.2, 3.2], [0, 3., 1], '%d', [0., 1.05], [0., 1.01, 0.2], '%2.1f', 'A' + '$_{6}$' + '-DNA T6-H3 20mM NH' + '$_{3}$', 1, 20)
plot_curve(ax[1], 'A6DNA_8.8_25C_40mM_Resi_6/', 'data_A6DNA_8.8_25C_40mM_Resi_6.csv', 0.329595522, 1.782285, [-0.2, 3.2], [0, 3., 1], '%d', [-0.1, 1.05], [0., 1.01, 0.2], '%2.1f', 'A' + '$_{6}$' + '-DNA T6-H3 40mM NH' + '$_{3}$', 1, 40)
plot_curve(ax[1], 'A6DNA_8.8_25C_100mM_Resi_6/', 'data_A6DNA_8.8_25C_100mM_Resi_6.csv', 0.331505686, 1.951047, [-0.1, 1.1], [0, 1.0, 0.5], '%2.1f', [-0.3, 1.1], [0., 1.01, 0.5], '%2.1f', 'A' + '$_{6}$' + '-DNA T6-H3 100mM NH' + '$_{3}$', 1, 100)
plot_curve(ax[1], 'A6DNA_8.8_25C_150mM_Resi_6/', 'data_A6DNA_8.8_25C_150mM_Resi_6.csv', 0.320788387, 1.929775, [-0.1, 1.1], [0, 1.0, 0.5], '%2.1f', [-0.3, 1.1], [0., 1.01, 0.5], '%2.1f', 'A' + '$_{6}$' + '-DNA T6-H3 150mM NH' + '$_{3}$', 1, 150)
#a6dna_t6_dg = plot_fit(ax[4], xpts, ypts, [0., 55.], [0., 50., 20.], '%d', [0.0, 0.04], [0.0, 0.041, 0.01], 1000., '%d', 'A' + '$_{6}$' + '-DNA T6-H3')
a6dna_t6_dg = plot_fit(ax[1], xpts, ypts, [0., 55.], [0., 50., 20.], '%d', [0.0, 0.45], [0.0, 0.401, 0.10], 1000., '%d', 'A' + '$_{6}$' + '-DNA T6-H3')

#for dummy in range(5):
#    ax[3, dummy].set_xlabel('Time (s)', fontsize=fs)
#for dummy in range(5):
#    ax[7, dummy].set_xlabel('k' + "$_{ex}$" + ' (s' + '$^{-1}$' + ')', fontsize=fs)

#for dummy in range(4):
#    ax[dummy, 0].set_ylabel('Norm. Volume', fontsize=fs)
#for dummy in range(4, 8):
#    ax[dummy, 0].set_ylabel('Norm. RSS', fontsize=fs)
#plt.show()
plt.tight_layout()
plt.savefig('figs8_c.pdf')

