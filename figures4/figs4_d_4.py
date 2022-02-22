# Plot NMR spectra
import matplotlib.pyplot as plt
import numpy as np
import nmrglue as ng
import math
import matplotlib

matplotlib.rcParams['axes.linewidth'] = 4.0
matplotlib.rcParams['font.family']='Sans-serif'
matplotlib.rcParams['font.sans-serif']=['Arial']

elw = 4.0
tick_length=8.0
fs = 37.5

def read_peaklist(filename, nucleus):
    ''' Read peak list of given nucleus '''
    f = open(filename, 'r')
    # headers
    line = f.readline()
    line = f.readline()
   
    # read data
    x = []
    y = []
    labels = []
    for line in f.readlines():
        line = line.strip('\n').split(' ')
        #print line
        line = [ele for ele in line if ele != '']
        line[1] = float(line[1])
        line[2] = float(line[2])
        line[0] = line[0].split("-")[0][:-1*int(len(nucleus))]
        labels.append(line[0])
        y.append(line[1])
        x.append(line[2])
    f.close()    
    return [x, y, labels]

def plot_spectrum(ax, dirname, contour_start, contour_num, contour_factor, xaxis_label, yaxis_label, xlim, xtick_lims, ylim, ytick_lims, peaklist, nucleus, xfmt, yfmt, plot_title):
    # define contour levels
    cl = contour_start * contour_factor ** np.arange(contour_num)

    # Read data
    dic, data = ng.pipe.read(dirname + 'test.ft2')
    
    # generate axes
    x_proton = ng.pipe.make_uc(dic, data, dim=1)
    x_proton_ppm = x_proton.ppm_scale()
    x_proton_ppm_start, x_proton_ppm_end = x_proton.ppm_limits()
    y_het = ng.pipe.make_uc(dic, data, dim=0)
    y_het_ppm = y_het.ppm_scale()
    y_het_ppm_start, y_het_ppm_end = y_het.ppm_limits()

    # Peak list
    [xlabel, ylabel, labels] = read_peaklist(dirname + peaklist, nucleus)
    for dummy in range(len(labels)):
       ax.text(xlabel[dummy], ylabel[dummy], labels[dummy], size=fs, color='k')    
    ax.errorbar(xlabel, ylabel, fmt='x', markersize=8, color='k')

    # Plotting
    ax.contour(data.transpose(), cl, colors='k', extent=(y_het_ppm_start, y_het_ppm_end, x_proton_ppm_start, x_proton_ppm_end))
    ax.set_xlabel(xaxis_label, fontsize=fs)
    ax.set_ylabel(yaxis_label, fontsize=fs)
    ax.set_ylim([ylim[0], ylim[1]])
    ax.set_xlim([xlim[0], xlim[1]])

    ax.set_xticks(np.arange(xtick_lims[1], xtick_lims[0]+0.01, xtick_lims[2]))
    ax.set_yticks(np.arange(ytick_lims[1], ytick_lims[0]+0.01, ytick_lims[2]))
    ax.set_xticklabels([xfmt%ele for ele in np.arange(xtick_lims[1], xtick_lims[0]+0.01, xtick_lims[2])], fontsize=fs)
    ax.set_yticklabels([yfmt%ele for ele in np.arange(ytick_lims[1], ytick_lims[0]+0.01, ytick_lims[2])], fontsize=fs)

    ax.tick_params(axis='x', which='major', direction='out', top=False, width=elw, length=tick_length, pad=6)
    ax.tick_params(axis='x', which='minor', direction='out', top=False, width=elw, length=0, pad=6)
    ax.tick_params(axis='y', which='major', direction='out', right=False, width=elw, length=tick_length, pad=6)
    ax.tick_params(axis='y', which='minor', direction='out', right=False, width=elw, length=0, pad=6)
  
    title = ax.set_title(plot_title, fontsize=fs) 
    title.set_position([0.5, 1.02])

fig, ax = plt.subplots(1, 6, figsize=(47.4, 8.9))

plot_spectrum(ax[0], 'A2DNAulb_6.8_25mM_25C_H2O_C8H8_HSQC_600/', 3.*math.pow(10, 6), 20, 1.4, "H6/8 (ppm)", "C6/8 (ppm)", [8.4, 7.0], [8.3, 7.2, 0.5], [144., 137.], [144., 137., 3.], 'A2DNAtlb_6.8_HSQC_25C_Aro.list', "C6", '%2.1f', '%d', "A" + "$_{2}$" + "-DNA(ulb)\npH 6.8 25 " + "$^{\circ}$" + "C")
plot_spectrum(ax[1], 'A2DNAs2lb_5.4_25mM_25C_H2O_C1H1_HSQC_700/', math.pow(10, 6), 20, 1.4, "H1' (ppm)", "C1' (ppm)", [6.1, 5.5], [6.1, 5.6, 0.4], [86.0, 83.5], [86.0, 84.0, 1.], 'A2DNA_T2lb_HSQC_Ali_5.4_25C.list', "C1'", '%2.1f', '%d', "A" + "$_{2}$" + "-DNA(s2lb)\npH 5.4 25 " + "$^{\circ}$" + "C")
plot_spectrum(ax[2], 'A2DNAs2lb_5.4_25mM_25C_H2O_C8H8_HSQC_700/', 5*math.pow(10, 5), 20, 1.4, "H6/8 (ppm)", "C6/8 (ppm)", [8.1, 6.9], [8.0, 7.0, 0.5], [140.0, 137], [140., 137., 1.], 'A2DNA_T2lb_HSQC_Aro_5.4_25C.list', "C8", '%2.1f', '%d', "A" + "$_{2}$" + "-DNA(s2lb)\npH 5.4 25 " + "$^{\circ}$" + "C")
plot_spectrum(ax[3], 'Scaf2CGTG6lb_5.4_25mM_30C_H2O_C1H1_HSQC_600/', 1.5*math.pow(10, 6), 20, 1.4, "H1' (ppm)", "C1' (ppm)", [6.2, 5.8], [6.2, 5.9, 0.2], [85.2, 84.0], [85.0, 84.0, 1.], 'scaf2-CGT_low_5.4_HSQC_30C_Ali.list', "C1'", '%2.1f', '%d', "scaf2_CGT(G6lb)" + "\npH 5.4 30 " + "$^{\circ}$" + "C")
plot_spectrum(ax[4], 'Scaf2CGTG6lb_5.4_25mM_30C_H2O_C8H8_HSQC_600/', 1.*math.pow(10, 6), 20, 1.4, "H8 (ppm)", "C8 (ppm)", [8.1, 7.6], [8.0, 7.7, 0.3], [139., 137.5], [139., 137., 1.], 'scaf2-CGT_low_5.4_HSQC_30C_Aro.list', "C6", '%2.1f', '%d', "scaf2_CGT(G6lb)" + "\npH 5.4 30 " + "$^{\circ}$" + "C")
plot_spectrum(ax[5], 'Scaf2TGCG6lb_5.4_25mM_30C_H2O_C1H1_HSQC_600/', 1.*math.pow(10, 6), 20, 1.4, "H1' (ppm)", "C1' (ppm)", [6.1, 5.6], [6.1, 5.7, 0.3], [85.2, 83.2], [85., 83., 1.], 'scaf2-TGC_low_5.4_HSQC_30C.list', "C1'", '%2.1f', '%d', "scaf2_TGC(G6lb)" + "\npH 5.4 30 " + "$^{\circ}$" + "C")



#plt.show()
plt.tight_layout()
plt.savefig('figs4_d_4.pdf')
