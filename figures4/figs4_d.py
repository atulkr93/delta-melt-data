# Plot NMR spectra
import matplotlib.pyplot as plt
import numpy as np
import nmrglue as ng
import math
import matplotlib
matplotlib.rcParams['font.family']='sans-serif'
matplotlib.rcParams['font.sans-serif']=['Arial']
matplotlib.rcParams['axes.linewidth'] = 4.0
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
       ax.text(xlabel[dummy]+0.05, ylabel[dummy], labels[dummy], size=fs, color='k')    
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

fig, ax = plt.subplots(1, 6, figsize=(48, 9))
plot_spectrum(ax[0], 'A6DNAulb_6.8_25mM_25C_H2O_C1H1_HSQC_700/', 3.*math.pow(10, 6), 20, 1.4, "H1' (ppm)", "C1' (ppm)", [6.5, 5.0], [6.5, 5.1, 0.5], [88.5, 83.5], [88.0, 84.0, 2.0], 'A6DNA_ulb_HSQC_Ali_6.8_25C.list', "C1'", '%2.1f', '%d', "A" + "$_{6}$" + "-DNA(ulb)\npH 6.8 25 " + "$^{\circ}$" + "C")
plot_spectrum(ax[1], 'A6DNAulb_6.8_25mM_25C_H2O_C8H8_HSQC_700/', 1.*math.pow(10, 6), 20, 1.4, "H6/8 (ppm)", "C6/8 (ppm)", [8.3, 6.9], [8.0, 7.0, 0.5], [144.5, 137.], [144., 137., 3.], 'A6DNA_ulb_HSQC_Aro_6.8_25C.list', "C6", '%2.1f', '%d', "A" + "$_{6}$" + "-DNA(ulb)\npH 6.8 25 " + "$^{\circ}$" + "C")
plot_spectrum(ax[2], 'A6DNAs1lb_5.4_25mM_25C_D2O_C1H1_HSQC_700/', 3.*math.pow(10, 6), 20, 1.4, "H1' (ppm)", "C1' (ppm)", [6.3, 5.1], [6.0, 5.3, 0.7], [87.2, 83.5], [87.0, 84.0, 1.0],  'A6lbDNA_HSQC_C1p.list', "C1'", '%2.1f', '%d', "A" + "$_{6}$" + "-DNA(s1lb)\npH 5.4 (D" + "$_{2}$" + "O) 25 " + "$^{\circ}$" + "C")
plot_spectrum(ax[3], 'A6DNAs1lb_5.4_25mM_25C_D2O_C8H8_HSQC_700/', 3.*math.pow(10, 5), 20, 1.4, "H6/8 (ppm)", "C6/8 (ppm)", [8.3, 6.9], [8.0, 7.0, 0.5], [144.4, 137.8], [144., 138., 3.0], 'A6lbDNA_HSQC_Aro.list', "C6", '%2.1f', '%d', "A" + "$_{6}$" + "-DNA(s1lb)\npH 5.4 (D" + "$_{2}$" + "O) 25 " + "$^{\circ}$" + "C")
plot_spectrum(ax[4], 'A6DNAs2lb_5.4_25mM_25C_D2O_C1H1_HSQC_700/', 3.*math.pow(10, 6), 20, 1.4, "H1' (ppm)", "C1' (ppm)", [6.5, 5.3], [6.4, 5.4, 0.5], [89., 83.], [89., 83., 2.0], 'T6lbDNA_HSQC_C1p.list', "C1'", '%2.1f', '%d', "A" + "$_{6}$" + "-DNA(s2lb)\npH 5.4 (D" + "$_{2}$" + "O) 25 " + "$^{\circ}$" + "C")
plot_spectrum(ax[5], 'A6DNAs2lb_5.4_25mM_25C_D2O_C8H8_HSQC_700/', 5.*math.pow(10, 5), 20, 1.4, "H6/8 (ppm)", "C6/8 (ppm)", [8.5, 6.9], [8.5, 7.0, 0.5], [144., 137.], [144., 138., 2.], 'T6lbDNA_HSQC_Aro.list', "C6", '%2.1f', '%d', "A" + "$_{6}$" + "-DNA(s2lb)\npH 5.4 (D" + "$_{2}$" + "O) 25 " + "$^{\circ}$" + "C")

#plt.show()
plt.tight_layout()
plt.savefig('figs4_d.pdf')
