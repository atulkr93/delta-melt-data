# Plot NMR spectra
import matplotlib.pyplot as plt
import numpy as np
import nmrglue as ng
import math
import matplotlib
import matplotlib.gridspec as gridspec

matplotlib.rcParams['axes.linewidth'] = 4.0
elw = 4.0
tick_length=8.0
fs = 24

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
        if len(nucleus) != 0:
            line[0] = line[0].split("-")[0][:-1*int(len(nucleus))]
        labels.append(line[0])
        y.append(line[1])
        x.append(line[2])
    f.close()    
    return [x, y, labels]

def plot_spectrum(ax, dirname, contour_start, contour_num, contour_factor, xaxis_label, yaxis_label, xlim, xtick_lims, ylim, ytick_lims, peaklist, nucleus, xfmt, yfmt, plot_title, color_plot, filename, yticklabel_status):
    # define contour levels
    cl = contour_start * contour_factor ** np.arange(contour_num)

    # Read data
    dic, data = ng.pipe.read(dirname + filename)
    
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
       ax.text(xlabel[dummy], ylabel[dummy], labels[dummy], size=fs, color=color_plot)    
    ax.errorbar(xlabel, ylabel, fmt='x', markersize=8, color=color_plot)

    # Plotting
    ax.contour(data.transpose(), cl, colors=color_plot, extent=(y_het_ppm_start, y_het_ppm_end, x_proton_ppm_start, x_proton_ppm_end))
    ax.set_xlabel(xaxis_label, fontsize=fs)
    ax.set_ylim([ylim[0], ylim[1]])
    ax.set_xlim([xlim[0], xlim[1]])
    ax.set_xticks(np.arange(xtick_lims[1], xtick_lims[0]+0.01, xtick_lims[2]))

    if yaxis_label[:3] == "C1'":
        ax.set_yticks(np.arange(ylim[1], ylim[0]+0.01, ytick_lims[2]))
    else:
        ax.set_yticks(np.arange(ytick_lims[1], ytick_lims[0]+0.01, ytick_lims[2]))
    ax.set_xticklabels([xfmt%ele for ele in np.arange(xtick_lims[1], xtick_lims[0]+0.01, xtick_lims[2])], fontsize=fs)
    if yticklabel_status == 1:
        ax.set_yticklabels([yfmt%ele for ele in np.arange(ytick_lims[1], ytick_lims[0]+0.01, ytick_lims[2])], fontsize=fs)
        ax.set_ylabel(yaxis_label, fontsize=fs)
    else:
        ax.set_yticklabels([])

    ax.tick_params(axis='x', which='major', direction='out', top=False, width=elw, length=tick_length, pad=6)
    ax.tick_params(axis='x', which='minor', direction='out', top=False, width=elw, length=0, pad=6)
    ax.tick_params(axis='y', which='major', direction='out', right=False, width=elw, length=tick_length, pad=6)
    ax.tick_params(axis='y', which='minor', direction='out', right=False, width=elw, length=0, pad=6)
  
    title = ax.set_title(plot_title, fontsize=fs) 
    title.set_position([0.5, 1.02])

#fig, ax = plt.subplots(figsize=(18, 18))
fig = plt.figure(figsize=(40, 20))
gs = gridspec.GridSpec(ncols=27, nrows=12)
#ax_aro = fig.add_subplot(gs[:9,9:18])
#ax_c2 = fig.add_subplot(gs[9:,9:18])
#ax_c1p = fig.add_subplot(gs[:6, 18:])
#ax_noe1 = fig.add_subplot(gs[6:9, 18:21])
#ax_noe2 = fig.add_subplot(gs[6:9, 21:24])
#ax_noe3 = fig.add_subplot(gs[6:9, 24:])
#ax_noe4 = fig.add_subplot(gs[9:, 18:21])
#ax_noe5 = fig.add_subplot(gs[9:, 21:24])
#ax_noe6 = fig.add_subplot(gs[9:, 24:])
#ax_imino1 = fig.add_subplot(gs[6:9, :9])
ax_imino2 = fig.add_subplot(gs[9:, :9])

#plot_spectrum(ax_aro, 'A6DNA-na_6.8_25mM_25C_H2O_13C_HSQC_600/', 6.*math.pow(10, 6), 20, 1.4, "H6/8 (ppm)", "C6/8 (ppm)", [8.4, 6.8], [8.4, 7.0, 0.5], [145, 137], [144., 138.0, 2.0], 'A6DNA_adia_HSQC_Aro.list', "C6", '%2.1f', '%d', "", 'k', 'test.ft2', 1)
#plot_spectrum(ax_c2, 'A6DNA-na_6.8_25mM_25C_H2O_13C_HSQC_600/', 6.*math.pow(10, 6), 20, 1.4, "H2 (ppm)", "C2 (ppm)", [8.1, 6.7], [8.0, 7.0, 0.5], [155.5, 152.5], [155., 153.0, 1.0], 'A6DNA_adia_HSQC_C2.list', "C2", '%2.1f', '%d', "", 'k', 'test.ft2', 1)
#plot_spectrum(ax_aro, 'A6DNAm3T9-na_6.8_25mM_25C_H2O_13C_HSQC_700/', 4.*math.pow(10, 6), 20, 1.4, "H6/8 (ppm)", "C6/8 (ppm)", [8.4, 6.8], [8.4, 7.0, 0.5], [145, 137], [144., 138.0, 2.0], 'A6_m3T9_adia_HSQC_Aro.list', "C6", '%2.1f', '%d', "", 'r', 'test.ft2', 1)
#plot_spectrum(ax_c2, 'A6DNAm3T9-na_6.8_25mM_25C_H2O_13C_HSQC_700/', 4.*math.pow(10, 6), 20, 1.4, "H2 (ppm)", "C2 (ppm)", [8.1, 6.7], [8.0, 7.0, 0.5], [155.5, 152.5], [155., 153.0, 1.0], 'A6_m3T9_adia_HSQC_C2.list', "C2", '%2.1f', '%d', "", 'r', 'test.ft2', 1)

plot_spectrum(ax_imino2, 'A6DNA-na_6.8_25mM_25C_H2O_Imino_HSQC_600/', 1.5*math.pow(10, 6), 20, 1.4, "H1 (ppm)", "N1 (ppm)", [13.2, 12.6], [13.2, 12.6, 0.3], [149, 146.3], [149., 147.0, 1.0], 'A6DNA_sofast_HMQC_Imino_G.list', "N1", '%2.1f', '%d', "", 'k', 'test.ft2', 1)
plot_spectrum(ax_imino2, 'A6DNAm3T9-na_6.8_25mM_25C_H2O_Imino_HSQC_600/', 1.*math.pow(10, 6), 20, 1.4, "H1 (ppm)", "N1 (ppm)", [13.2, 12.6], [13.2, 12.6, 0.3], [149, 146.3], [149., 147.0, 1.0], 'A6DNAm3T9_sofast_HMQC_Imino_G.list', "N3", '%2.1f', '%d', "", 'r', 'test.ft2', 1)


plt.savefig('figs7_f_2.pdf')
#plt.show()
