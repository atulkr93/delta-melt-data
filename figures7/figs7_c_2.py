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
ax_noe4 = fig.add_subplot(gs[9:, 18:21])
ax_noe5 = fig.add_subplot(gs[9:, 21:24])

plot_spectrum(ax_noe4, 'A6DNAm3T9-na_6.8_25mM_25C_H2O_NOESY_HSQC_600/', 5.*math.pow(10, 6), 20, 1.4, "H2/H6 (ppm)", "Methyl H (ppm)", [7.3, 7.075], [7.3, 7.1, 0.2], [2.9, 2.6], [2.9, 2.6, 0.2], 'A6DNAm3T9_150ms_NOESY_h2h62.list', "", '%2.1f', '%2.1f', "", 'r', 'test.ft2', 1)
plot_spectrum(ax_noe5, 'A6DNAm3T9-na_6.8_25mM_25C_H2O_NOESY_HSQC_600/', 5.*math.pow(10, 6), 20, 1.4, "H4a (ppm)", "Methyl H (ppm)", [7.1, 6.9], [7.1, 6.9, 0.1], [2.9, 2.6], [2.9, 2.6, 0.2], 'A6DNAm3T9_150ms_NOESY_c15h4a.list', "", '%2.1f', '%2.1f', "", 'r', 'test.ft2', 0)


plt.savefig('figs7_c_2.pdf')
#plt.show()
