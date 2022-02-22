# Plot NMR spectra
import matplotlib.pyplot as plt
import numpy as np
import nmrglue as ng
import math
import matplotlib

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

def plot_1d_spectrum(ax, dirname, filename, xlim, xticks, xfmt_string, title, ylim, xlabel):
    dic, data = ng.pipe.read(dirname + filename)
    uc = ng.pipe.make_uc(dic, data)
    ax.plot(uc.ppm_scale(), data, 'k', linewidth=elw)
    ax.set_xlim(xlim)    
    ax.set_xticks(np.arange(xticks[0], xticks[1]-0.1, -1.0*xticks[2]))
    ax.set_xticklabels([xfmt_string%ele for ele in np.arange(xticks[0], xticks[1]-0.1, -1.0*xticks[2])], fontsize=fs)
    ax.set_yticks([])
    ax.tick_params(axis='x', which='major', direction='out', top=False, width=elw, length=tick_length, pad=6)
    ax.tick_params(axis='x', which='minor', direction='out', top=False, width=elw, length=0, pad=6)
    ax.tick_params(axis='y', which='major', direction='out', right=False, width=elw, length=tick_length, pad=6)
    ax.tick_params(axis='y', which='minor', direction='out', right=False, width=elw, length=0, pad=6)
    title = ax.set_title(title, fontsize=fs)
    title.set_position([0.5, 1.02])
    ax.set_ylim(ylim)
    ax.set_xlabel(xlabel, fontsize=fs)

fig, ax = plt.subplots(1, 5, figsize=(48, 8))
#fig, ax = plt.subplots(1, 2)
plot_1d_spectrum(ax[0], 'waj_data/1H1D_25C_030415.fid/', 'test.ft2', [14.3, 12.5], [14.0, 13.0, 0.5], "%2.1f", "A" + "$_{6}$" + "-DNA pH 6.8 25 " + "$^{\circ}$" + "C", [-10000., 130000], "H1/H3 (ppm)")
plot_1d_spectrum(ax[1], 'GDNA_1D/', 'test.ft2', [14.0, 12.0], [14.0, 12.0, 1.0], "%d", "GDNA pH 8.8 " + "$^{\circ}$" + "C", [-2000000, 38000000], "H1/H3 (ppm)")
ax[2].remove()
ax[3].remove()
ax[4].remove()
#plt.show()
plt.savefig('figs4_b.pdf')