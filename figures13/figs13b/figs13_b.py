# Plot NMR spectra
import matplotlib as mpl
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

def plot_1d_spectrum(ax, dirname, filename, xlim, xticks, xfmt_string, title, ylim, xlabel, color_plot, multiplier):
    dic, data = ng.pipe.read(dirname + filename)
    uc = ng.pipe.make_uc(dic, data)
    ax.plot(uc.ppm_scale(), data*multiplier, color_plot, linewidth=elw, label=title)
    ax.set_xlim(xlim)    
    ax.set_xticks(np.arange(xticks[0], xticks[1]-0.1, -1.0*xticks[2]))
    ax.set_yticks([])
    ax.tick_params(axis='x', which='major', direction='out', top=False, width=elw, length=tick_length, pad=6)
    ax.tick_params(axis='x', which='minor', direction='out', top=False, width=elw, length=0, pad=6)
    ax.tick_params(axis='y', which='major', direction='out', right=False, width=elw, length=tick_length, pad=6)
    ax.tick_params(axis='y', which='minor', direction='out', right=False, width=elw, length=0, pad=6)
    #title = ax.set_title(title, fontsize=fs)
    #title.set_position([0.5, 1.02])
    ax.set_ylim(ylim)
    ax.set_xticklabels([]) 

def clear_legends(ax): 
    handles, labels = ax.get_legend_handles_labels()
    handles = [h[0] if type(h) is not mpl.lines.Line2D else h for h in handles]
    leg = ax.legend(handles, labels, title='', fancybox=False, ncol=1, handlelength=0, frameon=False, columnspacing=0.0, markerscale=0.0, handletextpad=0.2, borderpad=0.0, handleheight=0.0, labelspacing=0.2, fontsize=fs)

    color_slps = ['k', 'r']
    counter = 0
    for text in leg.get_texts():
        text.set_color(color_slps[counter])
        counter = counter + 1
    for l in leg.get_lines():
        l.set_linestyle('None')

#fig, ax = plt.subplots(1, 5, figsize=(48, 8))
fig, ax = plt.subplots(4, 4, figsize=(24, 22))
plot_1d_spectrum(ax[0, 0], 'atul-20200627-scaf2AGA_GC_pH5.4HSNB/4/', 'test.ft2', [15.0, 9.0], [15.0, 9.0, 2.0], "%2.1f", "scaf2_AGA" + "$^{GC}$", np.array([-0.25, 3]) * math.pow(10, 7), "H1/H3 (ppm)", 'k', 1.0)
plot_1d_spectrum(ax[0, 0], 'atul-20200629-scaf2AGA_m1GC_pH5.4HSNB/4/', 'test.ft2', [15.0, 9.0], [15.0, 9.0, 2.0], "%2.1f", "scaf2_AGA" + "$^{m1GC}$", np.array([-0.25, 3]) * math.pow(10, 7), "H1/H3 (ppm)", 'r', 0.7)
clear_legends(ax[0, 0])

plot_1d_spectrum(ax[0, 1], 'atul-20200627-scaf2AGC_GC_pH5.4HSNB/4/', 'test.ft2', [15.0, 9.0], [15.0, 9.0, 2.0], "%2.1f", "scaf2_AGC" + "$^{GC}$", np.array([-0.25, 3.2]) * math.pow(10, 7), "H1/H3 (ppm)", 'k', 1.0)
plot_1d_spectrum(ax[0, 1], 'atul-20200629-scaf2AGC_m1GC_pH5.4HSNB/4/', 'test.ft2', [15.0, 9.0], [15.0, 9.0, 2.0], "%2.1f", "scaf2_AGC" + "$^{m1GC}$", np.array([-0.25, 3.2]) * math.pow(10, 7), "H1/H3 (ppm)", 'r', 1.8)
clear_legends(ax[0, 1])

plot_1d_spectrum(ax[0, 2], 'atul-20200627-scaf2AGG_GC_pH5.4HSNB/4/', 'test.ft2', [15.0, 9.0], [15.0, 9.0, 2.0], "%2.1f", "scaf2_AGG" + "$^{GC}$", np.array([-0.25, 2]) * math.pow(10, 7), "H1/H3 (ppm)", 'k', 1.0)
plot_1d_spectrum(ax[0, 2], 'atul-20200629-scaf2AGG_m1GC_pH5.4HSNB/4/', 'test.ft2', [15.0, 9.0], [15.0, 9.0, 2.0], "%2.1f", "scaf2_AGG" + "$^{m1GC}$", np.array([-0.25, 2]) * math.pow(10, 7), "H1/H3 (ppm)", 'r', 1.0)
clear_legends(ax[0, 2])

plot_1d_spectrum(ax[0, 3], 'atul-20200627-scaf2AGT_GC_pH5.4HSNB/4/', 'test.ft2', [15.0, 9.0], [15.0, 9.0, 2.0], "%2.1f", "scaf2_AGT" + "$^{GC}$", np.array([-0.25, 3]) * math.pow(10, 7), "H1/H3 (ppm)", 'k', 1.0)
plot_1d_spectrum(ax[0, 3], 'atul-20200629-scaf2AGT_m1GC_pH5.4HSNB/4/', 'test.ft2', [15.0, 9.0], [15.0, 9.0, 2.0], "%2.1f", "scaf2_AGT" + "$^{m1GC}$", np.array([-0.25, 3]) * math.pow(10, 7), "H1/H3 (ppm)", 'r', 1.3)
clear_legends(ax[0, 3])

plot_1d_spectrum(ax[1, 0], 'atul-20200627-scaf2CGA_GC_pH5.4HSNB/4/', 'test.ft2', [15.0, 9.0], [15.0, 9.0, 2.0], "%2.1f", "scaf2_CGA" + "$^{GC}$", np.array([-0.25, 3.2]) * math.pow(10, 7), "H1/H3 (ppm)", 'k', 1.0)
plot_1d_spectrum(ax[1, 0], 'atul-20200627-scaf2CGA_m1GC_pH5.4HSNB/4/', 'test.ft2', [15.0, 9.0], [15.0, 9.0, 2.0], "%2.1f", "scaf2_CGA" + "$^{m1GC}$", np.array([-0.25, 3.2]) * math.pow(10, 7), "H1/H3 (ppm)", 'r', 1.3)
clear_legends(ax[1, 0])

plot_1d_spectrum(ax[1, 1], 'atul-20200627-scaf2CGC_GC_pH5.4HSNB/4/', 'test.ft2', [15.0, 9.0], [15.0, 9.0, 2.0], "%2.1f", "scaf2_CGC" + "$^{GC}$", np.array([-0.25, 5]) * math.pow(10, 7), "H1/H3 (ppm)", 'k', 1.0)
plot_1d_spectrum(ax[1, 1], 'atul-20200627-scaf2CGC_m1GC_pH5.4HSNB/4/', 'test.ft2', [15.0, 9.0], [15.0, 9.0, 2.0], "%2.1f", "scaf2_CGC" + "$^{m1GC}$", np.array([-0.25, 5]) * math.pow(10, 7), "H1/H3 (ppm)", 'r', 1.1)
clear_legends(ax[1, 1])

plot_1d_spectrum(ax[1, 2], 'atul-20200627-scaf2CGG_GC_pH5.4HSNB/4/', 'test.ft2', [15.0, 9.0], [15.0, 9.0, 2.0], "%2.1f", "scaf2_CGG" + "$^{GC}$", np.array([-0.25, 4.]) * math.pow(10, 7), "H1/H3 (ppm)", 'k', 1.0)
plot_1d_spectrum(ax[1, 2], 'atul-20200629-scaf2CGG_m1GC_pH5.4HSNB/4/', 'test.ft2', [15.0, 9.0], [15.0, 9.0, 2.0], "%2.1f", "scaf2_CGG" + "$^{m1GC}$", np.array([-0.25, 4.]) * math.pow(10, 7), "H1/H3 (ppm)", 'r', 1.2)
clear_legends(ax[1, 2])

plot_1d_spectrum(ax[1, 3], 'atul-20200627-scaf2CGT_GC_pH5.4HSNB/4/', 'test.ft2', [15.0, 9.0], [15.0, 9.0, 2.0], "%2.1f", "scaf2_CGT" + "$^{GC}$", np.array([-0.25, 3.2]) * math.pow(10, 7), "H1/H3 (ppm)", 'k', 1.0)
plot_1d_spectrum(ax[1, 3], 'atul-20200629-scaf2CGT_m1GC_pH5.4HSNB/4/', 'test.ft2', [15.0, 9.0], [15.0, 9.0, 2.0], "%2.1f", "scaf2_CGT" + "$^{m1GC}$", np.array([-0.25, 3.2]) * math.pow(10, 7), "H1/H3 (ppm)", 'r', 1.2)
clear_legends(ax[1, 3])

plot_1d_spectrum(ax[2, 0], 'atul-20200627-scaf2GGA_GC_pH5.4HSNB/4/', 'test.ft2', [15.0, 9.0], [15.0, 9.0, 2.0], "%2.1f", "scaf2_GGA" + "$^{GC}$", np.array([-0.25, 4.2]) * math.pow(10, 7), "H1/H3 (ppm)", 'k', 1.0)
plot_1d_spectrum(ax[2, 0], 'atul-20200629-scaf2GGA_m1GC_pH5.4HSNB/4/', 'test.ft2', [15.0, 9.0], [15.0, 9.0, 2.0], "%2.1f", "scaf2_GGA" + "$^{m1GC}$", np.array([-0.25, 4.2]) * math.pow(10, 7), "H1/H3 (ppm)", 'r', 1.4)
clear_legends(ax[2, 0])

plot_1d_spectrum(ax[2, 1], 'atul-20200627-scaf2GGC_GC_pH5.4HSNB/4/', 'test.ft2', [15.0, 9.0], [15.0, 9.0, 2.0], "%2.1f", "scaf2_GGC" + "$^{GC}$", np.array([-0.25, 4.2]) * math.pow(10, 7), "H1/H3 (ppm)", 'k', 1.0)
plot_1d_spectrum(ax[2, 1], 'atul-20200629-scaf2GGC_m1GC_pH5.4HSNB/4/', 'test.ft2', [15.0, 9.0], [15.0, 9.0, 2.0], "%2.1f", "scaf2_GGC" + "$^{m1GC}$", np.array([-0.25, 4.2]) * math.pow(10, 7), "H1/H3 (ppm)", 'r', 1.1)
clear_legends(ax[2, 1])

plot_1d_spectrum(ax[2, 2], 'atul-20200627-scaf2GGG_GC_pH5.4HSNB/4/', 'test.ft2', [15.0, 9.0], [15.0, 9.0, 2.0], "%2.1f", "scaf2_GGG" + "$^{GC}$", np.array([-0.25, 3.9]) * math.pow(10, 7), "H1/H3 (ppm)", 'k', 1.0)
plot_1d_spectrum(ax[2, 2], 'atul-20200629-scaf2GGG_m1GC_pH5.4HSNB/4/', 'test.ft2', [15.0, 9.0], [15.0, 9.0, 2.0], "%2.1f", "scaf2_GGG" + "$^{m1GC}$", np.array([-0.25, 3.9]) * math.pow(10, 7), "H1/H3 (ppm)", 'r', 1.1)
clear_legends(ax[2, 2])

plot_1d_spectrum(ax[2, 3], 'atul-20200627-scaf2GGT_GC_pH5.4HSNB/4/', 'test.ft2', [15.0, 9.0], [15.0, 9.0, 2.0], "%2.1f", "scaf2_GGT" + "$^{GC}$", np.array([-0.25, 2.8]) * math.pow(10, 7), "H1/H3 (ppm)", 'k', 1.0)
plot_1d_spectrum(ax[2, 3], 'atul-20200629-scaf2GGT_m1GC_pH5.4HSNB/4/', 'test.ft2', [15.0, 9.0], [15.0, 9.0, 2.0], "%2.1f", "scaf2_GGT" + "$^{m1GC}$", np.array([-0.25, 2.8]) * math.pow(10, 7), "H1/H3 (ppm)", 'r', 1.1)
clear_legends(ax[2, 3])

plot_1d_spectrum(ax[3, 0], 'atul-20200627-scaf2TGA_GC_pH5.4HSNB/4/', 'test.ft2', [15.0, 9.0], [15.0, 9.0, 2.0], "%2.1f", "scaf2_TGA" + "$^{GC}$", np.array([-0.25, 2.5]) * math.pow(10, 7), "H1/H3 (ppm)", 'k', 1.0)
plot_1d_spectrum(ax[3, 0], 'atul-20200629-scaf2TGA_m1GC_pH5.4HSNB/4/', 'test.ft2', [15.0, 9.0], [15.0, 9.0, 2.0], "%2.1f", "scaf2_TGA" + "$^{m1GC}$", np.array([-0.25, 2.5]) * math.pow(10, 7), "H1/H3 (ppm)", 'r', 1.1)
clear_legends(ax[3, 0])

plot_1d_spectrum(ax[3, 1], 'atul-20200627-scaf2TGC_GC_pH5.4HSNB/4/', 'test.ft2', [15.0, 9.0], [15.0, 9.0, 2.0], "%2.1f", "scaf2_TGC" + "$^{GC}$", np.array([-0.25, 4.7]) * math.pow(10, 7), "H1/H3 (ppm)", 'k', 1.0)
plot_1d_spectrum(ax[3, 1], 'atul-20200629-scaf2TGC_m1GC_pH5.4HSNB/4/', 'test.ft2', [15.0, 9.0], [15.0, 9.0, 2.0], "%2.1f", "scaf2_TGC" + "$^{m1GC}$", np.array([-0.25, 4.7]) * math.pow(10, 7), "H1/H3 (ppm)", 'r', 1.2)
clear_legends(ax[3, 1])

plot_1d_spectrum(ax[3, 2], 'atul-20200627-scaf2TGG_GC_pH5.4HSNB/4/', 'test.ft2', [15.0, 9.0], [15.0, 9.0, 2.0], "%2.1f", "scaf2_TGG" + "$^{GC}$", np.array([-0.25, 4.3]) * math.pow(10, 7), "H1/H3 (ppm)", 'k', 1.0)
plot_1d_spectrum(ax[3, 2], 'atul-20200629-scaf2TGG_m1GC_pH5.4HSNB/4/', 'test.ft2', [15.0, 9.0], [15.0, 9.0, 2.0], "%2.1f", "scaf2_TGG" + "$^{m1GC}$", np.array([-0.25, 4.3]) * math.pow(10, 7), "H1/H3 (ppm)", 'r', 1.2)
clear_legends(ax[3, 2])

plot_1d_spectrum(ax[3, 3], 'atul-20200627-scaf2TGT_GC_pH5.4HSNB/4/', 'test.ft2', [15.0, 9.0], [15.0, 9.0, 2.0], "%2.1f", "scaf2_TGT" + "$^{GC}$", np.array([-0.25, 2.]) * math.pow(10, 7), "H1/H3 (ppm)", 'k', 1.0)
plot_1d_spectrum(ax[3, 3], 'atul-20200627-scaf2TGT_m1GC_pH5.4HSNB/4/', 'test.ft2', [15.0, 9.0], [15.0, 9.0, 2.0], "%2.1f", "scaf2_TGT" + "$^{m1GC}$", np.array([-0.25, 2.]) * math.pow(10, 7), "H1/H3 (ppm)", 'r', 1.0)
clear_legends(ax[3, 3])

for dummy in range(4):
    ax[3, dummy].set_xlabel('H1/3 (ppm)', fontsize=fs)
    ax[3, dummy].set_xticklabels(['%d'%ele for ele in np.arange(15, 8, -2)], fontsize=fs)

#plt.show()
plt.savefig('figs13_b.pdf')
