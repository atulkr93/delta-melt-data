# Plot NMR spectra
import matplotlib.pyplot as plt
import numpy as np
import nmrglue as ng
import math
import matplotlib as mpl
import matplotlib

matplotlib.rcParams['axes.linewidth'] = 2.0
matplotlib.rcParams['font.family']='Sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
matplotlib.rcParams['axes.linewidth'] = 2.0
elw = 2.0
tick_length=8.0
fs = 45
fs2 = 37.0

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

def plot_1d_spectrum(ax, dirname, filename, xlim, xticks, xfmt_string, title, ylim, xlabel, plot_color, scaling, plot_status):
    dic, data = ng.pipe.read(dirname + filename)
    uc = ng.pipe.make_uc(dic, data)
    data = data / scaling
    if plot_status == 1:
        ax.plot(uc.ppm_scale(), data, 'k', linewidth=elw, color=plot_color, label=title)
    else:
        ax.plot(uc.ppm_scale(), data, 'k', linewidth=0, color=plot_color, label=title)
    ax.set_xlim(xlim)    
    ax.set_xticks(np.arange(xticks[0], xticks[1]-0.1, -1.0*xticks[2]))
    ax.set_xticklabels([xfmt_string%ele for ele in np.arange(xticks[0], xticks[1]-0.1, -1.0*xticks[2])], fontsize=fs)
    ax.set_yticks([])
    ax.tick_params(axis='x', which='major', direction='out', top=False, width=elw, length=tick_length, pad=6)
    ax.tick_params(axis='x', which='minor', direction='out', top=False, width=elw, length=0, pad=6)
    ax.tick_params(axis='y', which='major', direction='out', right=False, width=elw, length=tick_length, pad=6)
    ax.tick_params(axis='y', which='minor', direction='out', right=False, width=elw, length=0, pad=6)
    #title = ax.set_title(title, fontsize=fs)
    #title.set_position([0.5, 1.02])
    ax.set_ylim(ylim)
    ax.set_xlabel(xlabel, fontsize=fs)

def deal_legends(ax, fs):

    handles, labels = ax.get_legend_handles_labels()
    handles = [h[0] if type(h) is not mpl.lines.Line2D else h for h in handles]
    leg = ax.legend(handles, labels, fancybox=False, ncol=1, handlelength=0, frameon=False, columnspacing=0.0, markerscale=0.0, handletextpad=0.2, borderpad=0.0, handleheight=0.0, labelspacing=0.2, fontsize=fs2, loc=1)

    counter = 0
    for text in leg.get_texts():
        text.set_ha('left')
        text.set_color(color_slps[counter])
        counter = counter + 1
    for l in leg.get_lines():
        l.set_linestyle('None')
    

#fig, ax = plt.subplots(1, 5, figsize=(48, 8))
fig, ax = plt.subplots(1, 4, figsize=(37, 8))
#fig, ax = plt.subplots(1, 4)
color_slps = ['#99519f', '#63cae3', 'r', 'k']
plot_1d_spectrum(ax[0], 'a6dna_data/waj_data_m1a/1H1D_25C_030415.fid/', 'test.ft2', [15, 9], [15.0, 9.0, 2.0], "%d", "A" + "$_{6}$" + "-DNA" + "$^{m1A16}$" + " pH 6.8", [-1000., 2*math.pow(10,4)], "H1/H3 (ppm)", '#99519f', 1, 1)
plot_1d_spectrum(ax[0], 'a6dna_data/a6dna_m1g/jane_pH5.4_25C/', 'test.ft2', [15, 9], [15.0, 9.0, 2.0], "%d", "A" + "$_{6}$" + "-DNA" + "$^{m1G10}$" + " pH 5.4", [-1000., 2*math.pow(10,4)], "H1/H3 (ppm)", '#63cae3', 10, 1)
plot_1d_spectrum(ax[0], 'a6dna_data/evgenia_AT6ul_m1dG10A16_pH5.2_09Jul11/H1D/', 'test.ft2', [15, 9], [15.0, 9.0, 2.0], "%d", "A" + "$_{6}$" + "-DNA" + "$^{m1G10,m1A16}$" + " pH 5.2", [-1000., 2*math.pow(10,4)], "H1/H3 (ppm)", 'r', 10, 1)
plot_1d_spectrum(ax[0], 'a6dna_data/honglue_data/pH6.8_25C/', 'test.ft2', [15, 9], [15.0, 9.0, 2.0], "%d", "25 " + "$^{\circ}$" + "C", [-1000., 2*math.pow(10,4)], "H1/H3 (ppm)", '#808080', 10000, 0)
deal_legends(ax[0], fs2)

plot_1d_spectrum(ax[1], 'a2dna_g6g20/Henry-061818-A2DNA-G6m1G-pH5.4-10C-150NaCl-450uM/H1D/', 'test.ft2', [15, 8.8], [15.0, 9.0, 2.0], "%d", "A" + "$_{2}$" + "-DNA" + "$^{m1G6}$", [-20000., 2*math.pow(10,5)], "H1/H3 (ppm)", '#99519f', 10, 1)
plot_1d_spectrum(ax[1], 'a2dna_g6g20/atul-20180312-dsA2DNAm1G20-10C/H1D/', 'test.ft2', [15, 8.8], [15.0, 9.0, 2.0], "%d", "A" + "$_{2}$" + "-DNA" + "$^{m1G20}$", [-20000., 2*math.pow(10,5)], "H1/H3 (ppm)", '#63cae3', 10, 1)
plot_1d_spectrum(ax[1], 'a2dna_g6g20/Henry-061818-A2DNA-G620m1G-pH5.4-10C-150NaCl/H1D/', 'test.ft2', [15, 8.8], [15.0, 9.0, 2.0], "%d", "A" + "$_{2}$" + "-DNA" + "$^{m1G6,20}$", [-20000., 2*math.pow(10,5)], "H1/H3 (ppm)", 'r', 10, 1)
plot_1d_spectrum(ax[1], 'a2dna_g6g20/atul-20171215-dsA2DNA-pH6.8NB-HS/H1D/', 'test.ft2', [15, 8.8], [15.0, 9.0, 2.0], "%d", "pH 5.4HS\n10 " + "$^{\circ}$" + "C", [-20000., 2*math.pow(10,5)], "H1/H3 (ppm)", '#808080', 100, 0)
deal_legends(ax[1], fs2)

plot_1d_spectrum(ax[2], 'a2dna_athg/m1a7/Henry-061818-A2DNA-m1A7-pH6.8-10C-150NaCl/SFIMINO/', 'test.ft2', [15, 8.8], [15.0, 9.0, 2.0], "%d", "A" + "$_{2}$" + "-DNA" + "$^{m1A7}$", [-20000., 2*math.pow(10,5)], "H1/H3 (ppm)", '#99519f', 100, 1)
plot_1d_spectrum(ax[2], 'a2dna_athg/atul-20171215-dsA2DNAm1A17-pH6.8NB-HS-10C/SOFAST_IMINO/', 'test.ft2', [15, 8.8], [15.0, 9.0, 2.0], "%d", "A" + "$_{2}$" + "-DNA" + "$^{m1A17}$", [-20000., 2*math.pow(10,5)], "H1/H3 (ppm)", '#63cae3', 100, 1)
plot_1d_spectrum(ax[2], 'a2dna_athg/m1a717/Henry-061818-A2DNA-m1A717-pH6.8-10C-150NaCl/H1D/', 'test.ft2', [15, 8.8], [15.0, 9.0, 2.0], "%d", "A" + "$_{2}$" + "-DNA" + "$^{m1A7,17}$", [-20000., 2*math.pow(10,5)], "H1/H3 (ppm)", 'r', 10, 1)
plot_1d_spectrum(ax[2], 'a2dna_g6g20/atul-20171215-dsA2DNA-pH6.8NB-HS/H1D/', 'test.ft2', [15, 8.8], [15.0, 9.0, 2.0], "%d", "pH 6.8HS\n10 " + "$^{\circ}$" + "C", [-20000., 2*math.pow(10,5)], "H1/H3 (ppm)", '#808080', 100, 0)
deal_legends(ax[2], fs2)

plot_1d_spectrum(ax[3], 'a2dna_athg2/m1a16/atul-20180312-dsA2DNAm1A16-10C/H1D/', 'test.ft2', [15, 8.8], [15.0, 9.0, 2.0], "%d", "A" + "$_{2}$" + "-DNA" + "$^{m1A16}$" + " 10 " + "$^{\circ}$" + "C", [-20000., 2*math.pow(10,5)], "H1/H3 (ppm)", '#99519f', 100, 1)
plot_1d_spectrum(ax[3], 'a2dna_athg2/m1a1617/atul-20180312-dsA2DNAm1A1617-1C/SFIMINO/', 'test.ft2', [15, 8.8], [15.0, 9.0, 2.0], "%d", "A" + "$_{2}$" + "-DNA" + "$^{m1A16,17}$" + " 1 " + "$^{\circ}$" + "C", [-20000., 2*math.pow(10,5)], "H1/H3 (ppm)", '#63cae3', 100, 1)
deal_legends(ax[3], fs2)
#plt.show()
#plt.tight_layout()
plt.savefig('figs9_b.pdf')
