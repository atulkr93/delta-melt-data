# Plot NMR spectra
import matplotlib.pyplot as plt
import numpy as np
import nmrglue as ng
import math
import matplotlib

matplotlib.rcParams['axes.linewidth'] = 2.0
matplotlib.rcParams['font.family'] = 'Sans-serif'
matplotlib.rcParams['font.sans-serif']=['Arial']

elw = 2.0
tick_length=8.0
fs = 45

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
    global counter

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
    #[xlabel, ylabel, labels] = read_peaklist(dirname + peaklist, nucleus)
    #for dummy in range(len(labels)):
    #   ax.text(xlabel[dummy], ylabel[dummy], labels[dummy], size=fs, color='k')    
    #ax.errorbar(xlabel, ylabel, fmt='x', markersize=8, color='k')

    # Plotting
    ax.contour(data.transpose(), cl, colors='r', extent=(y_het_ppm_start, y_het_ppm_end, x_proton_ppm_start, x_proton_ppm_end))
    ax.set_ylim([ylim[0], ylim[1]])
    ax.set_xlim([xlim[0], xlim[1]])

    ax.set_xticks(np.arange(xtick_lims[1], xtick_lims[0]+0.01, xtick_lims[2]))
    ax.set_yticks(np.arange(ytick_lims[1], ytick_lims[0]+0.01, ytick_lims[2]))
    ax.set_xticklabels([xfmt%ele for ele in np.arange(xtick_lims[1], xtick_lims[0]+0.01, xtick_lims[2])], fontsize=fs)
    ax.set_xlabel(xaxis_label, fontsize=fs)
    if counter == 0:
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
    counter = counter + 1

fig, ax = plt.subplots(1, 4, figsize=(37, 8))
#fig, ax = plt.subplots(1, 4)
counter = 0
plot_spectrum(ax[0], 'a6dna_data/evgenia_AT6ul_m1dG10A16_pH5.2_09Jul11/11/', 1.1*math.pow(10, 6), 20, 1.4, "H6/8 (ppm)", "C6/8 (ppm)", [8.5, 6.8], [8.0, 7.0, 1], [148.5, 137.0], [148., 138., 2.], 'test.list', "N3", '%d', '%d', "A" + "$_{6}$" + "-DNA" + "$^{m1G10,m1A16}$" + "\npH 5.2 25 " + "$^{\circ}$" + "C")
plot_spectrum(ax[1], 'a2dna_g6g20/Henry-061818-A2DNA-G620m1G-pH5.4-10C-150NaCl/10/', 1*math.pow(10, 6), 20, 1.4, "H6/8 (ppm)", "C6/8 (ppm)", [8.5, 6.8], [8.0, 7.0, 1], [148.5, 137.0], [148., 138., 2.], 'test.list', "N3", '%d', '%d', "A" + "$_{2}$" + "-DNA" + "$^{m1G6,20}$" + "\npH 5.4HS 10 " + "$^{\circ}$" + "C")
plot_spectrum(ax[2], 'a2dna_athg/m1a717/Henry-061818-A2DNA-m1A717-pH6.8-10C-150NaCl/Caro/', 1.2*math.pow(10, 6), 20, 1.4, "H6/8 (ppm)", "C6/8 (ppm)", [8.5, 6.8], [8.0, 7.0, 1], [148.5, 137.0], [148., 138., 2.], 'test.list', "N3", '%d', '%d', "A" + "$_{2}$" + "-DNA" + "$^{m1A7,17}$" + "\npH 6.8HS 10 " + "$^{\circ}$" + "C")
plot_spectrum(ax[3], 'a2dna_athg2/m1a1617/Henry-061818-A2DNA-m1A1617-pH6.8-10C-150NaCl/10/', 1.0*math.pow(10, 6), 20, 1.4, "H6/8 (ppm)", "C6/8 (ppm)", [8.5, 6.8], [8.0, 7.0, 1], [148.5, 137.0], [148., 138., 2.], 'test.list', "N3", '%d', '%d', "A" + "$_{2}$" + "-DNA" + "$^{m1A16,17}$" + "\npH 6.8HS 10 " + "$^{\circ}$" + "C")

#plt.show()
#plt.tight_layout()
plt.savefig('figs9_c.pdf')
