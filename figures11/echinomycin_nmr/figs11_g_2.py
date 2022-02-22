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

def plot_spectrum(ax, dirname, contour_start, contour_num, contour_factor, xaxis_label, yaxis_label, xlim, xtick_lims, ylim, ytick_lims, peaklist, nucleus, xfmt, yfmt, plot_title, spectra_name, color_plot, xshift):
    # define contour levels
    cl = contour_start * contour_factor ** np.arange(contour_num)

    # Read data
    dic, data = ng.pipe.read(dirname + spectra_name)
    #data = ng.process.proc_base.add(data, xshift)

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
    #   ax.text(xlabel[dummy]+0.05, ylabel[dummy], labels[dummy], size=fs, color='k')    
    #ax.errorbar(xlabel, ylabel, fmt='x', markersize=8, color='k')

    # Plotting
    ax.contour(data.transpose(), cl, colors=color_plot, extent=(y_het_ppm_start, y_het_ppm_end, x_proton_ppm_start, x_proton_ppm_end))
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

fig, ax = plt.subplots()
plot_spectrum(ax, 'dna-ca-echin/Henry-080115-CA-echio-pH5-25/Henry-072715-CA-Echin-pH5-T25-HSQC2/4/', 1.9*math.pow(10, 6), 20, 1.4, "H6/8 (ppm)", "C6/8 (ppm)", [8.3, 6.7], [8.1, 7.0, 0.5], [146, 138], [146, 148, 2.0], 'test', "C6", '%2.1f', '%d', "E12DNA-HG", 'out.ft2', '#ff0000', 0.0)
plot_spectrum(ax, 'dna-ca-echin/Henry-050215-DNA-CA-echino-ph6.8/5/', 8*math.pow(10, 5), 20, 1.4, "H6/8 (ppm)", "C6/8 (ppm)", [8.3, 6.7], [8.1, 7.0, 0.5], [146, 138], [146, 138., 2.0], 'test', "C6", '%2.1f', '%d', "E12DNA-HG", 'out.ft2', '#8e4e9e', 0.0)
plot_spectrum(ax, 'dna-ca-echin/Henry-091716-CA-deazaG-Echin-pH5.3-25C/102/', 2.3*math.pow(10, 6), 20, 1.4, "H6/8 (ppm)", "C6/8 (ppm)", [8.3, 6.7], [8.1, 7.0, 0.5], [146, 138], [146, 138., 2.0], 'test', "C6", '%2.1f', '%d', "E12DNA-HG", 'out.ft2', '#30c2e3', 0.058)
plot_spectrum(ax, 'dna-ca-echin/Henry-110418-CA-m1G-Echin-pH6.8-25C-150NaCl/102/', 5*math.pow(10, 6), 20, 1.4, "H6/8 (ppm)", "C6/8 (ppm)", [8.3, 6.7], [8.1, 7.0, 0.5], [146, 138], [146, 138., 2.0], 'test', "C6", '%2.1f', '%d', "E12DNA-HG", 'out.ft2', 'g', -0.049)



#plt.show()
#plt.tight_layout()
plt.savefig('figs11_g_2.pdf')
