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
       ax.text(xlabel[dummy]+0.1, ylabel[dummy], labels[dummy], size=fs, color='k')    
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

fig, ax = plt.subplots(1, 6, figsize=(50.4, 9))
#plot_spectrum(ax[0], 'A6DNAulb_6.8_25mM_25C_H2O_C1H1_HSQC_700/', 3.*math.pow(10, 6), 20, 1.4, "H1' (ppm)", "C1' (ppm)", [6.5, 5.0], [6.5, 5.0, 0.5], [88.5, 83.5], [88.0, 84.0, 2.0], 'A6DNA_ulb_HSQC_Ali_6.8_25C.list', "C1'", '%2.1f', '%d', "A" + "$_{6}$" + "-DNA(ulb)\npH 6.8 25 " + "$^{\circ}$" + "C")
#plot_spectrum(ax[1], 'A6DNAulb_6.8_25mM_25C_H2O_C8H8_HSQC_700/', 1.*math.pow(10, 6), 20, 1.4, "H6/8 (ppm)", "C6/8 (ppm)", [8.3, 6.9], [8.0, 7.0, 0.5], [144.5, 137.], [144., 137., 3.], 'A6DNA_ulb_HSQC_Aro_6.8_25C.list', "C6", '%2.1f', '%d', "A" + "$_{6}$" + "-DNA(ulb)\npH 6.8 25 " + "$^{\circ}$" + "C")
#plot_spectrum(ax[2], 'A6DNAs1lb_5.4_25mM_25C_D2O_C1H1_HSQC_700/', 3.*math.pow(10, 6), 20, 1.4, "H1' (ppm)", "C1' (ppm)", [6.3, 5.1], [6.0, 5.0, 0.5], [87.2, 83.5], [87.0, 84.0, 1.0],  'A6lbDNA_HSQC_C1p.list', "C1'", '%2.1f', '%d', "A" + "$_{6}$" + "-DNA(s1lb)\npH 5.4 (D" + "$_{2}$" + "O) 25 " + "$^{\circ}$" + "C")
#plot_spectrum(ax[3], 'A6DNAs1lb_5.4_25mM_25C_D2O_C8H8_HSQC_700/', 3.*math.pow(10, 5), 20, 1.4, "H6/8 (ppm)", "C6/8 (ppm)", [8.3, 6.9], [8.0, 7.0, 0.5], [144.4, 137.8], [144., 138., 3.0], 'A6lbDNA_HSQC_Aro.list', "C6", '%2.1f', '%d', "A" + "$_{6}$" + "-DNA(s1lb)\npH 5.4 (D" + "$_{2}$" + "O) 25 " + "$^{\circ}$" + "C")
#plot_spectrum(ax[4], 'A6DNAs2lb_5.4_25mM_25C_D2O_C1H1_HSQC_700/', 3.*math.pow(10, 6), 20, 1.4, "H1' (ppm)", "C1' (ppm)", [6.5, 5.3], [6.5, 5.5, 0.5], [89., 83.], [89., 83., 2.0], 'T6lbDNA_HSQC_C1p.list', "C1'", '%2.1f', '%d', "A" + "$_{6}$" + "-DNA(s2lb)\npH 5.4 (D" + "$_{2}$" + "O) 25 " + "$^{\circ}$" + "C")
#plot_spectrum(ax[0], 'A6DNAs2lb_5.4_25mM_25C_D2O_C8H8_HSQC_700/', 5.*math.pow(10, 5), 20, 1.4, "H6/8 (ppm)", "C6/8 (ppm)", [8.5, 6.9], [8.5, 7.0, 0.5], [144., 137.], [144., 138., 2.], 'T6lbDNA_HSQC_Aro.list', "C6", '%2.1f', '%d', "A" + "$_{6}$" + "-DNA(s2lb)\npH 5.4 (D" + "$_{2}$" + "O) 25 " + "$^{\circ}$" + "C")
#plot_spectrum(ax[1], 'A6DNAs2lb_5.4_25mM_10C_H2O_NH_HSQC_700/', 5.*math.pow(10, 6), 20, 1.4, "H3 (ppm)", "N3 (ppm)", [14.25, 13.65], [14.2, 13.7, 0.2], [160.6, 158.], [160., 158., 1.], 'A6DNA_T6lb_NHHSQC_Imino_5.4_10C_N3.list', "N3", '%2.1f', '%d', "A" + "$_{6}$" + "-DNA(s2lb)\npH 5.4 10 " + "$^{\circ}$" + "C")
#plot_spectrum(ax[2], 'A6DNAs2lb_5.4_25mM_10C_H2O_NH_HSQC_700/', 5.*math.pow(10, 6), 20, 1.4, "H1 (ppm)", "N1 (ppm)", [13.2, 12.75], [13.2, 12.8, 0.2], [148., 146.], [148., 146., 1.], 'A6DNA_T6lb_NHHSQC_Imino_5.4_10C_N1.list', "N1", '%2.1f', '%d', "A" + "$_{6}$" + "-DNA(s2lb)\npH 5.4 10 " + "$^{\circ}$" + "C")
#plot_spectrum(ax[3], 'A6DNAA16lb_5.4_150mM_25C_H2O_C1H1_HSQC_600/', 1.*math.pow(10, 6), 20, 1.4, "H1' (ppm)", "C1' (ppm)", [5.9, 5.5], [5.9, 5.5, 0.2], [85., 84.], [85., 84., 1.], 'A6DNAA16lb_high_5.4_HSQC_25C_Ali.list', "C1'", '%2.1f', '%d', "A" + "$_{6}$" + "-DNA(A16lb)\npH 5.4HS 25 " + "$^{\circ}$" + "C")
#plot_spectrum(ax[4], 'A6DNAA16lb_5.4_150mM_25C_H2O_C8H8_HSQC_600/', 1.*math.pow(10, 6), 20, 1.4, "H8 (ppm)", "C8 (ppm)", [8.4, 8.0], [8.4, 8.0, 0.2], [142., 141.], [142., 141., 1.0], 'A6DNAA16lb_high_5.4_HSQC_25C_Aro.list', "C8", '%2.1f', '%d', "A" + "$_{6}$" + "-DNA(A16lb)\npH 5.4HS 25 " + "$^{\circ}$" + "C")



#plot_spectrum(ax[0], 'A6DNAm1A16s2lb_5.4_150mM_25C_H2O_C1H1_HSQC_600/', 3.*math.pow(10, 6), 20, 1.4, "H1' (ppm)", "C1' (ppm)", [6.5, 5.3], [6.5, 5.3, 0.4], [88.5, 83.5], [87.0, 85.0, 1.0], 'A6DNAm1A16T6lb_high_5.4_HSQC_25C_Ali.list', "C1'", '%2.1f', '%d', "A" + "$_{6}$" + "-DNA" + "$^{m1A16}$" + "(s2lb)\npH 5.4HS 25 " + "$^{\circ}$" + "C")
#plot_spectrum(ax[1], 'A6DNAm1A16s2lb_5.4_150mM_25C_H2O_C8H8_HSQC_600/', 2.*math.pow(10, 6), 20, 1.4, "H6/8 (ppm)", "C6/8 (ppm)", [8.38, 7], [8.2, 7., 0.4], [144, 136.6], [144., 137., 3.], 'A6DNAm1A16T6lb_high_5.4_HSQC_25C_Aro.list', "C6", '%2.1f', '%d', "A" + "$_{6}$" + "-DNA" + "$^{m1A16}$" + "(s2lb)\npH 5.4HS 25 " + "$^{\circ}$" + "C")
#plot_spectrum(ax[0], 'A6DNAm1A16s2lb_5.4_25mM_25C_D2O_C1H1_HSQC_700/', 8.*math.pow(10, 6), 20, 1.4, "H1' (ppm)", "C1' (ppm)", [6.5, 5.3], [6.5, 5.3, 0.4], [88.5, 83.5], [87.0, 85.0, 1.0], 'A6DNAm1A16_T6lb_HSQC_C1p.list', "C1'", '%2.1f', '%d', "A" + "$_{6}$" + "-DNA" + "$^{m1A16}$" + "(s2lb)\npH 5.4 (D" + "$_{2}$" + "O) 25 " + "$^{\circ}$" + "C")
#plot_spectrum(ax[1], 'A6DNAm1A16s2lb_5.4_25mM_25C_D2O_C8H8_HSQC_700/', 6.*math.pow(10, 5), 20, 1.4, "H6/8 (ppm)", "C6/8 (ppm)", [8.38, 7], [8.2, 7., 0.4], [144, 136.6], [144., 137., 3.], 'A6DNAm1A16_T6lb_HSQC_Aro.list', "C6", '%2.1f', '%d', "A" + "$_{6}$" + "-DNA" + "$^{m1A16}$" + "(s2lb)\npH 5.4 (D" + "$_{2}$" + "O) 25 " + "$^{\circ}$" + "C")


plot_spectrum(ax[0], 'A6DNAm1A16s2lb_5.4_25mM_10C_H2O_NH_HSQC_700/', 5.*math.pow(10, 6), 20, 1.4, "H3 (ppm)", "N3 (ppm)", [14.3, 11.5], [14., 12., 1.], [160.5, 155.7], [160., 156., 2.], 'A6DNAm1A16_sofast_HMQC_Imino_10C_N3.list', "N3", '%d', '%d', "A" + "$_{6}$" + "-DNA" + "$^{m1A16}$" + "(s2lb)\npH 5.4 10 " + "$^{\circ}$" + "C")
plot_spectrum(ax[1], 'A6DNAm1A16s2lb_5.4_25mM_10C_H2O_NH_HSQC_700/', 5.*math.pow(10, 6), 20, 1.4, "H1 (ppm)", "N1 (ppm)", [13.2, 12.6], [13.2, 12.7, 0.4], [148.5, 146.0], [148, 146, 1.], 'A6DNAm1A16_sofast_HMQC_Imino_10C_N1.list', "N1", '%2.1f', '%d', "A" + "$_{6}$" + "-DNA" + "$^{m1A16}$" + "(s2lb)\npH 5.4 10 " + "$^{\circ}$" + "C")
plot_spectrum(ax[2], 'A6DNAm1G10s1lb_5.4_25mM_25C_D2O_C1H1_HSQC_700/', 1.*math.pow(10, 6), 20, 1.4, "H1' (ppm)", "C1' (ppm)", [6.3, 5.5], [6.1, 5.6, 0.4], [87.5, 83.5], [87.0, 84.0, 1.0], 'A6DNAm1G10_A6lb_HSQC_C1p.list', "C1'", '%2.1f', '%d', "A" + "$_{6}$" + "-DNA" + "$^{m1G10}$" + "(s1lb)\npH 5.4 (D" + "$_{2}$" + "O) 25 " + "$^{\circ}$" + "C")
plot_spectrum(ax[3], 'A6DNAm1G10s1lb_5.4_25mM_25C_D2O_C8H8_HSQC_700/', 9.*math.pow(10, 5), 20, 1.4, "H6/8 (ppm)", "C6/8 (ppm)", [8.4, 6.9], [8., 7., 0.5], [146., 138.], [146., 138., 4.], 'A6DNAm1G10_A6lb_HSQC_Aro.list', "C6", '%2.1f', '%d', "A" + "$_{6}$" + "-DNA" + "$^{m1G10}$" + "(s1lb)\npH 5.4 (D" + "$_{2}$" + "O) 25 " + "$^{\circ}$" + "C")
plot_spectrum(ax[4], 'hz-120115-lblrA16dC15A6DNA-10C/Caro/', 1.1*math.pow(10, 5), 20, 1.4, "H6/8 (ppm)", "C6/8 (ppm)", [8.4, 7.2], [8.4, 7.3, 0.4], [143.5, 140.5], [143., 141., 1.], 'rA16dC15dsA6DNA_pH5.4_10C.list', "C6", '%2.1f', '%d', "A" + "$_{6}$" + "-DNA" + "$^{rA16}$" + "(rA16" + "$^{C8}$" + ",C15lb)\npH 5.4 10 " + "$^{\circ}$" + "C")
plot_spectrum(ax[5], 'ijk-021615-dl-a-rG10/Caro/Caro/', 2.35*math.pow(10, 5), 20, 1.4, "H6/8 (ppm)", "C6/8 (ppm)", [8.2, 6.9], [8.0, 7.0, 0.5], [144., 138], [144., 138., 3.], 'S1lb_rA16dsA6DNA_pH5.4_25C_Caro.list', "C6", '%2.1f', '%d', "A" + "$_{6}$" + "-DNA" + "$^{rG10}$" + "(s1lb)\npH 5.4 25 " + "$^{\circ}$" + "C")




#plt.show()
plt.tight_layout()
plt.savefig('figs4_d_3.pdf')
