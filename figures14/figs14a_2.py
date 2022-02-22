from collections import OrderedDict
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import matplotlib

matplotlib.rcParams['axes.linewidth'] = 2.0
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = 'Arial'
fs = 24

colors_plot = {'a':'k', 'c':'g', 'g':'brown', 't':'cyan'}
comp_nts = {'a':'t', 't':'a', 'g':'c', 'c':'g'}
bf = pd.read_csv('../figure9/cary_uv/m1a_25c_new.csv')
contexts_exclude = ['aag', 'gac', 'tag', 'cac']
bf = bf[~bf['context'].isin(contexts_exclude)]


def plot_contexts(bf, ax, middle_letter, contexts_exclude, yticks, fmt, ylab):
    counter = 1

    sequence_identifiers = []
    sequence_identifiers_new = []
   
    for start_letter in ['a', 'c', 'g', 't']:
        #sequence_identifier = "5'-" + start_letter + middle_letter + "_" + "-3'"
        #sequence_identifier_new = "5'-" + "_" + middle_letter + start_letter + "-3'"
        sequence_identifier = start_letter + middle_letter + "_" 
        sequence_identifier_new = "_" + middle_letter + start_letter
        sequence_identifiers.append(sequence_identifier)
        sequence_identifiers_new.append(sequence_identifier_new)

        energies = []
        energies_error = []
        energies_color = []
        labels = []
        energies_new = []
        energies_error_new = []
        energies_color_new = []
        labels_new = []
 
        for end_letter in ['a', 'c', 'g', 't']:
            sequence = start_letter + middle_letter + end_letter
            sequence_new = end_letter + middle_letter + start_letter
            if sequence not in contexts_exclude:
                bf_new = bf[bf['context'] == sequence]
                energies.append(bf_new['deldelG_25C'].values[0])
                energies_error.append(bf_new['deldelG_25C_error'].values[0])
                energies_color.append(colors_plot[end_letter])
                labels.append(end_letter)
                print sequence, energies[-1], energies_error[-1]
            if sequence_new not in contexts_exclude:
                bf_new = bf[bf['context'] == sequence_new]
                energies_new.append(bf_new['deldelG_25C'].values[0])
                energies_error_new.append(bf_new['deldelG_25C_error'].values[0])
                energies_color_new.append(colors_plot[end_letter])
                print sequence, energies_new[-1], energies_error_new[-1]
                labels_new.append(end_letter)
        ms = 12
        lw = 4 
        cs = 4
        # Plot
        for dummy in range(len(energies)):
            ax[0].plot(counter, energies[dummy], color=energies_color[dummy], label=labels[dummy]) 
            ax[0].errorbar(counter, energies[dummy], yerr=energies_error[dummy], fmt='o', color=energies_color[dummy], elinewidth=lw, capthick=lw, markeredgecolor=energies_color[dummy], markersize=ms, zorder=dummy*50, capsize=cs) 
        for dummy in range(len(energies_new)):
            ax[1].plot(counter, energies_new[dummy], color=energies_color_new[dummy], label=labels_new[dummy]) 
            ax[1].errorbar(counter, energies_new[dummy], yerr=energies_error_new[dummy], fmt='o', color=energies_color_new[dummy], elinewidth=lw, capthick=lw, markeredgecolor=energies_color_new[dummy], markersize=ms, zorder=dummy*50, capsize=cs)
        counter = counter + 1
        energies = []
        energies_error = []
        energies_color = []
        energies_new = []
        energies_error_new = []
        energies_color_new = []
        ax[0].set_xlim([0.0, 5.0])
        ax[1].set_xlim([0.0, 5.0])
    print sequence_identifiers
    print sequence_identifiers_new
    sequence_identifiers = [ele.upper() for ele in sequence_identifiers]
    sequence_identifiers_new = [ele.upper() for ele in sequence_identifiers_new]
    ax[0].set_xticks(range(1, 5))
    ax[0].set_xticklabels(sequence_identifiers, fontsize=fs)
    ax[1].set_xticks(range(1, 5))
    ax[1].set_xticklabels(sequence_identifiers_new, fontsize=fs)
    ax[0].tick_params(axis='x', direction='out', which='minor', length=0.0, width=2.0)
    ax[0].tick_params(axis='x', direction='out', which='major', length=6.0, width=2.0, top=False)
    ax[0].tick_params(axis='y', direction='out', which='minor', length=0.0, width=2.0)
    ax[0].tick_params(axis='y', direction='out', which='major', length=6.0, width=2.0, right=False)
    ax[1].tick_params(axis='x', direction='out', which='minor', length=0.0, width=2.0)
    ax[1].tick_params(axis='x', direction='out', which='major', length=6.0, width=2.0, top=False)
    ax[1].tick_params(axis='y', direction='out', which='minor', length=0.0, width=2.0)
    ax[1].tick_params(axis='y', direction='out', which='major', length=6.0, width=2.0, right=False)
    ax[0].set_yticks(yticks) 
    ax[0].set_yticklabels([fmt%ele for ele in yticks], fontsize=fs)
    ax[1].set_yticks(yticks) 
    ax[1].set_yticklabels([fmt%ele for ele in yticks], fontsize=fs)

    # Legend for 0
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    #ax[0].legend(by_label.values(), by_label.keys())
    labels = by_label.values()
    handles = by_label.keys()
    leg = ax[0].legend(by_label.values(), by_label.keys(), title='', fancybox=False, ncol=4, handlelength=0, frameon=False, columnspacing=0.0, markerscale=0.0, handletextpad=0.2, borderpad=0.0, handleheight=0.0, labelspacing=0.2, fontsize=fs, )
    leg2 = ax[1].legend(by_label.values(), by_label.keys(), title='', fancybox=False, ncol=4, handlelength=0, frameon=False, columnspacing=0.0, markerscale=0.0, handletextpad=0.2, borderpad=0.0, handleheight=0.0, labelspacing=0.2, fontsize=fs, )
    counter = 0 
    for text in leg.get_texts():
        text.set_color(colors_plot[text.get_text()])
        text.set_text(text.get_text().upper())
        counter = counter + 1
    for text in leg2.get_texts():
        text.set_color(colors_plot[text.get_text()])
        text.set_text(text.get_text().upper())
        counter = counter + 1

    for l in leg.get_lines():
        l.set_linestyle('None')
    for l in leg2.get_lines():
        l.set_linestyle('None')
    #leg.set_title('$\omega_{1}/2\pi$' + ' ' + '$(Hz)$', prop={'size':fs})
    ax[0].set_ylabel(ylab, fontsize=fs)
    ax[0].set_xlabel('Sequence Context', fontsize=fs)
    ax[1].set_xlabel('Sequence Context', fontsize=fs)

fig, ax = plt.subplots(1, 2, figsize=(14, 6))
plot_contexts(bf, ax, 'a', contexts_exclude, np.arange(2.0, 6.1, 2.0), '%d', 'A-T Hoogsteen\n' + '$\Delta\Delta G_{delta-melt, 25^{\\circ}C}^{\\circ}$'+ ' (kcal/mol)')

plt.savefig('figs14a_2.pdf')
#plt.show()

