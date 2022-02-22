from matplotlib import rcParams
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from uncertainties import unumpy 
import sys

fs=20
wd=2
lw=4
rcParams['axes.linewidth'] = 2.0
rcParams['patch.linewidth'] = 2.0
rcParams['font.family']='Sans-serif'
rcParams['font.sans-serif']=['Arial']

# Define sequence contexts from CHARTT
c_contexts = []

# exclude artefact signatures
signatures_exclude = ['SBS27','SBS43','SBS45','SBS46','SBS47','SBS48','SBS49','SBS50','SBS51','SBS52','SBS53','SBS54','SBS55','SBS56','SBS57','SBS58','SBS59','SBS60']

# Read chartt data
bf = pd.read_csv('m1a_37c_new_probs.csv')

# get order of interest
for dummy_outer in ['A', 'C', 'G', 'T']:
    for dummy_inner in ['A', 'C', 'G', 'T']:
        sequence_context = dummy_outer + 'T' + dummy_inner
        # exclude AGG and TGG as UV is non-two state
        if sequence_context not in ['CTT', 'GTC', 'CTA', 'GTG']:
            c_contexts.append(sequence_context) 
        
# get signatures
f = open('corr_athg_cosmic.csv', 'w')
f.write('signature,pearson_r\n')
bf2 = pd.read_csv('COSMIC_Mutational_Signatures_v3.1_sbs.csv')

# get signature names
signatures = [ele for ele in list(bf2) if ele[:3] == "SBS"]
signatures = [ele for ele in signatures if ele not in signatures_exclude]
print len(signatures)

filenames_net = []
counter = 1

# Loop over signatures
for dummy_sig in signatures:
    #################################################
    # First parse the C>T part of the signature
    sig_wts = []
    contexts_exclude = []
    c_contexts_new = []
    flag = 0
    for ele in c_contexts:
        bf2_subset = bf2.loc[((bf2['Type'] == "T>C") & (bf2['Subtype'] == ele))]
        if bf2_subset[dummy_sig].iloc[0] == ' -   ': 
            contexts_exclude.append(ele)
            flag = 1
        else:
            sig_wts.append(float(bf2_subset[dummy_sig].iloc[0]))
            c_contexts_new.append(ele)
        #print ele, float(bf2_subset[dummy_sig])

    # Get a corresponding chartt order
    # context_new already has exactly what we want
    prob_vals = []
    prob_vals_error = []
    for context in c_contexts_new:
        bf_subset = bf[bf['sequence'] == context]
        prob_vals.append(bf_subset['prob'].iloc[0])
        prob_vals_error.append(bf_subset['prob_error'].iloc[0]) 
    prob_vals_net = unumpy.uarray(prob_vals, prob_vals_error)
    
    sig_wts = np.array(sig_wts)
    sig_wts = sig_wts / np.sum(sig_wts)
    print dummy_sig
    for dummy in range(len(sig_wts)):
        print c_contexts_new[dummy], prob_vals_net[dummy], sig_wts[dummy]
    print

    #if flag == 1:
    #    break

    plt.figure(counter)
    counter = counter + 1
    plt.tick_params(axis='x', direction='out', top=False, length=lw, width=wd)
    plt.tick_params(axis='y', direction='out', right=False, length=lw, width=wd, labelsize=fs)
    plt.bar(range(len(c_contexts_new)), unumpy.nominal_values(prob_vals_net), width=0.4, color='r', linewidth=3, label='delta-melt', align='edge')
    plt.bar(np.array(range(len(c_contexts_new)))+0.4, sig_wts, width=0.4, color='k', linewidth=3, label=dummy_sig, align='edge')
    plt.xticks(np.array(range(len(c_contexts_new)))+0.4, [str(ele) for ele in c_contexts_new], fontsize=fs, rotation='vertical')
    plt.xlim([-0.2, len(c_contexts_new)-1+0.4+0.6])
    plt.legend(loc=0, fontsize=fs)

    ##plt.show()
    
    # Get pearson r
    print np.corrcoef(unumpy.nominal_values(prob_vals_net), sig_wts)[0,1]
    print np.corrcoef(sig_wts, unumpy.nominal_values(prob_vals_net))[0,1]
    r = np.corrcoef(unumpy.nominal_values(prob_vals_net), sig_wts)[0,1]

    title = plt.title("T>C " + str(dummy_sig) + " vs. AT Hoogsteen, r=" + '%0.2f'%r, fontsize=fs)
    title.set_position([0.5, 1.02])

    plt.savefig(dummy_sig + "_TtoC.pdf")
    filenames_net.append(dummy_sig + "_TtoC.pdf")
    f.write(dummy_sig+"_TtoC.pdf," + str(r) + "\n")

    #################################################
    # First parse the T>G part of the signature
    sig_wts = []
    contexts_exclude = []
    c_contexts_new = []
    flag = 0
    for ele in c_contexts:
        bf2_subset = bf2.loc[((bf2['Type'] == "T>G") & (bf2['Subtype'] == ele))]
        if bf2_subset[dummy_sig].iloc[0] == ' -   ': 
            contexts_exclude.append(ele)
            flag = 1
        else:
            sig_wts.append(float(bf2_subset[dummy_sig].iloc[0]))
            c_contexts_new.append(ele)
        #print ele, float(bf2_subset[dummy_sig])

    # Get a corresponding chartt order
    # context_new already has exactly what we want
    prob_vals = []
    prob_vals_error = []
    for context in c_contexts_new:
        bf_subset = bf[bf['sequence'] == context]
        prob_vals.append(bf_subset['prob'].iloc[0])
        prob_vals_error.append(bf_subset['prob_error'].iloc[0]) 
    prob_vals_net = unumpy.uarray(prob_vals, prob_vals_error)
    
    sig_wts = np.array(sig_wts)

    sig_wts = sig_wts / np.sum(sig_wts)
    for dummy in range(len(sig_wts)):
        print c_contexts_new[dummy], prob_vals_net[dummy], sig_wts[dummy], sig_wts[dummy]/np.sum(sig_wts[dummy]), np.sum(sig_wts[dummy])
    print
    print dummy_sig

    #if flag == 1:
    #    break

    plt.figure(counter)
    counter = counter + 1
    plt.tick_params(axis='x', direction='out', top=False, length=lw, width=wd)
    plt.tick_params(axis='y', direction='out', right=False, length=lw, width=wd, labelsize=fs)
    plt.bar(range(len(c_contexts_new)), unumpy.nominal_values(prob_vals_net), width=0.4, color='r', linewidth=3, label='delta-melt', align='edge')
    plt.bar(np.array(range(len(c_contexts_new)))+0.4, sig_wts, width=0.4, color='k', linewidth=3, label=dummy_sig, align='edge')
    plt.xticks(np.array(range(len(c_contexts_new)))+0.4, [str(ele) for ele in c_contexts_new], fontsize=fs, rotation='vertical')
    plt.xlim([-0.2, len(c_contexts_new)-1+0.4+0.6])
    plt.legend(loc=0, fontsize=fs)
    ##plt.show()
    
    # Get pearson r
    print np.corrcoef(unumpy.nominal_values(prob_vals_net), sig_wts)[0,1]
    print np.corrcoef(sig_wts, unumpy.nominal_values(prob_vals_net))[0,1]
    r = np.corrcoef(unumpy.nominal_values(prob_vals_net), sig_wts)[0,1]

    title = plt.title("T>G " + str(dummy_sig) + " vs. AT Hoogsteen, r=" + '%0.2f'%r, fontsize=fs)
    title.set_position([0.5, 1.02])
    plt.savefig(dummy_sig + "_TtoG.pdf")
    filenames_net.append(dummy_sig + "_TtoG.pdf")
    f.write(dummy_sig+"_TtoG.pdf," + str(r) + "\n")

    #################################################
    # First parse the C>G part of the signature
    sig_wts = []
    contexts_exclude = []
    c_contexts_new = []
    flag = 0
    for ele in c_contexts:
        bf2_subset = bf2.loc[((bf2['Type'] == "T>A") & (bf2['Subtype'] == ele))]
        if bf2_subset[dummy_sig].iloc[0] == ' -   ': 
            contexts_exclude.append(ele)
            flag = 1
        else:
            sig_wts.append(float(bf2_subset[dummy_sig].iloc[0]))
            c_contexts_new.append(ele)
        #print ele, float(bf2_subset[dummy_sig])

    # Get a corresponding chartt order
    # context_new already has exactly what we want
    prob_vals = []
    prob_vals_error = []
    for context in c_contexts_new:
        bf_subset = bf[bf['sequence'] == context]
        prob_vals.append(bf_subset['prob'].iloc[0])
        prob_vals_error.append(bf_subset['prob_error'].iloc[0]) 
    prob_vals_net = unumpy.uarray(prob_vals, prob_vals_error)
    
    sig_wts = np.array(sig_wts)

    sig_wts = sig_wts / np.sum(sig_wts)
    for dummy in range(len(sig_wts)):
        print c_contexts_new[dummy], prob_vals_net[dummy], sig_wts[dummy], sig_wts[dummy]/np.sum(sig_wts[dummy]), np.sum(sig_wts[dummy])
    print
    print dummy_sig

    #if flag == 1:
    #    break

    plt.figure(counter)
    plt.tick_params(axis='x', direction='out', top=False, length=lw, width=wd)
    plt.tick_params(axis='y', direction='out', right=False, length=lw, width=wd, labelsize=fs)
    counter = counter + 1
    plt.bar(range(len(c_contexts_new)), unumpy.nominal_values(prob_vals_net), width=0.4, color='r', linewidth=3, label='delta-melt', align='edge')
    plt.bar(np.array(range(len(c_contexts_new)))+0.4, sig_wts, width=0.4, color='k', linewidth=3, label=dummy_sig, align='edge')
    plt.xticks(np.array(range(len(c_contexts_new)))+0.4, [str(ele) for ele in c_contexts_new], fontsize=fs, rotation='vertical')
    plt.xlim([-0.2, len(c_contexts_new)-1+0.4+0.6])
    plt.legend(loc=0, fontsize=fs)
    ##plt.show()
    
    # Get pearson r
    print np.corrcoef(unumpy.nominal_values(prob_vals_net), sig_wts)[0,1]
    print np.corrcoef(sig_wts, unumpy.nominal_values(prob_vals_net))[0,1]
    r = np.corrcoef(unumpy.nominal_values(prob_vals_net), sig_wts)[0,1]

    title = plt.title("T>A " + str(dummy_sig) + " vs. AT Hoogsteen, r=" + '%0.2f'%r, fontsize=fs)
    title.set_position([0.5, 1.02])
    plt.savefig(dummy_sig + "_TtoA.pdf")
    filenames_net.append(dummy_sig + "_TtoA.pdf")
    f.write(dummy_sig+"_TtoA.pdf," + str(r) + "\n")


f.close()
command = "pdfunite "   
command2 = "rm "
for filename in filenames_net:
    command = command + str(filename) + " " 
    command2 = command2 + str(filename) + " " 
command = command + "out.pdf"
os.system(command)
os.system(command2)
