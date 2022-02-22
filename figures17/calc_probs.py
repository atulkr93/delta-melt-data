import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from uncertainties import unumpy
import sys

# Define nt complementarity
ctot_contexts = []
corr_nts = {'a':'t', 't':'a', 'g':'c', 'c':'g'}

bf = pd.read_csv('m1a_37c_new.csv')
for dummy_context in list(bf['context'].values):
    print dummy_context
    ctot_contexts.append(corr_nts[dummy_context[2]].upper() + corr_nts[dummy_context[1]].upper() + corr_nts[dummy_context[0]].upper())
   
# convert energies to probabilities
T = 273.16 + 37.0
eng = unumpy.uarray(bf['deldelG_37C'].values, bf['deldelG_37C'].values)
print eng
weights = unumpy.exp(-1*eng/(8.314*T/4184))
weights = weights / np.sum(weights)
print np.sum(unumpy.nominal_values(weights))

fout = open('m1a_37c_new_probs.csv', 'w')
fout.write('sequence,prob,prob_error\n')
for dummy in range(len(weights)):
    fout.write(str(ctot_contexts[dummy]) + "," + str(unumpy.nominal_values(weights[dummy])) + "," + str(unumpy.std_devs(weights[dummy])) + "\n")
fout.close()
