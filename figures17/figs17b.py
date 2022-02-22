import os

# Get the free energies for formation of 
# AT Hoogsteen bps at 37C using delta-melt
os.system('python athg_37c.py')

# Convert these free energies to probabilities
# for HG bp formation, assuming  all sequence
# contexts are equally abundant
os.system('python calc_probs.py')

# Then correlate these probabilities to cosmic
# cancer signaturs, plots in out.pdf
os.system('python corr_athg_cosmic.py')
