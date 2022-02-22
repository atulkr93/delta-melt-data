# Calculate delH and delS from Chen and Russu, 2004
# given delGs and delHs from Fig. 6 in the paper
import pandas as pd
from uncertainties import unumpy

# Fig. 6A
bf1 = pd.read_csv('TBP_base_opening_extraction.csv')
# Fig. 6B
bf2 = pd.read_csv('TBP_base_opening_dH_extraction.csv')

delg = unumpy.uarray(bf1['dG(kcal/mol)'].values, bf1['dG_err(kcal/mol)'].values)
delh = unumpy.uarray(bf2['dH(kcal/mol)'].values, bf2['dH_err(kcal/mol)'].values)
dels = (delh - delg)  / 288.16
print delg
print delh
print dels

bf1.insert(len(list(bf1)), 'delh(kcal/mol)', unumpy.nominal_values(delh))
bf1.insert(len(list(bf1)), 'delh_err(kcal/mol)', unumpy.std_devs(delh))
bf1.insert(len(list(bf1)), 'dels(kcal/mol)', unumpy.nominal_values(dels))
bf1.insert(len(list(bf1)), 'dels_err(kcal/mol)', unumpy.std_devs(dels))
bf1.to_csv('tbp_opening_delh_dels.csv')

