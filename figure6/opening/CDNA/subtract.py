import sys
import pandas as pd
from uncertainties import ufloat
from uncertainties import unumpy

bf1 = pd.read_csv('wt/test.csv')
delh1 = ufloat(bf1['delH'].iloc[-2], bf1['delH'].iloc[-1])
dels1 = ufloat(bf1['delS'].iloc[-2], bf1['delS'].iloc[-1])
f = open('chartt.csv', 'w')
f.write('resi,ddh(kcal/mol),ddh_error(kcal/mol),dds(cal/mol/k),dds_error(cal/mol/k)\n')

for dummy in ['m3T3', 'm3T5', 'm3T16', 'm3T17', 'm3T18', 'm3T19', 'm3T21']:
    filename2 = dummy + '/' + 'test.csv'
    bf2 = pd.read_csv(filename2)
    delh2 = ufloat(bf2['delH'].iloc[-2], bf2['delH'].iloc[-1])
    dels2 = ufloat(bf2['delS'].iloc[-2], bf2['delS'].iloc[-1])
    ddh = delh2 - delh1
    dds = (dels2 - dels1)
    print dummy
    print 'delH mutant - wt = ', unumpy.nominal_values(ddh), unumpy.std_devs(ddh) 
    print 'delS mutant - wt = ', unumpy.nominal_values(dds), unumpy.std_devs(dds)
    print
    line = str(int(dummy[3:])) + "," + str(unumpy.nominal_values(ddh)) + "," + str(unumpy.std_devs(ddh)) + "," + str(unumpy.nominal_values(dds)) + "," + str(unumpy.std_devs(dds)) + "\n" 
    f.write(line)

f.close() 
