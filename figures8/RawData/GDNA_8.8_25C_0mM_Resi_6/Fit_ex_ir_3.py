#!/usr/bin/python

import os
import sys
import pandas as pd
import nmrglue as ng
import pylab as pl
import numpy as np
from scipy.optimize import least_squares
from scipy.optimize import curve_fit
from uncertainties import *

def fun_ex_ir(x,kex,R1n):
    return (1 - (E*kex/(R1w-R1n)*(np.exp(-R1n*x)-np.exp(-R1w*x))))


# Input parameters
Sample = 'GDNA'
Condition = '8.8_25C_0mM'
PeakNum, Resi = [1,6]
#Output = '%s_%s_Resi_%s_test'%(Sample,Condition,str(Resi))
Output = 'test.csv'
trunc = 0
p0 = [10.,10.]
E = 1.877588
R1w = 0.308
R1n = 2.89
kex = 0.01

#os.system("mkdir -p %s"%Output)

#vdlistpath = "vdlist"
## Check vdlist exists
#if not os.path.isfile(vdlistpath):
#  print "No vdlist found."
#  sys.exit(-1)

## Read in delays
#FILE = open(vdlistpath, "rU")
#delays = np.array([float(x.strip()) for x in FILE])
#FILE.close()


#intensity = []


#for i in range(len(delays)):
    #tabpath = "fit_tab/T1_1DFit-%03d.tab"%(i+1)
    #pc,pf,rec = ng.fileio.pipe.read_table(tabpath)
    
    #height = rec[PeakNum-1][9]
    #vol =  rec[PeakNum-1][11]
    #print vol
    
    #intensity.append(float(vol))
    #intensity.append(float(height))

#intensity = np.array(intensity)
#intensity = intensity / intensity[0]

#if trunc != 0:  
#    delays = delays[:-trunc]
#    intensity = intensity[:-trunc]

bf = pd.read_csv(str(sys.argv[1]))
delays = np.array(bf['delay(s)'].values)
intensity = np.array(bf['vol_rel'].values)

popt, pcov = curve_fit(fun_ex_ir,delays,intensity,bounds=([0.,0.],[p0[0],p0[1]]), absolute_sigma=False)
kex = ufloat(popt[0],np.sqrt(np.diag(pcov))[0])
R1n = ufloat(popt[1],np.sqrt(np.diag(pcov))[1])
tau_ex = 1./kex

print("kex value: " + str(kex) + " (s-1)")
print("tau_ex value: " + str(tau_ex) + " (s)")
print("R1n value: " + str(R1n) + " (s-1)")

pl.plot(delays,intensity,'ok',markersize=8,label='raw data')
taulist = np.arange(0.0,np.max(delays),0.001)
pl.plot(taulist, 1 - E*popt[0]/(R1w-popt[1])*(np.exp(-popt[1]*taulist)-np.exp(-R1w*taulist)), '-r', linewidth=2, label='fitting')

pl.xlabel("Delay (s)")
pl.ylabel("Normed Volume")

#pl.legend(loc=4)
#pl.savefig("%s/profile_%s.pdf"%(Output,Output))

pl.show()

#output_df1 = pd.DataFrame([],columns=['delay(s)','vol_rel'])
#output_df1['delay(s)'] = pd.Series(delays)
#output_df1['vol_rel'] = pd.Series(intensity)
#output_df1.to_csv('%s/data_%s.csv'%(Output,Output),index=False)

#output2 = [['param','fit_val','fit_err']]
#output2.append(['kex(s-1)',kex.n,kex.s])
#output2.append(['tau_ex(s)',tau_ex.n,tau_ex.s])
#output2.append(['R1n(s-1)',R1n.n,R1n.s])

#output_df2 = pd.DataFrame(output2[1:],columns=output2[0])
#output_df2.to_csv('%s/fitparam_%s.csv'%(Output,Output),index=False,float_fmt='%.3f')



