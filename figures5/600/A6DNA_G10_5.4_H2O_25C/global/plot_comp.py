import os
import math
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def draw_subplot(dirname, ax):
    #rd_data_net = pd.read_csv(dirname + "2state/Local/Data_gc1p_1.csv")
    #params = pd.read_csv(dirname + "2state/LocalFits_gc1p.csv")
    rd_data_net = pd.read_csv(dirname + "2state/Local/Data_gc8_1.csv")
    params = pd.read_csv(dirname + "2state/LocalFits_gc8.csv")
    color_slp = {200.0: 'green', 300.0:'yellow', 400.0:'orange', 900.0:'red', 1000.0:'cyan'}
    color_slp = color_slp.values()

    # begin plotting
    counter = 0
    for dummy_slp in np.unique(rd_data_net['SLP'].values):
        print dummy_slp

        data = rd_data_net.loc[rd_data_net['SLP']==dummy_slp]
        ax.plot(data['Offset'], data['R2eff'], color='k', label=str(dummy_slp), linewidth=0.0, marker='o', markersize=3)
        ax.errorbar(data['Offset'], data['R2eff'], yerr=data['R2eff err'], fmt='o', color=color_slp[counter])

        counter = counter + 1
        ax.plot(data['Offset'], data[' Sim R2eff'], color=color_slp[counter], label=str(dummy_slp))
    r = params['kexAB'].iloc[0] / (params['dwB'].iloc[0] * params['lf'].iloc[0] * 2 * math.pi)

    #title = "$k_{ex}$=" + str(int(params['kexAB'].iloc[0])) + " " + "$p_b$=" + " " + "%.3f"%params['pB'].iloc[0] + " " + "R=" + "%.1f"%r
    #title = "$R_{1}$=" + '%.3f'%params['R1'].iloc[0] + " " + "$R_2$=" + " " + "%.3f"%params['R2'].iloc[0] 
    #title = "$R_{1}$=" + '%.3f'%params['R1'].iloc[0] + " " + "$dw$=" + " " + "%.3f"%params['dwB'].iloc[0] 
    title = "$k_{ex}$=" + str(int(params['kexAB'].iloc[0])) + " " + "$p_b$=" + " " + "%.3f"%(params['pB'].iloc[0]*100.0) + " " + "X2=" + '%3.2f'%params['RedChiSq']
    ax.set_title(title, fontsize=8)
    #ax.set_xticks([ele for ele in list(range(-3000, 3001, 1000))])
    #ax.set_xticklabels([str(int((ele+0.0)/1000.0)) for ele in range(-3000, 3001, 1000)])
    #ax.xlim([-1000.0,1000.0])
    #ax.set_xlabel('$\Omega/2\pi$' + '(kHz)')
    #ax.plot([0, 0], [15.0, 100.0], color='k')
    #ax.ylim([15.0, 100.0])
    ax.set_xticks([]) 
   
fig, axarr = plt.subplots(3, 6)
dirname_list = [0.1, 0.20, 0.30, 0.50, 0.6, 0.70, 0.8, 0.90, 1.0, 1.3, 1.6, 2.0, 2.5, 3.0, 3.5, 5.0, 8.0, 10.0] 
horiz_counter = 0
vert_counter = 0
for dirnum in dirname_list:
    if horiz_counter == 6:
        vert_counter = vert_counter + 1
        horiz_counter = 0
        plt.legend()
    dirname = 'pb_' + str(dirnum) + 'pct/'
    print dirname, vert_counter, horiz_counter
    draw_subplot(dirname, axarr[vert_counter, horiz_counter])
    horiz_counter = horiz_counter + 1

#draw_subplot(plt, 'base/', 'k')
#draw_subplot(plt, 'base_nocomplex/', 'r')

#draw_subplot(plt, 'pb_2_kex_1000/', 'g')

# Get largest difference in r2+rex
#bf1 = pd.read_csv('ls_sim/sim-r1p.csv')
#bf2 = pd.read_csv('bm_sim_gs/sim-r1p.csv')
#print "max eror in r2eff = ", np.amax(bf1['r2eff'].values - bf2['r2eff'].values)
#print "max eror in r1p = ", np.amax(bf1['r1p'].values - bf2['r1p'].values)

#    axarr[4, 0].set_xticks([str(ele) for ele in range(4001)])
#    axarr[dummy_row, dummy_column].set_xticklabels([str(int((ele+0.0)/1000.0)) for ele in range(-4000, 4001, 1000)])
plt.show()
