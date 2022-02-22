import os
import math
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

lf = 70.94

def draw_subplot(dirname, ax):
    global lf
    #rd_data_net = pd.read_csv(dirname + "2state/fit-result_hpGBCm6A_C2_37C.csv")
    #params = pd.read_csv(dirname + "2state/fitparams-result_hpGBCm6A_C2_37C.csv")
    rd_data_net = pd.read_csv(dirname + "2state/fit-result_methtlCEST_37C.csv")
    params = pd.read_csv(dirname + "2state/fitparams-result_methtlCEST_37C.csv")
    color_slp = {10.0: 'green', 25.0: 'black', 35.0:'k', 400.0:'orange', 800.0:'red', 800.0:'cyan'}
    # begin plotting
    counter = 0
    color_slp = color_slp.values()
    print color_slp
    for dummy_slp in rd_data_net['slp(hz)'].unique():
        data = rd_data_net.loc[rd_data_net['slp(hz)']==dummy_slp]
        ax.errorbar(data['offset(hz)'], data['norm_intensity'], yerr=data['norm_intensity_error'], color=color_slp[counter], marker='o', markersize=3, label=str(dummy_slp), linewidth=0.0)
        ax.plot(data['offset(hz)'], data['fit_norm_intensity'], color=color_slp[counter], linewidth=2)
        counter = counter + 1
        #ax.plot(data['offset'], data['r2eff'], color=color_slp[dummy_slp], label=str(dummy_slp))
    print params['fitval'].iloc[4], params['fitval'].iloc[2]
    r = params['fitval'].iloc[4] / (params['fitval'].iloc[2] * lf * 2 * math.pi)

    title = "$k_{ex}$=" + str(int(params['fitval'].iloc[4])) + " " + "$p_b$=" + " " + "%.2f"%(params['fitval'].iloc[0]*100.) + " " + "X2=" + '%3.2f'%params['fitval'].iloc[9]
    #title = "$k_{ex}$=" + str(int(params['fitval'].iloc[4])) + " " + "$p_b$=" + " " + "%.2f"%(params['fitval'].iloc[0]*100.) + " " + "dw=" + '%.1f'%params['fitval'].iloc[2] + ",R=" + '%.1f'%r
    #title = "$R_{1}$=" + '%.3f'%params['fitval'].iloc[7] + " " + "$R_2$=" + " " + "%.3f"%(params['fitval'].iloc[8])
    #title = "$R_{1}$=" + '%.3f'%params['fitval'].iloc[7] + " " + "$dw$=" + " " + "%.3f"%(params['fitval'].iloc[2]) + " " + "$R_2$=" + " " + "%.3f"%(params['fitval'].iloc[8])
    ax.set_title(title)
    #ax.set_xticks([ele for ele in list(np.arange(-10, 11, 5))])
    #ax.set_xticklabels([str(int(ele)) for ele in list(np.arange(-10,11,5))])
    ax.set_xlim([-1500,1500])
    #ax.set_xlabel('$\Omega/2\pi$' + '(kHz)')
    #ax.set_ylim([-0.05, 0.52])
    ax.set_xticks([0, 0])    

fig, axarr = plt.subplots(4, 4)
dirname_list = [0.01, 0.05, 0.1, 0.20, 0.30, 0.40, 0.50, 0.6, 0.75, 1.0, 2.0, 3.0, 5.0, 8.0] 

horiz_counter = 0
vert_counter = 0
for dirnum in dirname_list:
    if horiz_counter == 4:
        vert_counter = vert_counter + 1
        horiz_counter = 0
        plt.legend()
    dirname = 'pb_' + str(dirnum) + 'pct/'
    print dirname, vert_counter, horiz_counter
    draw_subplot(dirname, axarr[vert_counter, horiz_counter])
    horiz_counter = horiz_counter + 1



#    axarr[4, 0].set_xticks([str(ele) for ele in range(4001)])
#    axarr[dummy_row, dummy_column].set_xticklabels([str(int((ele+0.0)/1000.0)) for ele in range(-4000, 4001, 1000)])
plt.show()
