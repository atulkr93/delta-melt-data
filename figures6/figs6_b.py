# Plot the degeneracy
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import matplotlib as mpl
import os

mpl.rcParams['axes.linewidth'] = 4.

pb_master = [0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.45, 0.5, 0.75, 1.0, 2.0, 3.0, 5.0, 8.0]
fs = 36

r1p_simulation = '''
##################################################################################
# Run the BMNS simulation routine:
# > python BMNS.py -sim [BM Simulation Input File] (Optional Output directory)
##################################################################################
# "Params" block is where simulation parameters are defined.
#   Parameters can be defined manually, or read-in from a BM-fit CSV file
# If you read in fit CSV you can then manually define some parameters,
#   this will overwrite the parameters read in from CSV.
#---------------------------------------------------------------------------------
# - 'Read' reads in a BM fit csv file from local directory
#     only uses kexAB, kexAC, kexBC for exchange rates, not indv rate constants
#   This can replace or complement the individually specified values
# - 'lf' is larmor frequency of spin (e.g. 150.8MHz for 13C, etc)
# - 'AlignMag' is auto/gs/avg for automatic, ground-state, or average alignment
# - 'pB/pC' are populations of state B/C
# - 'kexAB/AC/BC' are exchange rates between the 3 states A,B, and C
# - 'R1/R2' are the intrinsic relax rates of states A/B/C
# Any of the above parameters can be commented out provided they are read-in
#  from a CSV file instead.
##################################################################################
+
Params
#Read pars.csv
lf %f
AlignMag Auto
pB %f 
pC 0.0
dwB %f
dwC 0.0
kexAB %f
kexAC 0.0
kexBC 0.0
R1 %f
R2 %f
R1b %f
R2b %f
R1c %f
R2c %f

##################################################################################
# "SLOFF" block defines spinlock powers and offsets to simulate with BM.
# Additionally, real data can be read in to be overlaid with simulated data.
# Additionally, simulated data can be error corrupted at the level of the
#     R1rho value to a certain percentage. Monte-Carlo iterations can be
#   defined for this error corruption.
#---------------------------------------------------------------------------------
# - 'Read' defines a .csv file that can be read in
#       that contains Offset(Hz), SLP(Hz) in columns
#       and these will be simulated. If commented out
#     then they will not be read in.
# - 'Data' defines a .csv file that contains real data.
#      Can directly read in data file generated by
#     the BM fitting routine.
#   Order is:
#    Col1: Corrected Offset(Hz)
#    Col2: SLP (Hz)
#    Col3: R1rho
#    Col4: R1rho error (must exist, even if all zeros)
#    Col5: R2eff
#    Col6: R2eff error (must exist, even if all zeros)
#   If not defined, then they will not be read in.
# - 'Error' noise corruption for simulated R1rho values.
#   e.g. 0.02 would error corrupt R1rho by 2%%
#   Generates error from MC error corruption, selecting
#   sigma from gaussian distribution around corrupted
#   R1rho value
# - 'MCNum' defines number of MC iterations for noise corruption.
# - 'on' Defines on-resonance R1rho values to be simulated
#  Add as many of these as you wish
#     Col1: 'on'
#   Col2: Lower SLP (Hz)
#   Col3: Upper SLP (Hz)
#   Col4: Number of onres SLPs to simulate between low/high
# - 'off' defines off-resonance R1rho values to be simulated
#   at a given SLP over a range of offsets.
#    Add as many 'off' rows as you need to generate more
#   more off-resonance points or spinlock powers
#     Col1: 'off'
#   Col2: SLP (Hz)
#   Col3: Lower Offset (Hz)
#   Col4: Upper Offset (Hz)
#   Col5: Number of offres SLPs to simulate between low/high
##################################################################################
+
SLOFF
#Read sloffs.csv
#Data data.csv
Error 0.0
MCNum 500
off %f %f %f 200
off %f %f %f 200
off %f %f %f 200
off %f %f %f 200
off %f %f %f 200

##################################################################################
# "Decay" block defines how each R1rho value is simulated by simulating decaying
#   magnetization as a function of time given parameters describing the chemical
#   exchange between 2/3 species.
# Decay is assumed to be monoexponential, and simulated R1rho values are given
#   by the monoexponential fit of decaying magnetization.
# Note: This assumption can be violated under some conditions, where decay
#       can be bi-exponential or other (not take in to account).
# Additionally, intensity values can be noise corrupted to give a noise
#   corrupted R1rho value originating from N-number of corrupted monoexponential
#   fits. This is approximating how we derive R1rho experimentally and its error.
#---------------------------------------------------------------------------------
# - 'vdlist' a number of delay points to simulate decaying magnetization over.
#     Col2: Lowest delay in seconds (usually 0)
#   Col3: Longest delay in seconds (>0.1 is good, but can use anything)
#   Col4: Number of delays between low and high
# - 'Read' defines a delay list to read in. This is any text file where each row
#   is a delay given in seconds (e.g. vdlist).
#   If commented out, it will not be read in. If given, it will be comined with
#   delay values simulated with the 'vdlist' command below.
# - 'PlotDec' can be 'yes' or 'no'. If 'yes', then it will plot the
#   simulated decay for each SLP/offset combination along with
#   the best-fit line for the R1rho value at that point.
#   WARNING: This can take a long time if you are simulating lots of data
# - 'Error' defines noise corruption value for simulated magnetization
#   at each time point. E.g. 0.02 would be 2%% noise corruption.
#   Error here translates to error in R1rho by simulating N=MCNum of
#     noise-corrupted monoexponential decays and fitting them and
#     calculating the error in R1rho from the distribution of fitted
#     R1rhos (error = sigma of gaussian distribution of fitted R1rhos)
# - 'MCNum' defines how many noise-corrupted decays to simulate and fit.
#     WARNING: This can take a long time if you are simulating a lot of data.
##################################################################################
+
Decay
vdlist 0.0 0.25 51
#Read delays
PlotDec no
Error 0.0
MCNum 500

##################################################################################
# "Plot" block lets you specify how to plot your simulated/real data.
#---------------------------------------------------------------------------------
# - 'Plot' can be 'line', 'symbol', or 'both'.
#   'Line' will plot a simulated line of R1rho values
#   'Symbol' will plot simulated R1rhos as symbol types defined below
#   'Both' with plot symbols over simulated lines
# - 'Line' defines the style of the line plot.
#   Col2: Type of line, see: 
#   http://matplotlib.org/examples/lines_bars_and_markers/line_styles_reference.html
#      -   -.  --  or  : 
#     Col3: Line-width, in pixels
# - 'Symbol' defines the style of the symbol plot.
#   Col2: Type of symbol, see: http://matplotlib.org/api/markers_api.html
#     Too many to list, but default is a filled circle: o
#   Col3: Size of symbol (pixels)
# - 'Overlay' defines how you plot data overlaid on simulation
# - 'OType' type of data to overlay, real or overlay.
# - 'OLine' line type for overlay
# - 'OSymbol' symbol type for overlay
# - 'Size' defines the plot width and height in inches
# - '(R1p/R2eff/On)_x/y' define the lower and upper limits of the respective axes
#   Comment out to let them be automatically defined.
#   Alternatively, set one or both values to 'None' to let the program
#   automatically define the limit of the lower/upper bounds, individually
#   e.g. 'R1p_x None 1000' would let lower x-axis limits be automatically
#   defined, but the upper limit would be set to 1000
# - 'Axis_FS' sets the axes numbers font sizes, X and Y axes, respectively
# - 'LabelFS' sets the X and Y axes labels font sizes
##################################################################################
+
Plot line
Line - 2
Symbol o 13
Overlay line
OType sim
OLine -- 2
OSymbol . 13
Size 10 8
#R1p_x None 1000
#R1p_y 0 100
#R2eff_x -1000 1000
#R2eff_y 0 100
On_x 0 None
#On_y 0 50
Axis_FS 32 32
Label_FS 32 32
Labels on
'''

cest_simulation = '''
##################################################################################
# Run the BMNS simulation routine:
# > python BMNS.py -sim [BM Simulation Input File] (Optional Output directory)
##################################################################################
# "Params" block is where simulation parameters are defined.
#   Parameters can be defined manually, or read-in from a BM-fit CSV file
# If you read in fit CSV you can then manually define some parameters,
#   this will overwrite the parameters read in from CSV.
#---------------------------------------------------------------------------------
# - 'Read' reads in a BM fit csv file from local directory
#     only uses kexAB, kexAC, kexBC for exchange rates, not indv rate constants
#   This can replace or complement the individually specified values
# - 'lf' is larmor frequency of spin (e.g. 150.8MHz for 13C, etc)
# - 'AlignMag' is auto/gs/avg for automatic, ground-state, or average alignment
# - 'pB/pC' are populations of state B/C
# - 'kexAB/AC/BC' are exchange rates between the 3 states A,B, and C
# - 'R1/R2' are the intrinsic relax rates of states A/B/C
# Any of the above parameters can be commented out provided they are read-in
#  from a CSV file instead.
##################################################################################
+
lf %f
pb %f
pc 0.0
dwb %f
dwc 0.0
kexAB %f
kexAC 0.0
kexBC 0.0
R1a %f
R2a %f
R1b %f
R2b %f
R1c %f
R2c %f
resn 1.0
T %f
mode %s
error_point 0.000
error_baseinten 0.00
inhomo %f
equil Y
ls N
J %f

##################################################################################
# "SLOFF" block defines spinlock powers and offsets to simulate with BM.
# Additionally, real data can be read in to be overlaid with simulated data.
# Additionally, simulated data can be error corrupted at the level of the
#     R1rho value to a certain percentage. Monte-Carlo iterations can be
#   defined for this error corruption.
#---------------------------------------------------------------------------------
# - 'Read' defines a .csv file that can be read in
#       that contains Offset(Hz), SLP(Hz) in columns
#       and these will be simulated. If commented out
#     then they will not be read in.
# - 'Data' defines a .csv file that contains real data.
#      Can directly read in data file generated by
#     the BM fitting routine.
#   Order is:
#    Col1: Corrected Offset(Hz)
#    Col2: SLP (Hz)
#    Col3: R1rho
#    Col4: R1rho error (must exist, even if all zeros)
#    Col5: R2eff
#    Col6: R2eff error (must exist, even if all zeros)
#   If not defined, then they will not be read in.
# - 'Error' noise corruption for simulated R1rho values.
#   e.g. 0.02 would error corrupt R1rho by 2%%
#   Generates error from MC error corruption, selecting
#   sigma from gaussian distribution around corrupted
#   R1rho value
# - 'MCNum' defines number of MC iterations for noise corruption.
# - 'on' Defines on-resonance R1rho values to be simulated
#  Add as many of these as you wish
#     Col1: 'on'
#   Col2: Lower SLP (Hz)
#   Col3: Upper SLP (Hz)
#   Col4: Number of onres SLPs to simulate between low/high
# - 'off' defines off-resonance R1rho values to be simulated
#   at a given SLP over a range of offsets.
#    Add as many 'off' rows as you need to generate more
#   more off-resonance points or spinlock powers
#     Col1: 'off'
#   Col2: SLP (Hz)
#   Col3: Lower Offset (Hz)
#   Col4: Upper Offset (Hz)
#   Col5: Number of offres SLPs to simulate between low/high
##################################################################################
+
off %f %f %f 250
off %f %f %f 250

#off 35.0 -4000.0 4000.0 300
#off 75.0 -4000.0 4000.0 300
#off 400.0 -4000.0 4000.0 300
#off 500.0 -4000.0 4000.0 300
#off 35.0 -1000.0 1000.0 200
#off 17.68 -1412.2 567.8 670
#off 27.90 -1422.2 577.8 510
#off 48.21 -1422.2 577.8 410


##################################################################################
# "Decay" block defines how each R1rho value is simulated by simulating decaying
#   magnetization as a function of time given parameters describing the chemical
#   exchange between 2/3 species.
# Decay is assumed to be monoexponential, and simulated R1rho values are given
#   by the monoexponential fit of decaying magnetization.
# Note: This assumption can be violated under some conditions, where decay
#       can be bi-exponential or other (not take in to account).
# Additionally, intensity values can be noise corrupted to give a noise
#   corrupted R1rho value originating from N-number of corrupted monoexponential
#   fits. This is approximating how we derive R1rho experimentally and its error.
#---------------------------------------------------------------------------------
# - 'vdlist' a number of delay points to simulate decaying magnetization over.
#     Col2: Lowest delay in seconds (usually 0)
#   Col3: Longest delay in seconds (>0.1 is good, but can use anything)
#   Col4: Number of delays between low and high
# - 'Read' defines a delay list to read in. This is any text file where each row
#   is a delay given in seconds (e.g. vdlist).
#   If commented out, it will not be read in. If given, it will be comined with
#   delay values simulated with the 'vdlist' command below.
# - 'PlotDec' can be 'yes' or 'no'. If 'yes', then it will plot the
#   simulated decay for each SLP/offset combination along with
#   the best-fit line for the R1rho value at that point.
#   WARNING: This can take a long time if you are simulating lots of data
# - 'Error' defines noise corruption value for simulated magnetization
#   at each time point. E.g. 0.02 would be 2%% noise corruption.
#   Error here translates to error in R1rho by simulating N=MCNum of
#     noise-corrupted monoexponential decays and fitting them and
#     calculating the error in R1rho from the distribution of fitted
#     R1rhos (error = sigma of gaussian distribution of fitted R1rhos)
# - 'MCNum' defines how many noise-corrupted decays to simulate and fit.
#     WARNING: This can take a long time if you are simulating a lot of data.
##################################################################################
+
Decay
vdlist 0.0 0.25 51
#Read delays
PlotDec no
Error 0.0
MCNum 500

##################################################################################
# "Plot" block lets you specify how to plot your simulated/real data.
#---------------------------------------------------------------------------------
# - 'Plot' can be 'line', 'symbol', or 'both'.
#   'Line' will plot a simulated line of R1rho values
#   'Symbol' will plot simulated R1rhos as symbol types defined below
#   'Both' with plot symbols over simulated lines
# - 'Line' defines the style of the line plot.
#   Col2: Type of line, see: 
#   http://matplotlib.org/examples/lines_bars_and_markers/line_styles_reference.html
#      -   -.  --  or  : 
#     Col3: Line-width, in pixels
# - 'Symbol' defines the style of the symbol plot.
#   Col2: Type of symbol, see: http://matplotlib.org/api/markers_api.html
#     Too many to list, but default is a filled circle: o
#   Col3: Size of symbol (pixels)
# - 'Overlay' defines how you plot data overlaid on simulation
# - 'OType' type of data to overlay, real or overlay.
# - 'OLine' line type for overlay
# - 'OSymbol' symbol type for overlay
# - 'Size' defines the plot width and height in inches
# - '(R1p/R2eff/On)_x/y' define the lower and upper limits of the respective axes
#   Comment out to let them be automatically defined.
#   Alternatively, set one or both values to 'None' to let the program
#   automatically define the limit of the lower/upper bounds, individually
#   e.g. 'R1p_x None 1000' would let lower x-axis limits be automatically
#   defined, but the upper limit would be set to 1000
# - 'Axis_FS' sets the axes numbers font sizes, X and Y axes, respectively
# - 'LabelFS' sets the X and Y axes labels font sizes
##################################################################################
+
Plot line
Line - 2
Symbol o 13
Overlay line
OType sim
OLine -- 2
OSymbol . 13
Size 10 8
#R1p_x None 1000
#R1p_y 0 100
#R2eff_x -1000 1000
#R2eff_y 0 100
On_x 0 None
#On_y 0 50
Axis_FS 32 32
Label_FS 32 32
Labels on
'''


def plot_r1p(ax, bf, trendlines, xaxis_lims, lim, step, ylims, fitpars, titleadd):
    ''' plot data frame inside given axis object '''
    color_slps = ['b', 'cyan', 'g', 'orange', 'brown']
    unique_slps = np.sort(np.unique(bf['SLP']))
    counter = 0
    for unique_slp in unique_slps:
        bf_subset = bf.loc[bf['SLP'] == unique_slp]
        ax.errorbar(bf_subset['Offset'], bf_subset['R2eff'], yerr=bf_subset['R2eff err'], fmt='o', color=color_slps[counter], elinewidth=2, markersize=5, label=str(int(unique_slp)))
        #ax.plot(bf_subset['Offset'], bf_subset[' Sim R2eff'], color='k', marker='o', linewidth=0.0) 

        if len(trendlines) != 0:
            trendlines_sim = trendlines.loc[trendlines['slp'] == unique_slp]
            ax.plot(trendlines_sim['offset'], trendlines_sim['r2eff'], color=color_slps[counter], linewidth=2, label='')
     
        ax.tick_params(axis='x', direction='out', width=2, top=False)
        ax.tick_params(axis='y', direction='out', width=2, right=False)
    
        ax.set_xlim(xaxis_lims)
        ax.set_xticks([int(ele) for ele in np.arange(-1.*lim, lim+1, step)])
        ax.set_xticklabels([str(int(float(ele)/1000.)) for ele in np.arange(-1.*lim, lim+1, step)], fontsize=fs)
         
        ax.set_ylim([ylims[0], ylims[1]])
        ax.set_yticks([int(ele) for ele in np.arange(ylims[0], ylims[1]+1, ylims[2])])
        ax.set_yticklabels([str(int(ele)) for ele in np.arange(ylims[0], ylims[1]+1, ylims[2])], fontsize=fs)

        counter = counter + 1

    ax.set_title(titleadd + "$r\chi^2$" + " " + "%.2f"%fitpars['RedChiSq'].iloc[0], fontsize=fs, y=1.02)

    counter = 0
    handles, labels = ax.get_legend_handles_labels()
    handles = [h[0] if type(h) is not mpl.lines.Line2D else h for h in handles]
    leg = ax.legend(handles, labels, title='$\omega_{1}/2\pi$' + ' ' + '$(Hz)$', fancybox=False, ncol=2, handlelength=0, frameon=False, columnspacing=0.0, markerscale=0.0, handletextpad=0.2, borderpad=0.0, handleheight=0.0, labelspacing=0.2, fontsize=fs)

    for text in leg.get_texts():
        text.set_color(color_slps[counter])
        counter = counter + 1
    for l in leg.get_lines():
        l.set_linestyle('None')
    leg.set_title('$\omega_{1}/2\pi$' + ' ' + '$(Hz)$', prop={'size':fs})

def get_simulation_r1p(slps, offsets, fitpars, sim_number):
    ''' run a simulation '''
    # Make the net number of slps = 5
    if len(slps) == 4:
        slps = list(slps)
        slps.append(slps[-1]+200.)
        slps = np.array(slps)
    # Write out simulation file
    f = open('BMNS_Simparams.txt', 'w')
    f.write(r1p_simulation%(fitpars['lf'].values[0], fitpars['pB'].values[0], fitpars['dwB'].values[0], fitpars['kexAB'].values[0], fitpars['R1'].values[0], fitpars['R2'].values[0], fitpars['R1b'].values[0], fitpars['R2b'].values[0], fitpars['R1c'].values[0], fitpars['R2c'].values[0], slps[0], offsets[0], offsets[1], slps[1], offsets[0], offsets[1], slps[2], offsets[0], offsets[1], slps[3], offsets[0], offsets[1], slps[4], offsets[0], offsets[1])) 
    f.close()

    # Run bm
    os.system('python /Users/ar359/rot2/rd/progs/BMNS-master/BMNS.py -sim BMNS_Simparams.txt output'+str(sim_number))
    os.system('rm BMNS_Simparams.txt')


def get_simulation_cest(slps, offsets, fitpars, alignmag, lf, trelax, inhomo, J):
    ''' run a cest simulation '''
    # Write out simulation file
    f = open('BMNS_Simparams.txt', 'w')
    
    fit_params = [str(ele) for ele in list(fitpars['param'])]

    f.write(cest_simulation%(lf, fitpars['fitval'].iloc[fit_params.index('pB')], fitpars['fitval'].iloc[fit_params.index('dwB')], fitpars['fitval'].iloc[fit_params.index('kexAB')], fitpars['fitval'].iloc[fit_params.index('R1')], fitpars['fitval'].iloc[fit_params.index('R2')], fitpars['fitval'].iloc[fit_params.index('R1')], fitpars['fitval'].iloc[fit_params.index('R2')], fitpars['fitval'].iloc[fit_params.index('R1')], fitpars['fitval'].iloc[fit_params.index('R2')], trelax, alignmag, inhomo, J, slps[0], offsets[0], offsets[1], slps[1], offsets[0], offsets[1]))   
 
    f.close()    

    os.system('python cest_sim_master_v10.py BMNS_Simparams.txt')
   
    filename = 'pb_' + '%4.3f'%fitpars['fitval'].iloc[fit_params.index('pB')] + '_kex_' + str(int(fitpars['fitval'].iloc[fit_params.index('kexAB')])) + '_T_' + '%3.2f'%trelax + '_r2b_' + '%2.1f'%fitpars['fitval'].iloc[fit_params.index('R2')] + '.csv'
    print "SEE", filename

    return filename


def plot_highsl_cest(ax, bf, trendlines, noes_trendlines, axes_lims, xstep, startx, ylims, location, fitpars):
    ''' plot highsl cest data '''
    color_slps = ['b', 'cyan', 'g', 'orange', 'brown']
    unique_slps  = np.sort(np.unique(bf['slp(hz)']))
    if len(trendlines) != 0:
        unique_slps2 = np.sort(np.unique(trendlines['slp(hz)']))
    if len(noes_trendlines) != 0:
        unique_slps3 = np.sort(np.unique(noes_trendlines['slp(hz)']))
    
    counter = 0
    for unique_slp in unique_slps:
        if len(noes_trendlines) != 0:
            noes_sim = noes_trendlines.loc[noes_trendlines['slp(hz)'] == unique_slps3[counter]]
            ax.plot(noes_sim['offset(hz)'], noes_sim['norm_intensity'], color=color_slps[counter], linewidth=2, linestyle="--", zorder=0, label='')

        bf_subset = bf.loc[bf['slp(hz)'] == unique_slp]
        ax.errorbar(bf_subset['offset(hz)'], bf_subset['norm_intensity'], yerr=bf_subset['norm_intensity_error'], fmt='o', color=color_slps[counter], zorder=10*(counter+1), elinewidth=2, markersize=5, label=str(int(math.ceil(unique_slp))))
        #ax.plot(bf_subset['Offset'], bf_subset[' Sim R2eff'], color='k', marker='o', linewidth=0.0) 

        if len(trendlines) != 0:
            trendlines_sim = trendlines.loc[trendlines['slp(hz)'] == unique_slps2[counter]]
            ax.plot(trendlines_sim['offset(hz)'], trendlines_sim['norm_intensity'], color=color_slps[counter], linewidth=2, zorder=10*(counter+1), label='')
 
        counter = counter + 1

        ax.tick_params(axis='x', direction='out', width=2, top=False)
        ax.tick_params(axis='y', direction='out', width=2, right=False)
    
        ax.set_xlim(axes_lims)
        ax.set_xticks([int(ele) for ele in np.arange(startx, axes_lims[1], xstep)])
        ax.set_xticklabels([int(float(ele)/1000.) for ele in np.arange(startx, axes_lims[1]+1, xstep)], fontsize=fs)

        ax.set_ylim([ylims[0], ylims[1]])
        ax.set_yticks([ele for ele in np.arange(ylims[2], ylims[1]+0.02, ylims[3])])
        ax.set_yticklabels(['%2.1f'%(ele) for ele in np.arange(ylims[2], ylims[1]+0.02, ylims[3])], fontsize=fs)
    #ax.set_ylabel('Norm. Intensity', fontsize=fs)
    ax.set_title("$r\chi^2$" + " " + "%.2f"%fitpars['fitval'].iloc[9], y=1.01, fontsize=fs)

    counter = 0
    handles, labels = ax.get_legend_handles_labels()
    handles = [h[0] if type(h) is not mpl.lines.Line2D else h for h in handles]
    leg = ax.legend(handles, labels, title='$\omega_{1}/2\pi$' + ' ' + '$(Hz)$', fancybox=False, ncol=2, handlelength=0, frameon=False, columnspacing=0.0, markerscale=0.0, handletextpad=0.5, borderpad=0.0, handleheight=0.0, labelspacing=0.2, loc=location, fontsize=fs)
    for text in leg.get_texts():
        text.set_color(color_slps[counter])
        counter = counter + 1
    for l in leg.get_lines():
        l.set_linestyle('None')
    leg.set_title('$\omega_{1}/2\pi$' + ' ' + '$(Hz)$', prop={'size':fs})


def plot_degeneracy(ax, dirname, param_name_r1p, param_name_hpcest, param_name_lpcest, title, ylims, pb):
    ''' Plot the degeneracy plot for HP/LP CEST and R1p '''
    global pb_master
    color_slps = []

    #######################################################################
    if param_name_r1p != "":
        # First, go for r1p
        pb_values = []
        r1p_rchi2_values = []

        # Loop over directories to get pB, red chi2
        for dummy_filename in os.listdir(dirname):
            # get pb directory
            if dummy_filename[:2] == 'pb':
                #print dummy_filename
                filename = dirname + dummy_filename + '/2state/LocalFits_' + str(param_name_r1p) + '.csv'
                bf = pd.read_csv(filename)
                pb_values.append(float(bf['pB'].iloc[0]))
                r1p_rchi2_values.append(float(bf['RedChiSq'].iloc[0]))

        # Convert values to arrays
        pb_values = np.array(pb_values)
        r1p_rchi2_values = np.array(r1p_rchi2_values)
        
        # sort them correctly
        r1p_rchi2_values_sorted = np.array([x for _,x in sorted(zip(pb_values, r1p_rchi2_values))])
        pb_values_sort = np.sort(pb_values)
        r1p_rchi2_values = r1p_rchi2_values_sorted
        pb_values = pb_values_sort
 
        print pb_values
        #print r1p_rchi2_values
        #print 

        # Normalize rchi2
        r1p_rchi2_values = r1p_rchi2_values/np.amin(r1p_rchi2_values)

        ax.errorbar(pb_values, r1p_rchi2_values, label="$R_{1\\rho}$", color='k', linewidth=4)
        color_slps.append('k')

    #######################################################################

    #######################################################################
    if param_name_hpcest != "":
        # Then go for high power CEST
        pb_values = []
        hpcest_rchi2_values = []

        # Loop over directories to get pB, red chi2
        for dummy_filename in os.listdir(dirname):
            # get pb directory
            if dummy_filename[:2] == 'pb':
                print dummy_filename
                filename = dirname + dummy_filename + '/2state/fitparams-' + str(param_name_hpcest) + '.csv'
                bf = pd.read_csv(filename)
                pb_values.append(float(bf['fitval'].iloc[0]))
                hpcest_rchi2_values.append(float(bf['fitval'].iloc[9]))

        # Convert values to arrays
        pb_values = np.array(pb_values)
        hpcest_rchi2_values = np.array(hpcest_rchi2_values)
        #print pb_values
        #print hpcest_rchi2_values
        #print 

        # sort them correctly
        hpcest_rchi2_values_sorted = np.array([x for _,x in sorted(zip(pb_values, hpcest_rchi2_values))])
        pb_values_sort = np.sort(pb_values)
        hpcest_rchi2_values = hpcest_rchi2_values_sorted
        pb_values = pb_values_sort
         
        # Normalize rchi2 values
        hpcest_rchi2_values = hpcest_rchi2_values/np.amin(hpcest_rchi2_values)
        ax.errorbar(pb_values, hpcest_rchi2_values, label="CEST", color='k', linewidth=4)
        color_slps.append('b')

    #######################################################################
    # Plot title
    ax.set_title(title, fontsize=fs, y=1.02)

    # plot pb
    ax.plot([pb/100, pb/100], ylims, linewidth=4, linestyle='--', color='k')
    print [pb/100, pb/100], ylims
 
    # Log scale
    ax.set_xscale('log')
    if (ylims[0] >= 0 and ylims[1] >= 0): 
        ax.set_yscale('log')

    # Legend
    counter = 0
    handles, labels = ax.get_legend_handles_labels()

    handles = [h[0] if type(h) is not mpl.lines.Line2D else h for h in handles]
    leg = ax.legend(handles, labels, title="", fancybox=False, ncol=1, handlelength=0, frameon=False, columnspacing=0.0, markerscale=0.0, handletextpad=0.2, borderpad=0.0, handleheight=0.0, labelspacing=0.2, fontsize=fs, loc=2)

    for text in leg.get_texts():
        text.set_color('k')
        counter = counter + 1
    for l in leg.get_lines():
        l.set_linestyle('None')
    #ax.legend(fontsize=fs)

    # Ticks
    ax.tick_params(axis='x', which='minor', direction='out', top=False, width=4, length=8)
    ax.tick_params(axis='y', which='minor', direction='out', right=False, width=4, length=8)
    ax.tick_params(axis='x', which='major', direction='out', top=False, width=4, length=8, pad=8)
    ax.tick_params(axis='y', which='major', direction='out', right=False, width=4, length=8, pad=8)
    ax.set_xticks([math.pow(10, ele) for ele in [-4., -3., -2., -1.]])
    ax.set_xticklabels(["$10^{%d}$"%ele for ele in [-4., -3., -2., -1.]], fontsize=fs)
    ax.set_yticks([math.pow(10, ele) for ele in [0., 1., 2.]])
    ax.set_yticklabels(["$10^{%d}$"%ele for ele in [0., 1., 2.]], fontsize=fs)
    # Axes limits
    print ylims
    ax.set_xlim(np.array([pb_values[0], pb_values[-1]]))
    ax.set_ylim(ylims)
      
#fig, ax = plt.subplots(4, 4, figsize=(20, 20))
fig, ax = plt.subplots(1, 7, figsize=(51.1, 8.4))

# degeneracy plots
plot_degeneracy(ax[0], '../figures5/700/A2DNA_G10_4.4_HS_H2O_25C/G10_C8/Fit/final/', 'G10_C8_Off', '', '', "A" + "$_{2}$" + "-DNA G10-C8\npH 4.4HS 25 " + "$^{\circ}$" + "C", [0., 350.], 2.586)
plot_degeneracy(ax[1], '../figures5/henry_data/pH5.3HS_25C/new/', 'G7N1_pH5.3HS_25C', '', '', "AcDNA G7-C8,C1',N1\n" + " pH 5.3HS 25 " + "$^{\circ}$" + "C", [0., 4.], 0.137)
plot_degeneracy(ax[2], '../figures5/700/A2DNA_G10_5.4_H2O_25C/global/new/', 'g10c8', '', '', "A" + "$_{2}$" + "-DNA G10-C8,C1'\npH 5.4 25 " + "$^{\circ}$" + "C", [0., 350.], 1.014)
plot_degeneracy(ax[3], '../figures5/600/A6DNA_G10_5.4_H2O_25C/global/final/', 'gc8', '', '', "A" + "$_{6}$" + "-DNA G10-C8,C1'\npH 5.4 25 " + "$^{\circ}$" + "C", [0., 350.], 0.717)
plot_degeneracy(ax[4], '../figures5/isaac_data/A6DNA_rG10/c15_c6_pH5.4_25C_700/', 'c15_c6', '', '', "A" + "$_{6}$" + "-DNA" + "$^{rG10}$" + " C15-C6\npH 5.4 25 " + "$^{\circ}$" + "C", [0., 350.], 0.577)
plot_degeneracy(ax[5], '../figures5/600/A6DNA_G10_5.4_HS_H2O_25C/global/final/', 'gc8', '', '', "A" + "$_{6}$" + "-DNA G10-C8,C1'\npH 5.4HS 25 " + "$^{\circ}$" + "C", [0., 350.], 0.251)
plot_degeneracy(ax[6], '../figures5/700/A6DNA_G10_6.8_H2O_25C/G_C1p/Fit/new/', 'G10_C1p_Off', '', '', "A" + "$_{6}$" + "-DNA G10-C1'\npH 6.8 25 " + "$^{\circ}$" + "C", [0., 3.], 0.063)





ax[0].set_ylabel("Normalized r" + "$\chi^2$", fontsize=fs)

for dummy in range(7):
    if dummy <= 1:
        ax[dummy].set_xlabel("$p_{i}$", fontsize=fs)
    else:
        ax[dummy].set_xlabel("$p_{i}$", fontsize=fs)

plt.tight_layout()
plt.savefig('figs6_b.pdf')
