# Gives points errors
# Also handles SLP inhomogeniety
# Corrected for proper alignment 
# Handles AUTO alignment
# Couplings in simulations also implemented

#!/usr/bin/python
from uncertainties import ufloat
import uncertainties.unumpy as unumpy
import pandas as pd
import numpy as np
import sys
from nmrglue import proc_base
from nmrglue import analysis
import math
import matplotlib.pyplot as plt
import sys

# Define BM Matrix for CEST
def MatrixBM_3state(k12, k21, k13, k31, k23, k32, delta1, delta2, delta3, w1, R1a, R2a, R1b, R2b, R1c, R2c, pa, pb, pc):
    ''' B-M matrix for 3-state CEST '''
    global params

    # Find range of SLPs to simulate inhomogeniety
    if params['inhomo'] != 0.0:
        sigma_slp = w1 * params['inhomo']
        slps_net = np.linspace(w1 - (2*sigma_slp), w1 + (2*sigma_slp), params['number_inhomo_slps'])
    else:
        slps_net = [w1]

    # Initialize net BM matrix variable
    BM_Mat = np.zeros((10, 10, params['number_inhomo_slps']))

    counter = 0
    for dummy_slp in slps_net: 
        # Reaction rates matrix
        K = np.array([[0.0, 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     ],
                      [0.0, -k12-k13, 0.0     , 0.0     , k21     , 0.0     , 0.0     , k31     , 0.0     , 0.0     ],
                      [0.0, 0.0     , -k12-k13, 0.0     , 0.0     , k21     , 0.0     , 0.0     , k31     , 0.0     ],
                      [0.0, 0.0     , 0.0     , -k12-k13, 0.0     , 0.0     , k21     , 0.0     , 0.0     , k31     ],
                      [0.0, k12     , 0.0     , 0.0     , -k21-k23, 0.0     , 0.0     , k32     , 0.0     , 0.0     ],
                      [0.0, 0.0     , k12     , 0.0     , 0.0     , -k21-k23, 0.0     , 0.0     , k32     , 0.0     ],
                      [0.0, 0.0     , 0.0     , k12     , 0.0     , 0.0     , -k21-k23, 0.0     , 0.0     , k32     ],
                      [0.0, k13     , 0.0     , 0.0     , k23     , 0.0     , 0.0     , -k31-k32, 0.0     , 0.0     ],
                      [0.0, 0.0     , k13     , 0.0     , 0.0     , k23     , 0.0     , 0.0     , -k31-k32, 0.0     ],
                      [0.0, 0.0     , 0.0     , k13     , 0.0     , 0.0     , k23     , 0.0     , 0.0     , -k31-k32]], dtype=float)
  
        # Relaxation rates matrix 
        R = np.array([[0.0   , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     ],
                      [0.0   , -R2a    , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     ],
                      [0.0   , 0.0     , -R2a    , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     ],
                      [R1a*pa, 0.0     , 0.0     , -R1a    , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     ],
                      [0.0   , 0.0     , 0.0     , 0.0     , -R2b    , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     ],
                      [0.0   , 0.0     , 0.0     , 0.0     , 0.0     , -R2b    , 0.0     , 0.0     , 0.0     , 0.0     ],
                      [R1b*pb, 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , -R1b    , 0.0     , 0.0     , 0.0     ],
                      [0.0   , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , -R2c    , 0.0     , 0.0     ],
                      [0.0   , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , -R2c    , 0.0     ],
                      [R1c*pc, 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , -R1c    ]], dtype=float)

        # Spin-lock matrix
        SL = np.array([[0.0, 0.0       , 0.0     , 0.0      , 0.0       , 0.0     , 0.0      , 0.0       , 0.0     , 0.0      ],
                       [0.0, 0.0       , 0.0     , dummy_slp, 0.0       , 0.0     , 0.0      , 0.0       , 0.0     , 0.0      ],
                       [0.0, 0.0       , 0.0     , 0.0      , 0.0       , 0.0     , 0.0      , 0.0       , 0.0     , 0.0      ],
                       [0.0, -dummy_slp, 0.0     , 0.0      , 0.0       , 0.0     , 0.0      , 0.0       , 0.0     , 0.0      ],
                       [0.0, 0.0       , 0.0     , 0.0      , 0.0       , 0.0     , dummy_slp, 0.0       , 0.0     , 0.0      ],
                       [0.0, 0.0       , 0.0     , 0.0      , 0.0       , 0.0     , 0.0      , 0.0       , 0.0     , 0.0      ],
                       [0.0, 0.0       , 0.0     , 0.0      , -dummy_slp, 0.0     , 0.0      , 0.0       , 0.0     , 0.0      ],
                       [0.0, 0.0       , 0.0     , 0.0      , 0.0       , 0.0     , 0.0      , 0.0       , 0.0     , dummy_slp],
                       [0.0, 0.0       , 0.0     , 0.0      , 0.0       , 0.0     , 0.0      , 0.0       , 0.0     , 0.0      ],
                       [0.0, 0.0       , 0.0     , 0.0      , 0.0       , 0.0     , 0.0      , -dummy_slp, 0.0     , 0.0      ]], dtype=float)

        # Offsets matrix
        OMEGA = np.array([[0.0, 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     ],
                          [0.0, 0.0     , -delta1 , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     ],
                          [0.0, delta1  , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     ],
                          [0.0, 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     ],
                          [0.0, 0.0     , 0.0     , 0.0     , 0.0     , -delta2 , 0.0     , 0.0     , 0.0     , 0.0     ],
                          [0.0, 0.0     , 0.0     , 0.0     , delta2  , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     ],
                          [0.0, 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     ],
                          [0.0, 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , -delta3 , 0.0     ],
                          [0.0, 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , delta3  , 0.0     , 0.0     ],
                          [0.0, 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     , 0.0     ]], dtype=float)

        # Store matrix slice in master variable
        BM_Mat[:, :, counter] = (K + R + SL + OMEGA) 

        counter = counter + 1

    # Return BM matrix
    return BM_Mat


def ls_3state(pa, pb, pc, k12, k21, k13, k31, k23, k32, wa, wb, wc, R2a, R2b, R2c, frq0, resn):
    ''' Do a 3-state lineshape simulation to get observed peak position '''
    # Pass the function the populations and individual rate constants 
    # populations = fractions or % of 1
    # rates = /s
    # w = rad/s
    # R2 = /s
    # frq0 = MHz
    # Resolution (resn) = Hz
 
    # Based on resolution, get acquisition time
    aq = 1/resn
  
    # Set sweep width, as 1.5x the greater of the two dWs. 
    # This sw is in Hz
    sw = np.amax(1.5*np.abs(np.array([wb, wc]))/(2*math.pi))
    
    # Now, determine number of points
    N = int(aq * 2 * sw)
    
    # Get dwell time 
    dt = aq / N
    
    # Set initial magnetization
    M0 = np.array([pa, pb, pc])
    
    # Title of plot
    # As the gyromagnetic ratio for Carbon is positive
    # w = -gamma * B
    # The - in the omega and in the matrix equation cancel out to 
    # give a +ve number
    R = np.zeros((3, 3), dtype=complex)
    R[0, 0] = -R2a + 1j*wa - k12 - k13
    R[0, 1] = k21
    R[0, 2] = k31
    R[1, 0] = k12
    R[1, 1] = -R2b + 1j*wb - k21 - k23
    R[1, 2] = k32
    R[2, 0] = k13
    R[2, 1] = k23
    R[2, 2] = -R2c + 1j*wc - k32 - k31
    
    # v is eigenvalue
    # G is eigenvector
    # G_1 is inverse of eigenvector matrix
    v,G = np.linalg.eig(R)
    G_1 = np.linalg.inv(G)
    
    # T is time domain data
    # dM/dT = A . M
    # Solution is M = M0 e^(At)
    # Then expand A according to eigen value equation
    T = np.linspace(0., aq, N)
    fid = np.zeros(T.shape, dtype=complex)
    for i,t in enumerate(T):
        A = np.dot(np.dot(G,np.diag(np.exp(v*t))), G_1)
        fid[i] = np.sum(np.dot(A, M0))
    
    data = proc_base.fft(fid).real
    
    # 1/dwell_time = 2 * sweep_width
    # we can go sweep width on either side
    # This interval is divided into the same # of points as time
    # domain data
    # Units of xf = Hz
    xf = np.arange(-1./(2*dt), 1./(2*dt), 1./dt/N)
    yf = data
    #print xf.shape, yf.shape 
    #plt.figure(1)
    #plt.plot(xf, yf, 'r-', lw=2, marker='o', markersize=5)
   

    # Find the peak width
    max_height = np.amax(yf)
    #print "SEE", max_height

    # Pick peaks
    #data = data[1:-1]
    #peaks = analysis.peakpick.pick(data, pthres=max_height/50.0, est_params=False, cluster=False, algorithm='downward')
    #print "peak fitting = ", peaks

    # Find peak maximum and convert to rad/s
    max_position = xf[np.argmax(yf)] * 2 * math.pi
    #print "COMP", np.argmax(yf)
    #print "peak position =", max_position, " Hz"
    #print "peak height =", max_height, " Hz"
    #plt.plot([max_position/(2*math.pi), max_position/(2*math.pi)], [0, max_height], color='k')

    #yf_side1 = np.absolute(yf-(max_height*0.5))
    #xf_side1 = xf[np.argmin(yf_side1)]
    #peak_width = (xf_side1 - xf[np.argmax(yf)])*frq0*2
    #print "Peak width = ", xf_side1, peak_width
    
    #plt.show()
 
    # return peak position in rad/s
    return [max_position, max_height]


def input_parser(input_filename):
    ''' parses input parameters file '''
    params_dict = {}

    f = open(input_filename, 'r')

    # Parse until you get to a +
    while 1:
        line = f.readline()
        if line[0] == '+':
            break

    # Now, keep reading until you get to a blank line
    # Parse the lines to get all relevant exchange parameters
    while 1:
        line = f.readline()
        line = line.strip('\n')
        if len(line) <= 1:
            break
        line = line.split(' ')
        if str(line[0]) != 'mode' and str(line[0]) != 'equil' and str(line[0]) != 'ls':
            params_dict[str(line[0])] = float(line[1])
        else:
            params_dict[str(line[0])] = str(line[1])

    # We need to parse AUTO alignment and decide whether it's equal 
    # to AVG or GS alignment here
    # So 1st, check whether alignment mode is set to AUTO 
    if params_dict['mode'] == "AUTO":
        # Compute exchange regime & decide alignment based on the same
        if (params_dict['dwb'] != 0.0):
            exch_regime = params_dict['kexAB']/(params_dict['dwb']*2*math.pi*params_dict['lf'])
        else:
            exch_regime = 0.1
        
        # Now, the exchange regime has been computed
        # decide alignment
        if exch_regime <= 1:
            params_dict['mode'] = 'GS'
            print "Auto decision, mode=GS"
        else:
            params_dict['mode'] = 'AVG'
            print "Auto decision, mode=AVG"

    # Sanity checking
    if params_dict['mode'] not in ['GS', 'AVG']:
        print "!!! Enter GS for GS alignment and AVG for avg alignment !!!" 
        sys.exit(0)
    if params_dict['equil'] not in ['Y', 'N']:
        print "!!! Please enter Y to allow equilibration and N for no equilibration !!!"
        sys.exit(0)
    if params_dict['ls'] not in ['Y', 'N']:
        print "!!! Wrong lineshape input. Enter Y for performing a lineshape simulation to get peak position, "
        print "and N to assume observed peak is GS peak !!!"
        sys.exit(0)

    if params_dict['inhomo'] == 0:
        params_dict['number_inhomo_slps'] = 1
    else:
        params_dict['number_inhomo_slps'] = 20

    # Convert into appropriate units
    # and define missing parameters
    params_dict['pa'] = 1 - params_dict['pb'] - params_dict['pc']
    params_dict['wa'] = 0.0
    params_dict['wb'] = (params_dict['dwb'] - params_dict['wa']) * params_dict['lf'] * 2 * math.pi
    params_dict['wc'] = (params_dict['dwc'] - params_dict['wa']) * params_dict['lf'] * 2 * math.pi
    params_dict['k12'] = params_dict['kexAB'] * params_dict['pb'] / (params_dict['pb'] + params_dict['pa'])
    params_dict['k21'] = params_dict['kexAB'] * params_dict['pa'] / (params_dict['pa'] + params_dict['pb'])
    params_dict['k13'] = params_dict['kexAC'] * params_dict['pc'] / (params_dict['pc'] + params_dict['pa'])
    params_dict['k31'] = params_dict['kexAC'] * params_dict['pa'] / (params_dict['pa'] + params_dict['pc'])
    if params_dict['pb'] + params_dict['pc'] > 0:
        params_dict['k23'] = params_dict['kexBC'] * params_dict['pc'] / (params_dict['pb'] + params_dict['pc'])
        params_dict['k32'] = params_dict['kexBC'] * params_dict['pb'] / (params_dict['pb'] + params_dict['pc'])
    else:
        params_dict['k23'] = 0.0
        params_dict['k32'] = 0.0

    # Now, key point - do the lineshape simulation
    # to get observed peak position
    if params_dict['ls'] == "Y":
        print "***LINESHAPE SIMULATION***"    
        [observed_peakpos, observed_peakheight] = ls_3state(params_dict['pa'], params_dict['pb'], params_dict['pc'], params_dict['k12'], params_dict['k21'], params_dict['k13'], params_dict['k31'], params_dict['k23'], params_dict['k32'], params_dict['wa'], params_dict['wb'], params_dict['wc'], params_dict['R2a'], params_dict['R2b'], params_dict['R2c'], params_dict['lf'], params_dict['resn'])
        print "Observed peakpos = ", '%4.3f'%(observed_peakpos/(2*math.pi)), " Hz, ", '%4.3f'%((observed_peakpos/(2*math.pi))/params_dict['lf']), " ppm"
    else:
        if params_dict['mode'] == "GS":
            observed_peakpos = 0.0
            print "GS Alignment"
        elif params_dict['mode'] == "AVG":
            observed_peakpos = (params_dict['pb'] * params_dict['wb']) + (params_dict['pc'] * params_dict['wc'])
            print "AVG Alignment"
 
    params_dict['obs_peakpos'] = observed_peakpos

    #for ele in params_dict.keys():
    #    print ele, params_dict[ele]
 
    # Keep reading all the comment lines
    while 1:
        line = f.readline()
        if line[0] == '+':
            break

    # Read until you get to blank line
    # Parse the spin-lock power offset combinations
    slps = np.array([])
    offset = np.array([])
    while 1:
        line = f.readline()
        if len(line) <= 1:
            break
        line = line.split(' ')
        lower_limit = float(line[2])
        upper_limit = float(line[3])
        no_points = int(line[4])
        spinlock_power = float(line[1])
        offset = np.append(offset, np.linspace(lower_limit, upper_limit, no_points))
        slps = np.append(slps, np.tile(np.array([spinlock_power]), no_points))

    f.close()
   
    # Add the slps and offsets to params_dict
    params_dict['offsets'] = offset
    params_dict['slps'] = slps

    # Key point - convert offset and slp into Hz
    # Also do the same for J coupling
    params_dict['offsets'] = params_dict['offsets'] * 2 * math.pi
    params_dict['slps'] = params_dict['slps'] * 2 * math.pi
    params_dict['J'] = params_dict['J'] * 2 * math.pi

    # Now, define the offsets, i.e., delta values for the 3 states 
    # Offset_peak = Offset + (obs_peakpos - w_peak) - this is wrong. 
    # What we are interested in is not the offset of the other peak, 
    # but the displacement of peak from offset field, when offset field is set at 0.0
    # Everything must be in rotating frame of the B1 field
    params_dict['delta1'] = (-1.*params_dict['offsets']) + (params_dict['wa']-params_dict['obs_peakpos'])
    params_dict['delta2'] = (-1.*params_dict['offsets']) + (params_dict['wb'] - params_dict['obs_peakpos'])
    params_dict['delta3'] = (-1.*params_dict['offsets']) + (params_dict['wc'] - params_dict['obs_peakpos'])
  
    # return parameter dict
    return params_dict


def matrix_exponential(bm_matrix, time_point, status_complexvals):
    ''' Matrix exponent using bm_matrix at given time point '''
    # Sanity check
    if np.isnan(bm_matrix).any() == True or np.isinf(bm_matrix).any() == True:
        print "There are nan elements in BM matrix!!!"
        sys.exit(0)
    
    # Get the eigen decomposition of the B-M matrix
    bm_matrix_evals, bm_matrix_evecs = np.linalg.eig(bm_matrix*time_point)

    if status_complexvals == 'off':
        if np.iscomplexobj(bm_matrix_evals):
            bm_matrix_evals = bm_matrix_evals.real

    # Compute matrix exponent using eigen decomposition
    mat_exp_value = np.dot(np.dot(bm_matrix_evecs, np.diag(np.exp(bm_matrix_evals))), np.linalg.inv(bm_matrix_evecs))
   
    # Return the real part of this matrix exponential matrix
    return mat_exp_value.real


def simulate_cest(params_dict, status_complex):
    # A pair of curves at wobs-piJ, wbox+piJ
    if params_dict['J'] != 0.0:
        print "non zero J", params_dict['J'] * 0.5, -0.5 * params_dict['J']
        intensity1 = simulate_cest_offset(params_dict, status_complex, params_dict['J']*0.5)
        intensity2 = simulate_cest_offset(params_dict, status_complex, params_dict['J']*-0.5)
        return (intensity1 + intensity2) * 0.5
    else:
        print "zero J"
        return simulate_cest_offset(params_dict, status_complex, 0.0)


def simulate_cest_offset(params_dict, status_complex, j_offset):
    ''' Simulate an acutal CEST profile given exchange parameters '''
    # Generate the weights for each SLP here
    # Inhomogeniety is same across all spins
    if params_dict['inhomo'] == 0.0:
        weights_slps = [1.0]
    else:
        w1_dummy = 10.0
        sigma_slp = w1_dummy * params_dict['inhomo']
        slps_net = np.linspace(w1_dummy - (2*sigma_slp), w1_dummy + (2 * sigma_slp), params_dict['number_inhomo_slps'])
        weights_slps = (np.exp((-1*np.square(slps_net-w1_dummy))/(2*sigma_slp*sigma_slp))/math.sqrt(2*math.pi*sigma_slp*sigma_slp))

    BM_mat = MatrixBM_3state(params_dict['k12'], params_dict['k21'], params_dict['k13'], params_dict['k31'], params_dict['k23'], params_dict['k32'], params_dict['delta1'][0]-j_offset, params_dict['delta2'][0]-j_offset, params_dict['delta3'][0]-j_offset, params_dict['slps'][0], params_dict['R1a'], params_dict['R2a'], params_dict['R1b'], params_dict['R2b'], params_dict['R1c'], params_dict['R2c'], params_dict['pa'], params_dict['pb'], params_dict['pc'])

    # Define starting magnetization
    # (unit_vector, GS_x, GS_y, GS_z, ES1_x, ES1_y, ES1_z, ES2_x, ES2_y, ES2_z
    if params_dict['equil'] == 'Y':
        print "*** COMPLETE EQUILIBRATIONi ***"
        M0_plusx = np.array([1.0, 0.0, 0.0, params_dict['pa'], 0.0, 0.0, params_dict['pb'], 0.0, 0.0, params_dict['pc']])
        M0_minusx = np.array([1.0, 0.0, 0.0, -1.0*params_dict['pa'], 0.0, 0.0, -1.0*params_dict['pb'], 0.0, 0.0, -1.0*params_dict['pc']])
    elif params_dict['equil'] == 'N':
        print "*** NO EQUILIBRATIONi ***"
        M0_plusx = np.array([1.0, 0.0, 0.0, params_dict['pa'], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        M0_minusx = np.array([1.0, 0.0, 0.0, -1.0*params_dict['pa'], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    # Now, simulate!
    # Initial time point
    # Take into account phase cycling
    time_point = 0.0
    intensities_initial = np.zeros(10)
    for dummy_inhomo in range(params_dict['number_inhomo_slps']):
        matrix_expo_value = matrix_exponential(BM_mat[:,:,dummy_inhomo], time_point, status_complex)
        intensities_plusx_initial = np.dot(matrix_expo_value, M0_plusx) 
        intensities_minusx_initial = np.dot(matrix_expo_value, M0_minusx) 
        intensities_initial = intensities_initial + ((intensities_plusx_initial-intensities_minusx_initial)*weights_slps[dummy_inhomo])

    # Get errors on intensity matrix
    #print 'base inten', params['error_baseinten']
    intensities_initial_pointerror = intensities_initial * params_dict['error_point']
    intensities_initial_error = intensities_initial * params_dict['error_baseinten']
    intensities_initial = intensities_initial 

    # Define variable for the magnetization of each state at each slp, offset value
    #intensities_net = np.zeros((10,params_dict['slps'].shape[0]))
    intensities_net = np.zeros((params_dict['slps'].shape[0]))
    intensities_net = unumpy.uarray(intensities_net, intensities_net)
 
    for dummy in range(len(params_dict['slps'])):
        #print 'SLP = ', params_dict['slps'][dummy]/(2*math.pi), ' Hz'
        #print 'Offset = ', params_dict['offsets'][dummy]/(2*math.pi), ' Hz'
        #print 'Offset GS = ', params_dict['delta1'][dummy]/(2*math.pi), ' Hz'
        #print 'Offset ES1 = ', params_dict['delta2'][dummy]/(2*math.pi), ' Hz'
        #print 'Offset ES2 = ', params_dict['delta3'][dummy]/(2*math.pi), ' Hz'
        ##print params_dict['wa'], params_dict['wb'], params_dict['wc'], params_dict['obs_peakpos']
        #print
    
        # Get the B-M matrix for this particular SLP and Offset combination
        BM_mat = MatrixBM_3state(params_dict['k12'], params_dict['k21'], params_dict['k13'], params_dict['k31'], params_dict['k23'], params_dict['k32'], params_dict['delta1'][dummy]-j_offset, params_dict['delta2'][dummy]-j_offset, params_dict['delta3'][dummy]-j_offset, params_dict['slps'][dummy], params_dict['R1a'], params_dict['R2a'], params_dict['R1b'], params_dict['R2b'], params_dict['R1c'], params_dict['R2c'], params_dict['pa'], params_dict['pb'], params_dict['pc'])

        # Now, simulate at time point of interest
        time_point = params_dict['T']
        #print params_dict['T']
        intensities = np.zeros(10)
        for dummy_inhomo in range(params_dict['number_inhomo_slps']):
            matrix_expo_value = matrix_exponential(BM_mat[:,:,dummy_inhomo], time_point, status_complex)
            intensities_plusx = np.dot(matrix_expo_value, M0_plusx) 
            intensities_minusx = np.dot(matrix_expo_value, M0_minusx) 
            intensities = intensities + ((intensities_plusx-intensities_minusx)*weights_slps[dummy_inhomo])

        # Now, normalize this intensity
        # Key point is that normalization will be different if exchange regime is different
        # Put alternatively, super slow exchange ==> GS_mag(t)/GS_mag(0)
        # At super fast exchange ==> (GS_mag(t) + ES_mag(t)) / (GS_mag(0) + ES_mag(0))
        # Things could be more complicated wherein GS <-> ES1 is fast while GS <-> ES2 is slow
        # One can accordingly define more complicated normalization schemes
        # Here it would be  ==> (GS_mag(t) + ES1_mag(t)) / (GS_mag(0) + ES1_mag(0))
        #intensities_normalized = np.divide(intensities, intensities_initial)
        #print params_dict['mode'], len(params_dict['mode']) 
        if params_dict['mode'] == 'GS':
            #print "GS Alignment"
            intensities_net[dummy] = ufloat(intensities[3]+np.random.normal(0.0, intensities_initial_pointerror[3],1), intensities_initial_error[3]) / intensities_initial[3]
        else:
            #print "AVG Alignment"
            denominator = (intensities_initial[3] + intensities_initial[6] + intensities_initial[9])
	    intensities_net[dummy] = ufloat((intensities[3]+intensities[6]+intensities[9]+np.random.normal(0.0, intensities_initial_pointerror[3],1)+np.random.normal(0.0, intensities_initial_pointerror[6],1)+np.random.normal(0.0,intensities_initial_pointerror[9],1)), intensities_initial_error[3]) / denominator 
 
    return intensities_net    

def plot_profiles(params, intensities):
    ''' Given the parameters, intensity from simulation, plot the CEST profile '''
    # Correct CEST output for fast exchange = run lineshape with pB at the end of the CEST period
    # Determine intensity at observed peak position
    # Plot intensity
    # For now, plot GS peak intensity only
    number_slps = np.unique(params['slps'])
    colors_plot = ['k', 'r', 'b', 'g', 'cyan', 'magenta', 'brown', 'yellow', 'teal', 'lightgreen']
    # Plot profile
    plt.figure(1)
    counter = 0
    for dummy_slp in number_slps:
        plt.errorbar((params['offsets'][np.where(params['slps'] == dummy_slp)]/(2*math.pi*params['lf'])), unumpy.nominal_values(intensities[np.where(params['slps'] == dummy_slp)]), yerr=unumpy.std_devs(intensities[np.where(params['slps'] == dummy_slp)]), color = colors_plot[counter], label='%4.2f'%float(dummy_slp/(2*math.pi)), linewidth=3, marker='o', markersize=6)
        #plt.plot((params['offsets'][np.where(params['slps'] == dummy_slp)]/(2*math.pi*params['lf']))+138.4, intensities[np.where(params['slps'] == dummy_slp)], color = colors_plot[counter], label=str(int(dummy_slp)), linewidth=3)
        #plt.plot((params['offsets'][np.where(params['slps'] == dummy_slp)]/(2*math.pi*params['lf']))+138.4, intensities[3][np.where(params['slps'] == dummy_slp)[0]], color = colors_plot[counter], label=str(int(dummy_slp)), linewidth=3)
        #plt.plot((params['offsets'][np.where(params['slps'] == dummy_slp)][::10]/(2*math.pi*params['lf']))+138.4, intensities[np.where(params['slps'] == dummy_slp)][::10], color = colors_plot[counter], label=str(int(dummy_slp)), linewidth=0.0, marker='o', markersize=8)

        counter = counter + 1
    #plt.xlim([((np.amin(params['offsets'])-(2*math.pi*0.8*params['lf']))/(2*math.pi*params['lf']))+138.4, ((np.amax(params['offsets'])+(2*math.pi*0.8*params['lf']))/(2*math.pi*params['lf']))+138.4])
    plt.xlim([((np.amin(params['offsets'])-(2*math.pi*0.8*params['lf']))/(2*math.pi*params['lf'])), ((np.amax(params['offsets'])+(2*math.pi*0.8*params['lf']))/(2*math.pi*params['lf']))])
    #plt.ylim([-0.05, 0.6])
    plt.legend(loc=4)
    plt.show()
 
  
def plot_profiles_comp(params, intensities_gs, bf_new):
    ''' Given the parameters, intensity from simulation, plot the CEST profile '''
    # Correct CEST output for fast exchange = run lineshape with pB at the end of the CEST period
    # Determine intensity at observed peak position
    # Plot intensity
    # For now, plot GS peak intensity only
    number_slps = np.unique(params['slps'])
    colors_plot = ['k', 'r', 'b', 'g', 'cyan', 'magenta', 'brown', 'yellow', 'teal', 'lightgreen']
    # Plot profile
    plt.figure(1)
    counter = 0
    #plt.plot(params['offsets']/(2*math.pi), unumpy.nominal_values(intensities_gs), color = 'r', label='Sim', linewidth=0.0, marker='o', markersize=5, markerfacecolor='r')
    #plt.plot(bf_new['offset(hz)'], bf_new['norm_intensity'], color = 'k', label='Exp', linewidth=0.0, marker='o', markersize=2, markerfacecolor='k')

    #plt.xlim([((np.amin(params['offsets'])-(2*math.pi*0.8*params['lf']))/(2*math.pi*params['lf'])), ((np.amax(params['offsets'])+(2*math.pi*0.8*params['lf']))/(2*math.pi*params['lf']))])
    #plt.ylim([-0.05, 0.6])
    plt.legend()
    plt.show()

def plot_profiles_split(params, intensity1, intensity2, label1, label2):
    ''' Plot 2 CEST profiles side by side '''
    number_slps = np.unique(params['slps'])
    colors_plot = ['k', 'r', 'b', 'g', 'cyan', 'magenta', 'brown', 'yellow', 'teal', 'lightgreen']
    # Plot profile
    fig, ax = plt.subplots(1,2)

    # Get y axis limits
    min_yaxis = np.amin(np.array([np.amin(unumpy.nominal_values(intensity1)), np.amin(unumpy.nominal_values(intensity2))]))
    max_yaxis = np.amin(np.array([np.amax(unumpy.nominal_values(intensity1)), np.amax(unumpy.nominal_values(intensity2))]))
    min_yaxis = min_yaxis - 0.05
    max_yaxis = max_yaxis + 0.05

    counter = 0
    for dummy_slp in number_slps:
        ax[0].errorbar(unumpy.nominal_values(intensity1[np.where(params['slps'] == dummy_slp)]), (params['offsets'][np.where(params['slps'] == dummy_slp)]/(2*math.pi*params['lf'])), xerr=unumpy.std_devs(intensity1[np.where(params['slps'] == dummy_slp)]), color = colors_plot[counter], label='%4.2f'%float(dummy_slp/(2*math.pi)), linewidth=3, marker='o', markersize=3)
        counter = counter + 1
    ax[0].set_title(label1)
    #ax[0].set_xlim([((np.amin(params['offsets'])-(2*math.pi*0.8*params['lf']))/(2*math.pi*params['lf'])), ((np.amax(params['offsets'])+(2*math.pi*0.8*params['lf']))/(2*math.pi*params['lf']))])
    ax[0].set_ylim([np.amin(params['offsets'])/(2*math.pi*params['lf']), np.amax(params['offsets']/(2*math.pi*params['lf']))])
    ax[0].set_xlim([-0.05, 0.52])
    ax[0].legend(loc=4)
    ax[0].set_xlabel('$I/I_{0}$')
    ax[0].set_ylabel('$\Omega/2\pi$' + ' (ppm)')

    counter = 0
    for dummy_slp in number_slps:
        ax[1].errorbar(unumpy.nominal_values(intensity2[np.where(params['slps'] == dummy_slp)]), (params['offsets'][np.where(params['slps'] == dummy_slp)]/(2*math.pi*params['lf'])), xerr=unumpy.std_devs(intensity2[np.where(params['slps'] == dummy_slp)]), color = colors_plot[counter], label='%4.2f'%float(dummy_slp/(2*math.pi)), linewidth=3, marker='o', markersize=3)
        counter = counter + 1
    ax[1].set_title(label2)
    #ax[1].set_xlim([((np.amin(params['offsets'])-(2*math.pi*0.8*params['lf']))/(2*math.pi*params['lf'])), ((np.amax(params['offsets'])+(2*math.pi*0.8*params['lf']))/(2*math.pi*params['lf']))])
    #ax[1].set_ylim([min_yaxis, max_yaxis])
    ax[1].set_ylim([np.amin(params['offsets'])/(2*math.pi*params['lf']), np.amax(params['offsets']/(2*math.pi*params['lf']))])
    ax[1].set_xlim([-0.05, 0.52])
    ax[1].legend(loc=4)
    ax[1].set_ylabel('$\Omega/2\pi$' + ' (ppm)')
    ax[1].set_xlabel('$I/I_{0}$')
    
    plt.show()





# parse the input parameters 
input_filename = str(sys.argv[1])
params = input_parser(input_filename)

# Do the actual simulations and get the intensity ratios
intensity_mag_complex = simulate_cest(params, 'on')
intensity_mag_nocomplex = simulate_cest(params, 'off')
plot_profiles(params, intensity_mag_complex)
#plot_profiles_split(params, intensity_mag_complex, intensity_mag_nocomplex, 'Complex', 'No Complex')

#plot_profiles(params, intensity_mag)
#output_filename = 'pb_' + '%5.4f'%params['pb'] + '_kex_' + str(int(params['kexAB'])) + '_T_' + '%2.2f'%params['T'] + '.csv'

# Alignment mode does not matter based on benchmark simulations
#params['R1b'] = 0.0
#intensity_mag_r2 = simulate_cest(params)
#params['mode'] = 'GS'
#intensity_mag_gs = simulate_cest(params)
#params['mode'] = 'AVG'
#intensity_mag_avg = simulate_cest(params)
#bf_new = pd.read_csv('data_1_edit.csv')
##plot_profiles_comp(params, intensity_mag, bf_new)
plot_profiles(params, intensity_mag)


output_dict = {}
output_dict['slp(hz)'] = params['slps']/(2*math.pi)
output_dict['offset(hz)'] = params['offsets']/(2*math.pi)
output_dict['norm_intensity'] = unumpy.nominal_values(intensity_mag_complex)
output_dict['norm_intensity_error'] = unumpy.std_devs(intensity_mag_complex)
output_dict['norm_intensity_nocomplex'] = unumpy.nominal_values(intensity_mag_nocomplex)
output_dict['norm_intensity_nocomplex_error'] = unumpy.std_devs(intensity_mag_nocomplex)
output_dict['trelax(s)'] = np.ones(intensity_mag_complex.shape)*params['T']

bf=pd.DataFrame(data=output_dict)

output_filename = 'pb_' + '%4.3f'%(params['pb'])  + '_kex_' + str(int(params['kexAB'])) + '_T_' + '%3.2f'%(params['T']) + '_r2b_' + '%3.1f'%(params['R2b']) + '.csv'

print output_filename
bf.to_csv(output_filename)
