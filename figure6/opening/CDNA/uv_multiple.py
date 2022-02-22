
# Curve fitting for melting curves
import matplotlib.pyplot as plt
import numpy as np
import sys
import math
from scipy.optimize import curve_fit
import os

def uv_duplex(t, mds, bds, mss, bss, delH, Tm):
    # Gas constant in kcal/mol
    R = 8.314 / 4184
    # Get the exponential factor
    exp_factor = np.exp(((1/Tm) - (1/(t+273.16))) * delH / R)
    # Get the fraction of double stranded species for a duplex
    f = ((1 + (4*exp_factor)) - np.sqrt(1 + (8*exp_factor)))/(4*exp_factor)
    # Use f to get absorbance of duplex/ss mixture
    absorbance = (((mds*t)+bds)*f) + (((mss*t)+bss)*(1-f))
    return absorbance

def uv_hairpin(t, mds, bds, mss, bss, delH, Tm):
    # Gas constant in kcal/mol
    R = 8.314 / 4184
    # Get the exponential factor
    exp_factor = np.exp(((1/Tm) - (1/(t+273.16))) * delH / R)
    # Get the fraction of hairpin that is intact
    f = exp_factor / (1 + exp_factor)
    # Use f to get absorbance of hairpin/unfolded species
    absorbance = (((mds*t)+bds)*f) + (((mss*t)+bss)*(1-f))
    return absorbance

def fit_data(temperature, absorbance, mode, Ct):
    # Call the appropriate fitting function
    R = 8.314 / 4184
    Ct = float(Ct)
    if mode == 'hairpin': 
        print "Hairpin"
        p0 = np.array([0.0001, 0.3, 0.0001, 0.4, -80, 340])
        [popt, pcov] = curve_fit(uv_hairpin, temperature, absorbance, p0)
        print_result(popt, pcov, ['mds', 'bds', 'mss', 'bss', 'delH', 'Tm'])
        delS = (popt[4] / popt[5]) 
        delG_15 = popt[4] - ((273.16+15)*delS)
        delG_25 = popt[4] - ((273.16+25)*delS)
        delG_37 = popt[4] - ((273.16+37)*delS)
        print "delS ", delS
        print "delG_15C ", delG_15
        print "delG_25C ", delG_25
        print "delG_37C ", delG_37
        plt.figure(1)
        plt.plot(temperature, absorbance, color='k')
        plt.plot(temperature, uv_hairpin(temperature, *popt), color='b')
        plt.show()

        labels = ['mds', 'bds', 'mss', 'bss', 'delH', 'Tm(K)', 'Tm(C)','delS', 'delG_25C', 'delG_37C', 'Ct']
        result = np.array(list(popt) + [popt[5]-273.16,delS, delG_25, delG_37, Ct])

    elif mode == 'duplex':
        print "Duplex"
        p0 = np.array([0.0001, 0.3, 0.0001, 0.4, -80, 300])
        [popt, pcov] = curve_fit(uv_duplex, temperature, absorbance, p0)
        print_result(popt, pcov, ['mds', 'bds', 'mss', 'bss', 'delH', 'Tm'])
        delS = (popt[4] / popt[5]) - (R * math.log(Ct * 0.000001 / 2))
        delG_10 = popt[4] - ((273.16+10)*delS)
        delG_15 = popt[4] - ((273.16+15)*delS)
        delG_20 = popt[4] - ((273.16+20)*delS)
        delG_25 = popt[4] - ((273.16+25)*delS)
        delG_30 = popt[4] - ((273.16+30)*delS)
        delG_37 = popt[4] - ((273.16+37)*delS)
        print "delS ", delS
        print "delG_10C ", delG_10
        print "delG_15C ", delG_15
        print "delG_20C ", delG_20
        print "delG_25C ", delG_25
        print "delG_30C ", delG_30
        print "delG_37C ", delG_37

        plt.figure(1)
        plt.plot(temperature, absorbance, color='k')
        plt.plot(temperature, uv_duplex(temperature, *popt), color='b')
        plt.show()
   
        labels = ['mds', 'bds', 'mss', 'bss', 'delH', 'Tm(K)', 'Tm(C)','delS', 'delG_25C', 'delG_37C', 'Ct']
        result = np.array(list(popt) + [popt[5]-273.16,delS, delG_25, delG_37, Ct])
 
    return [labels, result]

def print_result(parameters, errors, labels):
    ''' Print the results of the fitting '''
    for dummy in range(len(labels)):
        print labels[dummy], parameters[dummy]

def remove_spaces(line):
    line_temp = []
    for ele in line:
        if ele != "":
            line_temp.append(ele)
    return line_temp

def read_data(filename, mode, Ct):
    # Read the UV data file
    f = open(filename, "r")
    temperature = []
    absorbance = []
    # Read the 1st two lines
    line = f.readline()
    line = f.readline()
    # Now read the rest of the lines
    for line in f.readlines():
        if line[0] == 'D' or len(line) <= 2:
            break
        line = line.strip("\n")
        line = line.split("\t")
        line = remove_spaces(line)
        line[0] = float(line[0])
        line[1] = float(line[1])
        temperature.append(line[0])
        absorbance.append(line[1])
    # Now call the functions
    temperature = np.array(temperature)
    absorbance = np.array(absorbance)
    [labels, result] = fit_data(temperature, absorbance, mode, Ct)
 
    # Close the file
    f.close()
    
    return [labels, result]

def caller():
    ''' Call the fitting function and op stats '''
    result_net = []
    labels_net = []
    for filename in os.listdir('.'):
        if filename[-3:] == "txt":
            print(filename)
            [labels, result] = read_data(filename, sys.argv[1], sys.argv[2])
            labels_net = labels
            if result_net == []:
	        result_net = np.array([list(result)])
	    else:
		result_net = np.concatenate((result_net, [list(result)]), axis=0)

    # Now, compute stats
    result_avg = np.average(result_net, axis=0)
    result_std = np.std(result_net, axis=0)

    print result_avg
    print result_std
    print result_net

    # Now, o/p to file
    f = open(str(sys.argv[3]) + ".csv", "w")

    # Write header
    header = labels[0] + ''.join(["," + labels[dummy] for dummy in range(1, len(labels))]) + '\n'
    print header
    f.write(header)

    # Now, write the data
    for dummy in range(np.shape(result_net)[0]):
        data_line = str(result_net[dummy][0]) + ''.join(["," + str(result_net[dummy][dummy_int]) for dummy_int in range(1, np.shape(result_net)[1])]) + '\n'
	f.write(data_line)
  
    # Now, write the stats
    avg_line = str(result_avg[0]) + ''.join(["," + str(result_avg[dummy_int]) for dummy_int in range(1, np.shape(result_net)[1])]) + '\n'
    f.write(avg_line)
    
    std_line = str(result_std[0]) + ''.join(["," + str(result_std[dummy_int]) for dummy_int in range(1, np.shape(result_net)[1])]) + '\n'
    f.write(std_line)

    f.close()
caller()
