#!/usr/bin/python
import sys
from nmrglue import proc_base
from numpy import *
import numpy as np
import math
import matplotlib.pyplot as plt
# Set excited state populations
# and kinetic rates
pesB = 0.08
pesC = 0.0
kexAB = 33.0
kexAC = 0.0
kexBC = 0.0

# Chemical shifts of ES and GS in ppm
wgsA = 0.0
wesB = 1.4
wesC = 0.0  # in ppm

# R2 of ES and GS in s-1
R2_gsA = 15.0  # in 1/s
R2_esB = 15.0  # in 1/s
R2_esC = 15.0

# Larmor frequency of spin
frq0 = 176.4  # 13C for 600MHz spectrometer

# Set the acquisition time and number of points
aq = 6.63
N = 1024 * 16
dt = aq / N
resolution = (1/(aq))/frq0
print "Resolution = ", (1/(aq))/frq0

# Calculate quantities
pgsA = 1-pesB-pesC
k12 = kexAB * pesB / (pesB + pgsA)
k21 = kexAB * pgsA / (pgsA + pesB)
k13 = kexAC * pesC / (pesC + pgsA)
k31 = kexAC * pgsA / (pgsA + pesC)
if pesB + pesC <= 0:
    k23 = k32 = 0.0
else:
    k23 = kexBC * pesC / (pesB + pesC)
    k32 = kexBC * pesB / (pesB + pesC)

plt.figure(1)

for dummy_new in [0]:
    # Set initial magnetization
    M0 = array([pgsA, pesB, pesC])

    # Convert dW into frequency units
    wgsA = wgsA * frq0 * 2 * math.pi
    wesB = wesB * frq0 * 2 * math.pi
    wesC = wesC * frq0 * 2 * math.pi


    # Title of plot
    # As the gyromagnetic ratio for Carbon is positive
    # w = -gamma * B
    # The - in the omega and in the matrix equation cancel out to 
    # give a +ve number
    R = zeros((3, 3), dtype=complex)
    R[0, 0] = -R2_gsA + 1j*wgsA - k12 - k13
    R[0, 1] = k21
    R[0, 2] = k31
    R[1, 0] = k12
    R[1, 1] = -R2_esB + 1j*wesB - k21 - k23
    R[1, 2] = k32
    R[2, 0] = k13
    R[2, 1] = k23
    R[2, 2] = -R2_esC + 1j*wesC - k32 - k31
   
    # v is eigenvalue
    # G is eigenvector
    # G_1 is inverse of eigenvector matrix
    v,G = linalg.eig(R)
    G_1 = linalg.inv(G)
   
    # T is time domain data
    # dM/dT = A . M
    # Solution is M = M0 e^(At)
    # Then expand A according to eigen value equation
    T = linspace(0., aq, N)
    fid = zeros(T.shape, dtype=complex)
    for i,t in enumerate(T):
        A = dot(dot(G,diag(exp(v*t))), G_1)
        fid[i] = sum(dot(A, M0))

    data = proc_base.fft(fid).real

    # 1/dwell_time = 2 * sweep_width
    # we can go sweep width on either side
    # This interval is divided into the same # of points as time
    # domain data
    # Division by fre0 is to convert Hz to ppm
    xf = arange(-1./(2*dt), 1./(2*dt), 1./dt/N)/frq0
    yf = data

    #plt.plot(xf, yf, 'r-', lw=2, marker='o', markersize=5)
    plt.plot(xf, yf+5, 'r-', lw=2, color='k')

    # Find the peak width
    max_height = np.amax(yf)
    yf_side1 = np.absolute(yf-(max_height*0.5))
    xf_side1 = xf[np.argmin(yf_side1)]
    print "Peak center = ", xf[np.argmax(yf)]
    peak_width = (xf_side1 - xf[np.argmax(yf)])*frq0*2
    print "Peak width = ", xf_side1, peak_width
    #plt.xlim([-2, 2])
    #plt.ylim([0., 80.])

plt.legend()
plt.show()
