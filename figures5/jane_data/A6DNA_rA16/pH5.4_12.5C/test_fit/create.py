import os
pb_dir = [0.01, 0.05, 0.1, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.75, 1.0, 2.0, 3.0, 5.0, 8.0, 10.0]
pb_dir = [0.35]
pb_decimals = [float(ele) / 100.0 for ele in pb_dir]

opfile = '''
##################################################################################
# Run the BMNS fitting program:
# > python BMNS.py -fit [BM Parameter Input File] [R1rho Data Directory] (Optional Output directory)
##################################################################################
# Define fitting setup.
# FitType: can be 'global' or 'local' or 'brute'
#          This is for global or local optimizations, not shared parameter fits.
#          'Brute' designates brute-force fixed calculations of the range of parameter
#                   space designated by lower/upper bounds on parameters.
#          - 'brutep' will generate plots at each increment point.
#             WARNING: This can take a LONG time.
#          'Bruteint' brute-forces parameter space by fitting intensities instead of
#                     R1p values
#
#          'Local' uses Levenberg-Marquardt semi-gradient descent/Gauss-Newton algorithm
#          - 'localint' fits intensities directly rather than R1p
#          'Global' uses the "Adaptive Memory Programming for Global Optimizations"
#                   algorithm, with the local 'L-BFGS-B' function, and polishes the
#                   global optimum with L-M.
# FitEqn: fit equation, "BM" for Bloch-McConnell or "Lag" for Laguerre 2-/3-state
# NumFits: is number of fit minima to find (ie. loop over fitting algorithm)
# RandomFitStart : can be 'Yes' or 'No'
#                  if 'Yes', randomly selects initial guess from parameter bounds
##################################################################################
+
FitType local
FitEqn BM
NumFits 1
RandomFitStart No

##################################################################################
# Define fit parameter data, data names, base freqs,
#  initial parameter guesses, and paramter lower and upper bounds. 
#
# Add '+' to read in an additional set of parameters with given 'Name XYZ'
#   The 'Name' must match a .csv data file in given directory of the same name.
#
# Rows for parameters are as follows:
#  [Par name] [initial value] [lower bounds] [upper bounds] ([optional brute force number])
#
# If both lower and upper bounds are not given, they will be set to large values.
# '!' designates a fixed parameter that will not change throughout the fit.
# '*' designates a shared parameter that will be fitted for all data sets
#     also containing the 'x' flag, in a shared manner.
# '@' designates linear brute-force over parameter range of low to upper bounds
# '$' designates log brute-force over parameter range of low to upper bounds
#
# If R1b/c or R2b/c are fixed to 0, they will be shared with R1 / R2
#  e.g. "R1b! = 0.0" will be interpreted as "R1b = R1"
# 
# lf = Larmor frequency (MHz) of the nucleus of interest
#      15N:   60.76302 (600) or  70.960783 (700)
#      13C: 150.784627 (600) or 176.090575 (700)
#
# (optional) rnddel = Fraction of data to be randomly deleted before fit
#                     e.g 'rnddel 0.1' would randomly delete 10pct of data
#
# Temp [Celsius or Kelvin] : Define temperature to calculate free energies
#
# AlignMag [Auto/Avg/GS]
#          Auto : calculates kex/dw and aligns mag depending on slow (gs) vs. fast (avg)
#          Avg : Aligns magnetization/projects along average effective field of GS/ESs
#          GS : Aligns magnetization along ground-state
#
# x-axis Lower Upper (Hz): Sets lower and upper x-axis limits for both plots
#   if not given, will automatically set them
#
# y-axis Lower Upper : Sets lower and upper y-axis limits for both plots
#   if not given, will automatically set them
#
# Trelax increment Tmax (seconds) : sets the increment delay and maximum relaxation
#  delay to simulate R1rho at.
#  Use caution with this flag, recommended that is remains commented out.
#  Array of delays is given as a linear spacing from 0 - Tmax in Tmax/Tinc number of points
#  If not defined, the program will calculate the best Tmax from the experimental
#   R1rho data.
##################################################################################

+
Name ra16_c8_12.5c.csv
lf 150.784186
Temp 25.0
AlignMag AUTO
#Trelax 0.0005 0.5
#x-axis -2000 2000
#y-axis 0 50
pB! %f 1e-4 0.5
pC! 0.0 1e-6 0.5
dwB 3.0 -8. 8.
dwC! 0.0 -80 80
kexAB 4000.0 10.0 10000.0
kexAC! 0.0 1.0 500000.0
kexBC! 0.0 1.0 500000.0
R1 1.5 0.1 8.
R2 16.5 0.1 120.
R1b! 0.0
R2b! 0.0 
R1c! 0.0
R2c! 0.0
'''

for dummy in range(len(pb_dir)):
    ele = pb_dir[dummy]
    os.mkdir('pb_' + str(ele) + 'pct')
    os.chdir('pb_' + str(ele) + 'pct')
    os.system('ln -s ../ra16_c8_12.5c.csv .')
    f = open('BMNS_params.txt', 'w')
    f.write(opfile%pb_decimals[dummy])
    f.close()

    os.system('nohup python /Users/ar359/rot2/rd/progs/BMNS-master/BMNS.py -fit BMNS_params.txt . 2state > log.txt &')
    os.chdir('..')
