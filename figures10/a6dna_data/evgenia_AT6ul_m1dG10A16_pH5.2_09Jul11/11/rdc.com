#!/bin/csh

# Basic 2D Processing:
#   Cosine-Bells are used in both dimensions.
#   Use of "ZF -auto" rounds to power of 2 and doubles.
#   Use of "FT -auto" chooses correct Transform mode. Other options (-alt, -neg) may be needed.
#   Imaginaries are deleted during phasing with "-di" in each dimension.
#   Phase corrections should be inserted by hand. ALWAYS START WITH NO PHASING OR BASELINE CORRECTIONS!
#   More details below the script


nmrPipe -in test.fid \
#| nmrPipe  -fn SOL     									\
| nmrPipe  -fn POLY -time      								\
| nmrPipe  -fn EM -lb -6									\
| nmrPipe  -fn SP -off 0.5 -end 1. -pow 2 -c 0.5       					\
| nmrPipe  -fn ZF -auto -zf 3   								\
| nmrPipe  -fn FT -auto 									\
| nmrPipe  -fn PS -p0 94.0 -p1 0.00 -di -verb    						\
| nmrPipe  -fn EXT -left -sw   								\
#| nmrPipe  -fn BASE -nw .5ppm -nl 0% 5ppm 6ppm 9ppm 10ppm 100%  				\
#| nmrPipe  -fn POLY -nw .5ppm -nl 4.9ppm 5.0ppm 6.75ppm 7.0ppm 7.5ppm -sx1 0% -sxn 100% -avg -ord 0   	\
#| nmrPipe  -fn POLY -auto  								\
| nmrPipe  -fn TP        									\
| nmrPipe  -fn EM -lb -6									\
| nmrPipe  -fn SP -off 0.5 -end 1. -pow 2 -c 0.5       					\
| nmrPipe  -fn ZF -auto -zf 3    								\
| nmrPipe  -fn FT -auto  									\
| nmrPipe  -fn PS -p0 89.0 -p1 0.00 -di -verb     						\
#| nmrPipe  -fn EXT -x1 131ppm -xn 155ppm -sw    						\
#| nmrPipe  -fn BASE -nw 1ppm -nl 0% 130ppm 133ppm 147ppm 156ppm 100%    			\
#| nmrPipe  -fn POLY -nw 1ppm -nl 130ppm 133ppm 147ppm 156ppm -sx1 0% -sxn 100% -avg -ord 2     \
| nmrPipe  -fn POLY -auto   								\
#| nmrPipe  -fn EXT -x1 6.5ppm -xn 9ppm -sw       						\
   -ov -out test.ft2



#########################################################################################
#			nmrPipe Command Descriptions					#
#########################################################################################

#   Water Supression Options
# -fn SOL -- use only when peaks are far (>1ppm) from the water peak
# -fn POLY -time -- softer suppression


#   Window Functions
#	General Options:
#		  -size xxx  : # pts in function (ie. truncates shape)
#		  -start xxx : start point of function (for composite window functions)
#		  -c xxx     : Scaling of first point (0 to 1). Use 0.5 when no 1st order phase correction
#		  -one       : See manual (man SP)
#		  -hdr       : See manual (man SP)
#		  -inv       : Applies Inverse window function (for reverse processing)
#	Standard Functions:
# -fn EM -- Exponential window function
#	-lb xxx  : Line Broadening value (negative value sharpens peaks)
# -fn GM -- Lorentz-to-Gauss window function
# 	-g1 xxx  : Line sharpening exponentail term (opposite of EM), default is 0 (for pure gaussian)
#	-g2 xxx  : Gaussian broadening term (1.25-4x g1)
#	-g3 xxx  : Position of Gaussian (0 to 1) on FID. Default is 0
# -fn GMB -- Lorentz-to-Gauss window function (as in Bruker processing)
#	-lb xxx  : Opposite of g1 above (ie neg value sharpens).  Default is 0
#	-gb xxx  : Similar to g2 above, but as a fraction of 1
# -fn SP -- Sine bell window function
#	-off xxx : Starting point of sine-bell Offset in pi radians (0 = Sine, 0.5 = Cosine)
#	-end xxx : Ending point of sine-bell in pi radians. Default is 1
#	-pow xxx : Exponent of sine-bell. Typical values are 1 and 2.  Default is 1.


#   Baseline Corrections
# -fn BASE -- linear baseline correction
# -fn POLY -auto -- polynomial baseline correction.  May cause artifacts
# -fn POLY -nw __ -nl __ -avg -ord __ -- manually controlled polynomial baseline correction.


#   Other functions
# -fn PS -- Phase correction (1st: H-dim, 2nd: C-dim)
# -fn ZF -- Zero filling (Always doubles unless -zf 0, No significant gains using >2x zero filling)
#	-auto   : Rounds to the nearest power of 2.
#	-size   : Specify desired complex size of zero filling
#	-zf xxx : # times doubling size, use -zf 2 for RDC measurement
# -fn TP -- Transpose.  Switches processing to next dimension
# -fn LP -- Linear prediction.  Only used for qualitative and multi-D experiments
# -fn EXT -left -sw -- discards right portion of spectra
# -fn EXT -x1 xxx -xn xxx -sw -- selects region of spectrum to keep

