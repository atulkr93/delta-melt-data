#!/bin/csh

#
# Basic 2D Phase-Sensitive Echo-AntiEcho Processing:
#   Cosine-Bells are used in both dimensions.
#   Use of "ZF -auto" doubles size, then rounds to power of 2.
#   Use of "FT -auto" chooses correct Transform mode.
#   Imaginaries are deleted with "-di" in each dimension.
#   Phase corrections should be inserted by hand.

# -fn EXT -left -sw -- discards right portion of spectra
# POLY -- water suppression,  -neg after FT to reverse matrix

nmrPipe -in ./test.fid \
| nmrPipe  -fn SP -off 0.50 -end 1.00 -pow 2 -c 0.5 	    	\
| nmrPipe  -fn ZF -size 8192                             	\
| nmrPipe  -fn FT -auto                               		\
| nmrPipe  -fn EXT -x1 11ppm -xn 14.2ppm -sw \
| nmrPipe  -fn PS -p0 158.3 -p1 30.0 -di -verb         		\
| nmrPipe  -fn POLY -auto \
   -ov -out out.ft2
