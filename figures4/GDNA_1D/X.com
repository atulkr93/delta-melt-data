#!/bin/csh
# Basic 1D Processing:
#   POLY -time helps with supression of water peak artifacts as needed
#   GM is the invers exp width (Hz) -g1 to the Gaussian broaden width (Hz) -g2 and
#   and location of Gauss Maximum 0.0 - 1.0 (usually 0.0) for Lorentz-to-Gauss line broadening. 
#   Cosine-Bell
#   Use of "ZF -auto" doubles size, then rounds to power of 2, "ZF -size" makes data set
#   the indicated size, use a fourier number 4096, 8192 ...
#   Use of "FT -auto" chooses correct Transform mode.
#   Imaginaries are deleted with "-di" in each dimension.
#   Phase corrections should be inserted by hand.
#   Second polynomial correction for baseline
#   Extract region of spectrum between upper bound -x1 and lower bound -xn  
#

nmrPipe -in test.fid \
#| nmrPipe  -fn POLY -time \
#| nmrPipe  -fn GM -g1 5 -g2 15 -g3 0.0 -c 1.0 \
| nmrPipe  -fn SP -off 0.5 -end 0.98 -pow 1.0 -c 0.5   \
| nmrPipe  -fn ZF -size 16834  \
| nmrPipe  -fn FT    \
| nmrPipe  -fn PS -p0 138.0 -p1 0.0 -di -verb \
#| nmrPipe  -fn POLY -auto \
| nmrPipe  -fn EXT -x1 12ppm -xn 14ppm -sw \
 -ov -out test.ft2
 
 
