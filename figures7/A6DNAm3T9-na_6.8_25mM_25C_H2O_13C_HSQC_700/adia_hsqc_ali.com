#!/bin/csh

#
# Basic 2D Phase-Sensitive Processing:
#   Cosine-Bells are used in both dimensions.
#   Use of "ZF -auto" doubles size, then rounds to power of 2.
#   Use of "FT -auto" chooses correct Transform mode.
#   Imaginaries are deleted with "-di" in each dimension.
#   Phase corrections should be inserted by hand.

nmrPipe -in test.fid \
| nmrPipe  -fn  SOL      \
#| nmrPipe  -fn  POLY -time \
| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 1 -c 0.5    \
| nmrPipe  -fn ZF -size 16384                               \
| nmrPipe  -fn FT -auto                               \
| nmrPipe  -fn PS -p0 -48.30 -p1 -124.00 -di -verb         \
| nmrPipe  -fn EXT -xn 5ppm -x1 6.5ppm -sw -verb    \
| nmrPipe  -fn TP                                     \
| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 1 -c 0.5    \
| nmrPipe  -fn ZF -size 8192                               \
| nmrPipe  -fn FT -auto                               \
| nmrPipe  -fn PS -p0 111.00 -p1 0.00 -di -verb         \
   -ov -out test2.ft2
#pipe2ucsf test.ft2 A6_m3T9_adia_HSQC_Ali.ucsf